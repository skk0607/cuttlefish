
#include "Thread_Pool.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "CdBG.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"

#include <iostream>


template <uint16_t k>
/**
 * @brief 线程池构造函数
 *
 * 创建一个线程池对象，并初始化相关成员变量。
 * 在创建线程池的时候已经创建了消费者进程
 * @param thread_count 线程数量
 * @param dBG 指向线程池使用的数据结构的指针
 * @param task_type 任务类型
 */
Thread_Pool<k>::Thread_Pool(const uint16_t thread_count, void* const dBG, const Task_Type task_type):
    thread_count(thread_count),
    dBG(dBG),
    task_type(task_type),
    task_status(new volatile Task_Status[thread_count])//使用 new 关键字动态分配了一个 Task_Status 类型的数组，数组的大小由 thread_count 决定。volatile 关键字用于告诉编译器，这个数组的值可能会在程序的控制之外被改变，因此编译器不应该对这个数组的访问进行优化。
{
    // Mark the status of the task for each thread as `pending`.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
      task_status[t_id] = Task_Status::pending; // 每个线程都尚未提供任务

    // Resize the parameters collections.
    switch(task_type)
    {
    case Task_Type::classification:
        classify_params.resize(thread_count);
        break;
    
    case Task_Type::output_plain:
    case Task_Type::output_gfa:
    case Task_Type::output_gfa_reduced:
        output_params.resize(thread_count);
        break;

    case Task_Type::compute_states_read_space:
    case Task_Type::extract_unipaths_read_space:
        read_dBG_compaction_params.resize(thread_count);
        break;

    default:
        std::cerr << "Unrecognized task type encountered in thread pool. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // Launch the threads.
    // task: 每个线程的构造函数的输入函数参数
    // this: 作为参数传递给 emplace_back，确保了新创建的线程将使用当前对象的
    // task 函数实例。
    // t_id: 线程的 ID,作为task函数的参数
    // 由于成员函数属于类的对象，调用成员函数时需要一个隐式的 this
    // 指针来访问该对象的成员变量和成员函数。
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        thread_pool.emplace_back(&Thread_Pool::task, this, t_id);
}


template <uint16_t k>
/**
 * @brief 线程池任务函数
 *
 * 根据指定的线程ID执行线程池中的任务。
 *
 * @param thread_id 线程ID
 */
void Thread_Pool<k>::task(const uint16_t thread_id)
{
    while(true)
    {
        // Busy-wait for some task.
        while(task_status[thread_id] == Task_Status::pending);

        // No more tasks to come in on the future.
        if(task_status[thread_id] == Task_Status::no_more)
            return;


        // Some task is available for the thread number `thread_id`.
        if(task_status[thread_id] == Task_Status::available)
        {
            switch(task_type)
            {
            case Task_Type::classification:
                {
                    const Classification_Task_Params& params = classify_params[thread_id];
                    static_cast<CdBG<k>*>(dBG)->process_substring(params.seq, params.seq_len, params.left_end, params.right_end);
                }
                break;

            case Task_Type::output_plain:
                {
                    const Output_Task_Params& params = output_params[thread_id];
                    static_cast<CdBG<k>*>(dBG)->output_plain_off_substring(params.thread_id, params.seq, params.seq_len, params.left_end, params.right_end);
                }
                break;

            case Task_Type::output_gfa:
            case Task_Type::output_gfa_reduced:
                {
                    const Output_Task_Params& params = output_params[thread_id];
                    static_cast<CdBG<k>*>(dBG)->output_gfa_off_substring(params.thread_id, params.seq, params.seq_len, params.left_end, params.right_end);
                }
                break;
            
            case Task_Type::compute_states_read_space:
                {//线程池里所有线程任务都相同
                    const Read_dBG_Compaction_Params& params = read_dBG_compaction_params[thread_id];
                    static_cast<Read_CdBG_Constructor<k> *>(dBG)->process_edges(
                        static_cast<Kmer_SPMC_Iterator<k + 1> *>(params.parser),
                        params.thread_id);
                    // 将 void* 类型的指针 dBG 转换为 Read_CdBG_Constructor<k>*
                    // 类型的指针。这种转换是显式的，编译器会检查转换是否安全。
                    // 在src/Read_CdBG_Constructor.cpp 调用函数时传入的是this指针,也就是转换是安全的
                }
                break;

            case Task_Type::extract_unipaths_read_space:
                {
                    //取对应线程的参数
                    const Read_dBG_Compaction_Params& params = read_dBG_compaction_params[thread_id];
                    static_cast<Read_CdBG_Extractor<k>*>(dBG)->
                        process_vertices(static_cast<Kmer_SPMC_Iterator<k>*>(params.parser), params.thread_id);
                }
                break;
            }

            //任务做完了就立刻将线程的状态改为pending
            free_thread(thread_id);
        }
    }
}


template <uint16_t k>
/**
 * @brief 获取空闲线程
 *
 * 从线程池中获取一个处于空闲状态的线程ID。
 *
 * @return 空闲线程的ID，若线程池中没有空闲线程，则返回-1
 */
uint16_t Thread_Pool<k>::get_idle_thread() const
{
    int32_t idle_thread_id = -1;
    uint16_t t_id = 0;

    while(idle_thread_id == -1)
        if(task_status[t_id] == Task_Status::pending)//如果线程池里的线程有处于空闲状态的线程，就返回该线程的ID
            idle_thread_id = t_id;
        else
            t_id = (t_id + 1) % thread_count;

    
    return idle_thread_id;
}


template <uint16_t k>
void Thread_Pool<k>::get_thread(const uint16_t thread_id) const
{
    while(task_status[thread_id] != Task_Status::pending);
}


template <uint16_t k>
void Thread_Pool<k>::assign_classification_task(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end)
{
    classify_params[thread_id] = Classification_Task_Params(seq, seq_len, left_end, right_end);

    assign_task(thread_id);
}


template <uint16_t k>
void Thread_Pool<k>::assign_output_task(const uint16_t thread_id, const char* const seq, const size_t seq_len, const size_t left_end, const size_t right_end)
{
    output_params[thread_id] = Output_Task_Params(thread_id, seq, seq_len, left_end, right_end);

    assign_task(thread_id);
}


template <uint16_t k>
/**
 * @brief 分配读取dBG压缩任务
 *
 * 将读取dBG压缩任务分配给线程池中的指定线程。
 *
 * @param parser 解析器指针
 * @param thread_id 线程ID
 */
void Thread_Pool<k>::assign_read_dBG_compaction_task(void* const parser, const uint16_t thread_id)
{
    //给对应线程分配函数参数
    read_dBG_compaction_params[thread_id] = Read_dBG_Compaction_Params(parser, thread_id);

    assign_task(thread_id);
}


template <uint16_t k>
/**
 * @brief 分配任务给指定线程
 *
 * 将任务分配给指定线程。如果指定线程不是空闲状态，则输出错误信息并退出程序。
 *
 * @param thread_id 线程ID
 */
void Thread_Pool<k>::assign_task(const uint16_t thread_id)
{//如果线程不是处于等待中就报错
    if(task_status[thread_id] != Task_Status::pending)
    {
        std::cerr << "Expected thread " << thread_id << " to be idle while assigning a job, but found it unexpecedly busy. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    //表示线程状态是已经分配任务,待处理
    task_status[thread_id] = Task_Status::available;
}


template <uint16_t k>
/**
 * @brief 释放线程
 *
 * 将指定线程标记为待处理状态，释放线程以便再次使用。
 *
 * @param thread_id 线程ID
 */
void Thread_Pool<k>::free_thread(const uint16_t thread_id)
{
    if(task_status[thread_id] != Task_Status::available)
    {
        std::cerr << "Expected thread " << thread_id << " to be busy while assigning a job, but found it unexpecedly free. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    task_status[thread_id] = Task_Status::pending;
}


template <uint16_t k>
/**
 * @brief 等待线程池中的所有任务完成
 *
 * 等待线程池中的所有任务执行完成。该函数会阻塞当前线程，直到所有线程的任务状态都为非可用状态。
 */
void Thread_Pool<k>::wait_completion() const
{
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
        while(task_status[t_id] == Task_Status::available);
}


template <uint16_t k>
/**
 * @brief 关闭线程池
 *
 * 等待所有线程完成执行后，关闭线程池中的所有线程。
 */
void Thread_Pool<k>::close()
{
    // Wait for all the threads to finish.
    wait_completion();

    // Close all the threads.
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        // Signal each thread to stop running.
        task_status[t_id] = Task_Status::no_more;
        
        if(!thread_pool[t_id].joinable())
        {
            std::cerr << "Early termination of a worker thread encountered. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        thread_pool[t_id].join();
    }


    delete[] task_status;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Thread_Pool)
