
#include "Read_CdBG_Constructor.hpp"
#include "Edge.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Thread_Pool.hpp"

#include <chrono>


template <uint16_t k>
Read_CdBG_Constructor<k>::Read_CdBG_Constructor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table)
{}


template <uint16_t k>
/**
 * @brief 计算 DFA 状态
 *
 * 根据给定的边缘数据库路径，计算 DFA（确定有限自动机）的状态。
 *
 * @param edge_db_path 边缘数据库路径
 */
void Read_CdBG_Constructor<k>::compute_DFA_states(const std::string& edge_db_path)
{
    // std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    const Kmer_Container<k + 1> edge_container(edge_db_path);  // Wrapper container for the edge-database.
    // 创建读取 edge 集合的解析器
    Kmer_SPMC_Iterator<k + 1> edge_parser(&edge_container, params.thread_count());  // Parser for the edges from the edge-database.
    edge_count_ = edge_container.size();//刚好是edge uniuqe的数量
    std::cout << "Total number of distinct edges: " << edge_count_ << ".\n";
    // 输入的桶文件路径是 data/output/ceil.cf_hb
    const std::string &buckets_file_path = params.buckets_file_path();
    std::cout << "输入的桶文件路径是 " << params.buckets_file_path() << ".\n";
    if(!buckets_file_path.empty() && file_exists(buckets_file_path))    // The serialized hash table buckets, saved from some earlier execution, exists.
    {//检查有没有已经存储好的桶文件，如果存在，则直接加载
        std::cout <<    "Found the hash table buckets at file " << buckets_file_path << ".\n"
                        "Loading the buckets.\n";
        hash_table.load_hash_buckets(buckets_file_path);
        std::cout << "Loaded the buckets into memory.\n";
    }
    else
    {
        // Construct a thread pool.
        const uint16_t thread_count = params.thread_count();
        // 创建 thread_count 数量线程的线程池, 任务类型是
        // compute_states_read_space
        // 根据不同任务创建不同线程池
        Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::compute_states_read_space);

        // Launch the reading (and parsing per demand) of the edges from
        // disk.从磁盘启动边缘读取(并按需解析)。
        edge_parser.launch_production();

        // Launch (multi-threaded) computation of the states.
        // 启动(多线程)状态计算。
        // 计算每个线程应该处理的 edge_count_ 的百分数
        const uint64_t thread_load_percentile = static_cast<uint64_t>(std::round((edge_count_ / 100.0) / params.thread_count()));
        progress_tracker.setup(edge_count_, thread_load_percentile, "Computing DFA states");
        distribute_states_computation(&edge_parser, thread_pool);

        // Wait for the edges to be depleted from the database.
        // 等待数据库中的边缘耗尽。
        edge_parser.seize_production();

        // Wait for the consumer threads to finish parsing and processing the edges.
        thread_pool.close();

        std::cout << "\nNumber of processed edges: " << edges_processed << "\n";


        // Save the hash table buckets, if a file path is provided.
        if(params.save_buckets())
        {
            hash_table.save_hash_buckets(buckets_file_path); 
            std::cout << "Saved the hash buckets at " << buckets_file_path << "\n";
        }
    }


    // std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    // double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    // std::cout << "Done computing the DFA states. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
/**
 * @brief 分发状态计算
 *
 * 使用给定的 Kmer_SPMC_Iterator 和 Thread_Pool 对象，分发状态计算任务到线程池中。
 *
 * @param edge_parser Kmer_SPMC_Iterator 对象的指针，用于解析边信息
 * @param thread_pool Thread_Pool 对象的引用，用于执行计算任务
 */
void Read_CdBG_Constructor<k>::distribute_states_computation(Kmer_SPMC_Iterator<k + 1>* const edge_parser, Thread_Pool<k>& thread_pool)
{
    // 获取线程池中的线程数量
    const uint16_t thread_count = params.thread_count();

    // 遍历所有线程
    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        // 获取空闲的线程ID
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        // 将读取压缩dBG的任务分配给空闲的线程
        // 每个线程都是同样的Kmer_SPMC_Iterator<k + 1>
        // 这里本质的工作是让对应线程的状态改为available
        thread_pool.assign_read_dBG_compaction_task(edge_parser, idle_thread_id);
    }
}


template <uint16_t k>
/**
 * @brief 处理边信息
 *
 * 根据给定的参数，处理 CdBG（Compressed De Bruijn Graph）的边信息。
 *
 * @tparam k Kmer 长度
 *
 * @param edge_parser 指向 Kmer_SPMC_Iterator<k + 1> 类型的常量指针，用于解析边信息
 * @param thread_id 线程 ID
 */
void Read_CdBG_Constructor<k>::process_edges(Kmer_SPMC_Iterator<k + 1>* const edge_parser, const uint16_t thread_id)
{
    if(params.path_cover())
        process_path_cover_edges(edge_parser, thread_id);
    else
        process_cdbg_edges(edge_parser, thread_id);
}


template <uint16_t k>
/**
 * @brief 处理 CdBG 图的边
 *
 * 遍历输入的边解析器，对每一条边进行处理。
 *
 * @param edge_parser 边解析器指针
 * @param thread_id 线程 ID
 */
void Read_CdBG_Constructor<k>::process_cdbg_edges(Kmer_SPMC_Iterator<k + 1>* const edge_parser, const uint16_t thread_id)
{
    // Data locations to be reused per each edge processed.
    // 每个处理的边都要重用的数据位置。
    Edge<k> e;  // For the edges to be processed one-by-one.
/*
    cuttlefish::edge_encoding_t e_front, e_back;    // Edges incident to the front and to the back of a vertex with a crossing loop.
    cuttlefish::edge_encoding_t e_u_old, e_u_new;   // Edges incident to some particular side of a vertex `u`, before and after the addition of a new edge.
    cuttlefish::edge_encoding_t e_v_old, e_v_new;   // Edges incident to some particular side of a vertex `v`, before and after the addition of a new edge.
*/
    printf("当前正在处理数据的消费者线程是%d\n",thread_id);
    uint64_t edge_count = 0;    // Number of edges processed by this thread. 该线程处理的边数。
    uint64_t progress = 0;  // Number of edges processed by the thread; is reset at reaching 1% of its approximate workload. 线程处理的边数;重置为其近似工作量的1%。

    //存在线程不是在 no_more
    while (edge_parser->tasks_expected(thread_id))
      //每次只读1个 kmer,如果当前线程读过的kmer = 所能读的最大kmer则 else
        if(edge_parser->value_at(thread_id, e.e()))
        {
            //这里基本是根据 边获取 prefix 和 suffix,获取prefix和suffix的canonical形式,然后存储在哈希表中
            e.configure(hash_table);    // A new edge (k + 1)-mer has been parsed; set information for its two endpoints.

            if(e.is_loop())
                if(e.u().side() != e.v().side())    // It is a crossing loop.
                {
                    while(!add_crossing_loop(e.u()));
/*
                    while(!add_crossing_loop(e.u(), e_front, e_back));
                    
                    propagate_discard(e.u(), e.u().side() == cuttlefish::side_t::front ? e_front : e_back);
                    propagate_discard(e.v(), e.v().side() == cuttlefish::side_t::front ? e_front : e_back);
*/
                }
                else    // A one-sided loop.
                {
                    while(!add_one_sided_loop(e.u()));
/*
                    while(!add_one_sided_loop(e.u(), e_u_old));

                    propagate_discard(e.u(), e_u_old);
*/
                }
            else// It connects two endpoints `u` and `v` of two distinct vertex.
            {
                // 这里类似原子操作,只要更新状态
                while(!add_incident_edge(e.u()));
                while(!add_incident_edge(e.v()));
/*
                while(!add_incident_edge(e.u(), e_u_old, e_u_new));
                while(!add_incident_edge(e.v(), e_v_old, e_v_new));

                if(e_u_new == cuttlefish::edge_encoding_t::N)
                    propagate_discard(e.u(), e.v(), e_u_old);

                if(e_v_new == cuttlefish::edge_encoding_t::N)
                    propagate_discard(e.v(), e.u(), e_v_old);
*/
            }

            edge_count++;
            if(progress_tracker.track_work(++progress))
                progress = 0;
        }

    
    lock.lock();
    edges_processed += edge_count;//共享变量加锁
    lock.unlock();
}


template <uint16_t k>
void Read_CdBG_Constructor<k>::process_path_cover_edges(Kmer_SPMC_Iterator<k + 1>* const edge_parser, const uint16_t thread_id)
{
    Edge<k> e;  // For the edges to be processed one-by-one; say this is between the vertices `u` and `v`.

    uint64_t edge_count = 0;    // Number of edges processed by this thread.
    uint64_t progress = 0;  // Number of edges processed by the thread; is reset at reaching 1% of its approximate workload.


    while(edge_parser->tasks_expected(thread_id))
        if(edge_parser->value_at(thread_id, e.e()))
        {
            e.configure(hash_table);    // A new edge (k + 1)-mer has been parsed; set information for its two endpoints.

            if(e.is_loop())
                continue;
            else    // It connects two endpoints `u` and `v` of two distinct vertex.
                add_path_cover_edge(e);

            edge_count++;
            if(progress_tracker.track_work(++progress))
                progress = 0;
        }

    
    lock.lock();
    edges_processed += edge_count;
    lock.unlock();
}


template <uint16_t k>
uint64_t Read_CdBG_Constructor<k>::vertex_count() const
{
    return hash_table.size();
}


template <uint16_t k>
uint64_t Read_CdBG_Constructor<k>::edge_count() const
{
    return edge_count_;
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Constructor)
