
#include "Read_CdBG_Extractor.hpp"
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "Character_Buffer.hpp"
#include "Thread_Pool.hpp"


template <uint16_t k>
/**
 * @brief Read_CdBG_Extractor 构造函数
 *
 * 使用给定的构建参数和 kmer 哈希表来初始化 Read_CdBG_Extractor 对象。
 *
 * @param params 构建参数对象
 * @param hash_table kmer 哈希表对象
 */
Read_CdBG_Extractor<k>::Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table):
    params(params),
    hash_table(hash_table)
{}


template <uint16_t k>
/**
 * @brief 提取最大单元路径
 *
 * 从给定的顶点数据库路径中读取数据，并提取出最大的单元路径，然后将结果写入到输出文件路径中。
 *
 * @param vertex_db_path 顶点数据库路径
 * @param output_file_path 输出文件路径
 */
void Read_CdBG_Extractor<k>::extract_maximal_unitigs(const std::string& vertex_db_path, const std::string& output_file_path)
{
    std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();


    // Construct a thread pool.
    const uint16_t thread_count = params.thread_count();
    // 创建消费者线程池,同时thread_count个消费者线程也开始工作
    Thread_Pool<k> thread_pool(thread_count, this, Thread_Pool<k>::Task_Type::extract_unipaths_read_space);

    // Launch the reading (and parsing per demand) of the vertices from disk. 这里仅仅是读取顶点数据集的一些信息,并没有读取具体的边
    const Kmer_Container<k> vertex_container(vertex_db_path);  // Wrapper container for the vertex-database.
    Kmer_SPMC_Iterator<k> vertex_parser(&vertex_container, params.thread_count());  // Parser for the vertices from the vertex-database.
    std::cout << "Number of distinct vertices: " << vertex_container.size() << ".\n";
    //生产者线程开始
    vertex_parser.launch_production();

    // Clear the output file and initialize the output sink.
    // 清空输出文件并初始化输出接收器。
    clear_file(output_file_path);
    init_output_sink(output_file_path);

    // Launch (multi-threaded) extraction of the maximal unitigs.
    // 启动(多线程)提取maximal unitigs。
    // 计算每个线程需要处理的顶点的百分比
    const uint64_t thread_load_percentile = static_cast<uint64_t>(std::round((vertex_count() / 100.0) / params.thread_count()));
    // 启动进度追踪器，设置进度条的百分比
    // total_work = vertex_count() * 2, 总工作量是vertex_count()*2
    progress_tracker.setup(vertex_count() * 2, thread_load_percentile,
                            params.path_cover() ? "Extracting maximal path cover" :  "Extracting maximal unitigs");
    distribute_unipaths_extraction(&vertex_parser, thread_pool);

    // Wait for the vertices to be depleted from the database.
    vertex_parser.seize_production();

    // Wait for the consumer threads to finish parsing and processing edges.
    thread_pool.close();

    // Close the output sink.
    close_output_sink();

    std::cout << "\nNumber of scanned vertices: " << vertices_scanned << ".\n";
    unipaths_meta_info_.print();


    std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();
    double elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start).count();
    std::cout << "Extracted the paths. Time taken = " << elapsed_seconds << " seconds.\n";
}


template <uint16_t k>
/**
 * @brief 分配单路径提取任务
 *
 * 根据给定的顶点解析器和线程池，将单路径提取任务分配给空闲线程执行。
 *
 * @param vertex_parser 顶点解析器指针
 * @param thread_pool 线程池引用
 */
void Read_CdBG_Extractor<k>::distribute_unipaths_extraction(Kmer_SPMC_Iterator<k>* const vertex_parser, Thread_Pool<k>& thread_pool)
{
    const uint16_t thread_count = params.thread_count();

    for(uint16_t t_id = 0; t_id < thread_count; ++t_id)
    {
        //在线程池里找到pending的线程，将任务分配给空闲线程
        const uint16_t idle_thread_id = thread_pool.get_idle_thread();
        thread_pool.assign_read_dBG_compaction_task(vertex_parser, idle_thread_id);
    }
}


template <uint16_t k>
/**
 * @brief 处理顶点的函数
 *
 * 根据给定的顶点解析器和线程ID，处理顶点并提取相关的最大单元。
 *
 * @param vertex_parser 顶点解析器指针
 * @param thread_id 线程ID
 */
void Read_CdBG_Extractor<k>::process_vertices(Kmer_SPMC_Iterator<k>* const vertex_parser, const uint16_t thread_id)
{
    // 每个线程进入相同函数,局部变量都是不同的,相当于创建副本
    // Data structures to be reused per each vertex scanned.
    // 每个扫描的顶点都重用数据结构。
    Kmer<k> v_hat;  // The vertex copy to be scanned one-by-one.逐一扫描顶点副本。
    // 用于构建`v_hat`的最大单位的暂存空间。
    Maximal_Unitig_Scratch<k> maximal_unitig;  // The scratch space to be used to construct the containing maximal unitig of `v_hat`.
    // 此线程扫描的顶点数。
    uint64_t vertex_count = 0;  // Number of vertices scanned by this thread.
    // 该线程提取的Maximal_Unitig的元信息。
    Unipaths_Meta_info<k> extracted_unipaths_info;  // Meta-information over the maximal unitigs extracted by this thread.
    // 线程扫描的顶点数;重置为其近似工作量的1%。
    uint64_t progress = 0;  // Number of vertices scanned by the thread; is reset at reaching 1% of its approximate workload.

    Character_Buffer<BUFF_SZ, sink_t> output_buffer(output_sink.sink());  // The output buffer for maximal unitigs.


    while(vertex_parser->tasks_expected(thread_id))
    //value_at: 读取一条kmer
        if(vertex_parser->value_at(thread_id, v_hat))
        {
            if(extract_maximal_unitig(v_hat, maximal_unitig))
            {
                //标记前向和后向的状态为已经输出
                mark_maximal_unitig(maximal_unitig);

                extracted_unipaths_info.add_maximal_unitig(maximal_unitig);
                // output_buffer += maximal_unitig.fasta_rec();
                maximal_unitig.add_fasta_rec_to_buffer(output_buffer);

                if(progress_tracker.track_work(progress += maximal_unitig.size()))
                    progress = 0;
            }

            vertex_count++;//每个线程各自统计顶点数
            if(progress_tracker.track_work(++progress))
                progress = 0;
        }


    // Aggregate the meta-information over the extracted maximal unitigs and the thread-executions.
    lock.lock();

    vertices_scanned += vertex_count;
    unipaths_meta_info_.aggregate(extracted_unipaths_info);

    lock.unlock();
}


template <uint16_t k>
/**
 * @brief 初始化输出接收器
 *
 * 使用给定的输出文件路径初始化输出接收器。
 *
 * @param output_file_path 输出文件路径
 */
void Read_CdBG_Extractor<k>::init_output_sink(const std::string& output_file_path)
{
    output_sink.init_sink(output_file_path);
}


template <uint16_t k>
void Read_CdBG_Extractor<k>::close_output_sink()
{
    output_sink.close_sink();
}


template <uint16_t k>
const Build_Params& Read_CdBG_Extractor<k>::get_params() const
{
    return params;
}


template <uint16_t k>
const Unipaths_Meta_info<k>& Read_CdBG_Extractor<k>::unipaths_meta_info() const
{
    return unipaths_meta_info_;
}


template <uint16_t k>
uint64_t Read_CdBG_Extractor<k>::vertex_count() const
{
    return hash_table.size();
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, Read_CdBG_Extractor)
