
#include "globals.hpp"
#include "kmer_Enumerator.hpp"
#include "Kmer_Container.hpp"
#include "kmer_Enumeration_Stats.hpp"


template <uint16_t k> const std::size_t kmer_Enumerator<k>::min_memory;

/**
 * @brief 	 【简介】
 * 枚举 CdBG 中的边，并返回统计信息。
 * ip_type:输入文件类型
 * seqs: 输入序列,此为已经拼接为字符串的版本
 * cutoff: 枚举的阈值
 * thread_count: 线程数
 * max_memory: 最大内存
 * stric_memory: 是否使用严格内存
 * bits_per_vertex: 每个顶点占用的比特数
 * working_dir: 工作目录
 * ouput_path: 输出路径
 * @return 	 【解释返回值】
 * @warning	 【不可重入、阻塞等属性特殊说明】
 * @note	 【重大修改】
 */
template <uint16_t k>
kmer_Enumeration_Stats<k> kmer_Enumerator<k>::enumerate(
    const KMC::InputFileType input_file_type, const std::vector<std::string>& seqs, const uint32_t cutoff, const uint16_t thread_count,
    const std::size_t max_memory, const bool strict_memory, const bool estimate_mem_usage, const double bits_per_kmer,
    const std::string& working_dir_path, const std::string& output_db_path)
{
    // FunnyProgress progress;

    const bool estimate_mem = (k > small_k_threshold && estimate_mem_usage);   // Histogram estimation is not supported for small enough k's yet.
    std::size_t memory = std::max(max_memory, min_memory);
    stage1_params
        .SetInputFileType(input_file_type)
        .SetInputFiles(seqs)
        .SetKmerLen(k)
        .SetNThreads(thread_count)
        .SetTmpPath(working_dir_path)
        .SetEstimateHistogramCfg(estimate_mem ? KMC::EstimateHistogramCfg::ESTIMATE_AND_COUNT_KMERS : KMC::EstimateHistogramCfg::DONT_ESTIMATE)
        // .SetPercentProgressObserver(&progress)
    ;

    if(strict_memory)
        stage1_params
            .SetMaxRamGB(memory)
            .SetSignatureLen(signature_len)
            .SetNBins(bin_count)
        ;

    stage1_results = kmc.RunStage1(stage1_params);

    /**
    这段代码通过一系列的std::max函数调用和条件判断来确定memory的值。

    如果estimate_mem为true，则：

    首先计算solid_kmer_count_approx(cutoff)，这可能是一个近似地计算k-mer数量的函数，其中cutoff是某种阈值或参数。
    然后调用memory_limit函数，将上面得到的近似k-mer数量和bits_per_kmer作为参数，这个函数可能返回一个基于k-mer数量和每个k-mer所需的位数来估计的内存限制。
    使用std::max将memory_limit的结果和max_memory进行比较，取两者中的较大值。
    如果estimate_mem为false，则直接取max_memory作为std::max的第一个参数。

    最后，将上述std::max的结果与min_memory再次使用std::max进行比较，确保memory的值不会低于min_memory。
     */
    memory = std::max(
        (estimate_mem ? std::max(memory_limit(solid_kmer_count_approx(cutoff), bits_per_kmer), max_memory) : max_memory),
        min_memory);
    stage2_params
        .SetCutoffMin(cutoff)
        .SetNThreads(thread_count)
        .SetStrictMemoryMode(strict_memory)
#ifndef CF_VALIDATION_MODE
        .SetCounterMax(counter_max)
#endif
        .SetOutputFileName(output_db_path)
    ;

    if(strict_memory)
        stage2_params.SetMaxRamGB(memory);

    stage2_results = kmc.RunStage2(stage2_params);
    //计算 边的数据库有多大
    const std::size_t db_size = Kmer_Container<k>::database_size(output_db_path);
    // 返回一个边集合的对象
    return kmer_Enumeration_Stats<k>(stage1_results, stage2_results, memory, db_size);
}


template <uint16_t k>
uint64_t kmer_Enumerator<k>::solid_kmer_count_approx(const uint16_t cutoff) const
{
    uint64_t solid_kmer_count = 0;
    for (std::size_t freq = cutoff; freq < stage1_results.estimatedHistogram.size(); ++freq)
        solid_kmer_count += stage1_results.estimatedHistogram[freq];

    return solid_kmer_count;
}


template <uint16_t k>
std::size_t kmer_Enumerator<k>::memory_limit(const uint64_t unique_kmer_count, const double bits_per_kmer) const
{
    const double memory_in_bits = bits_per_kmer * unique_kmer_count;
    const double memory_in_bytes = memory_in_bits / 8.0;
    std::size_t memory_in_gb = static_cast<std::size_t>(memory_in_bytes / (1024 * 1024 * 1024));

    return memory_in_gb;
}



// Template instantiations for the required instances. 
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_ALL, kmer_Enumerator)
