
#include "dBG_Info.hpp"
#include "Read_CdBG_Constructor.hpp"
#include "Read_CdBG_Extractor.hpp"
#include "CdBG.hpp"
#include "Unipaths_Meta_info.hpp"
#include "Build_Params.hpp"
#include "utility.hpp"

#include <iomanip>
#include <fstream>
#include <iostream>


template <uint16_t k>
dBG_Info<k>::dBG_Info(const std::string& file_path):
    file_path_(file_path)
{
    if(file_exists(file_path))
        load_from_file();
}


template <uint16_t k>
std::string dBG_Info<k>::file_path() const
{
    return file_path_;
}


template <uint16_t k>
void dBG_Info<k>::load_from_file()
{
    std::ifstream input(file_path_.c_str());
        
    input >> dBg_info;

    if(input.fail())
    {
        std::cerr << "Error loading JSON object from file " << file_path_ << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    input.close();
}


template <uint16_t k>
/**
 * @brief 添加基础信息
 *
 * 将给定的 Read_CdBG_Constructor 对象的基本信息添加到 dBG_Info 对象中。
 *
 * @tparam k 模板参数
 *
 * @param cdbg_constructor Read_CdBG_Constructor 对象，提供基础信息
 */
void dBG_Info<k>::add_basic_info(const Read_CdBG_Constructor<k>& cdbg_constructor)
{
    dBg_info[basic_field]["vertex count"] = cdbg_constructor.vertex_count();
    dBg_info[basic_field]["edge count"] = cdbg_constructor.edge_count();
    std::cout << "basic info added: edge count "
              << cdbg_constructor.edge_count() << ",vertex count"<< cdbg_constructor.vertex_count() << std::endl;
}


template <uint16_t k>
void dBG_Info<k>::add_basic_info(const CdBG<k>& cdbg)
{
    dBg_info[basic_field]["vertex count"] = cdbg.vertex_count();
}


template <uint16_t k>
void dBG_Info<k>::add_short_seqs_info(const std::vector<std::pair<std::string, std::size_t>>& short_seqs)
{
    dBg_info[short_seqs_field] = short_seqs;
}


template <uint16_t k>
/**
 * @brief 添加单元路径信息
 *
 * 将给定的单元路径元信息添加到当前对象的 dBg_info 字典中。
 *
 * @param unipaths_info 单元路径元信息对象
 */
void dBG_Info<k>::add_unipaths_info(const Unipaths_Meta_info<k>& unipaths_info)
{
    // 添加最大单元格数量信息
    dBg_info[contigs_field]["maximal unitig count"] = unipaths_info.unipath_count();
    // 添加最大单元格中的顶点数量信息
    dBg_info[contigs_field]["vertex count in the maximal unitigs"] = unipaths_info.kmer_count();
    // 添加最短最大单元格长度信息
    dBg_info[contigs_field]["shortest maximal unitig length"] = unipaths_info.min_len();
    // 添加最长最大单元格长度信息
    dBg_info[contigs_field]["longest maximal unitig length"] = unipaths_info.max_len();
    // 添加最大单元格长度的总和信息
    dBg_info[contigs_field]["sum maximal unitig length"] = unipaths_info.sum_len();
    // 添加最大单元格长度的平均值信息
    dBg_info[contigs_field]["avg. maximal unitig length"] = unipaths_info.avg_len();
    // 添加注释信息，说明长度是以碱基为单位的
    dBg_info[contigs_field]["_comment"] = "lengths are in bases";
}


template <uint16_t k>
void dBG_Info<k>::add_unipaths_info(const Read_CdBG_Extractor<k>& cdbg_extractor)
{
    const Unipaths_Meta_info<k>& unipaths_info = cdbg_extractor.unipaths_meta_info();
    add_unipaths_info(unipaths_info);

    dBg_info[dcc_field]["DCC count"] = unipaths_info.dcc_count();
    if(unipaths_info.dcc_count() > 0)
    {
        dBg_info[dcc_field]["vertex count in the DCCs"] = unipaths_info.dcc_kmer_count();
        dBg_info[dcc_field]["sum DCC length (in bases)"] = unipaths_info.dcc_sum_len();
    }
}


template <uint16_t k>
void dBG_Info<k>::add_unipaths_info(const CdBG<k>& cdbg)
{
    const Unipaths_Meta_info<k>& unipaths_info = cdbg.unipaths_meta_info();
    add_unipaths_info(unipaths_info);
}


template <uint16_t k>
/**
 * @brief 添加构建参数
 *
 * 将给定的构建参数添加到dBG信息中。
 *
 * @tparam k k值类型
 *
 * @param params 构建参数对象
 */
void dBG_Info<k>::add_build_params(const Build_Params& params)
{
  // 嵌套字典, dg_info[params_field][key] = value
  // 这里将输入序列 都被拼接成1个字符串,序列之间放置来分隔符
    dBg_info[params_field]["input"] = concat_strings(params.sequence_input().seqs());
    dBg_info[params_field]["k"] = params.k();
    dBg_info[params_field]["output prefix"] = params.output_prefix();
}


template <uint16_t k>
void dBG_Info<k>::dump_info() const
{
    std::ofstream output(file_path_.c_str());
    output << std::setw(4) << dBg_info << "\n"; // Pretty-print the JSON wrapper with overloaded `std::setw`.

    if(output.fail())
    {
        std::cerr << "Error writing to the information file " << file_path_ << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    output.close();

    std::cout << "\nStructural information for the de Bruijn graph is written to " << file_path_ << ".\n";
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE, dBG_Info)
