
#include "Data_Logistics.hpp"
#include "utility.hpp"


Data_Logistics::Data_Logistics(const Build_Params& build_params):
    params(build_params)
{}


/**
 * @brief 获取输入路径集合
 *
 * 返回包含输入路径的字符串向量。
 *
 * @return 输入路径的字符串向量
 */
const std::vector<std::string> Data_Logistics::input_paths_collection() const
{
    return params.sequence_input().seqs();
}


const std::string Data_Logistics::working_dir_path() const
{
    return dirname(params.output_prefix());
}


const std::string Data_Logistics::edge_db_path() const
{
#ifdef CF_DEVELOP_MODE
    if(!params.edge_db_path().empty())
        return params.edge_db_path();
#endif

    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::edges_ext;
}

/**
 * @brief 获取顶点数据库路径
 *
 * 返回顶点数据库的路径。在开发模式下，如果参数中指定了顶点数据库路径，则返回该路径；
 * 否则，返回工作目录路径加上输出文件名前缀，并添加顶点数据库的扩展名(cf_V)
 *
 * @return 顶点数据库路径
 */
const std::string Data_Logistics::vertex_db_path() const
{
#ifdef CF_DEVELOP_MODE
    if(!params.vertex_db_path().empty())
        return params.vertex_db_path();
#endif
    
    return params.working_dir_path() + filename(params.output_prefix()) + cuttlefish::file_ext::vertices_ext;
}


const std::string Data_Logistics::output_file_path() const
{
    return params.output_file_path();
}
