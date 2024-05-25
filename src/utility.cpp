
#include "utility.hpp"

#include <cstring>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <sstream>
#include <iterator>
#include <filesystem>
#include <fstream>
#include <cstdio>


std::string get_random_string(const size_t len, const char* const alphabet)
{
    std::string str;
    str.reserve(len);

    const unsigned int seed = static_cast<unsigned int>(std::time(NULL));
    std::srand(seed);
    for (size_t i = 0; i < len; ++i)
        str += alphabet[(std::rand() % (sizeof(alphabet) - 1))];

    return str;
}


bool is_prefix(const std::string& s, const std::string& pref)
{
    if(s.length() < pref.length())
        return false;

    size_t idx = 0;
    for(; idx < pref.length() && s[idx] == pref[idx]; ++idx);

    return idx == pref.length();
}


bool file_exists(const std::string& file_path)
{
    return std::filesystem::exists(file_path);
}


bool dir_exists(const std::string& dir_path)
{
    return std::filesystem::is_directory(dir_path);
}


std::size_t file_size(const std::string& file_path)
{
    std::error_code ec;
    const uintmax_t size = std::filesystem::file_size(file_path, ec);
    return ec ? 0 : static_cast<std::size_t>(size);
}


bool file_prefix_exists(const std::string& path, const std::string& prefix)
{
    for(const auto& entry: std::filesystem::directory_iterator(path))
        if(is_prefix(filename(entry.path()), prefix))
            return true;

    return false;
}


std::string remove_whitespaces(const char* s)
{
    std::string str;
    str.reserve(strlen(s));

    for(const char* p = s; *p; ++p)
        if(!std::isspace(*p))
            str += *p;

    return str;
}

/**
 * @brief 拼接字符串
 *
 * 使用指定的分隔符将字符串向量中的字符串拼接成一个字符串，并返回结果。
 * 返回一个由s中所有字符串拼接而成的字符串，其中字符串之间用delimiter分隔。
 * @param s 字符串向量
 * @param delimiter 分隔符
 *
 * @return 拼接后的字符串
 */
const std::string concat_strings(const std::vector<std::string>& s, const std::string& delimiter)
{
  std::ostringstream concat_stream;
  // 使用std::copy函数将s中的字符串复制到concat_stream中。这里使用了std::ostream_iterator作为目标迭代器，它会将每个元素（在这里是std::string）插入到concat_stream中，并在每个元素之间插入delimiter。
  std::copy(
      s.begin(), s.end(),
      std::ostream_iterator<std::string>(concat_stream, delimiter.c_str()));
  // 从concat_stream中提取字符串并保存到concat_str中。
  std::string concat_str(concat_stream.str());
  // 删除concat_str末尾的delimiter，因为在步骤4中，最后一个字符串后面也加上了delimiter。
  concat_str.erase(concat_str.size() - delimiter.size(), delimiter.size());
  return concat_str;
}


/**
 * @brief 删除文件
 *
 * 删除指定路径下的文件。
 * std::filesystem 是 C++17 引入的一个库，用于处理文件系统和目录路径。
 * std::filesystem::remove 函数会尝试删除 file_path 指定的文件或目录。如果删除成功，它会返回 true；如果删除失败（例如，文件不存在或没有删除权限），它会返回 false。
 * @param file_path 文件路径
 *
 * @return 如果删除成功，返回 true；否则返回 false
 */
bool remove_file(const std::string& file_path)
{
    return std::filesystem::remove(file_path);
}


/**
 * @brief 清空文件内容
 *
 * 打开指定路径的文件，并将其内容清空。
 *
 * @param file_path 文件路径
 */
void clear_file(const std::string& file_path)
{
    // 打开文件，以输出和截断模式
    /**
    std::ofstream::out: 表示以输出模式打开文件。这是写入文件的默认模式。
    std::ofstream::trunc: 如果文件已经存在，则将其长度截断为0，即删除其所有内容。如果文件不存在，则创建它。
    */
    std::ofstream file(file_path.c_str(), std::ofstream::out | std::ofstream::trunc);

    // 判断文件打开是否失败
    if(file.fail())
    {
        // 输出错误信息并退出程序
        std::cerr << "Error opening file " << file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    // 关闭文件
    file.close();
}


const std::string filename(const std::string& file_path)
{
    return std::filesystem::path(file_path).filename().string();
}


const std::string dirname(const std::string& file_path)
{
    const std::string path = std::filesystem::path(file_path).remove_filename().string();
    return path.empty() ? "." : path;
}


void move_file(const std::string& from_path, const std::string& to_path)
{
    std::filesystem::copy(from_path, to_path);
    std::filesystem::remove(from_path);
}


/**
 * @brief 获取进程峰值内存使用量
 *
 * 读取进程状态文件，获取进程的峰值内存使用量（以字节为单位）。
 *
 * @return 进程的峰值内存使用量（以字节为单位），若读取文件失败则返回0
 */
std::size_t process_peak_memory()
{
    
    constexpr const char* process_file = "/proc/self/status";
    constexpr const char* peak_mem_field = "VmHWM:";
    const std::size_t field_len = std::strlen(peak_mem_field);

    std::FILE* fp = std::fopen(process_file, "r");
    if(fp == NULL)
    {
        std::cerr << "Error opening the process information file.\n";
        return 0;
    }

    char line[1024];
    std::size_t peak_mem = 0;
    while(std::fgets(line, sizeof(line) - 1, fp))
        if(std::strncmp(line, peak_mem_field, field_len) == 0)
        {
          // 将后面的字符串转换为无符号长整数，并存储在 peak_mem 中。
          peak_mem = std::strtoul(line + field_len, NULL, 0);
          break;
        }

    
    if(std::ferror(fp))
    {
        std::cerr << "Error reading the process information file.\n";
        return 0;
    }


    return peak_mem * 1024;
}
