
#include "Kmer_Container.hpp"
#include "Kmer_SPMC_Iterator.hpp"
#include "utility.hpp"
#include "globals.hpp"


template <uint16_t k>
/**
 * @brief Kmer_Container 构造函数
 *
 * 使用给定的 KMC 数据库文件路径创建 Kmer_Container 对象。
 * 在对应kmc_file_path路径里读取一些参数,然后存储到kmer_database_info里去,也就是该类的成员变量
 * @param kmc_file_path KMC 数据库文件路径
 */
Kmer_Container<k>::Kmer_Container(const std::string& kmc_file_path):
    kmc_file_path(kmc_file_path)
{
    CKMC_DB kmer_database;
    // 打开文件`.kmc_pre` & `.kmc_suf '，并将KMC DB参数读取到RAM。
   //std::cout<<"kmer_database 读取前:"<< kmer_database <<std::endl;
    if(!kmer_database.read_parameters(kmc_file_path))
    {
        std::cout << "Error opening KMC database files with prefix " << kmc_file_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    if(!kmer_database.Info(kmer_database_info))
    {
        std::cout << "Error reading from KMC database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    kmer_database.Close();


    const uint16_t kmer_len = kmer_length();
    if(kmer_len != k)
    {
        std::cerr << "Expected k value " << k << ", but is provided with a " << kmer_len << "-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
const std::string& Kmer_Container<k>::container_location() const
{
    return kmc_file_path;
}


template <uint16_t k>
uint32_t Kmer_Container<k>::kmer_length() const
{
    return kmer_database_info.kmer_length;
}


template <uint16_t k>
uint64_t Kmer_Container<k>::size() const
{
   return kmer_database_info.total_kmers;
}


template <uint16_t k>
uint64_t Kmer_Container<k>::size(const std::string& kmc_db_path)
{
    const Kmer_Container<k> kmer_container(kmc_db_path);
    return kmer_container.size();
}


template <uint16_t k>
bool Kmer_Container<k>::exists(const std::string& kmc_db_path)
{
    const std::string kmc_pref_file(kmc_db_path + ".kmc_pre");
    const std::string kmc_suff_file(kmc_db_path + ".kmc_suf");

    return file_exists(kmc_pref_file) && file_exists(kmc_suff_file);
}


template <uint16_t k>
/**
 * @brief 从指定路径中移除 KMC 数据库文件
 *
 * 从给定的路径中移除 KMC 数据库文件，包括前缀和后缀文件。
 *
 * @param kmc_db_path KMC 数据库文件的路径
 */
void Kmer_Container<k>::remove(const std::string& kmc_db_path)
{
    const std::string kmc_pref_file(kmc_db_path + ".kmc_pre");
    const std::string kmc_suff_file(kmc_db_path + ".kmc_suf");
    //如果存在1个删除失败，则返回错误
    if(!remove_file(kmc_pref_file) || !remove_file(kmc_suff_file))
    {
        std::cerr << "Error removing the KMC database file from path prefix " << kmc_db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
/**
 * @brief 获取 KMC 数据库大小
 *
 * 根据给定的 KMC 数据库前缀，计算 KMC 数据库的大小（以字节为单位）。
 *
 * @param kmc_db_prefix KMC 数据库前缀
 *
 * @return 返回 KMC 数据库的大小（以字节为单位）
 */
std::size_t Kmer_Container<k>::database_size(const std::string& kmc_db_prefix)
{
    const std::string kmc_pref_file(kmc_db_prefix + ".kmc_pre");
    const std::string kmc_suff_file(kmc_db_prefix + ".kmc_suf");

    const std::size_t pref_sz = file_size(kmc_pref_file);
    const std::size_t suff_sz = file_size(kmc_suff_file);

    if(!pref_sz || !suff_sz)
    {
        std::cerr << "Error computing size of KMC database at " << kmc_db_prefix << ". Possibly missing file(s). Aborting.\n";
        std::exit(EXIT_FAILURE);
    }


    return pref_sz + suff_sz;
}


// template <uint16_t k>
// typename Kmer_Container<k>::iterator Kmer_Container<k>::end() const
// {
//     return iterator(this, false);
// }

// template <uint16_t k>
// typename Kmer_Container<k>::buf_iterator Kmer_Container<k>::buf_begin() const
// {
//     return buf_iterator(this, true, false);
// }


// template <uint16_t k>
// typename Kmer_Container<k>::buf_iterator Kmer_Container<k>::buf_end() const
// {
//     return buf_iterator(this, false, true);
// }


template <uint16_t k>
/**
 * @brief 返回 SPMC 迭代器的起始位置
 *
 * 为给定的消费者数量返回 SPMC（Single Producer Multiple Consumer）迭代器的起始位置。
 *
 * @param consumer_count 消费者数量
 *
 * @return SPMC 迭代器的起始位置
 */
typename Kmer_Container<k>::spmc_iterator Kmer_Container<k>::spmc_begin(const size_t consumer_count) const
{
    return spmc_iterator(this, consumer_count);   
}


template <uint16_t k>
typename Kmer_Container<k>::spmc_iterator Kmer_Container<k>::spmc_end(const size_t consumer_count) const
{
    return spmc_iterator(this, consumer_count, false, true);
}



// Template instantiations for the required instances.
ENUMERATE(INSTANCE_COUNT, INSTANTIATE_ALL, Kmer_Container)
