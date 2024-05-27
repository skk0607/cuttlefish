
#ifndef KMER_SPMC_ITERATOR_HPP
#define KMER_SPMC_ITERATOR_HPP



#include "Kmer.hpp"
#include "Kmer_Container.hpp"
#include "kmc_api/kmc_file.h"

#include <cstdint>
#include <cstddef>
#include <memory>
#include <vector>
#include <string>
#include <thread>


// Data required by the consumers to correctly parse raw binary k-mers.
/**
 * @brief 消费者数据结构体
 *
 * 该结构体用于存储消费者数据，包括k-mer的后缀缓冲区、当前缓冲区中的k-mer数量、已解析的k-mer数量、k-mer的前缀缓冲区以及指向要开始解析的k-mer前缀的迭代器。
 *
 * 结构体使用L1缓存行大小对齐，以避免false-sharing。
 */
struct alignas(L1_CACHE_LINE_SIZE)
    Consumer_Data
{
    uint8_t* suff_buf{nullptr}; // Buffer for the raw binary suffixes of the k-mers.
    uint64_t kmers_available;   // Number of k-mers present in the current buffer.
    uint64_t kmers_parsed;      // Number of k-mers parsed from the current buffers.
    std::vector<std::pair<uint64_t, uint64_t>> pref_buf;    // Buffer for the raw binary prefixes of the k-mers, in the form: <prefix, #corresponding_suffix>
    // 指向开始解析k-mers的前缀的指针。
    std::vector<std::pair<uint64_t, uint64_t>>::iterator pref_it;   // Pointer to the prefix to start parsing k-mers from.
    // uint64_t pad_[1];           // Padding to avoid false-sharing.
};

// An "iterator" class to iterate over a k-mer database on disk, where a single producer thread
// (sequentially) reads the raw binary representations of the k-mers from disk, and a number of
// different consumer threads fetch (and parse) the raw binary k-mers.
// Note: in a technical sense, it's not an iterator.
template <uint16_t k>
class Kmer_SPMC_Iterator
{
    typedef Kmer_SPMC_Iterator iterator;


private:

    const Kmer_Container<k>* const kmer_container;  // The associated k-mer container over which to iterate.
    CKMC_DB kmer_database; // The k-mer database object.

    const uint64_t kmer_count;  // Number of k-mers present in the underlying database.
    const size_t consumer_count;  // Total number of consumer threads of the iterator.

    uint64_t kmers_read;    // Number of raw k-mers read (off disk) by the iterator.

    std::unique_ptr<std::thread> reader{nullptr};   // The thread doing the actual disk-read of the binary data, i.e. the producer thread.

    static constexpr size_t BUF_SZ_PER_CONSUMER = (1 << 24);   // Size of the consumer-specific buffers (in bytes): 16 MB.

    std::vector<Consumer_Data> consumer;   // Parsing data required for each consumer.

    // Status of the tasks for each consumer thread.
    enum class Task_Status: uint8_t
    {
        pending,    // k-mers yet to be provided;
        available,  // k-mers are available and waiting to be parsed and processed;
        no_more,    // no k-mers will be provided anymore.
    };

    volatile Task_Status* task_status{nullptr}; // Collection of the task statuses of the consumers.


    // Opens the k-mer database file with the path prefix `db_path`.
    void open_kmer_database(const std::string& db_path);

    // Closes the k-mer database file.
    void close_kmer_database();

    // Reads raw binary k-mer representations from the underlying k-mer database, and
    // makes those available for consumer threads. Reading continues until the database
    // has been depleted.
    void read_raw_kmers();

    // Returns the id (number) of an idle consumer thread.
    size_t get_idle_consumer() const;


public:

    // Constructs an iterator for the provided container `kmer_container`, on either
    // its beginning or its ending position, based on `at_begin` and `at_end`. The
    // iterator is to support `consumer_count` number of different consumers.
    Kmer_SPMC_Iterator(const Kmer_Container<k>* kmer_container, size_t consumer_count, bool at_begin = true, bool at_end = false);

    // Copy constructs an iterator from another one `other`.
    // Note: this should be prohibited, like the `operator=`. But the BBHash code
    // requires this to be implemented.
    Kmer_SPMC_Iterator(const iterator& other);

    // Destructs the iterator.
    ~Kmer_SPMC_Iterator();

    // Prohibits assignment-copying. This is a complex object with a background disk-reader
    // thread and a KMC database object in some *arbitrary* state. These should not be allowed
    // to be copied. Besides, there is a number of constant fields for the iterator.
    iterator& operator=(const iterator& rhs) = delete;

    // Tries to fetch and parse the next k-mer for the consumer with id `consumer_id` into `kmer`.
    // Returns `true` iff it's successful, i.e. k-mers were remaining for this consumer.
    bool value_at(size_t consumer_id, Kmer<k>& kmer);

    // Returns `true` iff this and `rhs` — both the iterators refer to the same container and
    // the same number of raw k-mers have been read (from disk) for both.
    bool operator==(const iterator& rhs) const;

    // Returns `true` iff the iterators, this and `rhs` — either they refer to different containers,
    // or a different number of raw k-mers have been read (from disk) for them.
    bool operator!=(const iterator& rhs) const;

    // Launches the background disk-read of raw binary k-mers.
    void launch_production();

    // Whether production has been launched yet — and more specifically, whether the production data
    // structures have been instantiated yet. This must be found true before attempting any sort of
    // access into the data structure. The only exception is the `launch_production` invokation.
    bool launched() const;

    // Waits for the disk-reads of the raw k-mers to be completed, and then waits for the consumers
    // to finish their ongoing tasks; then signals them that no more data are to be provided, and
    // also closes the k-mer database. 
    void seize_production();

    // Returns `true` iff tasks might be provided to the consumer with id `consumer_id` in future.
    bool tasks_expected(size_t consumer_id) const;

    // Returns `true` iff a task is available for the consumer with id `consumer_id`.
    bool task_available(size_t consumer_id) const;

    // Returns the memory (in bytes) used by the iterator.
    std::size_t memory() const;

    // Returns the memory (in bytes) to be used by an iterator supporting `consumer_count` consumers.
    static std::size_t memory(std::size_t consumer_count);

    // Dummy methods.
    const iterator& operator++() { return *this; }
    Kmer<k> operator*() { return Kmer<k>(); }
};


template <uint16_t k>
/**
 * @brief Kmer_SPMC_Iterator 构造函数
 *
 * 构造一个 Kmer_SPMC_Iterator 对象，用于迭代 Kmer_Container 中的 k-mer 序列。
 *
 * @param kmer_container 指向 Kmer_Container 的常量指针
 * @param consumer_count 消费者数量
 * @param at_begin 是否从迭代器开始位置构造
 * @param at_end 是否从迭代器结束位置构造
 */
inline Kmer_SPMC_Iterator<k>::Kmer_SPMC_Iterator(const Kmer_Container<k>* const kmer_container, const size_t consumer_count, const bool at_begin, const bool at_end):
    kmer_container(kmer_container),
    kmer_count{kmer_container->size()},
    consumer_count{consumer_count},
    kmers_read{at_end ? kmer_count : 0}
{
  // ^: 不同为1相同为0
  // 检查at_begin和at_end是否只有一个为真
  if (!(at_begin ^ at_end)) {
    std::cerr << "Invalid position provided for SPMC k-mer iterator "
                 "construction. Aborting.\n";
    std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::Kmer_SPMC_Iterator(const iterator& other):
    kmer_container(other.kmer_container),
    kmer_count{other.kmer_count},
    consumer_count{other.consumer_count},
    kmers_read{other.kmers_read}
{}


template <uint16_t k>
inline Kmer_SPMC_Iterator<k>::~Kmer_SPMC_Iterator()
{
    if(task_status != nullptr)
    {
        delete[] task_status;

        for(size_t id = 0; id < consumer_count; ++id)
           delete[] consumer[id].suff_buf;

        std::cerr << "\nCompleted a pass over the k-mer database.\n";
    }
}


template <uint16_t k>
/**
 * @brief 打开 k-mer 数据库
 *
 * 打开指定的 k-mer 数据库文件，并进行初始化操作。
 *
 * @param db_path 数据库文件路径
 */
inline void Kmer_SPMC_Iterator<k>::open_kmer_database(const std::string& db_path)
{
    if(!kmer_database.open_for_cuttlefish_listing(db_path))
    {
        std::cerr << "Error opening k-mer database with prefix " << db_path << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
inline void Kmer_SPMC_Iterator<k>::close_kmer_database()
{
    if(!kmer_database.Close())
    {
        std::cerr << "Error closing k-mer database. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
}


template <uint16_t k>
/**
 * @brief 启动生产函数
 *
 * 如果已经启动，则直接返回。
 * 初始化缓冲区以及解析数据结构。
 *
 * @tparam k Kmer 的长度
 */
inline void Kmer_SPMC_Iterator<k>::launch_production()
{
    //生产者线程是否已经初始化
    //实际就是检查是否已经创建了reader指针
    if(launched())
        return;

    // Initialize the buffers and the parsing data structures.
    // 初始化缓冲区和解析数据结构。
    // 每个消费者线程创建对应状态,构成数组
    // 线程的状态读取需要从内存读取,不允许优化
    task_status = new volatile Task_Status[consumer_count];
    // 将consumer vector的容量设置为consumer_count
    // 里面为空值
    consumer.resize(consumer_count);
    // 初始化consumer 的 消费者数据
    for(size_t id = 0; id < consumer_count; ++id)
    {
        auto& consumer_state = consumer[id];
        consumer_state.suff_buf = new uint8_t[BUF_SZ_PER_CONSUMER];
        consumer_state.kmers_available = 0;
        consumer_state.kmers_parsed = 0;
        consumer_state.pref_buf.clear();
        consumer_state.pref_it = consumer_state.pref_buf.begin();
        task_status[id] = Task_Status::pending;//等待读取kmer
    }

    // Open the underlying k-mer database.
    // 打开对应的 pre 和 suf 文件
    open_kmer_database(kmer_container->container_location());

    // Launch the background disk-reader thread.
    // 启动后台磁盘读取线程。
    // new std::thread([this]() { read_raw_kmers(); })
    // 创建了一个新的线程对象，该线程对象将执行上述的lambda表达式，即调用read_raw_kmers()函数。
    // reset:std::shared_ptr 的所有权转移到新的指针上，同时释放原有指针所管理的资源。
    // reader指针指向 new thread,同时释放原来指向的线程对象
    // 这里创建线程已经开始读取了
    reader.reset(
        new std::thread([this]()
            {
                read_raw_kmers();
            }
        )
    );
}

template <uint16_t k>
/**
 * @brief 判断是否已启动迭代器
 *
 * 判断迭代器是否已启动。如果迭代器已启动，则返回 true；否则返回 false。
 * return 生产者线程 != nullptr
 * @return 如果迭代器已启动，则返回 true；否则返回 false
 */
inline bool Kmer_SPMC_Iterator<k>::launched() const {
  return reader != nullptr;
}

template <uint16_t k>
/**
 * @brief 读取原始 k-mer
 *
 * 从 k-mer 数据库中读取原始 k-mer 数据，并将其存储到消费者的缓冲区中。
 *
 * @tparam k k-mer 的长度
 */
inline void Kmer_SPMC_Iterator<k>::read_raw_kmers()
{
    // 检查是否到达文件末尾
    while(!kmer_database.Eof())
    {
        // 找到空闲线程
        const size_t consumer_id = get_idle_consumer();
        Consumer_Data& consumer_state = consumer[consumer_id];
        // 这里返回成功读取的kmer个数
        consumer_state.kmers_available = kmer_database.read_raw_suffixes(consumer_state.suff_buf, consumer_state.pref_buf, BUF_SZ_PER_CONSUMER);
       // printf("生产者线程给予的当前的线程是%lu,读取的kmer个数是%lu\n", consumer_id, consumer_state.kmers_available);
        consumer_state.pref_it = consumer_state.pref_buf.begin();
//1462169
        if(!consumer_state.kmers_available)
        {
            std::cerr << "Error reading the suffix file. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        kmers_read += consumer_state.kmers_available;

        consumer_state.kmers_parsed = 0;
        task_status[consumer_id] = Task_Status::available;
    }
}


template <uint16_t k>
/**
 * @brief 获取空闲消费者
 *
 * 返回一个空闲的消费者的索引。如果所有消费者都在忙，则执行忙等待直到找到一个空闲的消费者。
 *
 * @return 空闲消费者的索引
 */
inline size_t Kmer_SPMC_Iterator<k>::get_idle_consumer() const
{
    size_t id{0};//定义一个 `size_t` 类型的局部变量 `id` 并初始化为0。这个变量将用于追踪消费者的ID。
    //从0开始找到一个空闲的消费者
    while(task_status[id] != Task_Status::pending)  // busy-wait
    {
        // id = (id + 1) % consumer_count;
        id++;
        if(id == consumer_count)
            id = 0;//如果到末尾就从头开始
    }

    return id;     
}


template <uint16_t k>
/**
 * @brief 生产者结束生产
 *
 * 等待磁盘读取完成，并等待消费者完成消费，同时通知消费者生产已经结束。
 * 关闭底层的 k-mer 数据库。
 */
inline void Kmer_SPMC_Iterator<k>::seize_production()
{
    // Wait for the disk-reads to be completed.
    if(!reader->joinable())//检查reader是否join
    {
        std::cerr << "Early termination encountered for the database reader thread. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
    // 等待磁盘读取完成
    // join() 方法的作用是阻塞当前线程（在示例中是主线程），直到目标线程（在示例中是 t1 线程）完成执行。join() 并不启动线程，它只是等待线程完成：
    reader->join();


    // Wait for the consumers to finish consumption, and signal them that the means of production have been seized.
    // 等待消费者完成消费，并向他们发出生产资料已被夺取的信号。
    for(size_t id = 0; id < consumer_count; ++id)
    {
        while(task_status[id] != Task_Status::pending); // busy-wait
        
        task_status[id] = Task_Status::no_more;//等待每个消费者消费完成，然后通知消费者生产已经结束
    }

    // Close the underlying k-mer database.
    close_kmer_database();
}


template <uint16_t k>
/**
 * @brief 获取指定消费者的 Kmer 值
 *
 * 根据给定的消费者 ID，从迭代器中获取 Kmer 值，并将其存储在提供的 Kmer 对象中。
 *
 * @param consumer_id 消费者 ID
 * @param kmer 用于存储获取的 Kmer 值的对象
 *
 * @return 如果成功获取到 Kmer 值，则返回 true；否则返回 false
 *
 * @note 尽可能延迟对 volatile 的访问，因为每次访问都会直接命中实际的内存位置。
 */
inline bool Kmer_SPMC_Iterator<k>::value_at(const size_t consumer_id, Kmer<k>& kmer)
{
    // TODO: try to delay this `volatile` access as much as possible, as each access directly hits the actual memory location.尽量延迟这种`volatile`访问，因为每次访问都会直接到达实际的内存位置。
    if(!task_available(consumer_id))//检查线程的状态是否为available
        return false;

    auto& ts = consumer[consumer_id];
    if(ts.kmers_parsed == ts.kmers_available)
    {
        task_status[consumer_id] = Task_Status::pending;
        return false;
    }
  //  printf("consumer_id: %lu\n", consumer_id);
    kmer_database.parse_kmer_buf<k>(ts.pref_it, ts.suff_buf, ts.kmers_parsed * kmer_database.suff_record_size(), kmer);
    ts.kmers_parsed++;

    return true;
}


template <uint16_t k>
/**
 * @brief 判断两个迭代器是否相等
 *
 * 判断当前迭代器与给定的迭代器 `rhs` 是否相等。如果两个迭代器的 `kmer_container` 和 `kmers_read` 成员变量都相等，则返回 `true`；否则返回 `false`。
 *
 * @param rhs 要比较的另一个迭代器
 *
 * @return 如果两个迭代器相等，则返回 `true`；否则返回 `false`
 */
inline bool Kmer_SPMC_Iterator<k>::operator==(const iterator& rhs) const
{
    return kmer_container == rhs.kmer_container && kmers_read == rhs.kmers_read;
}


template <uint16_t k>
inline bool Kmer_SPMC_Iterator<k>::operator!=(const iterator& rhs) const
{
    return !operator==(rhs);
}


template <uint16_t k>
/**
 * @brief 判断是否还有任务需要处理
 *
 * 根据给定的消费者ID，判断该消费者是否还有任务需要处理。
 *
 * @param consumer_id 消费者ID
 *
 * @return 如果还有任务需要处理，返回true；否则返回false
 */
inline bool Kmer_SPMC_Iterator<k>::tasks_expected(const size_t consumer_id) const
{
    return task_status[consumer_id] != Task_Status::no_more;
}


template <uint16_t k>
/**
 * @brief 判断任务是否可用
 *
 * 根据给定的消费者ID，判断该消费者是否有一个可用的任务。
 *
 * @param consumer_id 消费者ID
 *
 * @return 如果任务可用，返回true；否则返回false
 */
inline bool Kmer_SPMC_Iterator<k>::task_available(const size_t consumer_id) const
{
    return task_status[consumer_id] == Task_Status::available;
}


template <uint16_t k>
inline std::size_t Kmer_SPMC_Iterator<k>::memory() const
{
    return CKMC_DB::pref_buf_memory() + (consumer_count * BUF_SZ_PER_CONSUMER);
}


template <uint16_t k>
inline std::size_t Kmer_SPMC_Iterator<k>::memory(const std::size_t consumer_count)
{
    return CKMC_DB::pref_buf_memory() + (consumer_count * BUF_SZ_PER_CONSUMER);
}



#endif
