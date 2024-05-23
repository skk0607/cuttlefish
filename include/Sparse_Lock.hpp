
#ifndef SPARSE_LOCK_HPP
#define SPARSE_LOCK_HPP



#include <cstddef>
#include <cstdint>
#include <vector>
#include <cmath>

// A collection of locks, of type `T_Lock`.
// Intended to be used when a set of sparsely distributed locks over some index range is required.
// //锁的集合，类型为' T_Lock '。用于需要一组分布在某个索引范围上的稀疏锁时。
template <typename T_Lock>
class Sparse_Lock
{
private:

    // Each lock is assigned a power-of-two number of entries to guard.

    const size_t num_entries;           // Number of entries to guard.
    const uint8_t lg_per_lock_range;    // Base-2 log of the number of entries assigned to each lock.
    const size_t per_lock_range;        // Number of contiguous entries (indices) that each lock is assigned to.
    const size_t num_locks;             // Number of locks in the collection.
    std::vector<T_Lock> lock_;          // The collection of locks.


    // Returns the ID of the lock that the index `idx` corresponds to.
    std::size_t lock_id(std::size_t idx) const;


public:

    // Constructs a sparse-lock collection consisting of `lock_count` locks, for `range_size` number of entries.
    Sparse_Lock(size_t range_size, size_t lock_count);

    // Acquires lock for the entry with index `idx`.
    void lock(size_t idx);

    // Releases lock for the entry with index `idx`.
    void unlock(size_t idx);

    // Acquires lock for the entry with index `curr_idx` iff the corresponding lock for the index `prev_idx`
    // is a different lock.
    void lock_if_different(std::size_t prev_idx, std::size_t curr_idx);

    // Releases lock for the entry with index `curr_idx` iff the corresponding lock for the index `prev_idx`
    // is a different lock.
    void unlock_if_different(std::size_t prev_idx, std::size_t curr_idx);
};


template <typename T_Lock>
/**
 * @brief Sparse_Lock 类的构造函数
 *
 * 初始化 Sparse_Lock 类的实例，并设置相关参数。
 *
 * @param range_size 范围大小
 * @param lock_count 锁的数量
 */
inline Sparse_Lock<T_Lock>::Sparse_Lock(const size_t range_size, const size_t lock_count):
    num_entries(range_size),
    // static_cast<uint8_t>(...) 将括号里的数值转为uint8_t
    lg_per_lock_range(static_cast<uint8_t>(std::floor(std::log2((num_entries + lock_count - 1) / lock_count)))),
    per_lock_range(static_cast<size_t>(1) << lg_per_lock_range),
    num_locks((num_entries + per_lock_range - 1) / per_lock_range),
    lock_(num_locks)
{}

template <typename T_Lock>
/**
 * @brief id值 >> 每个锁负责的碱基数量 即为索引
 *
 * 根据给定的索引值计算并返回对应的锁的ID。
 * id值 >> 每个锁负责的碱基数量 即为索引
 * @param idx 索引值
 *
 * @return 锁的ID
 */
inline std::size_t Sparse_Lock<T_Lock>::lock_id(const std::size_t idx) const {
  return idx >> lg_per_lock_range;
}

template <typename T_Lock>
/**
 * @brief 锁定指定索引的锁
 *
 * 根据给定的索引，锁定对应的锁。
 *
 * @param idx 索引值
 */
inline void Sparse_Lock<T_Lock>::lock(const size_t idx)
{
    lock_[lock_id(idx)].lock();
}


template <typename T_Lock>
inline void Sparse_Lock<T_Lock>::unlock(const size_t idx)
{
    lock_[lock_id(idx)].unlock();
}


template <typename T_Lock>
inline void Sparse_Lock<T_Lock>::lock_if_different(const std::size_t prev_idx, const std::size_t curr_idx)
{
    if(lock_id(curr_idx) != lock_id(prev_idx))
        lock_[lock_id(curr_idx)].lock();
}


template <typename T_Lock>
inline void Sparse_Lock<T_Lock>::unlock_if_different(const std::size_t prev_idx, const std::size_t curr_idx)
{
    if(lock_id(curr_idx) != lock_id(prev_idx))
        lock_[lock_id(curr_idx)].unlock();
}



#endif
