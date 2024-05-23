
#ifndef SPIN_LOCK_HPP
#define SPIN_LOCK_HPP



#include <atomic>


// A lightweight lock-free mutex class.
// It is based on `std::atomic_flag`, which is guaranteed to be a lock-free atomic construct .
// Reference: https://en.cppreference.com/w/cpp/atomic/atomic_flag
class Spin_Lock
{
private:
  // The atomic flag.
  // 宏定义: ATOMIC_FLAG_INIT = 0
  std::atomic_flag lock_ = ATOMIC_FLAG_INIT;
  // lock_ 是一个 std::atomic_flag
  // 类型的变量，用于实现锁定和解锁操作。ATOMIC_FLAG_INIT 是其初始化器，确保
  // atomic_flag 的初始状态为未设置。

public:

    // Acquires the lock for mutually-exlcusive access to it.
    void lock();

    // Releases the lock, giving up the exclusive access to it.
    void unlock();
};


/**
 * @brief 锁定自旋锁
 *
 * 尝试获取自旋锁，如果锁已经被其他线程持有，则不断尝试直到获取成功。
 * 使用 `memory_order_acquire` 确保在当前线程中的读写操作不会被重新排序到加载变量 `lock_` 之前，
 * 保证 `lock` 调用之后的内存访问指令顺序不变。
 */
inline void Spin_Lock::lock()
{
    // Due to the memory access order `memory_order_acquire`, no reads or writes in the current thread can be
    // reordered before this load of the variable `lock_` (enforced by the compiler and the processor) —
    // ensuring that memory-access instructions after a `lock` invokation stays after it.

    while(lock_.test_and_set(std::memory_order_acquire))
        ;// while(lock_.test(std::memory_order_relaxed));   // C++20 optimization to avoid the redundant stores from the spinning `test_and_set`.
}


/**
 * @brief 解锁Spin_Lock
 *
 * 解锁Spin_Lock，释放锁定的资源。
 *
 * 使用`memory_order_release`内存访问顺序，确保在当前线程中，此存储操作（`lock_.clear`）之后不会发生任何读写操作（由编译器和处理器强制执行），
 * 从而确保在`unlock`调用之前的内存访问指令都保持在它之前。
 */
inline void Spin_Lock::unlock()
{
    // Due to the memory access order `memory_order_release`, no reads or writes in the current thread can be
    // reordered after this store of the variable `lock_` (enforced by the compiler and the processor) —
    // ensuring that memory-access instructions before an `unlock` invokation stays before it.
    
    lock_.clear(std::memory_order_release);
}



#endif
