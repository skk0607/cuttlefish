
#ifndef CHARACTER_BUFFER_HPP
#define CHARACTER_BUFFER_HPP



#include "Spin_Lock.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "FASTA_Record.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <fstream>
#include <iostream>


// A buffer class to contain contiguous characters. The buffer is to have a maximum
// capacity of `CAPACITY` (although it is non-binding when a string with length 
// larger than that is added), and it flushes to a sink of type `T_sink_` when it
// overflows or is destructed. Writing to the provided sink (in the constructor)
// is thread-safe.
//包含连续字符的buffer类该缓冲区的最大容量为` capacity `(尽管当添加长度大于该长度的字符串时，它是非绑定的)，并且当缓冲区溢出或被销毁时，它会刷新到类型为`T_sink_`的sink。写入提供的接收器(在构造函数中)是线程安全的。
template <std::size_t CAPACITY, typename T_sink_>
class Character_Buffer
{
private:

    std::vector<char> buffer;   // The character buffer.
    T_sink_& sink;  // Reference to the sink to flush the buffer content to.


    // Ensures that `buffer` has enough space for additional `append_size`
    // number of bytes, using flush and allocation as necessary.
    void ensure_space(std::size_t append_size);

    // Flushes the buffer content to the sink, and clears the buffer.
    void flush();


public:

    // Constructs a character buffer object that would flush its content to `sink`.
    Character_Buffer(T_sink_& sink);

    // Appends the content of `str` to the buffer. Flushes are possible.
    template <typename T_container_>
    void operator+=(const T_container_& str);

    // Appends the content of the FASTA record `fasta_rec` to the buffer. Flushes
    // are possible.
    template <typename T_container_>
    void operator+=(const FASTA_Record<T_container_>& fasta_rec);

    // Appends the content of the FASTA record `fasta_cycle` to the buffer, that is
    // supposed to be a cycle in a de Bruijn graph `G(·, k)`. The cyclic FASTA
    // sequence is rotated around its index `pivot` — the entire sequence is
    // right-rotated so that the `pivot`-index character is at index 0 finally. A
    // line-break is added at the end of the sequence, since the user might not be
    // able to provide it with the "to be rotated" sequence.
    template <uint16_t k, typename T_container_>
    void rotate_append_cycle(const FASTA_Record<T_container_>& fasta_cycle, std::size_t pivot);

    // Destructs the buffer object, flushing it if content are present.
    ~Character_Buffer();
};


// Helper class to actually flush the content of the `Character_Buffer` class to its
// sink of type `T_sink`.
// It's used to circumvent the C++ constraint that partial specialization of a
// a member function is not possible without partially specializing the entire
// class. We need to specialize the actual flushing mechanism to support various
// types of sinks, e.g. `std::ofstream`, `spdlog::logger` etc.
// Since the sole purpose of the class is to support the `Character_Buffer` class
// circumvent some contraint, everything is encapsulated in its specializations
// as private, with `Character_Buffer` as friend.
template <typename T_sink_>
class Character_Buffer_Flusher
{};


template <>
class Character_Buffer_Flusher<std::ofstream>
{
    template <std::size_t, typename> friend class Character_Buffer;

private:

    // Mutual-exclusion lock to control multi-threaded access to otherwise not thread-
    // safe sinks (e.g. `std::ofstream`). Note that, the lock is per sink-type, not per
    // actual sink — which is a limitation.
    static Spin_Lock lock;


    // Writes the content of the vector `buf` to the sink `sink`.
    static void write(std::vector<char>& buf, std::ofstream& sink);
};


template <>
class Character_Buffer_Flusher<Async_Logger_Wrapper>
{
    template <std::size_t, typename> friend class Character_Buffer;

private:

    // Writes the content of the vector `buf` to the sink `sink`. Note that the vector
    // `buf` is modified in the process — a null-terminator (`\0`) is appended at the
    // end — which is expected to be not problematic under the assumption that the
    // buffer is cleared after the write (i.e. flush).
    static void write(std::vector<char>& buf, const Async_Logger_Wrapper& sink);
};


template <std::size_t CAPACITY, typename T_sink_>
inline Character_Buffer<CAPACITY, T_sink_>::Character_Buffer(T_sink_& sink):
    sink(sink)
{
    buffer.reserve(CAPACITY);
}


template <std::size_t CAPACITY, typename T_sink_>
template <typename T_container_>
inline void Character_Buffer<CAPACITY, T_sink_>::operator+=(const T_container_& str)
{
    ensure_space(str.size());

    // `std::memcpy` at the end of `buffer` does not update the size of the vector `buffer`.
    buffer.insert(buffer.end(), str.begin(), str.end());
}


template <std::size_t CAPACITY, typename T_sink_>
template <typename T_container_>
/**
 * @brief 向字符缓冲区追加FASTA记录
 *
 * 将给定的FASTA记录追加到字符缓冲区中。
 *
 * @param fasta_rec 要追加的FASTA记录
 */
inline void Character_Buffer<CAPACITY, T_sink_>::operator+=(const FASTA_Record<T_container_>& fasta_rec)
{
    // 确保有足够的空间存储头信息、换行符、序列和换行符
    ensure_space(fasta_rec.header_size() + 1 + fasta_rec.seq_size() + 1);   // Two extra bytes for the line-breaks.

    // 追加头信息
    fasta_rec.append_header(buffer);    // Append the header.

    // 添加换行符
    buffer.emplace_back('\n');  // Break line.

    // 追加序列
    fasta_rec.append_seq(buffer);   // Append the sequence.

    // 添加换行符
    buffer.emplace_back('\n');  // Break line.
}


template <std::size_t CAPACITY, typename T_sink_>
template <uint16_t k, typename T_container_>
inline void Character_Buffer<CAPACITY, T_sink_>::rotate_append_cycle(const FASTA_Record<T_container_>& fasta_rec, const std::size_t pivot)
{
    ensure_space(fasta_rec.header_size() + 1 + fasta_rec.seq_size() + 1);   // Two extra bytes for two line-breaks.

    fasta_rec.append_header(buffer);    // Append the header.
    buffer.emplace_back('\n');  // Break line.
    fasta_rec.template append_rotated_cycle<k>(buffer, pivot);  // Append the sequence right-rotated around index `pivot`.
    buffer.emplace_back('\n');  // End the sequence.
}


template <std::size_t CAPACITY, typename T_sink_>
inline void Character_Buffer<CAPACITY, T_sink_>::ensure_space(const std::size_t append_size)
{
    if(buffer.size() + append_size >= CAPACITY) // Using `>=` since for async logging, a `\0` is inserted at the end of `buffer`.
    {
        flush();
        
        if(append_size >= CAPACITY)
        {
            // std::cerr <<    "A single output string overflows the string-buffer capacity.\n"
            //                 "Output string length: " << str.size() << ", string-buffer capacity: " << CAPACITY << ".\n"
            //                 "Please consider increasing the buffer capacity parameter in build for future use.\n";
            
            buffer.reserve(append_size);
        }
    }
}


template <std::size_t CAPACITY, typename T_sink_>
inline void Character_Buffer<CAPACITY, T_sink_>::flush()
{
    Character_Buffer_Flusher<T_sink_>::write(buffer, sink);

    buffer.clear();
}


template <std::size_t CAPACITY, typename T_sink_>
/**
 * @brief 销毁 Character_Buffer 对象
 *
 * 在销毁 Character_Buffer 对象时，如果缓冲区不为空，则调用 flush() 函数清空缓冲区。
 */
inline Character_Buffer<CAPACITY, T_sink_>::~Character_Buffer()
{
    if(!buffer.empty())
        flush();
}


/**
 * @brief 将字符缓冲区写入输出流
 *
 * 将给定的字符缓冲区写入指定的输出流中。
 *
 * @param buf 字符缓冲区
 * @param output 输出流
 */
inline void Character_Buffer_Flusher<std::ofstream>::write(std::vector<char>& buf, std::ofstream& output)
{
    lock.lock();

    output.write(buf.data(), buf.size());

    if(output.fail())
    {
        std::cerr << "Error writing the output. Aborting.\n";
        std::exit(EXIT_FAILURE);
    }

    lock.unlock();
}


inline void Character_Buffer_Flusher<Async_Logger_Wrapper>::write(std::vector<char>& buf, const Async_Logger_Wrapper& sink)
{
    buf.emplace_back('\0');

    sink.write(buf.data());
}


#endif
