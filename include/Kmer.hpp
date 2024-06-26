
#ifndef KMER_HPP
#define KMER_HPP



#include "DNA_Utility.hpp"
#include "Kmer_Utility.hpp"
#include "utility.hpp"
#include "kmc_api/kmc_file.h"
#include "xxHash/xxh3.h"

#include <cstdint>
#include <cstddef>
#include <cstring>
#include <string>
#include <algorithm>
#include <cassert>

/*
定义这个宏表明了我们的意图，即de Bruijn图顶点只使用奇数个k值。
因此，从(k + 1)-mers(边的顶点)中提取k-mers只会在k为奇数时发生。
*/
// Defining this macro states our intent that only odd k-values will be used for de Bruijn graph vertices.
// Hence, extraction of k-mers from (k + 1)-mers — vertices from edges — will only happen when k is odd.
#define ODD_K


template <uint16_t k>
class Kmer
{
    // Make k-mers friend for (k + 1)-mer, so that de Bruijn graph vertices, i.e. k-mers,
    // may access private information (the raw data) from edges, i.e. (k + 1)-mers.
    friend class Kmer<k - 1>;

    // Minimizers can be represented using 32-bit integers.
    typedef uint32_t minimizer_t;

private:

    // Number of 64-bit integers required to compactly represent the underlying k-mer with 2-bits/base encoding.
    static constexpr uint16_t NUM_INTS = (k + 31) / 32;

    // Bitmask used to clear the most significant DNA base character, i.e. the first base of the k-mer which is at the bits `2k-1 : 2k-2`.
    static constexpr const uint64_t CLEAR_MSN_MASK = ~(uint64_t(0b11) << (2 * ((k - 1) % 32)));

    // The underlying k-mer represented with 2-bit encoding, as a collection of
    // 64-bit integers. A k-mer `n_{k - 1} ... n_1 n_0` is stored in the array
    // `kmer_data` such that, `kmer_data[0]` stores the suffix `n_63 ... n_0`,
    // then `kmer_data[1]` stores `n_127 ... n_64`, and so on. That is, the
    // suffix is aligned with a byte boundary.
    // TODO: reverse this store-order of the data — i.e. `n_{k - 1}` as the
    // least significant base (so stored in `kmer_data[0]`) and `n_0` as the
    // most significant base (so stored in `kmer_data[0]`). This would optimize
    // at least the following:
    //      i) `from_KMC_data` () — aligning with the KMC data alignment and
    //      thus `memcpy` instead of bit-twiddling; ii) `operator<` — `memcmp`
    //      instead of highest-to-lowest index looping comparison.
    // 实际用于存储 kmer的数据,即uint64_t的数组~
    uint64_t kmer_data[NUM_INTS];


    // Left-shifts the collection of the bits at the `kmer_data` array by one base (2-bits).
    void left_shift();

    // Right-shifts the collection of the bits at the `kmer_data` array by one base (2-bits).
    void right_shift();

    // Left-shifts the collection of the bits at the `kmer_data` array by `B` bases (2B-bits).
    template <uint16_t B> void left_shift();


public:

    // Default constructs the k-mer with a 0-value, equivalent to "AA...A".
    Kmer();

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer(const char* label, std::size_t kmer_idx);

    // Constructs a k-mer from the provided characters at
    // `label[0, ..., k - 1]`.
    Kmer(const char* label);

    // Constructs a k-mer from the provided string `label`.
    Kmer(const std::string& label);

    // Constructs a k-mer from the provided characters at
    // `label[kmer_idx,...,kmer_idx + k - 1]`.
    Kmer(const std::string& label, std::size_t kmer_idx);

    // Constructs a k-mer from `kmer_api`, a k-mer object built from KMC.
    Kmer(const CKmerAPI& kmer_api);

    // Copy constructs the k-mer from another k-mer `rhs`.
    Kmer(const Kmer<k>& rhs);

    // Copy assignment operator.
    Kmer<k>& operator=(const Kmer<k>& rhs);

    // Returns a 64-bit hash value for the k-mer.
    uint64_t to_u64(uint64_t seed=0) const;

    // Gets the k-mer from the KMC api object `kmer_api`.
    void from_CKmerAPI(const CKmerAPI& kmer_api);

    // Gets the k-mer from its KMC raw-binary representation.
    void from_KMC_data(const uint64_t* kmc_data);

    // Gets the k-mer that is a prefix of the provided
    // (k + 1)-mer `k_plus_1_mer`.
    void from_prefix(const Kmer<k + 1>& k_plus_1_mer);

    // Gets the k-mer that is a suffix of the provided
    // (k + 1)-mer `k_plus_1_mer`.
    void from_suffix(const Kmer<k + 1>& k_plus_1_mer);

    // Returns the reverese complement of the k-mer.
    Kmer<k> reverse_complement() const;

    // Gets the k-mer that is the reverse complement of
    // the provided k-mer `other`.
    void as_reverse_complement(const Kmer<k>& kmer);

    // Returns true iff the bitwise encoding of this k-mer is lesser to the
    // encoding of the other k-mer `rhs`.
    bool operator<(const Kmer<k>& rhs) const;

    // Returns `true` iff the bitwise encoding of this k-mer is larger to
    // the encoding of the other k-mer `rhs`.
    bool operator>(const Kmer<k>& rhs) const;

    // Returns true iff this k-mer is identical to the other k-mer `rhs`.
    bool operator==(const Kmer<k>& rhs) const;

    // Returns true iff this k-mer is not identical to the other k-mer `rhs`.
    bool operator!=(const Kmer<k>& rhs) const;

    // Returns the `DNA::Base` (2-bit) encoding of the character at the front,
    // i.e. at the first index of the literal representation. For a k-mer
    // `n_{k - 1} ... n_1 n_0`, this is the base `n_{k - 1}`.
    DNA::Base front() const;

    // Returns the `DNA::Base` (2-bit) encoding of the character at the back,
    // i.e. at the last index of the literal representation. For a k-mer
    // `n_{k - 1} ... n_1 n_0`, this is the base `n_0`.
    DNA::Base back() const;

    // Returns `true` iff the k-mer is in the forward direction relative to
    // the other k-mer `kmer_hat`.
    bool in_forward(const Kmer<k>& kmer_hat) const;

    // Transforms this k-mer by chopping off the first base and
    // appending the next base character `next_base` to the end,
    //  i.e. rolls the k-mer by one base. Also sets the passed
    // reverse complement `rev_compl` of the k-mer accordingly.
    void roll_to_next_kmer(char next_base, Kmer<k>& rev_compl);

    // Transforms this k-mer by chopping off the first base and
    // appending the next base `base` to the end, i.e.  rolls
    // the k-mer by one base. Also sets the passed reverse
    // complement `rev_compl` of the k-mer accordingly.
    void roll_to_next_kmer(DNA::Base base, Kmer<k>& rev_compl);

    // Transforms this k-mer by chopping off the first base and
    // appending the next base coded with the edge encoding `edge`,
    // i.e. rolls the k-mer by one base. Also sets the passed
    // reverse complement `rev_compl` of the k-mer accordingly.
    void roll_to_next_kmer(DNA::Extended_Base edge, Kmer<k>& rev_compl);

    // Transforms this k-mer by chopping off the first base and
    // appending the base coded with the edge encoding `edge` to
    // the end, i.e. rolls the k-mer to the "right" by one base.
    void roll_forward(DNA::Extended_Base edge);

    // Transforms this k-mer by chopping off the last base and
    // appending the base coded with the edge encoding `edge` to
    // the beginning, i.e. rolls the k-mer to the "left" by one base.
    void roll_backward(DNA::Extended_Base edge);

    // Returns the canonical version of the k-mer, comparing it to its
    // reverse complement `rev_compl`.
    Kmer canonical(const Kmer<k>& rev_compl) const;

    // Returns the canonical version of the k-mer.
    Kmer canonical() const;

    // Given a k-mer `kmer` and its reverse complement `rev_compl`,
    // returns a pointer to one of these, which represents the
    // canonical form.
    static const Kmer<k>* canonical(const Kmer<k>& kmer, const Kmer<k>& rev_compl);

    // Returns the string label of the k-mer.
    std::string string_label() const;

    // Gets the string label of the k-mer into the container `label`.
    template <typename T_container_>
    void get_label(T_container_& label) const;

    // Implicitly converts the k-mer to a `std::string`.
    operator std::string() const;

    // Returns a randomly generated k-mer.
    static Kmer<k> random_kmer();

    // Prints the literal representation of the K-mer `kmer` to the
    // stream `ostream`.
    template <uint16_t K>
    friend std::ostream& operator<<(std::ostream& out, const Kmer<K>& kmer);

    // Returns the lexicographic l-minimizer for the k-mer.
    template <uint8_t l>
    minimizer_t minimizer() const;

    // Returns the l-minimizer for the k-mer where the vector `order`
    // determines the minimizer-ordering of the l-mers, i.e. the order
    // of the l-mer `i` is `order[i]`.
    template <uint8_t l>
    minimizer_t minimizer(const std::vector<uint32_t>& order) const;

    // Accumulates the counts of the l-mers of the k-mer into `count`.
    template <uint8_t l>
    void count_lmers(std::vector<uint64_t>& count) const;
};


template <uint16_t k>
inline void Kmer<k>::left_shift()
{
    left_shift<1>();
}


template <uint16_t k>
/**
 * @brief 将 Kmer 向右移动两位
 *
 * 将 Kmer 向右移动两位，通过循环遍历 kmer_data 数组，并将每个元素右移两位。
 * 同时，将下一个元素的最低两位左移 62 位，然后与前一个元素进行或运算。
 * 相当于整个数组向右移2位
 * @tparam k Kmer 的长度
 */
inline void Kmer<k>::right_shift()
{
    constexpr uint64_t mask_LSN = 0b11;

    for(uint16_t idx = 0; idx < NUM_INTS - 1; ++idx)
    //idx >> 2位,然后idx+1 取最后2位 << 62位
        kmer_data[idx] = (kmer_data[idx] >> 2) | ((kmer_data[idx + 1] & mask_LSN) << 62);

    kmer_data[NUM_INTS - 1] >>= 2;
}


template <uint16_t k>
template <uint16_t B>
inline void Kmer<k>::left_shift()
{
    static_assert(B < 32, "invalid bit-shift amount");

    if constexpr(B != 0)
    {
        constexpr uint16_t num_bit_shift = 2 * B;
        // if B = 1, mask_MSNs = 11 << 62
        constexpr uint64_t mask_MSNs = ((static_cast<uint64_t>(1) << num_bit_shift) - 1) << (64 - num_bit_shift);
        // idx的末尾62 bit 和 idx-1 的最高 2bit拼接
        for(uint16_t idx = NUM_INTS - 1; idx > 0; --idx)
            kmer_data[idx] = (kmer_data[idx] << num_bit_shift) | ((kmer_data[idx - 1] & mask_MSNs) >> (64 - num_bit_shift));

        kmer_data[0] <<= num_bit_shift;
    }
}


template <uint16_t k>
/**
 * @brief 将 Kmer 转换为 64 位无符号整数
 *
 * 使用给定的种子值，将 Kmer 对象转换为 64 位无符号整数。
 *
 * @param seed 种子值
 *
 * @return 转换后的 64 位无符号整数
 */
inline uint64_t Kmer<k>::to_u64(const uint64_t seed) const
{
    //k值向上取整自4的倍数
    constexpr uint16_t NUM_BYTES = (k + 3) / 4;
    /**
在C++中，这段代码似乎是为了将一个长度为k的k-mer（即一个由k个核苷酸组成的序列）转换为一个64位无符号整数（uint64_t）。转换过程中使用了XXH3_64bits_withSeed函数，它可能是一个来自某个哈希库（如xxHash）的哈希函数。
NUM_BYTES计算的是k-mer的字节长度。k-mer通常由DNA序列的核苷酸（A, C, G,
T）组成，每个核苷酸可以编码为2比特信息（如果使用二进制编码）。因此，理论上一个k-mer需要的比特数是2
* k。但是，在实际处理中，我们通常使用字节（8比特）作为数据的基本单位。 (k + 3) /
4这个表达式是用来计算将k-mer转换为字节所需要的最少字节数。由于每个字节包含8比特，所以每个k-mer需要的字节数应该是2
* k除以8。但是，由于我们不能有“部分字节”，所以我们需要向上取整来确保有足够的空间来存储整个k-mer。 (k + 3) / 4这个计算等同于(k + 7) / 8（因为除以4等于除以8然后乘以2，加上3等于加上7），后者是常见的向上取整到最接近的字节数的技巧。加7然后除以8确保了当k是8的倍数时，结果不会向下取整。然后，结果乘以2（因为k是按比特计算的，而我们需要字节），所以实际上是(k + 7) / 4。但在这个代码中，用了(k + 3) / 4，这实际上是一个简化，因为当k不是4的倍数时，结果仍然是正确的（向上取整到最近的4的倍数），但是当k是4的倍数时，它不会增加额外的字节。这可能是因为k-mer数据已经以某种方式确保了足够的存储，或者这个简化在上下文中是可以接受的。 简而言之，(k + 3) / 4是计算存储k-mer所需的最小字节数的近似值，并向上取整到最近的4的倍数。这种近似在某些情况下可能是可以接受的，尽管更精确的计算可能是(k + 7) / 8然后乘以2（转换为字节）。但需要注意的是，这里的解释假设k-mer的编码是以比特为基础的，并且每个核苷酸编码为2比特。实际实现可能有所不同，所以具体细节可能需要查看Kmer类的其余部分或相关文档。
     */
    return XXH3_64bits_withSeed(kmer_data, NUM_BYTES, seed);
}


template <uint16_t k>
inline Kmer<k>::Kmer():
    kmer_data() // Value-initializes the data array, i.e. zeroes it out.
{}


template <uint16_t k>
inline Kmer<k>::Kmer(const char* const label, const std::size_t kmer_idx):
    Kmer(label + kmer_idx)
{}


template <uint16_t k>
__attribute__((optimize("unroll-loops")))
inline Kmer<k>::Kmer(const char* const label)
{
    assert(std::strlen(label) == k);

    constexpr uint16_t packed_word_count = k / 32;

    // Get the fully packed words' binary representations.
    for(uint16_t data_idx = 0; data_idx < packed_word_count; ++data_idx)
        kmer_data[data_idx] = Kmer_Utility::encode<32>((label + k) - (data_idx << 5) - 32);

    // Get the partially packed (highest index) word's binary representation.
    if constexpr((k & 31) > 0)
        kmer_data[NUM_INTS - 1] = Kmer_Utility::encode<k & 31>(label);
}


template <uint16_t k>
inline Kmer<k>::Kmer(const std::string& label):
    Kmer(label.c_str())
{}


template <uint16_t k>
inline Kmer<k>::Kmer(const std::string& label, const std::size_t kmer_idx):
    Kmer(label.c_str(), kmer_idx)
{}


template <uint16_t k>
inline Kmer<k>::Kmer(const CKmerAPI& kmer_api)
{
    from_CKmerAPI(kmer_api);
}


template <uint16_t k>
inline Kmer<k>::Kmer(const Kmer<k>& rhs)
{
    std::memcpy(kmer_data, rhs.kmer_data, NUM_INTS * sizeof(uint64_t));
}


template <uint16_t k>
inline Kmer<k>& Kmer<k>::operator=(const Kmer<k>& rhs)
{
    std::memcpy(kmer_data, rhs.kmer_data, NUM_INTS * sizeof(uint64_t));

    return *this;
}


template <uint16_t k>
inline void Kmer<k>::from_CKmerAPI(const CKmerAPI& kmer_api)
{
    kmer_api.to_u64<NUM_INTS>(kmer_data);
}


template <uint16_t k>
inline void Kmer<k>::from_KMC_data(const uint64_t* const kmc_data)
{
    // The endianness of the k-mer data array in the KMC database is in the opposite
    // order of the one that Cuttlefish uses. So the fetch gathers it in reverse.
    // Roughly, the KMC database stores a literal k-mer `b_{k - 1} ... b_1 b_0` such
    // that the prefix aligns with a byte boundary, i.e. `b_{k - 1} ... b_{k - 1 - 63}`
    // is stored as one 64-bit collection at its `kmer_data[0]` entry, and the rest
    // follows this alignment. This is why, the 64-bits corresponding to the substring
    // `b_63 ... b_0` might be shared between the two maximum indices of `kmer_data`;
    // so can be all the next disjoint 64-bit substrings. The remainder substring
    // (of < 64-bits) can be found from just the 0'th entry of `kmer_data`.
    
    constexpr uint8_t byte_alignment = (k % 4 != 0 ? 4 - (k % 4) : 0);
    constexpr uint32 offset = 62 - ((k - 1 + byte_alignment) & 31) * 2;
    //std::memcpy(kmer_data, kmc_data, NUM_INTS * sizeof(kmer_data[0]));
    
    if(offset)
    {
        for(int32 i = NUM_INTS - 1; i >= 1; --i)
        {
            kmer_data[NUM_INTS - 1 - i] = kmc_data[i] >> offset;
            kmer_data[NUM_INTS - 1 - i] += kmc_data[i - 1] << (64 - offset);
        }

        kmer_data[NUM_INTS - 1] = kmc_data[0] >> offset;
    }
    else
        for (int32 i = NUM_INTS - 1; i >= 0; --i)
            kmer_data[NUM_INTS - 1 - i] = kmc_data[i];
            
}


template <uint16_t k>
inline void Kmer<k>::from_prefix(const Kmer<k + 1>& k_plus_1_mer)
{
    // Note: `Kmer<k>` and `Kmer<k + 1>` always have the same number of words (i.e. `NUM_INTS`) for odd k-values.
    // The only time that they have different numbers of words is when `k` is a multiple of 32. In such cases,
    // a (k + 1)-mer contains just one base in its highest index word, and a k-mer's words are fully packed.
    //先将 k+1 mer 的前 k 位复制到 k-mer->kmer_data 中
    std::memcpy(kmer_data, k_plus_1_mer.kmer_data, NUM_INTS * sizeof(uint64_t));
    right_shift();  // Clear the LSN of the (k + 1)-mer, from the k-mer.

    #ifndef ODD_K   // The following `if` conditional can only be `true` when `k` is a multiple of 32.
    constexpr uint16_t kp1_NUM_INTS = ((k + 1) + 31) / 32;
    if(kp1_NUM_INTS != NUM_INTS)    // Fetch the only base not copied from the (k + 1)-mer as the MSN for this k-mer.
        kmer_data[NUM_INTS - 1] |= (k_plus_1_mer.kmer_data[kp1_NUM_INTS - 1] << 62);
    #endif
}


template <uint16_t k>
inline void Kmer<k>::from_suffix(const Kmer<k + 1>& k_plus_1_mer)
{
    std::memcpy(kmer_data, k_plus_1_mer.kmer_data, NUM_INTS * sizeof(uint64_t));

    #ifndef ODD_K   // The following `if` conditional can only be `true` when `k` is a multiple of 32.
    constexpr uint16_t kp1_NUM_INTS = ((k + 1) + 31) / 32;
    if(kp1_NUM_INTS != NUM_INTS)    // The only base not copied from the (k + 1)-mer isn't required to be fetched — it will be cleared out anyways.
        return;
    #endif

    kmer_data[NUM_INTS - 1] &= Kmer<k + 1>::CLEAR_MSN_MASK; // Clear the MSN of the (k + 1)-mer from this k-mer.
}


template <uint16_t k>
/**
 * @brief 返回 Kmer 的反向互补序列
 *
 * 生成当前 Kmer 对象的反向互补序列，并返回一个新的 Kmer 对象。
 *
 * @return 返回一个包含反向互补序列的 Kmer 对象
 */
inline Kmer<k> Kmer<k>::reverse_complement() const
{
    Kmer<k> rev_compl;
    rev_compl.as_reverse_complement(*this);

    return rev_compl;
}


template <uint16_t k>
/**
 * @brief 获取另一个Kmer的反向互补序列
 *
 * 将给定的Kmer对象的反向互补序列赋值给当前Kmer对象。
 * 在处理过程中，每次处理一个字节而不是64位整数。
 *
 * @param other 另一个Kmer对象
 */
inline void Kmer<k>::as_reverse_complement(const Kmer<k>& other)
{
    // Working with bytes instead of 64-bit words at a time.

    uint8_t* const rev_compl = reinterpret_cast<uint8_t*>(kmer_data);
    const uint8_t* const data = reinterpret_cast<const uint8_t*>(other.kmer_data);


    // Get the reverse complement for the fully packed bytes.

    constexpr uint16_t packed_byte_count = k / 4;
    for(uint16_t byte_idx = 0; byte_idx < packed_byte_count; ++byte_idx)
        rev_compl[packed_byte_count - 1 - byte_idx] = Kmer_Utility::reverse_complement(data[byte_idx]);


    // Get the reverse complement for the only byte that might be partially packed (possible for the highest-indexed byte only).

    constexpr uint16_t rem_base_count = k % 4;
    if constexpr(rem_base_count == 0)
        return;
    
    rev_compl[packed_byte_count] = 0;
    left_shift<rem_base_count>();

    for(int i = 0; i < rem_base_count; ++i)
        rev_compl[0] |= (DNA_Utility::complement(DNA::Base((data[packed_byte_count] & (0b11 << (2 * i))) >> (2 * i)))
                                        << (2 * (rem_base_count - 1 - i)));
}


template <uint16_t k>
inline bool Kmer<k>::operator<(const Kmer<k>& rhs) const
{
    for(int16_t idx = NUM_INTS - 1; idx >= 0; --idx)
        if(kmer_data[idx] != rhs.kmer_data[idx])
            return kmer_data[idx] < rhs.kmer_data[idx];

    return false;
}


template <uint16_t k>
inline bool Kmer<k>::operator>(const Kmer<k>& rhs) const
{
    for(int16_t idx = NUM_INTS - 1; idx >= 0; --idx)
        if(kmer_data[idx] != rhs.kmer_data[idx])
            return kmer_data[idx] > rhs.kmer_data[idx];

    return false;
}


template <uint16_t k>
inline bool Kmer<k>::operator==(const Kmer<k>& rhs) const
{
    return std::memcmp(kmer_data, rhs.kmer_data, NUM_INTS * sizeof(uint64_t)) == 0;
}


template <uint16_t k>
inline bool Kmer<k>::operator!=(const Kmer<k>& rhs) const
{
    return !operator==(rhs);
}


template <uint16_t k>
/**
 * @brief 获取Kmer对象的第一个碱基
 *
 * 获取Kmer对象的第一个碱基，即最高位的碱基。
 *
 * @return DNA::Base Kmer对象的第一个碱基
 */
inline DNA::Base Kmer<k>::front() const
{
    // Relative index of the most significant nucleotide in it's 64-bit word.
    constexpr uint16_t rel_idx_MSN = 2 * ((k - 1) % 32);

    // Mask to extract the most significant nucleotide.
    constexpr uint64_t mask_MSN = (static_cast<uint64_t>(0b11) << rel_idx_MSN);

    return DNA::Base((kmer_data[NUM_INTS - 1] & mask_MSN) >> rel_idx_MSN);
}


template <uint16_t k>
inline DNA::Base Kmer<k>::back() const
{
    // Mask to extract the least significant nucleotide.
    constexpr uint64_t mask_LSN = static_cast<uint64_t>(0b11);

    return DNA::Base(kmer_data[0] & mask_LSN);
}


template <uint16_t k>
inline bool Kmer<k>::in_forward(const Kmer<k>& kmer_hat) const
{
    return operator==(kmer_hat);
}


template <uint16_t k>
/**
 * @brief 滚动到下一个 k-mer
 *
 * 根据给定的下一个碱基，将当前 k-mer 滚动到下一个 k-mer，并计算其反向互补序列。
 *
 * @param next_base 下一个碱基
 * @param rev_compl 反向互补序列
 */
inline void Kmer<k>::roll_to_next_kmer(const char next_base, Kmer<k>& rev_compl)
{
    const DNA::Base mapped_base = DNA_Utility::map_base(next_base);

    roll_to_next_kmer(mapped_base, rev_compl);
}


template <uint16_t k>
/**
 * @brief 将当前Kmer滚动到下一个Kmer
 *
 * 将当前Kmer滚动到下一个Kmer，并更新其反向互补序列。
 *
 * @param base DNA碱基
 * @param rev_compl 反向互补序列的引用
 */
inline void Kmer<k>::roll_to_next_kmer(const DNA::Base base, Kmer<k>& rev_compl)
{
    // Logically, since a left shift moves the MSN out of the length `k` boundary, the clearing of the base
    // may seem redundant. But, the `to_u64` hashing method implementation works with bytes — not clearing
    // out this base breaks the consistency of the hashing.
    kmer_data[NUM_INTS - 1] &= CLEAR_MSN_MASK;
    left_shift();
    //整体左移2bit,然后 | base
    kmer_data[0] |= base;
    //整体右移2bit,然后 | base
    rev_compl.right_shift();
    rev_compl.kmer_data[NUM_INTS - 1] |= (static_cast<uint64_t>(DNA_Utility::complement(base)) << (2 * ((k - 1) & 31)));
}


template <uint16_t k>
inline void Kmer<k>::roll_to_next_kmer(const DNA::Extended_Base edge, Kmer<k>& rev_compl)
{
    const DNA::Base mapped_base = DNA_Utility::map_base(edge);

    roll_to_next_kmer(mapped_base, rev_compl);
}


template <uint16_t k>
inline void Kmer<k>::roll_forward(const DNA::Extended_Base edge)
{
    const DNA::Base mapped_base = DNA_Utility::map_base(edge);

    kmer_data[NUM_INTS - 1] &= CLEAR_MSN_MASK;
    left_shift<1>();
    kmer_data[0] |= static_cast<uint64_t>(mapped_base);
}


template <uint16_t k>
inline void Kmer<k>::roll_backward(const DNA::Extended_Base edge)
{
    // Relative index of the most significant nucleotide in it's 64-bit word.
    constexpr uint16_t rel_idx_MSN = 2 * ((k - 1) % 32);

    const DNA::Base mapped_base = DNA_Utility::map_base(edge);

    right_shift();
    kmer_data[NUM_INTS - 1] |= (static_cast<uint64_t>(mapped_base) << rel_idx_MSN);
}


template <uint16_t k>
inline Kmer<k> Kmer<k>::canonical(const Kmer<k>& rev_compl) const
{
    return std::min(*this, rev_compl);
}


template <uint16_t k>
inline Kmer<k> Kmer<k>::canonical() const
{
    return canonical(reverse_complement());
}


template <uint16_t k>
inline const Kmer<k>* Kmer<k>::canonical(const Kmer<k>& kmer, const Kmer<k>& rev_compl)
{
    return kmer < rev_compl ? &kmer : &rev_compl;
}


template <uint16_t k>
inline std::string Kmer<k>::string_label() const
{
    Kmer<k> kmer(*this);
    char* label = new char[k + 1];

    for(uint16_t idx = 0; idx < k; ++idx)
    {
        switch(kmer.kmer_data[0] & 0b11)
        {
        case DNA::A:
            label[idx] = 'A';
            break;
        
        case DNA::C:
            label[idx] = 'C';
            break;

        case DNA::G:
            label[idx] = 'G';
            break;

        case DNA::T:
            label[idx] = 'T';
            break;

        default:
            label[idx] = 'N';
        }


        kmer.right_shift();
    }


    std::reverse(label, label + k);
    std::string str_label(label, label + k);

    
    delete[] label;
    
    return str_label;
}


template <uint16_t k>
template <typename T_container_>
/**
 * @brief 获取 Kmer 的标签
 *
 * 从 Kmer 对象中获取标签，并将其存储到给定的容器中。
 *
 * @param label 存储标签的容器引用
 */
inline void Kmer<k>::get_label(T_container_& label) const
{
    label.resize(k);//将大小调整为k
    /**
    由于DNA序列的每个碱基通常使用2位表示（即4种碱基：A, T, C,
    G），所以一个64位的uint64_t可以存储32个碱基的信息。packed_word_count表示可以完全使用64位单词来存储的碱基对的数量。
     */
    constexpr uint16_t packed_word_count = k / 32;

    // Get the fully packed words' representations.
    for(uint16_t data_idx = 0; data_idx < packed_word_count; ++data_idx)
        for(uint16_t bit_pair_idx = 0; bit_pair_idx < 32; ++bit_pair_idx)
            label[(k - 1) - ((data_idx << 5) + bit_pair_idx)] =
                DNA_Utility::map_char(static_cast<DNA::Base>((kmer_data[data_idx] & (0b11ULL << (2 * bit_pair_idx))) >> (2 * bit_pair_idx)));

    // Get the partially packed (highest index) word's representation.
    for(uint16_t bit_pair_idx = 0; bit_pair_idx < (k & 31); ++bit_pair_idx)
        label[(k - 1) - (((NUM_INTS - 1) << 5) + bit_pair_idx)] =
            DNA_Utility::map_char(static_cast<DNA::Base>((kmer_data[NUM_INTS - 1] & (0b11ULL << (2 * bit_pair_idx))) >> (2 * bit_pair_idx)));
}


template <uint16_t k>
inline Kmer<k>::operator std::string() const
{
    std::string label;
    get_label(label);

    return label;
}


template <uint16_t k>
inline Kmer<k> Kmer<k>::random_kmer()
{
    return Kmer<k>(get_random_string(k, "ACGT"));
}


template <uint16_t k>
std::ostream& operator<<(std::ostream& out, const Kmer<k>& kmer)
{
    out << kmer.string_label();

    return out;
}


template <uint16_t k>
template <uint8_t l>
inline typename Kmer<k>::minimizer_t Kmer<k>::minimizer() const
{
    // static_assert(l <= k);

    // TODO: SIMD?

    minimizer_t lmer = kmer_data[0] & ((1ULL << (2 * l)) - 1);
    minimizer_t minmzr = lmer;

    for(uint16_t idx = l; idx < k; ++idx)
    {
        const uint16_t word_idx = (idx >> 5);
        const uint16_t base_idx = (idx & 31);
        lmer =  (lmer >> 2) |
                (((kmer_data[word_idx] & (0b11ULL << (2 * base_idx))) >> (2 * base_idx)) << (2 * (l - 1)));

        if(minmzr > lmer)
            minmzr = lmer;
    }


    return minmzr;
}


template <uint16_t k>
template <uint8_t l>
inline typename Kmer<k>::minimizer_t Kmer<k>::minimizer(const std::vector<uint32_t>& order) const
{
    // static_assert(l <= k);

    // TODO: SIMD?

    minimizer_t lmer = kmer_data[0] & ((1ULL << (2 * l)) - 1);
    minimizer_t minmzr = lmer;

    for(uint16_t idx = l; idx < k; ++idx)
    {
        const uint16_t word_idx = (idx >> 5);
        const uint16_t base_idx = (idx & 31);
        lmer =  (lmer >> 2) |
                (((kmer_data[word_idx] & (0b11ULL << (2 * base_idx))) >> (2 * base_idx)) << (2 * (l - 1)));

        if(order[minmzr] > order[lmer])
            minmzr = lmer;
    }


    return minmzr;
}


template <uint16_t k>
template <uint8_t l>
inline void Kmer<k>::count_lmers(std::vector<uint64_t>& count) const
{
    // static_assert(l <= k);

    std::size_t lmer = kmer_data[0] & ((1ULL << (2 * l)) - 1);
    count[lmer]++;

    for(uint16_t idx = l; idx < k; ++idx)
    {
        const uint16_t word_idx = (idx >> 5);
        const uint16_t base_idx = (idx & 31);

        lmer =  (lmer >> 2) |
                (((kmer_data[word_idx] & (0b11ULL << (2 * base_idx))) >> (2 * base_idx)) << (2 * (l - 1)));

        count[lmer]++;
    }
}



#endif
