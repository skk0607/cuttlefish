
#ifndef DIRECTED_VERTEX_HPP
#define DIRECTED_VERTEX_HPP



#include "Kmer.hpp"
#include "globals.hpp"
#include "Kmer_Hash_Table.hpp"

#include <cstdint>
#include <iostream>


// A class denoting an instance of a vertex. It's "directed" in the sense that
// the k-mer observed for the vertex is in a particular orientation — although a
// vertex `v` has an unambiguous canonical k-mer `v_hat`, the vertex can be
// observed in two different k-mer forms: `v_hat` and `{v_hat}_bar` — the class
// keeps track of the particular k-mer form observed for the vertex instance.
// 表示顶点实例的类。从某种意义上说，它是“有向的”，因为顶点观察到的k-mer具有特定的方向——尽管顶点`v`具有明确的标准k-mer `v_hat`，但该顶点可以观察到两种不同的k-mer形式:`v_hat`和`{v_hat}_bar`——该类跟踪顶点实例观察到的特定k-mer形式。
template <uint16_t k>
class Directed_Vertex
{
private:

    Kmer<k> kmer_;  // The observed k-mer for the vertex.
    Kmer<k> kmer_bar_;  // Reverse complement of the k-mer observed for the vertex.
    const Kmer<k>* kmer_hat_ptr;    // Pointer to the canonical form of the k-mer associated to the vertex.
    uint64_t h; // Hash value of the vertex, i.e. hash of the canonical k-mer.

    // Initialize the data of the class once the observed k-mer `kmer_` is set.
    void init(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);


public:

    // Constructs an empty vertex.
    Directed_Vertex()
    {}

    // Constructs a vertex observed for the k-mer `kmer`. Gets the hash value of the vertex using
    // the hash table `hash`.
    Directed_Vertex(const Kmer<k>& kmer, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Copy constructs the vertex from `rhs`.
    Directed_Vertex(const Directed_Vertex<k>& rhs);

    // Assigns the vertex `rhs` to this one, and returns a constant reference to this object.
    const Directed_Vertex<k>& operator=(const Directed_Vertex<k>& rhs);

    // Returns `true` iff the k-mer observed for the vertex is in its canonical form.
    bool in_canonical_form() const;

    // Configures the vertex with the k-mer `v`, and uses the hash table `hash` to get the
    // hash value of the vertex.
    void from_kmer(const Kmer<k>& v, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Configures the vertex with the source (i.e. prefix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Configures the vertex with the sink (i.e. suffix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the observed k-mer for the vertex.
    const Kmer<k>& kmer() const;

    // Returns the reverse complement of the observed k-mer for the vertex.
    const Kmer<k>& kmer_bar() const;

    // Returns the canonical form of the vertex.
    const Kmer<k>& canonical() const;

    // Returns the hash value of the vertex.
    uint64_t hash() const;

    // Transforms this vertex to another by chopping off the first base from the associated
    // observed k-mer, and appending the nucleobase `b` to the end, i.e. effecitively
    // rolling the associated k-mer by one base "forward". The hash table `hash` is used
    // to get the hash value of the new vertex.
    void roll_forward(cuttlefish::base_t b, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the side of the vertex which is to be the incidence side of some bidirected
    // edge instance if this vertex instance were to be the source vertex (i.e. prefix k-mer)
    // of that edge.
    cuttlefish::side_t exit_side() const;

    // Returns the side of the vertex which is to be the incidence side of some bidirected
    // edge instance if this vertex instance were to be the sink vertex (i.e. suffix k-mer)
    // of that edge.
    cuttlefish::side_t entrance_side() const;

    // Returns `true` iff this vertex and the vertex `v` are the same vertex, without the
    // directionality.
    bool is_same_vertex(const Directed_Vertex<k>& v) const;
};


template <uint16_t k>
/**
 * @brief 初始化有向顶点
 *
 * 使用给定的 Kmer 哈希表初始化有向顶点对象。
 * 这里是比较Kmer和Kmer的反向互补,然后选择1个较小的存入当前成员变量h中
 * @param hash Kmer 哈希表引用
 */
inline void Directed_Vertex<k>::init(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_bar_.as_reverse_complement(kmer_);
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);
    // 计算canonical的hash值,然后存入对象的属性中
    // TODO: 是否存入了hash表?
    h = hash(*kmer_hat_ptr);
   // std::cout<<kmer_.string_label()<<"的hash值是"<<h<<std::endl;
}


template <uint16_t k>
/**
 * @brief 构造有向顶点
 *
 * 使用给定的 Kmer 和 Kmer 哈希表构造有向顶点对象。
 *
 * @param kmer Kmer 对象
 * @param hash Kmer 哈希表引用
 */
inline Directed_Vertex<k>::Directed_Vertex(const Kmer<k>& kmer, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash):
    kmer_(kmer)
{
    init(hash);
}


template <uint16_t k>
inline Directed_Vertex<k>::Directed_Vertex(const Directed_Vertex<k>& rhs):
    kmer_(rhs.kmer_),
    kmer_bar_(rhs.kmer_bar_),
    kmer_hat_ptr(rhs.kmer_hat_ptr == &rhs.kmer_ ? &kmer_ : &kmer_bar_),
    h(rhs.h)
{}


template <uint16_t k>
inline const Directed_Vertex<k>& Directed_Vertex<k>::operator=(const Directed_Vertex<k>& rhs)
{
    kmer_ = rhs.kmer_;
    kmer_bar_ = rhs.kmer_bar_;
    kmer_hat_ptr = (rhs.kmer_hat_ptr == &rhs.kmer_ ? &kmer_ : &kmer_bar_);
    h = rhs.h;

    return *this;
}


template <uint16_t k>
inline bool Directed_Vertex<k>::in_canonical_form() const
{
    return &kmer_ == kmer_hat_ptr;
}


template <uint16_t k>
inline void Directed_Vertex<k>::from_kmer(const Kmer<k>& v, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_ = v;
    init(hash);
}


template <uint16_t k>
/**
 * @brief 根据前缀生成有向顶点的表示
 *
 * 使用给定的长度为 k+1 的 Kmer 对象和 Kmer 哈希表，生成有向顶点的表示。
 *
 * @param e 长度为 k+1 的 Kmer 对象
 * @param hash Kmer 哈希表
 */
inline void Directed_Vertex<k>::from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_.from_prefix(e);//根据边获取前缀存储在成员变量kmer_中
    init(hash);//利用from_prefix获取的成员变量信息
}


template <uint16_t k>
inline void Directed_Vertex<k>::from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    kmer_.from_suffix(e);
    init(hash);
}


template <uint16_t k>
/**
 * @brief 获取 kmer 值
 *
 * 返回当前有向顶点中的 kmer 值。
 * 实际就是返回kmer_
 * @return 返回一个 const Kmer<k>& 类型的 kmer 值引用
 */
inline const Kmer<k>& Directed_Vertex<k>::kmer() const
{
    return kmer_;
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::kmer_bar() const
{
    return kmer_bar_;
}


template <uint16_t k>
inline const Kmer<k>& Directed_Vertex<k>::canonical() const
{
    return *kmer_hat_ptr;
}


template <uint16_t k>
inline uint64_t Directed_Vertex<k>::hash() const
{
    return h;
}


template <uint16_t k>
/**
 * @brief 向前滚动顶点
 *
 * 向前滚动当前顶点，更新其 k-mer 信息，并根据给定的哈希表重新计算哈希值。
 * 
 * @param b 下一个碱基
 * @param hash Kmer 哈希表引用
 */
inline void Directed_Vertex<k>::roll_forward(const cuttlefish::base_t b, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    //此时kmer_和kmer_bar都移动了2bit为了获取下一个kmer
    kmer_.roll_to_next_kmer(b, kmer_bar_);
    //比较获得 cannonical 的指针
    kmer_hat_ptr = Kmer<k>::canonical(kmer_, kmer_bar_);
    //计算hash值
    h = hash(*kmer_hat_ptr);
}


template <uint16_t k>
/**
 * @brief 返回有向顶点的出口边类型
 *
 * 判断当前有向顶点的kmer指针是否等于kmer_hat_ptr指针，根据判断结果返回出口边类型。
 * 如果等于，则返回cuttlefish::side_t::back，否则返回cuttlefish::side_t::front。
 *
 * @return 出口边类型，为cuttlefish::side_t枚举类型
 */
inline cuttlefish::side_t Directed_Vertex<k>::exit_side() const
{
    return &kmer_ == kmer_hat_ptr ? cuttlefish::side_t::back : cuttlefish::side_t::front;
}


template <uint16_t k>
/**
 * @brief 获取有向顶点的入口侧
 *
 * 根据当前顶点的 kmer 值和 kmer_hat_ptr 的指向关系，判断顶点的入口侧是前方还是后方，并返回对应的枚举值。
 *
 * @return 入口侧的枚举值
 */
inline cuttlefish::side_t Directed_Vertex<k>::entrance_side() const
{
    return &kmer_ == kmer_hat_ptr ? cuttlefish::side_t::front : cuttlefish::side_t::back;
}


template <uint16_t k>
inline bool Directed_Vertex<k>::is_same_vertex(const Directed_Vertex<k>& v) const
{
    return hash() == v.hash();
}



#endif
