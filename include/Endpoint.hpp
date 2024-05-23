
#ifndef ENDPOINT_HPP
#define ENDPOINT_HPP




#include "Kmer.hpp"
#include "Directed_Vertex.hpp"
#include "globals.hpp"

#include <cstdint>


template <uint16_t k, uint8_t BITS_PER_KEY> class Kmer_Hash_Table;


// A class denoting an endpoint of a bidirected edge instance.
template <uint16_t k>
class Endpoint
{
private:

    Directed_Vertex<k> v;   // The endpoint vertex.
    cuttlefish::side_t s;   // The side of `v` to which the edge instance is incident to.
    // 与此端点相关的边缘实例的`DNA::Extended_Base`编码。
    cuttlefish::edge_encoding_t e;  // The `DNA::Extended_Base` encoding of the edge instance incident to this endpoint.


    // Constructs an `Endpoint` object that appears in the form `kmer` in an edge instance, and
    // is the source (i.e. prefix) of that edge iff `is_source` is true — which decides the edge
    // incidence side to the corresponding vertex. Also gets the hash value of the vertex using
    // the hash table `hash`. The sole application of this constructor is get a specific side of
    // a vertex where the edge incidence information is to be discarded, hence no edge-encoding
    // is provided with the constructor, although the class has such a member.
    Endpoint(const Kmer<k>& kmer, bool is_source, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the side of the associated vertex to which the edge instance corresponding to this
    // endpoint is incident to, if this endpoint is the source endpoint of the edge.
    cuttlefish::side_t exit_side() const;

    // Returns the side of the associated vertex to which the edge instance corresponding to this
    // endpoint is incident to, if this endpoint is the sink endpoint of the edge.
    cuttlefish::side_t entrance_side() const;

    // Returns the `DNA::Extended_Base` encoding of the edge `e` corresponding to this endpoint,
    // given the endpoint is the source endpoint of the edge.
    cuttlefish::edge_encoding_t exit_edge(const Kmer<k + 1>& e) const;

    // Returns the `DNA::Extended_Base` encoding of the edge `e` corresponding to this endpoint,
    // given the endpoint is the sink endpoint of the edge.
    cuttlefish::edge_encoding_t entrance_edge(const Kmer<k + 1>& e) const;


public:

    // Constructs an empty endpoint.
    Endpoint()
    {}

    // Configures the endpoint with the source (i.e. prefix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Configures the endpoint with the sink (i.e. suffix) k-mer of the edge (k + 1)-mer `e`;
    // and uses the hash table `hash` to get the hash value of the vertex.
    void from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns the neighboring endpoint of this endpoint that's connected with an edge encoded
    // with the code `e`, from the point-of-view of this endpoint. Uses the hash table `hash`
    // to get the hash value of the corresponding neighbor vertex.
    Endpoint neighbor_endpoint(cuttlefish::edge_encoding_t e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash) const;

    // Returns the canonical form of the associated vertex.
    const Kmer<k>& canonical() const;

    // Returns the side of the endpoint to which the corresponding edge is incident to.
    cuttlefish::side_t side() const;

    // Returns the `DNA::Extended_Base` encoding of the corresponding edge
    // incident to the endpoint.
    // 其实是返回e
    cuttlefish::edge_encoding_t edge() const;

    // Returns the hash value of the vertex associated to this endpoint.
    uint64_t hash() const;
};


template <uint16_t k>
/**
 * @brief Endpoint 类构造函数
 *
 * 用于创建一个 Endpoint 对象，根据给定的 Kmer 对象、是否为源节点以及哈希表进行初始化。
 *
 * @param kmer Kmer 对象，用于初始化 Endpoint 对象的 Kmer 成员
 * @param is_source 是否为源节点，决定 Endpoint 对象的 s 成员是入口侧还是出口侧
 * @param hash Kmer_Hash_Table 对象，用于初始化 Endpoint 对象的 v 成员
 */
inline Endpoint<k>::Endpoint(const Kmer<k>& kmer, const bool is_source, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash):
    v(kmer, hash),
    s(is_source ? exit_side() : entrance_side())
{}


template <uint16_t k>
/**
 * @brief 从前缀生成端点
 *
 * 使用给定的前缀和哈希表生成端点。
 *
 * @param e 前缀 Kmer 对象
 * @param hash Kmer 哈希表对象
 */
inline void Endpoint<k>::from_prefix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
  // 利用边获取前缀,然后计算前缀的canonical form,根据canonical_form计算hash值
  v.from_prefix(e, hash);
  // 判断存储的v(Directed_Vertex)的  Kmer<k> kmer_ == kmer_hat_ptr
  s = exit_side();
  // 根据 s == 1(back),返回 e的最后一个碱基编码，否则返回其互补碱基编码。
  this->e = exit_edge(e);
}


template <uint16_t k>
/**
 * @brief 从后缀生成 Endpoint 对象
 *
 * 使用给定的 Kmer 对象和 Kmer 哈希表，从其后缀生成 Endpoint 对象。
 *
 * @param e Kmer 对象，用于生成 Endpoint 对象的后缀
 * @param hash Kmer 哈希表，用于辅助生成 Endpoint 对象
 */
inline void Endpoint<k>::from_suffix(const Kmer<k + 1>& e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    v.from_suffix(e, hash);

    s = entrance_side();
    this->e = entrance_edge(e);
}

template <uint16_t k>
/**
 * @brief 获取退出方向
 *
 * 返回当前端点的退出方向。
 * 判断存储的v(Directed_Vertex)的  Kmer<k> kmer_ == kmer_hat_ptr
 * 返回的 bool 类型
 * 和论文描述一致,如果这个端点是边的起点,且 kmer形式与 cannonical 形式一致,则返回 back,也就是边从当前端点的back出来
 * @return 退出方向，以 cuttlefish::side_t 枚举类型表示
 */
inline cuttlefish::side_t Endpoint<k>::exit_side() const {
  return v.exit_side();
}

template <uint16_t k>
/**
 * @brief 获取入口side
 *
 * 返回 &kmer_ == kmer_hat_ptr 如果是 返回front = 0,否则back = 1
 *
 * @return 入口边的类型
 */
inline cuttlefish::side_t Endpoint<k>::entrance_side() const {
  return v.entrance_side();
}

template <uint16_t k>
/**
 * @brief 获取退出边的编码
 *
 * 根据给定的 Kmer<k + 1> 类型的边 e，返回其对应的退出边的编码。
 * 如果 s == back,返回 e的最后一个碱基编码，否则返回其互补碱基编码。
 * @param e Kmer<k + 1> 类型的边
 *
 * @return 退出边的编码
 */
inline cuttlefish::edge_encoding_t Endpoint<k>::exit_edge(const Kmer<k + 1>& e) const
{
    return DNA_Utility::map_extended_base(s == cuttlefish::side_t::back ?
                                            e.back() : DNA_Utility::complement(e.back()));
}


template <uint16_t k>
inline cuttlefish::edge_encoding_t Endpoint<k>::entrance_edge(const Kmer<k + 1>& e) const
{
    return DNA_Utility::map_extended_base(s == cuttlefish::side_t::front ?
                                            e.front() : DNA_Utility::complement(e.front()));
}


template <uint16_t k>
inline Endpoint<k> Endpoint<k>::neighbor_endpoint(const cuttlefish::edge_encoding_t e, const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash) const
{
    Kmer<k> kmer(canonical());

    if(s == cuttlefish::side_t::back)
    {
        kmer.roll_forward(e);
        return Endpoint<k>(kmer, false, hash);
    }
    
    kmer.roll_backward(e);
    return Endpoint<k>(kmer, true, hash);
}


template <uint16_t k>
inline const Kmer<k>& Endpoint<k>::canonical() const
{
    return v.canonical();
}


template <uint16_t k>
/**
 * @brief 获取端点的方向
 *
 * 返回端点的side。
 *
 * @return 端点的side
 */
inline cuttlefish::side_t Endpoint<k>::side() const
{
    return s;
}


template <uint16_t k>
/**
 * @brief 返回端点的边编码
 *
 * 这是一个内联函数，用于返回 Endpoint 模板类中 k 类型参数的端点的边编码。
 *
 * @return 端点的边编码
 */
inline cuttlefish::edge_encoding_t Endpoint<k>::edge() const
{
    return e;
}


template <uint16_t k>
/**
 * @brief 计算端点的哈希值
 *
 * 根据端点的值计算哈希值。
 *
 * @return 端点的哈希值
 */
inline uint64_t Endpoint<k>::hash() const
{
    return v.hash();
}



#endif
