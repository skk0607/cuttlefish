
#ifndef EDGE_HPP
#define EDGE_HPP




#include "globals.hpp"
#include "Endpoint.hpp"

#include <cstdint>


template <uint16_t k, uint8_t BITS_PER_KEY> class Kmer_Hash_Table;


// Class for an instance of a bidirected edge.
// NB: for some (k + 1)-mer `e`, `e` and `e_bar` denote the same bidirected edge
// `e_hat`; but these being different (k + 1)-mers, they are treated as
// different *instances* of the same edge. Semantically, the underlying edges
// are the same. This edge instance is in the tuple form `(u, s_\hat{u}, v,
// s_\hat{v})`.
// 双向边实例的类
//注:对于某些(k + 1)-mer ' e '， ' e '和' e_bar '表示同一条双向边' e_hat ';
//但由于这些是不同的(k + 1)-mers，它们被视为同一条边的不同*实例*。从语义上讲，底层的边是相同的。这个边缘实例的元组形式是`(u, s_\hat{u}， v, s_\hat{v})`。
template <uint16_t k>
class Edge
{
private:

    Kmer<k + 1> e_; // The edge (k + 1)-mer (need not be in canonical form).边(k + 1)-mer(不一定是canonical形)。
    Endpoint<k> u_; // One endpoint k-mer of this edge instance — source of the `e` form.该边缘实例的一个端点k-mer - `e`形式的source。
    Endpoint<k> v_; // One endpoint k-mer of this edge instance — sink of the `e` form. 该边缘实例的一个端点k-mer - e形式的sink。


public:

    // Returns a mutable reference to the edge (k + 1)-mer.
    Kmer<k + 1>& e();

    // Returns the source endpoint `u` of the edge instance.
    const Endpoint<k>& u() const;

    // Returns the sink endpoint `v` of the edge instance.
    const Endpoint<k>& v() const;

    // Configures the edge data, i.e. sets the relevant information of the
    // edge from the underlying (k + 1)-mer. Uses the hash table `hash` to
    // get the hash values of the endpoint vertices. Must be used whenever
    // the edge (k + 1)-mer (updatable using `e()`) is modified.
    // 配置边数据，即设置边缘实例的相关信息。使用哈希表` hash `来获取端点顶点的哈希值。必须在使用` e() `更新边缘(k + 1)-mer时使用。
    void configure(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash);

    // Returns `true` iff the edge is a loop.
    bool is_loop() const;
};


template <uint16_t k>
/**
 * @brief 获取 Edge<k> 对象的 e_ 成员
 *
 * 返回一个对 Edge<k> 对象的 e_ 成员的引用，该成员是 Kmer<k + 1> 类型的对象。
 *
 * @return Kmer<k + 1>& Edge<k> 对象的 e_ 成员的引用
 */
inline Kmer<k + 1>& Edge<k>::e()
{
    return e_;
}


template <uint16_t k>
inline const Endpoint<k>& Edge<k>::u() const
{
    return u_;
}


template <uint16_t k>
inline const Endpoint<k>& Edge<k>::v() const
{
    return v_;
}


template <uint16_t k>
/**
 * @brief 配置边对象
 *
 * 使用给定的 Kmer 哈希表来配置边对象。根据边的哈希值，从哈希表中获取前缀和后缀信息，
 * 并将其分别设置到边的起始顶点和结束顶点中。
 *
 * @param hash Kmer 哈希表引用
 */
inline void Edge<k>::configure(const Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash)
{
    u_.from_prefix(e_, hash),
    v_.from_suffix(e_, hash);
}


template <uint16_t k>
inline bool Edge<k>::is_loop() const
{
    return u_.canonical() == v_.canonical();
}



#endif
