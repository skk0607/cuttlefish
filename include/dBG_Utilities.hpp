
#ifndef DBG_UTILITIES_HPP
#define DBG_UTILITIES_HPP



#include "DNA_Utility.hpp"
#include "globals.hpp"

#include <algorithm>


// =============================================================================
namespace cuttlefish
{
    // Returns `true` iff the edge encoding `e` is fuzzy, i.e. a unique encoding
    // is not known for the corresponding edge(s).
    bool is_fuzzy_edge(const edge_encoding_t e);

    // Returns the opposite (or complement) side of the vertex-side `s`.
    side_t opposite_side(const side_t s);


    // Replaces the sequence `seq` in-place with its reverse complement.
    template <typename T_container_> void reverse_complement(T_container_& seq);
}


/**
 * @brief 判断是否为模糊边
 *
 * 根据给定的边编码判断是否为模糊边。
 * 只要不是ACGT就返回True，否则返回False。
 * @param e 边编码
 *
 * @return 如果为模糊边则返回true，否则返回false
 */
inline bool cuttlefish::is_fuzzy_edge(const edge_encoding_t e)
{
    return e == edge_encoding_t::N || e == edge_encoding_t::E;
}


/**
 * @brief 获取对立的侧面
 *
 * 根据给定的侧面，返回其对立面的侧面。
 *
 * @param s 侧面类型
 *
 * @return 对立面的侧面类型
 */
inline cuttlefish::side_t cuttlefish::opposite_side(const side_t s)
{
    return s == side_t::back ? side_t::front : side_t::back;
}


template <typename T_container_>
/**
 * @brief 反转并补全序列
 *
 * 将给定的序列进行反转，并对每个碱基进行补全操作。
 *
 * @param seq 序列容器引用
 *
 * @note 序列必须非空
 */
inline void cuttlefish::reverse_complement(T_container_& seq)
{
    assert(!seq.empty());

    auto fwd = seq.begin();
    auto bwd = seq.end() - 1;

    for(; fwd < bwd; ++fwd, --bwd)
    {
        std::swap(*fwd, *bwd);
        
        *fwd = DNA_Utility::complement(*fwd),
        *bwd = DNA_Utility::complement(*bwd);
    }

    if(fwd == bwd)
        *fwd = DNA_Utility::complement(*fwd);
}



#endif
