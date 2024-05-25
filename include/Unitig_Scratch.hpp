
#ifndef UNITIG_SCRATCH_HPP
#define UNITIG_SCRATCH_HPP



#include "Directed_Vertex.hpp"
#include "dBG_Utilities.hpp"

#include <cstddef>
#include <cstdint>
#include <vector>


// =============================================================================
// A class to keep scratch data (i.e. working space) for unitigs.
// 一个为单元保存临时数据(即工作空间)的类。
template <uint16_t k>
class Unitig_Scratch
{
private:

    // 100K (soft limit) unitig vertices can be retained in memory, at most, before reallocations.
    static constexpr std::size_t BUFF_SZ = 100 * 1024UL;

    Directed_Vertex<k> anchor;      // The anchor vertex of the unitig traversal.
    // (虽然也可以从 "定向 "顶点推断，但扩展的侧边将由客户端代码处理）。
    Directed_Vertex<k> endpoint_;   // The current end of the unitig through which farther extensions can be done.
                                    // (The side for the extension is to be handled by the client code, although can
                                    // also be inferred from the "directed" vertex.)
    // unitigs中按字典序排列的最小顶点。
    Directed_Vertex<k> min_vertex_; // The lexicographically minimum vertex in the unitig.
    std::size_t vertex_idx;         // Index of the vertex in the path being traversed.
    std::size_t min_v_idx;          // Index of the lexicographically minimum vertex in the path.
    // unitigs的实际标签。
    std::vector<char> label_;       // Literal label of the unitig.
    // unitigs连续顶点的hash值
    std::vector<uint64_t> hash_;    // Hashes of the constituent vertices of the unitig.
    bool is_cycle_;                 // Whether the unitig is cyclical or not.


    // Clears the scratch data.
    void clear();

public:

    // Constructs an empty unitig scratch.
    Unitig_Scratch();

    // Initializes the unitig scratch with the vertex `v`.
    void init(const Directed_Vertex<k>& v);

    // Extends the unitig scratch with the vertex `v`, and its literal form
    // with the symbol `b`. Returns `true` iff adding `v` to the unitig does
    // not render itself a cycle.
    bool extend(const Directed_Vertex<k>& v, char b);

    // Reverse complements the unitig.
    void reverse_complement();

    // Returns the literal label of the unitig.
    const std::vector<char>& label() const;

    // Returns the hash collection of the unitig vertices.
    const std::vector<uint64_t>& hash() const;

    // Returns the current extension-end vertex of the unitig.
    const Directed_Vertex<k>& endpoint() const;

    // Returns the count of vertices in this unitig.
    std::size_t size() const;

    // Returns `true` iff unitig is a cycle.
    bool is_cycle() const;

    // Returns the lexicographically minimum vertex in the unitig.
    const Directed_Vertex<k>& min_vertex() const;

    // Returns the index of the lexicographically minimum vertex in the unitig.
    std::size_t min_vertex_idx() const;
};


template <uint16_t k>
/**
 * @brief 清空 Unitig_Scratch 对象
 *
 * 将 Unitig_Scratch 对象中的 label_ 和 hash_ 成员变量清空。
 */
inline void Unitig_Scratch<k>::clear()
{
    label_.clear();//这里面实现都是vector
    hash_.clear();
}


template <uint16_t k>
/**
 * @brief 初始化 Unitig_Scratch 对象
 *
 * 使用给定的 Directed_Vertex 对象初始化 Unitig_Scratch 对象。
 *
 * @param v Directed_Vertex 对象
 */
inline void Unitig_Scratch<k>::init(const Directed_Vertex<k>& v)
{
    clear();//清空label_和hash_的vector
    //初始化时,anchor = endpoint_ = min_vertex_
    min_vertex_ = endpoint_ = anchor = v;
    //顶点索引以及最小顶点索引都初始化为0
    min_v_idx = vertex_idx = 0;
    //这里的label_开始是空的
    endpoint_.kmer().get_label(label_);
    //hash_先放入当前顶点的hash值
    hash_.emplace_back(endpoint_.hash());
    //循环标志位设置为false
    is_cycle_ = false;
}


template <uint16_t k>
/**
 * @brief 扩展unitigs
 *
 * 在给定的unitigs中，根据给定的有向顶点和字符进行扩展。
 *
 * @param v 有向顶点对象
 * @param b 扩展的字符
 *
 * @return 如果扩展成功，返回 true；如果扩展到锚点形成循环，返回 false
 */
inline bool Unitig_Scratch<k>::extend(const Directed_Vertex<k>& v, const char b)
{
    //到达原来扩展的锚点
    if(v.is_same_vertex(anchor))
    {
        is_cycle_ = true;
        return false;
    }


    endpoint_ = v;//当前的点endpoint_ 
    vertex_idx++;//顶点索引加1
    // 获取最小字典序的顶点以及对应索引
    if(min_vertex_.canonical() > endpoint_.canonical())
        min_vertex_ = endpoint_,
        min_v_idx = vertex_idx;

    label_.emplace_back(b);//push进vector
    hash_.emplace_back(endpoint_.hash());

    return true;
}


template <uint16_t k>
/**
 * @brief 反转并补全
 *
 * 将 Unitig_Scratch 对象中的标签序列进行反转并补全操作。
 *
 * @note 该函数会改变对象的内部状态。
 */
inline void Unitig_Scratch<k>::reverse_complement()
{
    cuttlefish::reverse_complement(label_);
    min_v_idx = (hash_.size() - 1 - min_v_idx);
}


template <uint16_t k>
inline const std::vector<char>& Unitig_Scratch<k>::label() const
{
    return label_;
}


template <uint16_t k>
inline const std::vector<uint64_t>& Unitig_Scratch<k>::hash() const
{
    return hash_;
}


template <uint16_t k>
/**
 * @brief 获取端点
 *
 * 返回 Unitig_Scratch 对象的端点。
 *
 * @return 指向 Directed_Vertex<k> 类型的常量引用，表示端点
 */
inline const Directed_Vertex<k>& Unitig_Scratch<k>::endpoint() const
{
    return endpoint_;
}


template <uint16_t k>
inline std::size_t Unitig_Scratch<k>::size() const
{
    return hash_.size();
}


template <uint16_t k>
inline bool Unitig_Scratch<k>::is_cycle() const
{
    return is_cycle_;
}


template <uint16_t k>
inline const Directed_Vertex<k>& Unitig_Scratch<k>::min_vertex() const
{
    return min_vertex_;
}


template <uint16_t k>
inline std::size_t Unitig_Scratch<k>::min_vertex_idx() const
{
    return min_v_idx;
}



#endif
