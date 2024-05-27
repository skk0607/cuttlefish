
#ifndef MAXIMAL_UNITIG_SCRATCH_HPP
#define MAXIMAL_UNITIG_SCRATCH_HPP



#include "Unitig_Scratch.hpp"
#include "FASTA_Record.hpp"
#include "Character_Buffer.hpp"

#include <cstdint>
#include <cstddef>
#include <vector>

//一个保存临时数据的类，用于从覆盖它并在交汇点重叠的两个组成单元中构建最大单元。也就是说，最大单元g被分成两个单元`u_b`和`u_f`，在某个顶点`v`， `u_b`和`u_f`分别连接到`v`的前部和后部。单元的构建方式是路径从`v`开始。因此，文字形式的最大单位是' \bar(u_f) \glue_k u_b '(或其反向补码)。
// =============================================================================
// A class to keep scratch data for building maximal unitigs from two of its
// constituent unitigs that cover it and overlap at a meeting-point vertex.
// That is, the maximal unitig is split into two unitigs `u_b` and `u_f`, at
// some vertex `v`—`u_b` and `u_f` are connected to the front and to the back
// of `v`, respectively. The unitigs are built such that the paths start from
// `v`. Thus, the maximal unitig in literal form is `\bar(u_f) \glue_k u_b`
// (or its reverse complement).
template <uint16_t k>
class Maximal_Unitig_Scratch
{
private:

    Unitig_Scratch<k> unitig_back;  // The unitig `u_b` (see note above class body).
    Unitig_Scratch<k> unitig_front; // The unitig `u_f` (see note above class body).

    uint64_t id_;   // The unique ID of the maximal unitig.

    Unitig_Scratch<k>* cycle;   // Pointer to either `u_b` or `u_f`, whichever contains the maximal unitig in case of it's a cycle.


    // Returns whether the maximal unitig is in canonical form.
    bool is_canonical() const;


public:

    // Constructs an empty scratch space for the unitig.
    Maximal_Unitig_Scratch();

    // Returns the unitig scratch `u_b` or `u_f`, based on `s` (see note above
    // class body).
    Unitig_Scratch<k>& unitig(const cuttlefish::side_t s);

    // Returns the unique ID of the maximal unitig.
    uint64_t id() const;

    // Returns whether the maximal unitig is linear, i.e. it is a linear path
    // and not a Detached Chordless Cycle (DCC) in the underlying graph.
    bool is_linear() const;

    // Returns the hashes of the vertices of the unitig at side `s`.
    const std::vector<uint64_t>& unitig_hash(cuttlefish::side_t s) const;

    // Returns the hashes of the vertices in the maximal unitig in case it's
    // a DCC.
    const std::vector<uint64_t>& cycle_hash() const;

    // Returns the count of vertices in the maximal unitig.
    std::size_t size() const;

    // Returns the signature vertex of the maximal unitig, which is the first
    // vertex in the canonical form of the unitig.
    const Directed_Vertex<k>& sign_vertex() const;

    // Marks the maximal unitig as linear, i.e not a DCC.
    void mark_linear();

    // Marks the maximal unitig as a DCC, and signals that the cycle has been
    // extracted in the unitig scratch at side `s`.
    void mark_cycle(cuttlefish::side_t s);

    // Signals the scratch that the unitig pieces `u_b` and `u_f` are in their
    // final forms and will not be modified anymore. So it restructures the
    // maximal unitig so as to put its label in canonical form and sets its
    // unique ID.
    void finalize();

    // Returns `true` iff the maximal unitig has been marked as a cycle.
    bool is_cycle() const;

    // Returns a FASTA record of the maximal unitig (in canonical form).
    // Applicable when the maximal unitig is linear.
    const FASTA_Record<std::vector<char>> fasta_rec() const;

    // Adds a corresponding FASTA record for the maximal unitig into `buffer`.
    template <std::size_t CAPACITY, typename T_sink_> void add_fasta_rec_to_buffer(Character_Buffer<CAPACITY, T_sink_>& buffer) const;
};


template <uint16_t k>
/**
 * @brief 获取最大单元型(Unitig)的草稿对象
 *
 * 根据给定的方向参数，返回最大单元型(Unitig)的草稿对象。
 *
 * @param s 方向
 * @return Unitig_Scratch<k>& 最大单元型(Unitig)的草稿对象
 */
inline Unitig_Scratch<k>& Maximal_Unitig_Scratch<k>::unitig(const cuttlefish::side_t s)
{
    return s == cuttlefish::side_t::back ? unitig_back : unitig_front;
}


template <uint16_t k>
/**
 * @brief 判断是否为规范形式
 *
 * 判断 Maximal_Unitig_Scratch 对象是否为规范形式。
 * 如果front 扩展的最后一个顶点的 反向互补比 back 扩展的最后1个顶点的反向互补小，则返回 true；否则返回 false。
 * @return 如果为规范形式，则返回 true；否则返回 false。
 */
inline bool Maximal_Unitig_Scratch<k>::is_canonical() const
{
    return unitig_front.endpoint().kmer_bar() < unitig_back.endpoint().kmer_bar();
}


template <uint16_t k>
/**
 * @brief 判断是否为线性结构
 *
 * 判断当前 Maximal_Unitig_Scratch 对象是否为线性结构。
 *
 * @return 如果为线性结构则返回 true，否则返回 false
 */
inline bool Maximal_Unitig_Scratch<k>::is_linear() const
{
    return cycle == nullptr;
}


template <uint16_t k>
inline uint64_t Maximal_Unitig_Scratch<k>::id() const
{
    return id_;
}


template <uint16_t k>
/**
 * @brief 标记为线性
 *
 * 将当前对象标记为线性。
 * 线性表示对象不是循环结构。
 */
inline void Maximal_Unitig_Scratch<k>::mark_linear()
{
    cycle = nullptr;
}


template <uint16_t k>
inline const std::vector<uint64_t>& Maximal_Unitig_Scratch<k>::unitig_hash(const cuttlefish::side_t s) const
{
    return (s == cuttlefish::side_t::back ? unitig_back.hash() : unitig_front.hash());
}


template <uint16_t k>
inline const std::vector<uint64_t>& Maximal_Unitig_Scratch<k>::cycle_hash() const
{
    return cycle->hash();
}


template <uint16_t k>
inline std::size_t Maximal_Unitig_Scratch<k>::size() const
{
    return is_linear() ?    (unitig_back.size() + unitig_front.size() - 1) :
                            cycle->size();
}


template <uint16_t k>
/**
 * @brief 获取带符号顶点
 *
 * 返回当前 Maximal_Unitig_Scratch 对象中带符号的顶点。
 * 如果是线性,且为canpical形式，则返回 unitig_front 的 endpoint;
 * 如果是线性,且不是canpical形式，则返回 unitig_back 的 endpoint;
 * 如果是循环，则返回 cycle 的最小顶点。
 * @return 带符号的顶点引用
 */
inline const Directed_Vertex<k>& Maximal_Unitig_Scratch<k>::sign_vertex() const
{
    return is_linear() ?   (is_canonical() ? unitig_front.endpoint() : unitig_back.endpoint()) :
                            cycle->min_vertex();
}


template <uint16_t k>
inline void Maximal_Unitig_Scratch<k>::mark_cycle(const cuttlefish::side_t s)
{
    cycle = &(s == cuttlefish::side_t::back ? unitig_back : unitig_front);
}


template <uint16_t k>
/**
 * @brief 完成最大单元格的初始化
 *
 * 根据当前对象的状态，完成最大单元格的初始化操作。
 * 如果对象是线性的，则根据是否规范来决定其id和反转互补操作；
 * 如果对象是非线性的，则设置id为最小顶点的哈希值，并根据需要执行反转互补操作。
 */
inline void Maximal_Unitig_Scratch<k>::finalize()
{
    if(is_linear())
    {
        if(is_canonical())
            id_ = unitig_front.endpoint().hash(),//只存储标记顶点的哈希值为id
            unitig_front.reverse_complement();//存储更小的unititgs? 标准化存储?
        else
            id_ = unitig_back.endpoint().hash(),
            unitig_back.reverse_complement();
    }
    else
    {
        id_ = cycle->min_vertex().hash();
        if(!cycle->min_vertex().in_canonical_form())
            cycle->reverse_complement();
    }
}


template <uint16_t k>
/**
 * @brief 判断是否为环状结构
 *
 * 判断当前 Maximal_Unitig_Scratch 对象是否为环状结构。
 *
 * @return 如果为环状结构则返回 true，否则返回 false
 */
inline bool Maximal_Unitig_Scratch<k>::is_cycle() const
{
    return !is_linear();
}


template <uint16_t k>
/**
 * @brief 返回 FASTA 记录
 *
 * 返回一个表示当前对象的 FASTA 记录。如果当前对象是规范的，则记录中的序列标签为前向单位序列的标签；否则，记录中的序列标签为反向单位序列的标签。
 *
 * @return 返回一个包含 ID、序列标签和 k 值的 FASTA 记录对象
 */
inline const FASTA_Record<std::vector<char>> Maximal_Unitig_Scratch<k>::fasta_rec() const
{
    // 如果当前对象是规范的，则返回正向序列的FASTA记录
    return is_canonical() ?
            FASTA_Record<std::vector<char>>(id(), unitig_front.label(), unitig_back.label(), 0, k) :
    // 否则返回反向序列的FASTA记录
            FASTA_Record<std::vector<char>>(id(), unitig_back.label(), unitig_front.label(), 0, k);
}


template <uint16_t k>
template <std::size_t CAPACITY, typename T_sink_>
/**
 * @brief 将 Fasta 记录添加到缓冲区中
 *
 * 如果当前对象是线性的，则将 Fasta 记录添加到给定的缓冲区中；否则，以循环的方式将 Fasta 记录添加到缓冲区中。
 *
 * @param buffer 字符缓冲区引用
 */
inline void Maximal_Unitig_Scratch<k>::add_fasta_rec_to_buffer(Character_Buffer<CAPACITY, T_sink_>& buffer) const
{
    // 如果是线性结构
    if(is_linear())
        // 将FASTA记录添加到缓冲区中
        buffer += fasta_rec();
    else
        // 旋转追加循环结构的FASTA记录到缓冲区中
        buffer.template rotate_append_cycle<k>(FASTA_Record<std::vector<char>>(id(), cycle->label()), cycle->min_vertex_idx());
}



#endif
