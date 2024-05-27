
#ifndef READ_CDBG_EXTRACTOR_HPP
#define READ_CDBG_EXTRACTOR_HPP



#include "globals.hpp"
#include "DNA_Utility.hpp"
#include "Kmer_Hash_Table.hpp"
#include "Directed_Vertex.hpp"
#include "Maximal_Unitig_Scratch.hpp"
#include "Build_Params.hpp"
#include "Spin_Lock.hpp"
#include "Async_Logger_Wrapper.hpp"
#include "Output_Sink.hpp"
#include "Unipaths_Meta_info.hpp"
#include "Progress_Tracker.hpp"

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>


// Forward declarations.
template <uint16_t k> class Kmer_SPMC_Iterator;
template <uint16_t k> class Thread_Pool;


// A class to extract the vertices from a compacted de Bruin graph — which are the maximal unitigs of some ordinary de Bruijn graph.
template <uint16_t k>
class Read_CdBG_Extractor
{
    friend class Thread_Pool<k>;

private:

    const Build_Params params;  // Required parameters (wrapped inside).
    Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table; // Hash table for the vertices (i.e. canonical k-mers) of the original (uncompacted) de Bruijn graph.

    // typedef std::ofstream sink_t;
    typedef Async_Logger_Wrapper sink_t;
    // maximal unitigs的sink,实际可以理解为1个日志写入器
    Output_Sink<sink_t> output_sink;    // Sink for the output maximal unitigs.

    // TODO: give these limits more thoughts, especially their exact impact on the memory usage.
    static constexpr std::size_t BUFF_SZ = 100 * 1024ULL;   // 100 KB (soft limit) worth of maximal unitig records (FASTA) can be retained in memory, at most, before flushing.

    mutable uint64_t vertices_scanned = 0;    // Total number of vertices scanned from the database.
    mutable Spin_Lock lock; // Mutual exclusion lock to access various unique resources by threads spawned off this class' methods.

    mutable uint64_t vertices_marked = 0;   // Total number of vertices marked as present in maximal unitigs; used for the extraction of detached chordless cycle(s), if any.
    
    Unipaths_Meta_info<k> unipaths_meta_info_;  // Meta-information over the extracted maximal unitigs.

    Progress_Tracker progress_tracker;  // Progress tracker for the maximal unitigs extraction task.


    // Distributes the maximal unitigs extraction task — disperses the graph vertices (i.e. k-mers)
    // parsed by the parser `vertex_parser` to the worker threads in the thread pool `thread_pool`,
    // for the unitpath-flanking vertices to be identified and the corresponding unipaths to be extracted.
    void distribute_unipaths_extraction(Kmer_SPMC_Iterator<k>* vertex_parser, Thread_Pool<k>& thread_pool);

    // Prcesses the vertices provided to the thread with id `thread_id` from the parser
    // `vertex_parser`, i.e. for each vertex `v` provided to that thread, attempts to
    // piece-wise construct its containing maximal unitig.
    void process_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Extracts the maximal unitig `p` that contains the vertex `v_hat`, and `maximal_unitig` is
    // used as the working scratch for the extraction, i.e. to build and store the two unitigs
    // connecting to the two sides of `v_hat`. Returns `true` iff the extraction is successful,
    // which happens when `p` is attempted for output-marking first by this thread.
    bool extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig);

    // Traverses a unitig starting from the vertex `v_hat`, exiting it through the side `s_v_hat`.
    // The DFA of `v_hat` is supposed to have the state `st_v`. `unitig` is used as the working
    // scratch to build the unitig. Returns `true` iff the unitig could have been traversed
    // maximally up-to its endpoint in the direction of the walk from `v_hat`, which is possible
    // iff no other thread output-marks it in the meantime.
    bool walk_unitig(const Kmer<k>& v_hat, State_Read_Space st_v, cuttlefish::side_t s_v_hat, Unitig_Scratch<k>& unitig);

    // Marks all the vertices which have their hashes present in `path_hashes` as outputted.
    void mark_path(const std::vector<uint64_t>& path_hashes);

    // Marks all the vertices in the constituent unitigs of `maximal_unitig` as outputted.
    void mark_maximal_unitig(const Maximal_Unitig_Scratch<k>& maximal_unitig);

    // Marks the vertex `v` as outputted. Returns `true` iff `v` has not been marked yet and the hash
    // table update is successful.
    bool mark_vertex(const Directed_Vertex<k>& v);

    // Initializes the output sink, corresponding to the file `output_file_path`.
    void init_output_sink(const std::string& output_file_path);

    // Closes the output sink.
    void close_output_sink();


public:

    // Constructs a vertex-extractor object for some compacted read de Bruijn graph, with the required
    // parameters wrapped inside `params`, and uses the Cuttlefish hash table `hash_table`.
    // 为一些压缩的read de Bruijn图构造一个顶点提取器对象，将所需的参数包装在`params`中，并使用乌贼哈希表`hash_table`。
    Read_CdBG_Extractor(const Build_Params& params, Kmer_Hash_Table<k, cuttlefish::BITS_PER_READ_KMER>& hash_table);

    // Extracts the maximal unitigs of the de Bruijn graph with the vertex set at path prefix `vertex_db_path`,
    // into the output file at `output_file_path`.
    void extract_maximal_unitigs(const std::string& vertex_db_path, const std::string& output_file_path);    

    // Returns the parameters collection for the compacted graph construction.
    const Build_Params& get_params() const;

    // Returns a wrapper over the meta-information of the extracted unitigs.
    const Unipaths_Meta_info<k>& unipaths_meta_info() const;

    // Returns the number of vertices in the underlying graph.
    uint64_t vertex_count() const;


// The following stuffs are not used anymore with the current algorithm.
/*
private:

    static constexpr std::size_t SEQ_SZ = 1 * 1024ULL * 1024ULL;    // 1 MB (soft limit) sized maximal unitig, at most, is constructed at a time.

    // Scans the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`
    // for potential unipath-flanking vertices, i.e. for each vertex `v` provided to that thread,
    // identifies whether it is a unipath-flanking vertex, and if it is, then piece-wise constructs
    // the corresponding unipath.
    void scan_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Extracts the maximal unitig `p` that is flanked by the vertex `v_hat` and connects to `v_hat`
    // through its side `s_v_hat`. Returns `true` iff the extraction is successful, which happens when
    // the maximal unitig is encountered and attempted for output-marking _first_, by some thread. If
    // the attempt is successful, then the maximal unitig is extracted in its canonical form into
    // `unipath` (it is overwritten), a unique ID for it is put in `id`, and the hashes of the vertices
    // constituting the path overwrites `path_hashes` (when the user-option is specified). If not,
    // `unipath` and `path_hashes` may contain partial form of the path, and `id` is unaltered.
    bool extract_maximal_unitig(const Kmer<k>& v_hat, cuttlefish::side_t s_v_hat, uint64_t& id, std::vector<char>& unipath, std::vector<uint64_t>& path_hashes);

    // Marks all the vertices that are present in the maximal unitigs of the graph with its vertex
    // set being present at the path prefix `vertex_db_path`.
    void mark_maximal_unitig_vertices(const std::string& vertex_db_path);

    // Scans the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`
    // for potential unipath-flanking vertices. If a vertex `v` is found to be a flanking one, then
    // piece-wise constructs the corresponding (partial) maximal unitig starting the traversal from
    // `v`, and marks the vertices along the way. Premature halts before traversing the entire unitig
    // `p` is possible, in cases when some other thread is concurrently constructing `p`, but from
    // the opposite flank — the halt happens at the threads' meeting-point.
    void mark_maximal_unitig_vertices(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Marks (partially) the vertices of the maximal unitig `p` that is flanked by the vertex `v_hat`
    // from one side and connects to `v_hat` through its side `s_v_hat`. `p` might not be marked
    // completely by this thread if some other thread is concurrently traversing `p`, but from the
    // opposite flank. However, together these two threads mark `p` completely. Also returns the
    // number of vertices marked in this execution.
    std::size_t mark_maximal_unitig(const Kmer<k>& v_hat, cuttlefish::side_t s_v_hat);

    // Extracts all the detached chordless cycles present in the graph with its vertex set being
    // present at the path prefix `vertex_db_path`, into the output file at `output_file_path`.
    void extract_detached_chordless_cycles(const std::string& vertex_db_path, const std::string& output_file_path);

    // Scans the vertices provided to the thread with id `thread_id` from the parser `vertex_parser`
    // for potential detached chordless cycles. If a vertex `v` is found to be not marked as present
    // in the earlier extracted maximal unitigs, then it implies — by definition from the Cuttlefish
    // algorithm — that `v` belongs to a chordless cycle that is detached completely from the rest of
    // the graph. The method piece-wise constructs the cycle, starting the traversal from `v`.
    void extract_detached_chordless_cycles(Kmer_SPMC_Iterator<k>* vertex_parser, uint16_t thread_id);

    // Extracts the detached chordless cycle `p` that contains the vertex `v_hat`. Returns `true` iff
    // the extraction is successful, which happens when the cycle is encountered and attempted for
    // output-marking _first_ by this thread. If the attempt is successful, then the cycle is extracted
    // in its literal form that starts with `v_hat`'s canonical representation, into `cycle` (it is
    // overwritten); also, a unique ID for it is put in `id`. If not, `cycle` may contain partial form
    // of the unitig, and `id` is unaltered. The index of the lexicographically lowest (canonical) k-mer
    // in `cycle` is recorded at `pivot`.
    bool extract_cycle(const Kmer<k>& v_hat, uint64_t& id, std::vector<char>& cycle, std::size_t& pivot);

    // Marks the two endpoint vertices of a maximal unitig `p` as outputted: the first vertex in the
    // canonical form of `p`, `sign_vertex`, and the last vertex in the form, `cosign_vertex`. Returns
    // `true` iff the vertices have not been marked yet and the corresponding hash table updates are
    // successful.
    bool mark_flanking_vertices(const Directed_Vertex<k>& sign_vertex, const Directed_Vertex<k>& cosign_vertex);

    // Note: The following methods are only applicable when the heuristic of information-discarding
    // from branching vertices to their neighbors has been implemented in the DFA states computation
    // phase. In the general case, these functions with their specified input parameters and their
    // intended output values can not possibly be defined. Given just the state of a vertex, unless
    // for flanking vertices that are branching, it's not possible to determine —
    //      i.   whether it is a maximal unitig flanking vertex;
    //      ii.  whether some specific side of it is a flanking side;
    //      iii. and which side of it may connect to the containing maximal unitig.
    // The vertex itself, along with the hash table containing the DFA states, are required in general
    // — so that the neighboring vertices can also be probed to answer these queries. Nevertheless,
    // given that the heuristic of propagation of information-discarding from branching vertices has
    // been implemented in the DFA states computation phase, the following method definitions can
    // correctly respond to the queries given just the state. This is because the heuristic transforms
    // the flanking non-branching vertices into branching ones. Thus, although their states are not
    // technically "correct" as per the theoretical model — we are throwing away more information from
    // the model than it already is doing — this does not affect the output for the purposes of maximal
    // unitigs extraction.

    // Returns `true` iff some vertex `v` with the provided state `state` is a flanking vertex for the
    // maximal unitig `p` containing it. If yes, then stores the side of `v` to `unipath_side` through
    // which `v` is connected to `p`. If `p` is trivial, i.e. `p = v`, then the returned side is `back`,
    // —  necessitated by an optimization for the extraction of unipaths in their canonical forms.
    // NB: this method is only applicable if the heuristic of information-propagation from branching vertices
    // has been implemented in the DFA states computation phase. See the detailed comment in the class body.
    static bool is_flanking_state(State_Read_Space state, cuttlefish::side_t& unipath_side);

    // Returns `true` iff the vertex-side `side` for a vertex with state `state` flanks the maximal unitig
    // containing the vertex.
    // NB: this method is only applicable if the heuristic of information-propagation from branching vertices
    // has been implemented in the DFA states computation phase. See the detailed comment in the class body.
    static bool is_flanking_side(State_Read_Space state, cuttlefish::side_t side);


public:

    // Extracts the chordless cycles from the de Bruijn graph that are completely disconnected from the
    // rest of the graph. The graph is to contain its vertex set at the path prefix `vertex_db_path`,
    // and the cycles are appeneded to the output file at `output_file_path`. `dbg_info` is used to
    // determine whether the compacted graph had been constructed earlier—in which case some data
    // structures are re-used from the earlier construction.
    void extract_detached_cycles(const std::string& vertex_db_path, const std::string& output_file_path, const dBG_Info<k>& dbg_info);

    // Returns `true` iff the de Bruijn graph has DCCs (Detached Chordless Cycles).
    bool has_dcc() const;

    // Returns the number of vertices present in maximal unitigs (excluding the DCCs).
    uint64_t unipaths_vertex_count() const;
*/
};


template <uint16_t k>
/**
 * @brief 标记顶点
 *
 * 在 CdBG 提取器中标记给定的有向顶点。
 *
 * @param v 有向顶点对象
 *
 * @return 如果顶点已被标记过，则返回 false；否则返回 true
 */
inline bool Read_CdBG_Extractor<k>::mark_vertex(const Directed_Vertex<k>& v)
{
    //根据标志顶点的hash值获取Entry_API
    Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER> bucket = hash_table[v.hash()];
    //这个桶就是用来找对应state的
    State_Read_Space& state = bucket.get_state();
    //如果两端已经确定了是非分支,则直接返回false
    if(state.is_outputted())
        return false;
    //标记已输出
    state.mark_outputted();
    // 凡是更新 state,都先用拷贝,然后更新
    return hash_table.update(bucket);
}


template <uint16_t k>
/**
 * @brief 标记路径
 *
 * 将给定的路径哈希值标记为已输出状态。
 *
 * @param path_hashes 路径哈希值向量
 */
inline void Read_CdBG_Extractor<k>::mark_path(const std::vector<uint64_t>& path_hashes)
{   //每个顶点都标记为已输出
    for(const uint64_t hash: path_hashes)
        hash_table.update(hash, State_Read_Space::mark_outputted);
}


template <uint16_t k>
/**
 * @brief 标记最大单元体
 *
 * 根据给定的最大单元体信息，标记相应的路径。
 *
 * @param maximal_unitig 最大单元体信息
 */
inline void Read_CdBG_Extractor<k>::mark_maximal_unitig(const Maximal_Unitig_Scratch<k>& maximal_unitig)
{
    // 如果最大单元体是线性的
    if(maximal_unitig.is_linear())
    {
        // 标记后向路径
        mark_path(maximal_unitig.unitig_hash(cuttlefish::side_t::back));
        // 标记前向路径
        mark_path(maximal_unitig.unitig_hash(cuttlefish::side_t::front));
    }
    else
        // 如果最大单元体是环形的，则标记循环哈希
        mark_path(maximal_unitig.cycle_hash());
}


template <uint16_t k>
/**
 * @brief 提取最大单元图
 *
 * 从给定的Kmer中提取最大单元图，并存储到指定的Maximal_Unitig_Scratch对象中。
 * return false的情况:
 * 1. 当前顶点两端side是否分支已经确定,但是当前到达的side却不是分支
 * 2. 扩展到达的顶点两端side是否分支已经确定,但是当前到达的side却不是分支
 * 3. maximal unitigs 的最小顶点在标记顶点时被标记为已输出
 * @param v_hat Kmer对象，表示要提取的Kmer 读取的kmer
 * @param maximal_unitig Maximal_Unitig_Scratch对象，用于存储提取的最大单元图
 *
 * @return 提取成功返回true，否则返回false
 */
inline bool Read_CdBG_Extractor<k>::extract_maximal_unitig(const Kmer<k>& v_hat, Maximal_Unitig_Scratch<k>& maximal_unitig)
{
    static constexpr cuttlefish::side_t back = cuttlefish::side_t::back;
    static constexpr cuttlefish::side_t front = cuttlefish::side_t::front;

    //返回 v_hat的state_
    State_Read_Space state = hash_table[v_hat].state(); // State of the vertex `v_hat`.
    // 包含maximal unitig的元素已经输出。
    // 就是表示该点的两个side的 state都已经确定了。
    if(state.is_outputted())    // The containing maximal unitig has already been outputted.
        return false;

    // 让cycle = nullptr
    maximal_unitig.mark_linear();
    // 获取v_hat的状态state,然后从back开始walk,使用 unitig_back作为容器
    if(!walk_unitig(v_hat, state, back, maximal_unitig.unitig(back)))
        return false;//为false的情况,到达的顶点两端side是否分支已经确定,但是当前到达的side却不是分支
    // unitig() 返回对应side的 Unitig_Scratch<k>
    // 上面都标记了cycle = nullptr,为什么还需要这样判断环?
    if(maximal_unitig.unitig(back).is_cycle())
        maximal_unitig.mark_cycle(back);//是环就把对应side的unitig类赋值给cycle
    else
        if(!walk_unitig(v_hat, state, front, maximal_unitig.unitig(front)))
            return false;

    // sign_vertex() 返回两端端点中更小的顶点
    // 如果是环就返回所有列表里最小的顶点
    // 没有标记所有顶点,而是只标记标志顶点
    if(!mark_vertex(maximal_unitig.sign_vertex()))
        return false;
    // 没懂为什么要让back或者front其中之一反向互补
    maximal_unitig.finalize();
    return true;
}


template <uint16_t k>
/**
 * @brief 遍历单位图（unitig）
 *
 * 遍历给定单位图，根据指定的顶点、状态和边，更新单位图的信息。
 *
 * @param v_hat 当前顶点的 Kmer 对象
 * @param st_v 当前顶点的状态
 * @param s_v_hat 当前顶点的边侧信息
 * @param unitig 单位图对象
 *
 * @return 遍历是否成功
 */
inline bool Read_CdBG_Extractor<k>::walk_unitig(const Kmer<k>& v_hat, const State_Read_Space st_v, const cuttlefish::side_t s_v_hat, Unitig_Scratch<k>& unitig)
{
    // Data structures to be reused per each vertex extension of the unitig.
    // 需要被重用的数据结构
    // 当前顶点`v_hat`的一侧，通过它来扩展单元，即退出`v_hat`。
    // 如果是back开始扩展就是back，如果是front开始扩展就是front
    cuttlefish::side_t s_v = s_v_hat;   // The side of the current vertex `v_hat` through which to extend the unitig, i.e. exit `v_hat`.
    // 当前顶点被添加到unitigs。
    Directed_Vertex<k> v(s_v == cuttlefish::side_t::back ? v_hat : v_hat.reverse_complement(), hash_table); // Current vertex being added to the unitig.
    // 顶点`v`的状态。
    State_Read_Space state = st_v;  // State of the vertex `v`.
    // 潜在的下一个边缘从“v”包括到unitig。
    cuttlefish::edge_encoding_t e_v;    // The potential next edge from `v` to include into the unitig.
    // 与边' e_v '和从' v '出发的' s_v '相对应的碱基有可能添加到单位g的字面形式中。
    cuttlefish::base_t b_ext;   // The nucleobase corresponding to the edge `e_v` and the exiting side `s_v` from `v` to potentially add to the literal form of the unitig.
    // 主要是label_和 hash_装入v的相关信息
    unitig.init(v); // Initialize the unitig with the current vertex.


    while(true)
    {
        //根据边获取的下一个碱基
        e_v = state.edge_at(s_v);
        //如果是N or E,说明是非分支节点
        if(cuttlefish::is_fuzzy_edge(e_v))  // Reached an endpoint.
            break;
        //只有不是N or E的时候，才进行下一步
        b_ext = (s_v == cuttlefish::side_t::back ? DNA_Utility::map_base(e_v) : DNA_Utility::complement(DNA_Utility::map_base(e_v)));
        // std::cout << "v的起始标签是" << v.kmer().string_label() << std::endl;
        // v 属于重复使用的数据结构
        v.roll_forward(b_ext, hash_table);  // Walk to the next vertex.
        // std::cout << "v的移动后的标签是" << v.kmer().string_label() << std::endl;
        // 同理 state
        state = hash_table[v.hash()].state();
        // 同理 s_v
        // v的kmer是否为 canoical kmer,如果是就是back，如果不是就是front
        s_v = v.entrance_side();
        // 到达的顶点两边side都已经确定是否为非分支
        if(state.is_outputted())
            return state.was_branching_side(s_v);   // If `s_v` was a branching side, then the walk just crossed to a different unitig;
                                                    // so this unitig is depleted. Otherwise, `s_v` must belong to this unitig. In that
                                                    // case, the unitig has already been outputted earlier.
        // 到达的side状态为N,return True
        if(state.is_branching_side(s_v))    // Crossed an endpoint and reached a different unitig.
            break;

        // Still within the unitig.
        // 如果扩展到anchor,也就是起始点就break
        if(!unitig.extend(v, DNA_Utility::map_char(b_ext)))
            break;  // The unitig is a DCC (Detached Chordless Cycle).
        // 到达顶点的另一端side
        s_v = cuttlefish::opposite_side(s_v);
    }


    return true;
}



#endif
