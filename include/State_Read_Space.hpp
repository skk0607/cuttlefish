
#ifndef STATE_READ_SPACE_HPP
#define STATE_READ_SPACE_HPP



#include "globals.hpp"
#include "DNA.hpp"

#include <cstdint>


template <uint8_t BITS_PER_KEY> class Kmer_Hash_Entry_API;

// 类在read de Bruijn图中自动机的状态空间中表示一个状态。
// Class for a state in the state-space of the automata in read de Bruijn graphs.
class State_Read_Space
{
    friend class Kmer_Hash_Entry_API<cuttlefish::BITS_PER_READ_KMER>;

private:

    cuttlefish::state_code_t code;  // Numeric code of the state.

    static constexpr uint8_t BITS_PER_SIDE = 3; // Number of bits required to `Extended_Base`-encode edges incident to a side.
    static constexpr uint8_t FRONT_IDX = BITS_PER_SIDE; // Starting index of the three bits encoding the front-incident edge.
    static constexpr uint8_t BACK_IDX = 0;  // Starting index of the three bits encoding the back-incident edge.

    // Bitmask to extract the edge-encoding from some side. Has to be shifted to appropriate index before extraction.
    static constexpr uint8_t SIDE_MASK = (1 << BITS_PER_SIDE) - 1;

    // Bitmask used to extract the 'Extended_Base`-encoding of the edge(s) incident to the front side of a vertex.
    static constexpr cuttlefish::state_code_t FRONT_MASK = SIDE_MASK << FRONT_IDX;

    // Bitmask used to extract the 'Extended_Base`-encoding of the edge(s) incident to the back side of a vertex.
    static constexpr cuttlefish::state_code_t BACK_MASK = SIDE_MASK << BACK_IDX;


    // Constructs a state that wraps the provided numeric value `code`.
    State_Read_Space(cuttlefish::state_code_t code);

    // Sets the back-encoding of the state to the `Extended_Base`-encoding `edge`.
    void set_back_encoding(cuttlefish::edge_encoding_t edge);

    // Sets the front-encoding of the state to the `Extended_Base`-encoding `edge`.
    void set_front_encoding(cuttlefish::edge_encoding_t edge);


public:

    // Constructs the state of a vertex having both its sides unvisited.
    constexpr State_Read_Space();

    // Returns the wrapped state-code value.
    cuttlefish::state_code_t get_state() const;

    // Returns `true` iff some vertex having this state has been outputted.
    bool is_outputted() const;

    // Returns the `Extended_Base`-encoding of the edge(s) incident to the side
    // `side` of a vertex having this state.
    // 返回具有此状态的顶点在给定侧边的Extended_Base编码
    cuttlefish::edge_encoding_t edge_at(cuttlefish::side_t side) const;

    // Returns `true` iff some vertex having this state is branching (i.e. has
    // multiple incident edges) at its side `side`, and hasn't been outputted yet.
    bool is_branching_side(cuttlefish::side_t side) const;

    // Returns `true` iff some vertex having this state is branching (i.e. has
    // multiple incident edges) at its side `side`, and has already been outputted.
    bool was_branching_side(cuttlefish::side_t side) const;

    // Updates the `Extended_Base` encoding of the side `side` of this state, with
    // `edge`.
    void update_edge_at(cuttlefish::side_t side, cuttlefish::edge_encoding_t edge);

    // Marks the state as already been outputted.
    void mark_outputted();

    // Returns `true` iff the underlying code is the same as that one of `rhs`.
    bool operator==(const State_Read_Space& rhs) const;

    // Returns the state for the vertices that have been marked as outputted.
    static const State_Read_Space& get_outputted_state();

    // For the given code `code` of some state `s`, returns the code of the
    // state `s_op` which is the corresponding state where the vertices having
    // the DFA state `s` in the underlying graph transition to when outputted. 
    static cuttlefish::state_code_t mark_outputted(cuttlefish::state_code_t code);
};


/**
 * @brief 构造State_Read_Space对象
 *
 * 构造一个State_Read_Space对象，并通过位运算设置code成员变量的值。
 * 将cuttlefish::edge_encoding_t::E类型的值左移FRONT_IDX位，并与另一个cuttlefish::edge_encoding_t::E类型的值进行按位或运算，
 * 将结果赋值给code成员变量。
 */
inline constexpr State_Read_Space::State_Read_Space():
    code{(static_cast<cuttlefish::state_code_t>(cuttlefish::edge_encoding_t::E) << FRONT_IDX) | static_cast<cuttlefish::state_code_t>(cuttlefish::edge_encoding_t::E)}
{}


inline State_Read_Space::State_Read_Space(const cuttlefish::state_code_t code):
    code(code)
{}


inline void State_Read_Space::set_back_encoding(cuttlefish::edge_encoding_t edge)
{
    code = (code & FRONT_MASK) | (static_cast<cuttlefish::state_code_t>(edge) << BACK_IDX);
}


/**
 * @brief 设置前沿编码
 *
 * 将给定的边缘编码设置为前沿编码，并更新状态码。
 * next_base << 3 | back_code
 * @param edge 边缘编码
 */
inline void State_Read_Space::set_front_encoding(cuttlefish::edge_encoding_t edge)
{
    code = (code & BACK_MASK) | (static_cast<cuttlefish::state_code_t>(edge) << FRONT_IDX);
}


/**
 * @brief 获取状态
 *
 * 返回当前状态机的状态码。
 *
 * @return 状态码
 */
inline cuttlefish::state_code_t State_Read_Space::get_state() const
{
    return code;
}


/**
 * @brief 获取边的编码
 *
 * 根据给定的边所在的方向，返回该边的编码。
 *
 * @param side 边所在的方向
 *
 * @return 边的编码
 */
inline cuttlefish::edge_encoding_t State_Read_Space::edge_at(const cuttlefish::side_t side) const
{
    return static_cast<cuttlefish::edge_encoding_t>(side == cuttlefish::side_t::front ? (code & FRONT_MASK) >> FRONT_IDX : (code & BACK_MASK) >> BACK_IDX);
}


/**
 * @brief 判断是否为分支侧
 *
 * 判断给定的侧是否为分支侧。
 *
 * @param side 侧
 *
 * @return 如果给定的侧为分支侧，则返回 true；否则返回 false
 */
inline bool State_Read_Space::is_branching_side(const cuttlefish::side_t side) const
{
    return edge_at(side) == cuttlefish::edge_encoding_t::N;
}


/**
 * @brief 判断是否为分支边
 *
 * 判断给定的边是否为分支边。
 *
 * @param side 边所在的侧
 *
 * @return 如果给定的边是分支边，则返回 true；否则返回 false
 */
inline bool State_Read_Space::was_branching_side(const cuttlefish::side_t side) const
{
    return edge_at(side) == cuttlefish::edge_encoding_t::OP_branching;
}


/**
 * @brief 更新边的编码信息
 *
 * 根据给定的边和边所在的侧面，更新 `State_Read_Space` 对象中对应边的编码信息。
 * side = front就设置front的编码,反之back的编码。
 * @param side 边的侧面
 * @param edge 边的编码
 */
inline void State_Read_Space::update_edge_at(const cuttlefish::side_t side, const cuttlefish::edge_encoding_t edge)
{
    side == cuttlefish::side_t::front ? set_front_encoding(edge) : set_back_encoding(edge);
}


/**
 * @brief 标记输出
 *
 * 将当前对象标记为已输出状态。根据前后两个方向的分支情况，设置对应的编码值。
 */
inline void State_Read_Space::mark_outputted()
{
    static constexpr cuttlefish::edge_encoding_t OP_non_branch = cuttlefish::edge_encoding_t::OP_non_branch;
    static constexpr cuttlefish::edge_encoding_t OP_branching = cuttlefish::edge_encoding_t::OP_branching;
    
    if(!is_outputted())
    {
        set_back_encoding(is_branching_side(cuttlefish::side_t::back) ? OP_branching : OP_non_branch);
        set_front_encoding(is_branching_side(cuttlefish::side_t::front) ? OP_branching : OP_non_branch);
    }
}


/**
 * @brief 判断是否已经输出
 *
 * 判断当前状态是否已经输出。
 *
 * @return 如果已经输出则返回 true，否则返回 false
 */
inline bool State_Read_Space::is_outputted() const
{
    static constexpr uint8_t OP_non_branch = static_cast<uint8_t>(cuttlefish::edge_encoding_t::OP_non_branch);
    static constexpr uint8_t OP_branching = static_cast<uint8_t>(cuttlefish::edge_encoding_t::OP_branching);

    return  code == ((OP_non_branch << FRONT_IDX)   |   (OP_non_branch << BACK_IDX))    ||
            code == ((OP_non_branch << FRONT_IDX)   |   (OP_branching << BACK_IDX))     ||
            code == ((OP_branching << FRONT_IDX)    |   (OP_non_branch << BACK_IDX))    ||
            code == ((OP_branching << FRONT_IDX)    |   (OP_branching << BACK_IDX)); 
}


/**
 * @brief 判断两个 State_Read_Space 对象是否相等
 *
 * 通过比较两个 State_Read_Space 对象的 code 成员变量是否相等，来判断它们是否相等。
 *
 * @param rhs 要比较的另一个 State_Read_Space 对象
 *
 * @return 如果两个对象的 code 成员变量相等，则返回 true；否则返回 false
 */
inline bool State_Read_Space::operator==(const State_Read_Space& rhs) const
{
    return code == rhs.code;
}


inline cuttlefish::state_code_t State_Read_Space::mark_outputted(const cuttlefish::state_code_t code)
{
    State_Read_Space state(code);
    state.mark_outputted();

    return state.get_state();
}



#endif
