
#ifndef APPLICATION_HPP
#define APPLICATION_HPP



#include "Validator.hpp"

#include <cstdint>


class Build_Params;
class Validation_Params;


// The top-level application class for the compaction algorithm.
/*
定义了一个模板类Application，它接受两个模板参数：一个无符号16位整数k和一个模板模板参数T_App。这个类主要用于表示一个应用程序的实例，并且这个实例在某种层次结构（top-down
hierarchy）中与其它应用程序实例相关联。

typename:表示T_App是一个类型参数,与之相对比就是
k,不是类型参数,声明时直接使用数字即可.
T_APP是 需要 template <uint16_t k> class kmer_Enumeration_Stats; 声明的参数
*/
template <uint16_t k, template<uint16_t> typename T_App>
class Application
{
private:

    // Pointer to an application instance of the next `Application` class in the top-down hierarchy (on `k`).
    Application<k - 2, T_App>* const app_next_level;

    // Pointer to a driver object that operates with the k-value `k`.
    T_App<k>* const app;

    // TODO: Make the validator member generic, like `T_App`.
    // Pointer to a `Validator` object that operates with the k-value `k`.
    Validator<k>* const validator;


public:

    // Constructs an `Application` instance with the provided build-parameters,
    // if the provided `k` parameter matches to the specialized template argument `k`.
    Application(const Build_Params& params);

    // Constructs an `Application` instance with the provided validation-parameters,
    // if the provided `k` parameter matches to the specialized template argument `k`.
    Application(const Validation_Params& params);

    ~Application();

    // Executes the compaction algorithm.
    // 本质是在cpp文件里调用app->execute()
    void execute() const;

    // Validates the result of the compaction algorithm.
    bool validate() const;
};


template <template<uint16_t> typename T_App>
class Application<1, T_App>
{
private:

    T_App<1>* const app;

    Validator<1>* const validator;


public:

    Application(const Build_Params& params);

    Application(const Validation_Params& params);

    ~Application();

    void execute() const;

    bool validate() const;
};



#endif
