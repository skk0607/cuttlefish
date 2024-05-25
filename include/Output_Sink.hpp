
#ifndef OUTPUT_SINK_HPP
#define OUTPUT_SINK_HPP



#include "Async_Logger_Wrapper.hpp"
#include "spdlog/spdlog.h"

#include <fstream>
#include <string>


// A basic sink wrapper with minimal functionality — open, get reference to the wrapped sink, and close.
template <typename T_sink_>
class Output_Sink
{};


template <>
class Output_Sink<std::ofstream>
{
private:

    std::ofstream output_;


public:

    void init_sink(const std::string& output_file_path)
    {
        output_ = std::ofstream(output_file_path);
    }

    std::ofstream& sink()
    {
        return output_;
    }

    void close_sink()
    {
        output_.close();
    }
};

/**
 * @brief 	  这表示我们正在为一个模板类提供一个特化版本。< >中的内容是空的，因为我们正在为特定的类型（即Async_Logger_Wrapper）提供特化，而不是为特定的模板参数值。
 */
template <>
class Output_Sink<Async_Logger_Wrapper>
{
private:

    Async_Logger_Wrapper output_;


public:

    void init_sink(const std::string& output_file_path)
    {
        output_.init_logger(output_file_path);
    }

    /**
     * @brief 获取异步日志包装器的输出目标
     *
     * 返回异步日志包装器的输出目标对象引用，以便进一步操作。
     *
     * @return 异步日志包装器的输出目标对象引用
     */
    Async_Logger_Wrapper& sink()
    {
        return output_;
    }

    /**
     * @brief 关闭日志输出
     *
     * 关闭日志输出流，停止日志的写入操作。
     */
    void close_sink()
    {
        output_.close_logger();
    }
};



#endif
