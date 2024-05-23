
#ifndef PROGRESS_TRACKER_HPP
#define PROGRESS_TRACKER_HPP



#include "Spin_Lock.hpp"

#include <cstdint>
#include <string>
#include <cmath>
#include <iostream>


// A basic class to track and display progress for some work.
class Progress_Tracker
{
private:

    uint64_t total_work_load;  //   Total amount of work to be done over time.
    uint64_t work_chunk_threshold; // Granularity of the provided work chunk sizes that triggers tracking updates.

    uint64_t total_work_done;   // Amount of work done until now.
    uint16_t percent_work_done; // Percentage of the completed workload.
    std::string log_message;  // Message to display at the logs.

    Spin_Lock lock;  // Lock to ensure multiple threads can access the tracker safely.

public:
  // Sets up the tracker for some task with total size `total_work_load`;
  // updates to the tracking are to be triggered when some work-chunk of size at
  // least `work_chunk_threshold` is provided to it. The log message to be
  // displayed over the course of tracking is `log_message`.
  // 为某些任务设置跟踪器，总大小为`total_work_load`;当向跟踪提供了大小至少为`work_chunk_threshold`的工作块时，将触发对跟踪的更新。在跟踪过程中显示的日志消息是`log_message`。
  void setup(uint64_t total_work_load, uint64_t work_chunk_threshold,
             const std::string &log_message);

  // Tracks progress made for a work-chunk of size `work_chunk_size`. If the
  // passed chunk-size is large enough, then it is tracked and `true` is
  // returned. All lesser sized chunk update requests are ignored and `false` is
  // returned. So, repeated invocation is suggested.
  bool track_work(uint64_t work_chunk_size);
};


/**
 * @brief 跟踪工作进度
 *
 * 根据给定的工作块大小，跟踪并更新工作进度。如果工作块大小达到或超过设定的阈值，则更新总完成工作量，
 * 计算新的工作百分比，并打印进度信息。
 *
 * @param work_chunk_size 工作块大小
 *
 * @return 如果工作块大小达到或超过阈值，则返回 true；否则返回 false
 */
inline bool Progress_Tracker::track_work(const uint64_t work_chunk_size)
{
    if(work_chunk_size >= work_chunk_threshold)
    {
        lock.lock();
        
        total_work_done += work_chunk_size;
        // 计算已完成工作量占总工作量的百分比。首先，将 total_work_done 乘以
        // 100.0，然后除以 total_work_load（总工作量）。接着，使用 std::round
        // 对结果进行四舍五入，并将结果转换为 uint16_t 类型。
        const uint16_t new_percent = static_cast<uint16_t>(std::round((total_work_done * 100.0) / total_work_load));
        if(percent_work_done < new_percent)
        {
          // 如果新的百分比大于 当前的进度百分比 percent_work_done，则更新
          // percent_work_done
          // 并打印新的进度信息到标准错误流（通常用于输出错误信息或诊断信息）。
          percent_work_done = new_percent;
          std::cerr << "\r[" << log_message << "]\t" << percent_work_done
                    << "%";
        }

        lock.unlock();

        return true;
    }

    return false;
}



#endif
