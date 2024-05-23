
#include "Seq_Input.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>


const std::vector<std::string> Seq_Input::empty_collection;


Seq_Input::Seq_Input(   const std::vector<std::string>& seqs,
                        const std::vector<std::string>& lists,
                        const std::vector<std::string>& dirs):
    seq_paths_(seqs),
    list_paths_(lists),
    dir_paths_(dirs)
{}


Seq_Input::Seq_Input(   const std::optional<std::vector<std::string>>& seqs,
                        const std::optional<std::vector<std::string>>& lists,
                        const std::optional<std::vector<std::string>>& dirs):
    Seq_Input(seqs.value_or(empty_collection), lists.value_or(empty_collection), dirs.value_or(empty_collection))
{}


const std::vector<std::string>& Seq_Input::seq_paths() const
{
    return seq_paths_;
}


const std::vector<std::string>& Seq_Input::list_paths() const
{
    return list_paths_;
}


const std::vector<std::string>& Seq_Input::dir_paths() const
{
    return dir_paths_;
}


/**
 * @brief 获取序列
 *
 * 从给定的原始序列路径、序列列表和序列目录中收集序列，并返回包含所有序列的字符串向量。
 *
 * @return 包含所有序列的字符串向量
 */
const std::vector<std::string> Seq_Input::seqs() const
{
    std::vector<std::string> seqs;

    // Collect sequences from the raw sequence paths provided.
    //  从提供的原始序列路径中收集序列。
    seqs.insert(seqs.end(), seq_paths_.begin(), seq_paths_.end());


    // Collect sequences from the provided sequence lists.
    for(const std::string& list_path: list_paths_)
    {
        std::ifstream input(list_path.c_str(), std::ifstream::in);
        if(input.fail())
        {
            std::cerr << "Error opening list file " << list_path << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        std::string seq_path;
        while(input >> seq_path)
            seqs.emplace_back(seq_path);

        input.close();
    }


    // Collect sequences from the provided sequence directories.
    for(const std::string& dir_path: dir_paths_)
        for(const auto& entry: std::filesystem::directory_iterator(dir_path))
            seqs.emplace_back(entry.path());


    return seqs;
}


bool Seq_Input::empty() const
{
    return seq_paths_.empty() && list_paths_.empty() && dir_paths_.empty();
}
