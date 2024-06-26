#include "version.hpp"

#include <algorithm>
#include <string>
#include <vector>
#include <iostream>
#include <optional>


#ifdef __cplusplus
extern "C" {
#endif
  int cf_build(int argc, char** argv);
  int cf_validate(int argc, char** argv);
#ifdef __cplusplus
}
#endif


std::string executable_version()
{
    return "cuttlefish " VERSION;
}


void display_help_message()
{
    std::cout << executable_version() << "\n";
    std::cout << "Supported commands: `build`, `help`, `version`.\n";
    
    std::cout << "Usage:\n";
    std::cout << "\tcuttlefish build [options]\n";
}


int main(int argc, char** argv)
{
#ifdef CF_DEVELOP_MODE
    std::cout << "Warning: Executing in Develop Mode.\n";
#endif

    if(argc < 2)
        display_help_message();
    else
    {
      std::string command(argv[1]);
      printf("%s\n", command.c_str());
      // std::transform 的作用是将 command 字符串中的所有字符转换为小写。
      // 这个函数使用了lambda表达式，lambda表达式可以捕获变量。
      std::transform(command.begin(), command.end(), command.begin(),
                     [](const char ch) { return std::tolower(ch); });

      if (command == "build")
        return cf_build(argc - 1, argv + 1);
      else if (command == "validate")
        return cf_validate(argc - 1, argv + 1);
      else if (command == "help")
        display_help_message();
      else if (command == "version")
        std::cout << executable_version() << "\n";
      else
        display_help_message();
    }

    return EXIT_SUCCESS;
}
