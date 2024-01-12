#include "Painter.h"

void CMD::drawProgressBar(int progress, int total, int barWidth) {
    float ratio = static_cast<float>(progress) / total;
    int completedWidth = static_cast<int>(ratio * barWidth);
    std::cout << "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < completedWidth)
            std::cout << "=";
        else if (i == completedWidth)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << static_cast<int>(ratio * 100.0) << "%\r";//using std::cout<<"\r" + std::cout.flush() to rewrite a line
    std::cout.flush();
}

void CMD::printHeader() {
    std::cout << "\033[33;1m";//明黄色
    std::cout << R"(  _______________________________________________ )" << "\n";
    std::cout << R"( |  _   _   ____    _   _   ___   _____   ____   |)" << "\n";//生成工具：https://tools.kalvinbg.cn/txt/ascii
    std::cout << R"( | | | | | |___ \  | \ | | |_ _| |_   _| / ___|  |)" << "\n";
    std::cout << R"( | | | | |   __) | |  \| |  | |    | |   \___ \  |)" << "\n";
    std::cout << R"( | | |_| |  / __/  | |\  |  | |    | |    ___) | |)" << "\n";
    std::cout << R"( |  \___/  |_____| |_| \_| |___|   |_|   |____/  |)" << "\n";
    std::cout << R"( |_______________________________________________|)" << "\n";
    std::cout << "\033[36m";//天蓝
    std::cout << "Finite Volume Method Solver (version 1.0), created by"
        << "\033[32m"//绿
        << " tgl\n"
        << "\033[36m"//天蓝
        << "2023-10-31\n";
    std::cout << "\033[35m";//暗紫
    std::cout << "--------------------------------------------------------\n";
    std::cout << "\033[0m";//RESET

}
