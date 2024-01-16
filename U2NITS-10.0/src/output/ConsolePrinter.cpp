#include "ConsolePrinter.h"
#include <iostream>
#include "../global/StringProcessor.h"
#include <sstream>

void ConsolePrinter::printHeader() {
	//
	std::cout << R"(  _______________________________________________ )" << "\n";
	std::cout << R"( |  _   _   ____    _   _   ___   _____   ____   |)" << "\n";//生成工具：https://tools.kalvinbg.cn/txt/ascii
	std::cout << R"( | | | | | |___ \  | \ | | |_ _| |_   _| / ___|  |)" << "\n";
	std::cout << R"( | | | | |   __) | |  \| |  | |    | |   \___ \  |)" << "\n";
	std::cout << R"( | | |_| |  / __/  | |\  |  | |    | |    ___) | |)" << "\n";
	std::cout << R"( |  \___/  |_____| |_| \_| |___|   |_|   |____/  |)" << "\n";
	std::cout << R"( |_______________________________________________|)" << "\n";
	std::cout << "Finite Volume Method Solver (version 10.0), created by"
		<< " tgl\n";
	std::cout << "------------------------------------------------------------\n";


}

void ConsolePrinter::drawProgressBar(double percent) {
	std::cout << "Progress:";
	if (percent > 100)percent = 100;
	if (percent < 0)percent = 0;
	int nTotal = 45;//总长度，单位char
	const int nBlock = int(nTotal * percent / 100.0);//块个数
	for (int i = 0; i < nBlock; i++) {
		std::cout << "█";
	}
	for (int i = 0; i < nTotal - nBlock; i++) {
		std::cout << " ";
	}
	std::cout << "| " << percent << "%";
}



#ifdef _WIN32

COORD ConsolePrinter::getCursorPosition() {
	CONSOLE_SCREEN_BUFFER_INFO consoleScreenBufferInfo;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &consoleScreenBufferInfo);
	return COORD{ consoleScreenBufferInfo.dwCursorPosition.X, consoleScreenBufferInfo.dwCursorPosition.Y };
}

void ConsolePrinter::restoreCursorPosition() {
	m_cursorPosition = getCursorPosition();
}

void ConsolePrinter::MoveToLastCursorPosition() {
	setCursorPosition(m_cursorPosition);
}

void ConsolePrinter::EraseToLastCursorPosition() {
	clearDisplay(m_cursorPosition, getCursorPosition());
	MoveToLastCursorPosition();
}

COORD ConsolePrinter::getScreenBufferSize() {
	CONSOLE_SCREEN_BUFFER_INFO consoleScreenBufferInfo;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &consoleScreenBufferInfo);
	return consoleScreenBufferInfo.dwSize;//返回窗口缓冲区大小

	/*
	
	typedef struct _CONSOLE_SCREEN_BUFFER_INFO {
	  COORD      dwSize; the size of the console screen buffer, in character columns and rows
	  COORD      dwCursorPosition; the column and row coordinates of the cursor in the console screen buffer.
	  WORD       wAttributes;
	  SMALL_RECT srWindow; coordinates of the upper-left and lower-right corners of the display window
	  COORD      dwMaximumWindowSize;
	} CONSOLE_SCREEN_BUFFER_INFO;
	
	*/
}

void ConsolePrinter::setCursorPosition(COORD cursorPosition) {
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), cursorPosition);
}

void ConsolePrinter::clearDisplay(COORD p1, COORD p2){
	// 控制台界面函数 https://blog.51cto.com/dlican/3745102
	// (0,0)位于缓冲区左上角

	COORD currentCursor = getCursorPosition();
	int width = getScreenBufferSize().X;
	int x1 = p1.X;
	int y1 = p1.Y;
	int z1 = width * y1 + x1;
	int x2 = p2.X;
	int y2 = p2.Y;
	int z2 = width * y2 + x2;
	if (z2 > z1) {
		LPDWORD pWritten;//[output]接收实际写入缓冲区的字符数
		pWritten = new DWORD[z2 - z1];
		// 该函数效率更高，可避免闪烁(与std::cout相比)
		FillConsoleOutputCharacter(GetStdHandle(STD_OUTPUT_HANDLE), ' ', z2 - z1, p1, pWritten);
		delete []pWritten;
	}
	//for (int i = 0; i < z2 - z1; i++) {// 若z1<=z2，则不会输出
	//	std::cout << ' ';
	//}
	setCursorPosition(currentCursor);


}
std::string ConsolePrinter::printSolveInfo(double calTime, double calStep, double calSpeed, double nFile, double t, double T, const double* residual_vector) {
	std::stringstream info;
	info << "\n"
		<< "  Calculate time: \t" << StringProcessor::timeFormat((int)calTime) << "\t(" << calTime << " s)\n"
		<< "  Calculate step: \t" << calStep << "\n"
		<< "  Calculate speed: \t" << calSpeed << "\t step/s\n"
		<< "  Output file num: \t" << nFile << "\n"
		<< "  Residual rho: \t" << std::scientific << residual_vector[0] << "\n" << std::defaultfloat
		<< "  Physical time: \t" << t << " s\t/" << T << " s\n"
		<< "Press ESC to end Computation\n";
	std::string str = info.str();
	std::cout << str;
	return str;
}
#endif
