#include "ConsolePrinter.h"
#include <iostream>

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


#ifdef _WIN32

COORD ConsolePrinter::getCursorPosition() {
	CONSOLE_SCREEN_BUFFER_INFO pBuffer;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &pBuffer);
	return COORD{ pBuffer.dwCursorPosition.X, pBuffer.dwCursorPosition.Y };
}

void ConsolePrinter::restoreCursorPosition() {
	cursorPosition = getCursorPosition();
}

void ConsolePrinter::MoveToLastCursorPosition() {
	moveCursorToPosition(cursorPosition);
}

void ConsolePrinter::EraseToLastCursorPosition() {
	clearDisplayBetweenTwoCursors(cursorPosition, getCursorPosition());
	MoveToLastCursorPosition();
}

COORD ConsolePrinter::getScrnInfo() {                   //获取控制台窗口缓冲区大小
	HANDLE hStd = GetStdHandle(STD_OUTPUT_HANDLE);      //获得标准输出设备句柄
	CONSOLE_SCREEN_BUFFER_INFO scBufInf;                //定义一个窗口缓冲区信息结构体
	GetConsoleScreenBufferInfo(hStd, &scBufInf);        //获取窗口缓冲区信息
	return scBufInf.dwSize;                             //返回窗口缓冲区大小
}

void ConsolePrinter::moveCursorToPosition(COORD pstn) {
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pstn);
}

void ConsolePrinter::clearDisplayBetweenTwoCursors(COORD p1, COORD p2){
	COORD currentCursor = getCursorPosition();

	int yDValue(p2.Y - p1.Y);
	int xMax = getScrnInfo().X; 
	moveCursorToPosition(p1);            
	for (int y(0); y <= yDValue; y++) {
		for (int x(p1.X); x <= xMax; x++)	{
			std::cout << ' ';
			int px = 0;
			if (x != xMax)
				px = x + 1;
			else
				px = 0;
			if (y == yDValue && px == p2.X)
				break;
		}
	}

	moveCursorToPosition(currentCursor);
}

void ConsolePrinter::drawProgressBar(double percent) {
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

#elif defined __linux__


COORD ConsolePrinter::getCursorPosition() { return COORD(); }
//获取控制台窗口缓冲区大小
COORD ConsolePrinter::getScrnInfo() { return COORD(); }
//设置控制台光标位置
void ConsolePrinter::setCursorPosition(COORD pstn) {}
//清除p1, p2之间的控制台内容，并将光标定位于p1
void ConsolePrinter::clearDisplay(COORD p1, COORD p2) {}
//绘制进度条。0<=percent<=100
void ConsolePrinter::drawProgressBar(int percent) {}


#endif
