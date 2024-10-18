#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include "ConsoleProxy.h"

using namespace GUI;

void ConsoleProxy::showCursor(bool show) {
	CONSOLE_CURSOR_INFO cci;
	cci.bVisible = show;
	cci.dwSize = 1;
	SetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &cci);
}

COORD ConsoleProxy::getCursorPosition() {
	CONSOLE_SCREEN_BUFFER_INFO consoleScreenBufferInfo;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &consoleScreenBufferInfo);
	return COORD{ consoleScreenBufferInfo.dwCursorPosition.X, consoleScreenBufferInfo.dwCursorPosition.Y };
}

void ConsoleProxy::setCursorPosition(COORD p) {
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), p);
}

void ConsoleProxy::clearDisplay(COORD p1, COORD p2) {
	// reference: [Console coding 控制台窗口图形界面编程](https://blog.51cto.com/dlican/3745102)
	COORD p_old = getCursorPosition();
	CONSOLE_SCREEN_BUFFER_INFO consoleScreenBufferInfo;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &consoleScreenBufferInfo);
	int width = consoleScreenBufferInfo.dwSize.X;
	int num_write = (width * p2.Y + p2.X) - (width * p1.Y + p1.X);
	if (num_write > 0) {
		LPDWORD pStr = new DWORD[num_write];
		// This function can avoid screen flicker (compared to std::cout)
		FillConsoleOutputCharacter(GetStdHandle(STD_OUTPUT_HANDLE), ' ', num_write, p1, pStr);
		delete[] pStr;
	}
	setCursorPosition(p_old);
}