#include <windows.h>
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H


class ConsolePrinter {
private:
	COORD cursorPosition;

public:
	static void printHeader();
	//绘制进度条。0<=percent<=100
	static void drawProgressBar(double percent);

	ConsolePrinter() { restoreCursorPosition(); }
	//存储当前控制台光标位置。x轴朝右，y轴朝下
	void restoreCursorPosition();
	void MoveToLastCursorPosition();
	void EraseToLastCursorPosition();

private:
	//清除p1, p2之间的控制台内容，并将光标定位于p1
	void clearDisplayBetweenTwoCursors(COORD p1, COORD p2);
	//设置控制台光标位置
	void moveCursorToPosition(COORD pstn);
	//获取控制台窗口缓冲区大小
	COORD getScrnInfo();
	//获取光标
	COORD getCursorPosition();

};


#endif