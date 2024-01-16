
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H
#include <windows.h>
#include <string>

class ConsolePrinter {
private:
	COORD m_cursorPosition;

public:

	ConsolePrinter() { restoreCursorPosition(); }
	//存储当前控制台光标位置。x轴朝右，y轴朝下
	void restoreCursorPosition();
	void MoveToLastCursorPosition();
	void EraseToLastCursorPosition();

public:
	static COORD getCursorPosition();
	static void setCursorPosition(COORD pstn);
	//窗口缓冲区大小
	static COORD getScreenBufferSize();
	//打印初始界面
	static void printHeader();
	//绘制进度条。0<=percent<=100
	static void drawProgressBar(double percent);
	//清除p1, p2之间的控制台内容，光标不变
	static void clearDisplay(COORD p1, COORD p2);
	//输出，且返回字符
	static std::string printSolveInfo(double calTime, double calStep, double calSpeed, double nFile, double t, double T, const double* residual_vector);
};


#endif