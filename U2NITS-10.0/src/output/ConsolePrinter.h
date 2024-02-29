
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H
#include <windows.h>
#include <string>

class ConsolePrinter {
private:
	COORD m_cursorPosition;

public:
	enum InfoType {
		type_nan_detected
	};

public:
	ConsolePrinter() { restoreCursorPosition(); }
	// 存储当前控制台光标位置。x轴朝右，y轴朝下
	void restoreCursorPosition();
	void MoveToLastCursorPosition();
	void EraseToLastCursorPosition();

public:
	static COORD getCursorPosition();
	static void setCursorPosition(COORD pstn);
	// 窗口缓冲区大小
	static COORD getScreenBufferSize();
	// 打印初始界面
	static void printHeader();
	// 原神，启动！
	static void printGenshinStart();
	// 绘制进度条。0<=percent<=100
	static void drawProgressBar(double percent);
	// 清除p1, p2之间的控制台内容，光标不变
	static void clearDisplay(COORD p1, COORD p2);
	// 组装求解信息
	static std::string assemblySolveInfo(double calTime, int calStep, int maxIteration, double calSpeed, double nFile, double t, double T, const double* residual_vector);
	// 打印求解信息
	static void printInfo(InfoType type);
};


#endif