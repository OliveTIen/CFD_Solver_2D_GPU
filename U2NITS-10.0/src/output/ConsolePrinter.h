
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H
#include <windows.h>
#include <string>
#include <iostream>

class ConsolePrinter {
private:
	COORD m_cursorPosition;
	std::string m_solveInfo;

public:
	enum InfoType {
		type_nan_detected
	};
	enum HeaderStyle {
		normal,
		simple
	};
	enum InfoStyle {};
	

public:
	ConsolePrinter() { restoreCursorPosition(); }
	// 存储当前控制台光标位置。x轴朝右，y轴朝下
	void restoreCursorPosition();
	void restoreClearStartPosition();
	void setClearStartPosition();
	void update();
	void clear();
	void printWelcomeInterface();
	void MoveToLastCursorPosition();
	void EraseToLastCursorPosition();

	void print(std::string str) { std::cout << str; }

public:
	static COORD getCursorPosition();
	static void setCursorPosition(COORD pstn);
	// 窗口缓冲区大小
	static COORD getScreenBufferSize();
	// 打印初始界面
	static void printHeader(HeaderStyle h);
	// ...
	static void printGenshinStart();
	static void drawProgressBar(double value, double maxValue, int barLength = 45);
	// 清除p1, p2之间的控制台内容，光标不变
	static void clearDisplay(COORD p1, COORD p2);
	// 组装求解信息
	void assemblySolveInfo(double calTime, int calStep, int maxIteration, double calSpeed, double nFile, double t, double T, const double* residual_vector);
	// return solveInfo
	std::string setSolveInfo(int startStep, int currentStep, int endStep, int numOfFile, double usedTime, double physicalTime, double maxPhysicalTime, const double* residual_vector, double CFL);
	std::string getSolveInfo() { return m_solveInfo; }
	// 打印求解信息
	static void printInfo(InfoType type);

private:
	// 绘制进度条。0<=percent<=100
	static void m_drawProgressBar(double percent, int barLength);


};


#endif