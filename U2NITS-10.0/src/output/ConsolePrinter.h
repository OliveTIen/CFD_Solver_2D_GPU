
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
	// �洢��ǰ����̨���λ�á�x�ᳯ�ң�y�ᳯ��
	void restoreCursorPosition();
	void restoreClearStartPosition();
	void setClearStartPosition();
	void update();
	void clear();
	void printWelcome_2();
	void MoveToLastCursorPosition();
	void EraseToLastCursorPosition();

	void print(std::string str) { std::cout << str; }

public:
	static COORD getCursorPosition();
	static void setCursorPosition(COORD pstn);
	// ���ڻ�������С
	static COORD getScreenBufferSize();
	// ��ӡ��ʼ����
	static void printHeader(HeaderStyle h);
	// ...
	static void printGenshinStart();
	// ���ƽ�������0<=percent<=100
	static void drawProgressBar(double percent);
	// ���p1, p2֮��Ŀ���̨���ݣ���겻��
	static void clearDisplay(COORD p1, COORD p2);
	// ��װ�����Ϣ
	void assemblySolveInfo(double calTime, int calStep, int maxIteration, double calSpeed, double nFile, double t, double T, const double* residual_vector);
	std::string getSolveInfo() { return m_solveInfo; }
	// ��ӡ�����Ϣ
	static void printInfo(InfoType type);

};


#endif