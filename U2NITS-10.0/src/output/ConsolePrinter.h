
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H
#include <windows.h>
#include <string>

class ConsolePrinter {
private:
	COORD m_cursorPosition;

public:

	ConsolePrinter() { restoreCursorPosition(); }
	//�洢��ǰ����̨���λ�á�x�ᳯ�ң�y�ᳯ��
	void restoreCursorPosition();
	void MoveToLastCursorPosition();
	void EraseToLastCursorPosition();

public:
	static COORD getCursorPosition();
	static void setCursorPosition(COORD pstn);
	//���ڻ�������С
	static COORD getScreenBufferSize();
	//��ӡ��ʼ����
	static void printHeader();
	//���ƽ�������0<=percent<=100
	static void drawProgressBar(double percent);
	//���p1, p2֮��Ŀ���̨���ݣ���겻��
	static void clearDisplay(COORD p1, COORD p2);
	//������ҷ����ַ�
	static std::string printSolveInfo(double calTime, double calStep, double calSpeed, double nFile, double t, double T, const double* residual_vector);
};


#endif