#include <windows.h>
#ifndef CONSOLE_OUTPUT_H
#define CONSOLE_OUTPUT_H


class ConsolePrinter {
private:
	COORD cursorPosition;

public:
	static void printHeader();
	//���ƽ�������0<=percent<=100
	static void drawProgressBar(double percent);

	ConsolePrinter() { restoreCursorPosition(); }
	//�洢��ǰ����̨���λ�á�x�ᳯ�ң�y�ᳯ��
	void restoreCursorPosition();
	void MoveToLastCursorPosition();
	void EraseToLastCursorPosition();

private:
	//���p1, p2֮��Ŀ���̨���ݣ�������궨λ��p1
	void clearDisplayBetweenTwoCursors(COORD p1, COORD p2);
	//���ÿ���̨���λ��
	void moveCursorToPosition(COORD pstn);
	//��ȡ����̨���ڻ�������С
	COORD getScrnInfo();
	//��ȡ���
	COORD getCursorPosition();

};


#endif