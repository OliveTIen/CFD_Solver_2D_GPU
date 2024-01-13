#include "ConsolePrinter.h"
#include <iostream>

void ConsolePrinter::printHeader() {
	//
	std::cout << R"(  _______________________________________________ )" << "\n";
	std::cout << R"( |  _   _   ____    _   _   ___   _____   ____   |)" << "\n";//���ɹ��ߣ�https://tools.kalvinbg.cn/txt/ascii
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

COORD ConsolePrinter::getScrnInfo() {                   //��ȡ����̨���ڻ�������С
	HANDLE hStd = GetStdHandle(STD_OUTPUT_HANDLE);      //��ñ�׼����豸���
	CONSOLE_SCREEN_BUFFER_INFO scBufInf;                //����һ�����ڻ�������Ϣ�ṹ��
	GetConsoleScreenBufferInfo(hStd, &scBufInf);        //��ȡ���ڻ�������Ϣ
	return scBufInf.dwSize;                             //���ش��ڻ�������С
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
	int nTotal = 45;//�ܳ��ȣ���λchar
	const int nBlock = int(nTotal * percent / 100.0);//�����
	for (int i = 0; i < nBlock; i++) {
		std::cout << "��";
	}
	for (int i = 0; i < nTotal - nBlock; i++) {
		std::cout << " ";
	}
	std::cout << "| " << percent << "%";
}

#elif defined __linux__


COORD ConsolePrinter::getCursorPosition() { return COORD(); }
//��ȡ����̨���ڻ�������С
COORD ConsolePrinter::getScrnInfo() { return COORD(); }
//���ÿ���̨���λ��
void ConsolePrinter::setCursorPosition(COORD pstn) {}
//���p1, p2֮��Ŀ���̨���ݣ�������궨λ��p1
void ConsolePrinter::clearDisplay(COORD p1, COORD p2) {}
//���ƽ�������0<=percent<=100
void ConsolePrinter::drawProgressBar(int percent) {}


#endif
