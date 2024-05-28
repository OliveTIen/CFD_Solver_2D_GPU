#include "ConsolePrinter.h"
#include <iostream>
#include "../global/StringProcessor.h"
#include <sstream>
#include "LogWriter.h"
#include <iomanip>

void ConsolePrinter::printHeader(HeaderStyle h) {
	switch (h) {
	case HeaderStyle::normal:

		std::cout << R"(  _______________________________________________ )" << "\n";
		std::cout << R"( |  _   _   ____    _   _   ___   _____   ____   |)" << "\n";//生成工具：https://tools.kalvinbg.cn/txt/ascii
		std::cout << R"( | | | | | |___ \  | \ | | |_ _| |_   _| / ___|  |)" << "\n";
		std::cout << R"( | | | | |   __) | |  \| |  | |    | |   \___ \  |)" << "\n";
		std::cout << R"( | | |_| |  / __/  | |\  |  | |    | |    ___) | |)" << "\n";
		std::cout << R"( |  \___/  |_____| |_| \_| |___|   |_|   |____/  |)" << "\n";
		std::cout << R"( |_______________________________________________|)" << "\n";
		std::cout << "Finite Volume Method Solver (version 10.0), created by"
			<< " LASD\n";
		std::cout << "------------------------------------------------------------\n";

		break;
	case HeaderStyle::simple:
		std::cout << "U2NITS, by LASD\n";
		std::cout << "------------------------------------------------------------\n";
		break;

	}




}

void ConsolePrinter::printGenshinStart() {
	// 原神，启动
	if (false) {
	// doh 生成工具：https://tools.kalvinbg.cn/txt/ascii
	std::cout << R"(	    GGGGGGGGGGGGG                                                      hhhhhhh               iiii                    )" << "\n";
	std::cout << R"(     GGG::::::::::::G                                                      h:::::h              i::::i                   )" << "\n";
	std::cout << R"(   GG:::::::::::::::G                                                      h:::::h               iiii                    )" << "\n";
	std::cout << R"(  G:::::GGGGGGGG::::G                                                      h:::::h                                       )" << "\n";
	std::cout << R"( G:::::G       GGGGGG    eeeeeeeeeeee    nnnn  nnnnnnnn        ssssssssss   h::::h hhhhh       iiiiiiinnnn  nnnnnnnn     )" << "\n";
	std::cout << R"(G:::::G                ee::::::::::::ee  n:::nn::::::::nn    ss::::::::::s  h::::hh:::::hhh    i:::::in:::nn::::::::nn   )" << "\n";
	std::cout << R"(G:::::G               e::::::eeeee:::::een::::::::::::::nn ss:::::::::::::s h::::::::::::::hh   i::::in::::::::::::::nn  )" << "\n";
	std::cout << R"(G:::::G    GGGGGGGGGGe::::::e     e:::::enn:::::::::::::::ns::::::ssss:::::sh:::::::hhh::::::h  i::::inn:::::::::::::::n )" << "\n";
	std::cout << R"(G:::::G    G::::::::Ge:::::::eeeee::::::e  n:::::nnnn:::::n s:::::s  ssssss h::::::h   h::::::h i::::i  n:::::nnnn:::::n )" << "\n";
	std::cout << R"(G:::::G    GGGGG::::Ge:::::::::::::::::e   n::::n    n::::n   s::::::s      h:::::h     h:::::h i::::i  n::::n    n::::n )" << "\n";
	std::cout << R"(G:::::G        G::::Ge::::::eeeeeeeeeee    n::::n    n::::n      s::::::s   h:::::h     h:::::h i::::i  n::::n    n::::n )" << "\n";
	std::cout << R"( G:::::G       G::::Ge:::::::e             n::::n    n::::nssssss   s:::::s h:::::h     h:::::h i::::i  n::::n    n::::n )" << "\n";
	std::cout << R"(  G:::::GGGGGGGG::::Ge::::::::e            n::::n    n::::ns:::::ssss::::::sh:::::h     h:::::hi::::::i n::::n    n::::n )" << "\n";
	std::cout << R"(   GG:::::::::::::::G e::::::::eeeeeeee    n::::n    n::::ns::::::::::::::s h:::::h     h:::::hi::::::i n::::n    n::::n )" << "\n";
	std::cout << R"(     GGG::::::GGG:::G  ee:::::::::::::e    n::::n    n::::n s:::::::::::ss  h:::::h     h:::::hi::::::i n::::n    n::::n )" << "\n";
	std::cout << R"(        GGGGGG   GGGG    eeeeeeeeeeeeee    nnnnnn    nnnnnn  sssssssssss    hhhhhhh     hhhhhhhiiiiiiii nnnnnn    nnnnnn )" << "\n";
	}

	// big https://tool.cccyun.cc/ascii_art
	std::cout << R"(  _____                _     _         _____                            _   )" << "\n";
	std::cout << R"( / ____|              | |   (_)       |_   _|                          | |  )" << "\n";
	std::cout << R"(| |  __  ___ _ __  ___| |__  _ _ __     | |  _ __ ___  _ __   __ _  ___| |_ )" << "\n";
	std::cout << R"(| | |_ |/ _ \ '_ \/ __| '_ \| | '_ \    | | | '_ ` _ \| '_ \ / _` |/ __| __|)" << "\n";
	std::cout << R"(| |__| |  __/ | | \__ \ | | | | | | |  _| |_| | | | | | |_) | (_| | (__| |_ )" << "\n";
	std::cout << R"( \_____|\___|_| |_|___/_| |_|_|_| |_| |_____|_| |_| |_| .__/ \__,_|\___|\__|)" << "\n";
	std::cout << R"(                                                      | |                   )" << "\n";
	std::cout << R"(                                                      |_|                   )" << "\n";


}

void ConsolePrinter::m_drawProgressBar(myfloat percent, int barLength) {
	std::cout << "Progress:";
	if (percent > 1)percent = 1;
	if (percent < 0)percent = 0;
	const int nBlock = int(barLength * percent);//块个数
	for (int i = 0; i < nBlock; i++) {
		std::cout << "█";
	}
	for (int i = 0; i < barLength - nBlock; i++) {
		std::cout << " ";
	}
	// https://en.cppreference.com/w/cpp/io/manip/setprecision
	const auto default_precision{ std::cout.precision() };
	std::cout << std::setprecision(2);
	std::cout << std::fixed;
	std::cout << "| " << percent * 100.0 << "%";
	std::cout << std::defaultfloat;
	std::cout << std::setprecision(default_precision);
}

void ConsolePrinter::drawProgressBar(myfloat value, myfloat maxValue, int barLength) {
	m_drawProgressBar(value / maxValue, barLength);
}



#ifdef _WIN32

COORD ConsolePrinter::getCursorPosition() {
	CONSOLE_SCREEN_BUFFER_INFO consoleScreenBufferInfo;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &consoleScreenBufferInfo);
	return COORD{ consoleScreenBufferInfo.dwCursorPosition.X, consoleScreenBufferInfo.dwCursorPosition.Y };
}

void ConsolePrinter::restoreCursorPosition() {
	m_cursorPosition = getCursorPosition();
}

void ConsolePrinter::restoreClearStartPosition() {
	restoreCursorPosition();
}

void ConsolePrinter::setClearStartPosition() {
	restoreCursorPosition();
}

void ConsolePrinter::update() {
	auto& p1 = m_cursorPosition;
	ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
	ConsolePrinter::setCursorPosition(p1);
}

void ConsolePrinter::clear() {
	update();
}

void ConsolePrinter::printWelcomeInterface() {
	LogWriter::logAndPrint("Using VS Profiler. Please remove '/Profile' option on release.\n", LogWriter::Warning);
	ConsolePrinter::printHeader(ConsolePrinter::HeaderStyle::simple);
}

void ConsolePrinter::MoveToLastCursorPosition() {
	setCursorPosition(m_cursorPosition);
}

void ConsolePrinter::EraseToLastCursorPosition() {
	clearDisplay(m_cursorPosition, getCursorPosition());
	MoveToLastCursorPosition();
}

COORD ConsolePrinter::getScreenBufferSize() {
	CONSOLE_SCREEN_BUFFER_INFO consoleScreenBufferInfo;
	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &consoleScreenBufferInfo);
	return consoleScreenBufferInfo.dwSize;//返回窗口缓冲区大小

	/*
	
	typedef struct _CONSOLE_SCREEN_BUFFER_INFO {
	  COORD      dwSize; the size of the console screen buffer, in character columns and rows
	  COORD      dwCursorPosition; the column and row coordinates of the cursor in the console screen buffer.
	  WORD       wAttributes;
	  SMALL_RECT srWindow; coordinates of the upper-left and lower-right corners of the display window
	  COORD      dwMaximumWindowSize;
	} CONSOLE_SCREEN_BUFFER_INFO;
	
	*/
}

void ConsolePrinter::setCursorPosition(COORD cursorPosition) {
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), cursorPosition);
}

void ConsolePrinter::clearDisplay(COORD p1, COORD p2){
	// 控制台界面函数 https://blog.51cto.com/dlican/3745102
	// (0,0)位于缓冲区左上角

	COORD currentCursor = getCursorPosition();
	int width = getScreenBufferSize().X;
	int x1 = p1.X;
	int y1 = p1.Y;
	int z1 = width * y1 + x1;
	int x2 = p2.X;
	int y2 = p2.Y;
	int z2 = width * y2 + x2;
	if (z2 > z1) {
		LPDWORD pWritten;//[output]接收实际写入缓冲区的字符数
		pWritten = new DWORD[z2 - z1];
		// 该函数效率更高，可避免闪烁(与std::cout相比)
		FillConsoleOutputCharacter(GetStdHandle(STD_OUTPUT_HANDLE), ' ', z2 - z1, p1, pWritten);
		delete []pWritten;
	}
	//for (int i = 0; i < z2 - z1; i++) {// 若z1<=z2，则不会输出
	//	std::cout << ' ';
	//}
	setCursorPosition(currentCursor);


}

std::string ConsolePrinter::setSolveInfo(int startStep, int currentStep, int endStep, int numOfFile, myfloat usedTime, myfloat physicalTime, myfloat maxPhysicalTime, const myfloat* residual_vector, myfloat CFL) {
	myfloat speed = getSpeed(currentStep, startStep, usedTime);
	myfloat remainTime = (endStep - currentStep) / speed;

	//std::stringstream info_log;
	//info_log << "speed: " << speed << " step/s\n";
	//LogWriter::log(info_log.str());

	std::stringstream info;
	info << "\n"
		<< "  Time used: \t" << StringProcessor::timeFormat((int)usedTime) << "\t(Remain: " << StringProcessor::timeFormat((int)remainTime) << " s)\n"
		<< "  Step: \t" << currentStep << " \t/" << endStep << "\n"
		<< "  Speed: \t" << speed << "\t step/s\n"
		<< "  Output file num: \t" << numOfFile << "\n"
		<< "  res_rho: \t" << std::scientific << residual_vector[0] << std::defaultfloat << "\tCFL: " << CFL << "\n"
		<< "  Physical time: \t" << physicalTime << " s\t/" << maxPhysicalTime << " s\n"
		<< "Press ESC to end Computation\n";
	m_solveInfo = info.str();
	return m_solveInfo;
}

void ConsolePrinter::printInfo(InfoType type) {
	switch (type) {
	case InfoType::type_nan_detected:
	{
		std::cout << "\nWarning: \"NaN\" detected. ";
		std::cout << "Possible Reason: \n"
			<< "  1. Last time, this program terminated abnormally, leading to broken autosave files.\n"
			<< "  2. Invalid boundary condition.\n"
			<< "  3. When you continue to compute, you use a different boundary condition.\n";
		std::cout << "Computation stopped.\n";
		break;
	}
	default:
		break;
	}
}

myfloat ConsolePrinter::getSpeed(int currentStep, int startStep, myfloat usedTime) {
	return (currentStep - startStep) / usedTime; 
}

void ConsolePrinter::logSpeed(int currentStep, int startStep, myfloat usedTime) {
	myfloat speed = getSpeed(currentStep, startStep, usedTime);
	std::stringstream info_log;
	info_log << "speed: " << speed << " step/s\n";
	LogWriter::log(info_log.str());
}
#endif
