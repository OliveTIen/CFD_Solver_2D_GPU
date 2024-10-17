#pragma once
#include <windows.h>
#include <string>
#include <iostream>
#include "typeDefine.h"

namespace GUI {

	class ConsoleProxy {
	public:
		static void showCursor(bool show = true);
		static COORD getCursorPosition();
		static void setCursorPosition(COORD p);
		static void clearDisplay(COORD p1, COORD p2);
	};
}


