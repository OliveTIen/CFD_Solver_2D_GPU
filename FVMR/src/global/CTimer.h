#ifndef CTIMER_H
#define CTIMER_H
#include <ctime>
class CTimer {
private:
	clock_t m_clock;

public:
	CTimer() { start(); }
	// 重设起始时间
	void start() {
		m_clock = clock();
	}
	double getTime() {
		return (double)(clock() - m_clock) / CLOCKS_PER_SEC;
	}
};

#endif // !CTIMER_H
