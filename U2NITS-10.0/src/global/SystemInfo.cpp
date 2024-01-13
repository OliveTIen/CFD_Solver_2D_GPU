#include "SystemInfo.h"
#include <time.h>

std::string SystemInfo::getCurrentDateTime() {
	// 1900/01/01 00:00:00
	time_t nowtime;	time(&nowtime); tm p; localtime_s(&p, &nowtime);
	std::string str; char buf[6];
	sprintf_s(buf, _countof(buf), "%04d", p.tm_year + 1900);	str += buf;	str += "/";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_mon + 1);	str += buf;	str += "/";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_mday);	str += buf;	str += " ";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_hour);	str += buf;	str += ":";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_min);	str += buf;	str += ":";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_sec);	str += buf;
	return str;
}

std::string SystemInfo::getCurrentDateTime_suitableForFileName() {
	// 1900-01-01-000000
	time_t nowtime;	time(&nowtime); tm p; localtime_s(&p, &nowtime);
	std::string str; char buf[6];
	sprintf_s(buf, _countof(buf), "%04d", p.tm_year + 1900);	str += buf;	str += "-";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_mon + 1);	str += buf;	str += "-";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_mday);	str += buf;	str += "-";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_hour);	str += buf;	str += "";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_min);	str += buf;	str += "";
	sprintf_s(buf, _countof(buf), "%02d", p.tm_sec);	str += buf;
	return str;
}
