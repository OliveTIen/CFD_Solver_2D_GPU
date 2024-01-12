#include "head.h"
using namespace std;


//静态成员变量初始化，如果要用函数初始化，可以放在未命名namespace中，使其仅在当前文件可见。注意不要写成类的函数
namespace {
	
}

rapidjson::Document global::config = rapidjson::Document();
std::string global::filename = "";
//std::string global::lastfilename = "";
std::string global::exePath = "";//exe所在路径，有斜杠(\\)
int global::autosaveFileNum = 3;
//int global::dimension = 1;
//bool global::_continue = 1;//continue?[1-yes]
//double global::x1 = 0;
//double global::x2 = 1;
//double global::CFL = 0.4;
//double global::T = 0.01;//计算总时间
//double global::gamma = 1.4;


//double global::ruvp_inlet[4];
//double global::ruvp_outlet[4];//静态成员变量自动初始化为0
//double global::ruvp_inf[4];//静态成员变量自动初始化为0
//int global::nt = 100;
//int global::nElement = 100;
//int global::flag_equation = _EULER;
int global::flag_reconstruct = _REC_constant;
//int global::step_per_output = 50;//多少步输出一次流场（针对二维）
//int global::flag_output = _OUT_u;
//bool global::output_var_ruvp[4];

void global::getFileName_UseUserInput() {
	std::string lastf;
	if (config.FindMember("lastfilename") != config.MemberEnd()) lastf = config["lastfilename"].GetString();//将config的值赋给lastf
	if (filename.length() == 0) {


		//输入提示
		std::cout << "Please input filename (without suffix):" << std::endl;
		if(lastf!="[NULL]")
			//若lastf!="[NULL]"则说明可以直接读取上次filename
			std::cout << "(Press a single \"Enter\" to use last filename [" << lastf << "])\n";
		else {
			//若lastf=="[NULL]"则说明是原本没有config文件，因此不能使用last filename
			//检测所有inp文件
		}
		//接受用户输入，作为文件名。若输入为空(即直接按enter)，则使用lastf
		char a = 0;
		std::string str;
		while (a != '\n') {
			a = getchar();
			str.push_back(a);
		}
		str.pop_back();//删除最后的\n
		if (str == "")str = lastf;
		filename = str;
	}
	else {
		std::cout << "Filename:" << filename << std::endl;
	}
	//写配置文件，目的是更新last filename
	global::writeConfig();
}

std::vector<std::string> global::splitString(std::string tLine) {
	tLine = tLine + ' ';//目的是防止最后一个单词没有被push_back进tWords
	std::vector<std::string> tWords;
	std::string tWord;
	bool isCurrentMeaningful = 0;
	bool isLastMeaningful = 0;
	const int length = (int)tLine.size();
	//std::cout << "tLine= " << tLine << std::endl;
	//std::cout << "size= " << length << std::endl; 成功返回tLine长度
	for (int j = 0; j < length; j++) {
		if (tLine[j] != ' ' && tLine[j] != '\t' && tLine[j] != ',')
			isCurrentMeaningful = 1;
		else
			isCurrentMeaningful = 0;
		if (isCurrentMeaningful)
			tWord.push_back(tLine[j]);
		else if (isLastMeaningful) {
			tWords.push_back(tWord);
			tWord.clear();
		}
		isLastMeaningful = isCurrentMeaningful;
	}
	return tWords;
}

void global::iniExePath() {
	exePath = getExePath(1);
}

std::string global::getExePath(int flag) {
	if (exePath != "")return exePath;

#ifdef _WIN32
	//获取路径并返回字符串
	const int maxPathLength = 260;
	char buffer[maxPathLength];
	_getcwd(buffer, maxPathLength);
	std::string str = buffer;
	if (flag == 1) {
		str += std::string("\\") ;
	}
	return str;


#elif defined __linux__
	char* buffer;
	buffer = getcwd(NULL, 0);
	std::string str = buffer;
	free(buffer);
	if (flag == 1) {
		str += std::string("\\");
	}
	return str;
#endif 
}

std::vector<std::string> global::ls(std::string path) {
	path += "*";
	_finddata64i32_t fileInfo;
	std::vector<std::string> files;
	intptr_t hFile = _findfirst(path.c_str(), &fileInfo);
	if (hFile == -1) {
		std::cout << "Error in global::ls()" << std::endl;
		return files;
	}
	do {
		files.push_back(fileInfo.name);
	} while (_findnext(hFile, &fileInfo) == 0);
	return files;
}

void global::createFolder(std::string foldername) {
	std::string path = global::exePath;
	std::vector<std::string> files = global::ls(path);
	for (int i = 0; i < files.size(); i++) {
		if (files[i] == foldername)return;//文件夹已经存在
	}
	foldername = "mkdir " + path + foldername;
	system(foldername.c_str());
}

std::vector<int> global::Words2Ints(std::vector<std::string> words) {
	std::vector<int> ints;
	for (int i = 0; i < words.size(); i++) {
		ints.push_back(std::stoi(words[i]));
	}
	return ints;
}

void global::writeLog(std::string s, bool app) {
	//须在filename初始化后使用，否则会生成".LOG"
	ofstream f;
	if(app)f.open(exePath + "output\\" + filename + ".LOG", ios::app);
	else f.open(exePath + "output\\" + filename + ".LOG");
	f << s;
	f.close();

}

void global::writeLogAndCout(std::string str) {
	writeLog(str, 1);
	std::cout << str;
}

std::string global::currentDateTime() {
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

void global::writeConfig() {
	//变量赋给config
	config["lastfilename"].SetString(filename.c_str(), (rapidjson::SizeType)filename.size());//filename->"lastfilename"

	//config写入json https://blog.csdn.net/yang_aq/article/details/116934216
	rapidjson::StringBuffer buffer;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> prettyWriter(buffer);//PrettyWriter的大括号格式化了更好看
	config.Accept(prettyWriter);
	std::string content = buffer.GetString();
	std::ofstream outfile(exePath + ".config\\config.json");
	if (outfile.is_open()) {
		outfile << content;
		outfile.close();
	}
	else std::cout << "Error: fail to open config.txt, in void global::writeConfig()\n";

}

void global::readConfig() {
	//读取配置文件
	//读取json文件，赋给config
	std::ifstream inf(exePath + ".config\\config.json");
	if (inf.is_open()) {
		std::string json_content((std::istreambuf_iterator<char>(inf)), std::istreambuf_iterator<char>()); //将文件的数据流转为std::string类型
		inf.close(); 
		config.Parse(json_content.c_str());
	}
	else {
		//读取失败，则使用默认。注意，此时不必WriteConfig，因为Config要在获取文件名之后Write
		writeLog("Prompt: fail to open config.json, use backup config. (global::readConfig)\n");
		useBackupConfig();
	}
	//config的值 等后面需要的时候再读取
	//

}

void global::useBackupConfig() {
	//ch_json赋给config
	//由于是新生成的文件，用[NULL]指示
	const char* ch_json = R"({
		"lastfilename":"[NULL]",
		"version":"1.0"
		})";
	config.Parse(ch_json);
}

void cmd::printHeader() {
	//
	cout << R"(  ______            _    _     ______    ____      ____     )" << "\n";
	cout << R"( |        |      | | \  / |   /         |    \    |    \    )" << "\n";
	cout << R"( |_____   |      | |  \/  |  |          |     |   |     |   )" << "\n";
	cout << R"( |         \    /  |      |  |          |-----    |-----    )" << "\n";
	cout << R"( |          \__/   |      |   \______   |         |         )" << "\n";
	cout << "------------------------------------------------------------\n";
	cout << "Finite Volume Method Solver (version 0.0.0), created by"
		<< " tgl\n"
		<< "2023.8.18\n";
	cout << "------------------------------------------------------------\n";

	////原本是彩色的，但在win7上不适配
	//cout << "\033[33;5m";//闪烁
	////cout << "yellow\n";
	//cout << "\033[33;1m";//明黄色
	//cout << R"(  ______            _    _     ______    ____      ____     )" << "\n";
	//cout << R"( |        |      | | \  / |   /         |    \    |    \    )" << "\n";
	//cout << R"( |_____   |      | |  \/  |  |          |     |   |     |   )" << "\n";
	//cout << R"( |         \    /  |      |  |          |-----    |-----    )" << "\n";
	//cout << R"( |          \__/   |      |   \______   |         |         )" << "\n";
	//cout << "\033[35m";//紫色
	//cout << "------------------------------------------------------------\n";
	//cout << "\033[36m";//天蓝
	//cout << "Finite Volume Method Solver (version 0.0.0), created by"
	//	<< "\033[32m"//绿
	//	<< " tgl\n"
	//	<< "\033[36m"//天蓝
	//	<< "2023.7.13\n";
	//cout << "\033[35m";//暗紫
	//cout << "------------------------------------------------------------\n";
	//cout << "\033[0m";//RESET
	
}

#ifdef _WIN32

COORD cmd::getCursorPosition()                           //通过WindowsAPI函数获取光标的位置
{
	CONSOLE_SCREEN_BUFFER_INFO pBuffer;

	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &pBuffer);
	//利用标准输出句柄获得光标坐标信息

	return COORD{ pBuffer.dwCursorPosition.X, pBuffer.dwCursorPosition.Y };
	//封装为表示坐标的COORD结构
}

COORD cmd::getScrnInfo()                                     //获取控制台窗口缓冲区大小
{
	HANDLE hStd = GetStdHandle(STD_OUTPUT_HANDLE);      //获得标准输出设备句柄
	CONSOLE_SCREEN_BUFFER_INFO scBufInf;                //定义一个窗口缓冲区信息结构体

	GetConsoleScreenBufferInfo(hStd, &scBufInf);        //获取窗口缓冲区信息

	return scBufInf.dwSize;                             //返回窗口缓冲区大小
}

void cmd::setCursorPosition(COORD pstn) {
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pstn);
	//通过标准输出句柄控制光标的位置
}

void cmd::clearDisplay(COORD firstPst, COORD lastPst)        //清除部分屏幕内容，从firstPst坐标到lastPst坐标之间的内容
{
	int yDValue(lastPst.Y - firstPst.Y);                    //记录首末位置纵坐标差值，控制迭代次数

	COORD size(getScrnInfo());                              //记录目前控制台缓冲区大小

	setCursorPosition(firstPst);                           //移动光标到首位置
	for (int y(0); y <= yDValue; y++)           //一层循环控制清除行数
	{
		for (int x(firstPst.X); x <= size.X; x++)           //二层循环避免重复清除
		{
			std::cout << ' ';                       //输出一个空格来覆盖原内容，达到清除效果
			int px;                                 //记录光标当前位置的横坐标
			if (x != size.X)
				px = x + 1;
			else
				px = 0;
			if (y == yDValue && px == lastPst.X)      //与光标末位置作对比，达到末位置即退出循环
				break;
		}
	}
	setCursorPosition(firstPst);
}

void cmd::drawProgressBar(double percent) {
	if (percent > 100)percent = 100;
	if (percent < 0)percent = 0;
	int nTotal = 45;//总长度，单位char
	const int nBlock = int(nTotal * percent / 100.0);//块个数
	for (int i = 0; i < nBlock; i++) {
		std::cout << "█";
	}
	for (int i = 0; i < nTotal - nBlock; i++) {
		std::cout << " ";
	}
	std::cout << "| " << percent << "%";
}

#elif defined __linux__


COORD cmd::getCursorPosition() { return COORD(); }
//获取控制台窗口缓冲区大小
COORD cmd::getScrnInfo() { return COORD(); }
//设置控制台光标位置
void cmd::setCursorPosition(COORD pstn) {}
//清除p1, p2之间的控制台内容，并将光标定位于p1
void cmd::clearDisplay(COORD p1, COORD p2) {}
//绘制进度条。0<=percent<=100
void cmd::drawProgressBar(int percent) {}


#endif
