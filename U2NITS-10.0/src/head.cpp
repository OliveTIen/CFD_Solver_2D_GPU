#include "head.h"
using namespace std;


//��̬��Ա������ʼ�������Ҫ�ú�����ʼ�������Է���δ����namespace�У�ʹ����ڵ�ǰ�ļ��ɼ���ע�ⲻҪд����ĺ���
namespace {
	
}

rapidjson::Document global::config = rapidjson::Document();
std::string global::filename = "";
//std::string global::lastfilename = "";
std::string global::exePath = "";//exe����·������б��(\\)
int global::autosaveFileNum = 3;
//int global::dimension = 1;
//bool global::_continue = 1;//continue?[1-yes]
//double global::x1 = 0;
//double global::x2 = 1;
//double global::CFL = 0.4;
//double global::T = 0.01;//������ʱ��
//double global::gamma = 1.4;


//double global::ruvp_inlet[4];
//double global::ruvp_outlet[4];//��̬��Ա�����Զ���ʼ��Ϊ0
//double global::ruvp_inf[4];//��̬��Ա�����Զ���ʼ��Ϊ0
//int global::nt = 100;
//int global::nElement = 100;
//int global::flag_equation = _EULER;
int global::flag_reconstruct = _REC_constant;
//int global::step_per_output = 50;//���ٲ����һ����������Զ�ά��
//int global::flag_output = _OUT_u;
//bool global::output_var_ruvp[4];

void global::getFileName_UseUserInput() {
	std::string lastf;
	if (config.FindMember("lastfilename") != config.MemberEnd()) lastf = config["lastfilename"].GetString();//��config��ֵ����lastf
	if (filename.length() == 0) {


		//������ʾ
		std::cout << "Please input filename (without suffix):" << std::endl;
		if(lastf!="[NULL]")
			//��lastf!="[NULL]"��˵������ֱ�Ӷ�ȡ�ϴ�filename
			std::cout << "(Press a single \"Enter\" to use last filename [" << lastf << "])\n";
		else {
			//��lastf=="[NULL]"��˵����ԭ��û��config�ļ�����˲���ʹ��last filename
			//�������inp�ļ�
		}
		//�����û����룬��Ϊ�ļ�����������Ϊ��(��ֱ�Ӱ�enter)����ʹ��lastf
		char a = 0;
		std::string str;
		while (a != '\n') {
			a = getchar();
			str.push_back(a);
		}
		str.pop_back();//ɾ������\n
		if (str == "")str = lastf;
		filename = str;
	}
	else {
		std::cout << "Filename:" << filename << std::endl;
	}
	//д�����ļ���Ŀ���Ǹ���last filename
	global::writeConfig();
}

std::vector<std::string> global::splitString(std::string tLine) {
	tLine = tLine + ' ';//Ŀ���Ƿ�ֹ���һ������û�б�push_back��tWords
	std::vector<std::string> tWords;
	std::string tWord;
	bool isCurrentMeaningful = 0;
	bool isLastMeaningful = 0;
	const int length = (int)tLine.size();
	//std::cout << "tLine= " << tLine << std::endl;
	//std::cout << "size= " << length << std::endl; �ɹ�����tLine����
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
	//��ȡ·���������ַ���
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
		if (files[i] == foldername)return;//�ļ����Ѿ�����
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
	//����filename��ʼ����ʹ�ã����������".LOG"
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
	//��������config
	config["lastfilename"].SetString(filename.c_str(), (rapidjson::SizeType)filename.size());//filename->"lastfilename"

	//configд��json https://blog.csdn.net/yang_aq/article/details/116934216
	rapidjson::StringBuffer buffer;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> prettyWriter(buffer);//PrettyWriter�Ĵ����Ÿ�ʽ���˸��ÿ�
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
	//��ȡ�����ļ�
	//��ȡjson�ļ�������config
	std::ifstream inf(exePath + ".config\\config.json");
	if (inf.is_open()) {
		std::string json_content((std::istreambuf_iterator<char>(inf)), std::istreambuf_iterator<char>()); //���ļ���������תΪstd::string����
		inf.close(); 
		config.Parse(json_content.c_str());
	}
	else {
		//��ȡʧ�ܣ���ʹ��Ĭ�ϡ�ע�⣬��ʱ����WriteConfig����ΪConfigҪ�ڻ�ȡ�ļ���֮��Write
		writeLog("Prompt: fail to open config.json, use backup config. (global::readConfig)\n");
		useBackupConfig();
	}
	//config��ֵ �Ⱥ�����Ҫ��ʱ���ٶ�ȡ
	//

}

void global::useBackupConfig() {
	//ch_json����config
	//�����������ɵ��ļ�����[NULL]ָʾ
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

	////ԭ���ǲ�ɫ�ģ�����win7�ϲ�����
	//cout << "\033[33;5m";//��˸
	////cout << "yellow\n";
	//cout << "\033[33;1m";//����ɫ
	//cout << R"(  ______            _    _     ______    ____      ____     )" << "\n";
	//cout << R"( |        |      | | \  / |   /         |    \    |    \    )" << "\n";
	//cout << R"( |_____   |      | |  \/  |  |          |     |   |     |   )" << "\n";
	//cout << R"( |         \    /  |      |  |          |-----    |-----    )" << "\n";
	//cout << R"( |          \__/   |      |   \______   |         |         )" << "\n";
	//cout << "\033[35m";//��ɫ
	//cout << "------------------------------------------------------------\n";
	//cout << "\033[36m";//����
	//cout << "Finite Volume Method Solver (version 0.0.0), created by"
	//	<< "\033[32m"//��
	//	<< " tgl\n"
	//	<< "\033[36m"//����
	//	<< "2023.7.13\n";
	//cout << "\033[35m";//����
	//cout << "------------------------------------------------------------\n";
	//cout << "\033[0m";//RESET
	
}

#ifdef _WIN32

COORD cmd::getCursorPosition()                           //ͨ��WindowsAPI������ȡ����λ��
{
	CONSOLE_SCREEN_BUFFER_INFO pBuffer;

	GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &pBuffer);
	//���ñ�׼��������ù��������Ϣ

	return COORD{ pBuffer.dwCursorPosition.X, pBuffer.dwCursorPosition.Y };
	//��װΪ��ʾ�����COORD�ṹ
}

COORD cmd::getScrnInfo()                                     //��ȡ����̨���ڻ�������С
{
	HANDLE hStd = GetStdHandle(STD_OUTPUT_HANDLE);      //��ñ�׼����豸���
	CONSOLE_SCREEN_BUFFER_INFO scBufInf;                //����һ�����ڻ�������Ϣ�ṹ��

	GetConsoleScreenBufferInfo(hStd, &scBufInf);        //��ȡ���ڻ�������Ϣ

	return scBufInf.dwSize;                             //���ش��ڻ�������С
}

void cmd::setCursorPosition(COORD pstn) {
	SetConsoleCursorPosition(GetStdHandle(STD_OUTPUT_HANDLE), pstn);
	//ͨ����׼���������ƹ���λ��
}

void cmd::clearDisplay(COORD firstPst, COORD lastPst)        //���������Ļ���ݣ���firstPst���굽lastPst����֮�������
{
	int yDValue(lastPst.Y - firstPst.Y);                    //��¼��ĩλ���������ֵ�����Ƶ�������

	COORD size(getScrnInfo());                              //��¼Ŀǰ����̨��������С

	setCursorPosition(firstPst);                           //�ƶ���굽��λ��
	for (int y(0); y <= yDValue; y++)           //һ��ѭ�������������
	{
		for (int x(firstPst.X); x <= size.X; x++)           //����ѭ�������ظ����
		{
			std::cout << ' ';                       //���һ���ո�������ԭ���ݣ��ﵽ���Ч��
			int px;                                 //��¼��굱ǰλ�õĺ�����
			if (x != size.X)
				px = x + 1;
			else
				px = 0;
			if (y == yDValue && px == lastPst.X)      //����ĩλ�����Աȣ��ﵽĩλ�ü��˳�ѭ��
				break;
		}
	}
	setCursorPosition(firstPst);
}

void cmd::drawProgressBar(double percent) {
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


COORD cmd::getCursorPosition() { return COORD(); }
//��ȡ����̨���ڻ�������С
COORD cmd::getScrnInfo() { return COORD(); }
//���ÿ���̨���λ��
void cmd::setCursorPosition(COORD pstn) {}
//���p1, p2֮��Ŀ���̨���ݣ�������궨λ��p1
void cmd::clearDisplay(COORD p1, COORD p2) {}
//���ƽ�������0<=percent<=100
void cmd::drawProgressBar(int percent) {}


#endif
