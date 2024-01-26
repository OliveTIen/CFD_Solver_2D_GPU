
#ifndef STRING_PROCESSOR_H
#define STRING_PROCESSOR_H

#include <vector>
#include <string>

class StringProcessor {
public:
	//// �ַ�������ת��
	// �ַ�������תint����
	static std::vector<int> stringVector_2_intVector(std::vector<std::string> words);

	//// �ַ���ת��
	// �ַ���ת�ַ������顣���ո���Ʊ���ָ�
	static std::vector<std::string> splitString(std::string tLine);

	//// �ַ����޸�
	// �ַ�����'='�滻Ϊ" = " Դ��CFD-bighw
	static std::string replaceCharInString(std::string word, char oldChar = '=', std::string newChar = " = ");
	// ɾ��flag���������ַ�
	static std::string deleteCommentByChar(std::string word, char flag = '%');

	//// �ַ�������
	// ʱ��ת�ַ��������룺�룬�����ʱ:��:��
	static std::string timeFormat(int dt);
	// double����ת�ַ���������ת��Ϊ{}��ʾ����ʽ
	static std::string doubleArray_2_string(double* U, int length);

};

#endif // !STRING_PROCESSOR_H
