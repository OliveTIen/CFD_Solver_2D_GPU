#ifndef _CONTINUEFILEREADER_H_
#define _CONTINUEFILEREADER_H_

#include <string>
#include <vector>

/// <summary>
///  ����FVM_2D::readContinueFile()����������ָ����ʼ���Ȳ������ܸ��ӣ�
/// Ŀǰ�Ȳ�æ��ɸ���
/// </summary>
class ContinueFileReader {
	// ���齫reader writer�ϲ� �Ž�fileio�ļ��� ��ѯһ�����������Ŀ¼��֯

public:
	static int readContinueFile_1();
	// ���Ux Uy Ŀǰ�ò���
	static int readContinueFile_2();

private:
	static void read_1() {};// δ���
	static void process(int maxNodeID, int maxElementID, std::vector<std::vector<int>>& edges_of_all_sets);

};

#endif // _CONTINUEFILEREADER_H_