#ifndef _CONTINUEFILEREADER_H_
#define _CONTINUEFILEREADER_H_

#include <string>
#include <vector>


/*
����FVM_2D::readContinueFile()����������ָ����ʼ���Ȳ������ܸ��ӣ�
Ŀǰ�Ȳ�æ��ɸ���
*/
class ContinueFileReader {
	// ���齫reader writer�ϲ� �Ž�fileio�ļ��� ��ѯһ�����������Ŀ¼��֯

public:
	static int readContinueFile_1();
	// ���Ux Uy Ŀǰ�ò���
	static int readContinueFile_2_unused_addUxUy();

	static int readContinueFile_1_1();

private:
	static void process(int maxNodeID, int maxElementID, std::vector<std::vector<int>>& edges_of_all_sets);

};

#endif // _CONTINUEFILEREADER_H_