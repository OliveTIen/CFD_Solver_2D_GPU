#ifndef TECPLOT_FILE_READER
#define TECPLOT_FILE_READER
#include <string>

/*
2024-05-27 �����
�ڼ���Ma0.77����ʱ����;����NaN��û��pause�����ļ�����recovery�ļ������⸲��
һ�ֽ���취�Ǵ�Euler����������ʼ����
��һ�ֳ���֮���Ǳ�д��Tecplot����ļ�����Ĵ��롣Tecplot�ǻ��ڸ��ģ���ȡ����Ҫת��������
����tecplot�ļ������õı����GPUID������ID
*/

class TecplotFileReader {
public:
	static TecplotFileReader* getInstance();

private:
	static TecplotFileReader* p_instance;

	TecplotFileReader() {};
};

#endif