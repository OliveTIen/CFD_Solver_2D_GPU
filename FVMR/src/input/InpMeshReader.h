#ifndef INP_MESH_READER
#define INP_MESH_READER
#include <string>
// �������ڰ�װ��ȡ����Ĳ���������FVM_2D�����ӷ�׳̶�
class InpMeshReader {

public:
	// ��ȡinp�������ݣ���ʼ���������
	static int readGmeshFile(std::string filepath);
	// [todo]
	static int readGmeshFile_2(std::string filepath);
	// [�ѷ���]
	static int readSU2File(std::string suffix);

public:
	// [�ѷ���]
	class SU2 {
	public:
		enum State {
			S_START,
			S_NDIME,
			S_NPOIN,
			S_NELEM,
			S_mark,
			S_END
		};
		enum ErrorCode {
			E_no_error,
			E_invalid_first_word,
			E_invalid_squeue,
			E_open_file_error,
			E_meet_number_in_the_start,
			E_invalid_word_num,
			E_invalid_word
		};

	public:
		static int readFile(std::string suffix);

	};
};
#endif // !GMESH_MESH_READER
