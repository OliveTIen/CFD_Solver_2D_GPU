#include "CInput.h"
#include "FieldInitializer.h"
#include "InpMeshReader.h"
#include "SU2MeshReader.h"
#include "TomlFileManager.h"
#include "../legacy/FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/GlobalPara.h"
#include "../global/SystemInfo.h"
#include "../output/LogWriter.h"
#include "../output/ConsolePrinter.h"
#include "ContinueFileReader.h"
#include "../global/CExit.h"
#include "../global/StringProcessor.h"
#include <sstream>
#include <map>
#include "../drivers/RoutineController.h"


/*
�������ֵ�Ƿ������б�
value: �����ֵ
category: ֵ��������
table: �б�������ֵ֪������
b_print: �Ƿ��������Ļ
*/
void check_if_value_belongs_to_table(int value, std::string category, std::map<int,std::string> table, bool b_print = false) {
	/*
	2024-05-01
	�������ֵ�Ƿ������б�
	�����ڣ��������Ӧ�����������򱨴��˳�
	�б�Ϊ��ʱ���൱�ڡ������ڡ�
	*/
	std::stringstream ss;
	auto iter = table.find(value);
	if (iter != table.end()) {
		ss << category << ": " << iter->second << "\n";
		LogWriter::log(ss.str());
		if (b_print)LogWriter::print(ss.str());
	}
	else {
		ss << "invalid " << category << ": " << value << "\n";
		LogWriter::logError(ss.str());
		exit(-1);
	}
	
}

void U2NITS::CInput::readConfig() {
	LogWriter::log(SystemInfo::getCurrentDateTime() + "\n");
	TomlFileManager::getInstance()->readTomlFile(FilePathManager::getInstance()->getTomlFilePath());

	check_if_value_belongs_to_table(RoutineController::getInstance()->getStrategy(), "time strategy",
		{ {0, "null"}, {1, "adaptive CFL"} }, true);
	check_if_value_belongs_to_table(GlobalPara::basic::dimension, "dimension", { {2, "2D"} });
	check_if_value_belongs_to_table(GlobalPara::physicsModel::equation, "equation", { {_EQ_euler, "Euler"}, {_EQ_NS, "NS"} }, true);
	check_if_value_belongs_to_table(GlobalPara::basic::useGPU, "platform", { {0,"CPU(host)"},{1,"GPU(device)"},{2,"half GPU"},{3,"development mode"} }, true);
}


// ��¼�߽����
void log_boundary_condition() {
	using namespace GlobalPara::boundaryCondition::_2D;
	myfloat* inlet_ruvp = inlet::ruvp;
	myfloat* outlet_ruvp = outlet::ruvp;
	myfloat* inf_ruvp = inf::ruvp;
	const int num_ruvp = 4;
	// �߽����
	std::string str;
	str += "BoundaryCondition:\n";
	str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(inlet_ruvp, num_ruvp)
		+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(outlet_ruvp, num_ruvp)
		+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(inf_ruvp, num_ruvp)
		+ "\n";
	LogWriter::log(str);
	// �����ŵ���μ�FieldInitialize����ΪҪ����initial_typeȷ����inf, inlet����outlet


}

void read_mesh_file_and_initialize() {
	/*
	��ͷ��ʼ����ȡ���񣬳�ʼ�������ͱ߽�
	*/
	// ��ȡ����
	std::string dir = FilePathManager::getInstance()->getInputDirectory();// Ŀ¼
	std::string basename = GlobalPara::basic::filename;// �ļ�������
	std::string type = GlobalPara::basic::meshFileType;// ��׺
	if (type == "su2") {
		std::string filename = basename + "." + type;
		LogWriter::log("mesh file: " + filename + "\n");
		SU2MeshReader::read_mesh_and_process(dir + filename, true);
	}
	else {
		LogWriter::logAndPrint("Invalid mesh file type: " + type + " @CInput::readField_1()\n");
		exit(-1);
	}
	// ��ʼ�������ͱ߽�
	FieldInitializer::getInstance()->setInitialAndBoundaryCondition();
}

void U2NITS::CInput::readField_1() {

	/*
	��ConfigҪ�����㣬��ֱ�Ӷ�ȡ�����ļ��������ͷ��ʼ
	*/
	if (GlobalPara::basic::_continue) {
		// ���Ҳ��������ļ������ͷ��ʼ
		if (ContinueFileReader::readContinueFile_1() == -1) {
			LogWriter::logAndPrint("failed to read continue file(pause_*.dat), will read mesh file.\n");
			read_mesh_file_and_initialize();
			GlobalPara::basic::_continue = false;// ��֤_continueΪfalse����ֹ����histWriter��д�ļ�ͷ
		}
	}
	else {
		read_mesh_file_and_initialize();
	}

	// ��־��¼�߽����
	log_boundary_condition();
}

void U2NITS::CInput::readField_2_unused() {
	/*
	��ȡ�����ļ���
	��һ��������Ԫ������һ�к󣬾Ϳ��Կ���solver�е�Ԫ���ڴ�ռ䡣ȱ�������޷���֤��������Ƿ���ȷ
	����������std::vectorһ��һ��push����ɺ󣬿���solver���ݵ��ڴ沢��ʼ����Ȼ��ɾ����ʱ��std::vector
	*/
	LogWriter::logAndPrintError("unimplemented error @CInput::readField_2_unused\n");
	exit(-1);
	
}
