#include "CInput.h"
#include "FieldInitializer.h"
#include "InpMeshReader.h"
#include "SU2MeshReader.h"
#include "TomlFileManager.h"
#include "../FVM_2D.h"
#include "../global/FilePathManager.h"
#include "../global/GlobalPara.h"
#include "../global/SystemInfo.h"
#include "../output/LogWriter.h"
#include "../output/ConsolePrinter.h"
#include "ContinueFileReader.h"
#include "../global/CExit.h"
#include "../global/StringProcessor.h"


void U2NITS::CInput::readConfig() {
	LogWriter::log(SystemInfo::getCurrentDateTime() + "\n");
	TomlFileManager::getInstance()->readTomlFile(FilePathManager::getInstance()->getTomlFilePath());


	if (GlobalPara::basic::dimension != 2) {
		LogWriter::logAndPrintError("invalid dimension type.\n");
		CExit::saveAndExit(2);
	}
	if (GlobalPara::physicsModel::equation != _EQ_euler) {
		LogWriter::logAndPrintError("Invalid equation type.\n");
		CExit::saveAndExit(_EQ_euler);
	}
	if (GlobalPara::basic::useGPU == 0) {
		LogWriter::logAndPrint("use CPU.\n");
	}
	else if(GlobalPara::basic::useGPU == 1){
		LogWriter::logAndPrint("use GPU.\n");
	}
	else if (GlobalPara::basic::useGPU == 2) {
		LogWriter::logAndPrint("use half GPU.\n");
	}
	else if (GlobalPara::basic::useGPU == 3) {
		LogWriter::logAndPrint("enter development mode.\n",LogWriter::Warning,LogWriter::Warning);
	}
	else {
		LogWriter::logAndPrintError("invalid useGPU value.\n");
		CExit::saveAndExit(1);
	}
}

void U2NITS::CInput::readField_1() {

	//// ��ʼ��
	// ���Զ�ȡ�����ļ�������ȡʧ�ܻ���_continue==false���ͷ��ʼ
	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = ContinueFileReader::readContinueFile_1();
		if (flag_readContinue == -1) {
			LogWriter::logAndPrint("Fail to read previous mesh(pause_*.dat). Will try to start from zero again.\n");
			startFromZero = true;
			// ��ֹ����histWriter��д�ļ�ͷ
			GlobalPara::basic::_continue = false;
		}
	}
	else {
		startFromZero = true;
	}
	// ��ͷ��ʼ����ȡ���񡢳�ʼ��
	if (startFromZero) {
		const std::string& type = GlobalPara::basic::meshFileType;
		if (type == "inp") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = InpMeshReader::readGmeshFile(dir + GlobalPara::basic::filename + ".inp");
			if (flag_readMesh == -1) {
				LogWriter::logAndPrintError("Fail to read mesh. Program will exit.\n");
				return;//�˳�
			}
			FieldInitializer::setInitialAndBoundaryCondition();
		}
		else if (type == "su2") {
			std::string dir = FilePathManager::getInstance()->getInputDirectory();
			int flag_readMesh = SU2MeshReader::readFile_2(dir + GlobalPara::basic::filename + ".su2", true);
			if (flag_readMesh == -1) {
				LogWriter::logAndPrint("Fail to read mesh. Program will exit.\n");
				return;//�˳�
			}
			FieldInitializer::setInitialAndBoundaryCondition();
		}
		else {
			LogWriter::logAndPrint("Invalid mesh file type: " + type + ". Program will exit.\n");
			return;
		}

	}

	// ��־��¼�߽����
	auto writeBoundaryCondition = [] (double* inlet_ruvp, double* outlet_ruvp, double* inf_ruvp, const int num_ruvp)->void{
		std::string str;
		str += "BoundaryCondition:\n";
		str += "inlet::ruvp\t" + StringProcessor::doubleArray_2_string(inlet_ruvp, num_ruvp)
			+ "\noutlet::ruvp\t" + StringProcessor::doubleArray_2_string(outlet_ruvp, num_ruvp)
			+ "\ninf::ruvp\t" + StringProcessor::doubleArray_2_string(inf_ruvp, num_ruvp)
			+ "\n";
		LogWriter::log(str);
	};
	using namespace GlobalPara::boundaryCondition::_2D;
	writeBoundaryCondition(inlet::ruvp, outlet::ruvp, inf::ruvp, 4);
}

void U2NITS::CInput::readField_2_unused() {
	/*
	��ȡ�����ļ���
	��һ��������Ԫ������һ�к󣬾Ϳ��Կ���solver�е�Ԫ���ڴ�ռ䡣ȱ�������޷���֤��������Ƿ���ȷ
	����������std::vectorһ��һ��push����ɺ󣬿���solver���ݵ��ڴ沢��ʼ����Ȼ��ɾ����ʱ��std::vector
	*/

	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = ContinueFileReader::readContinueFile_2();
		if (flag_readContinue != 0) {
			LogWriter::logAndPrint("Fail to read previous mesh(pause_*.dat).\n");
			startFromZero = true;
			GlobalPara::basic::_continue = false;// ��ֹ����histWriter��д�ļ�ͷ
		}
	}
	else {
		startFromZero = true;
	}
	// ��ͷ��ʼ����ȡ���񡢳�ʼ��
	if (startFromZero) {
		const std::string& type = GlobalPara::basic::meshFileType;
		std::string dir = FilePathManager::getInstance()->getInputDirectory();
		int flag_readMesh = 502;
		if (type == "inp") {
			flag_readMesh = InpMeshReader::readGmeshFile_2(dir + GlobalPara::basic::filename + ".inp");
		}
		else if (type == "su2") {
			flag_readMesh = SU2MeshReader::readFile_2(dir + GlobalPara::basic::filename + ".su2", true);
		}
		else {
			LogWriter::logAndPrintError("Invalid mesh file type: " + type + "\n");
			flag_readMesh = 404;
		}
		if (flag_readMesh != 0) {
			LogWriter::logAndPrintError("Fail to read mesh. Program will exit.\n");
			CExit::saveAndExit(-3);//�˳�
		}
		FieldInitializer::setInitialAndBoundaryCondition();

	}
}
