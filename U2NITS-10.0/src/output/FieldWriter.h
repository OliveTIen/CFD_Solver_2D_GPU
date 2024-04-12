#ifndef FIELD_WRITER_H
#define FIELD_WRITER_H
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "../global/GlobalPara.h"
#include "../Node_2D.h"
#include "../Element_2D.h"
#include "../gpu/datatype/NodeSoA.h"
#include "../gpu/datatype/ElementSoA.h"
#include "../gpu/datatype/FieldSoA.h"
/*
һ��û�����ݵ��࣬���ȫ���Ǿ�̬��Ա��������ô���������ռ��Ч����ͬ��
*/

class FieldWriter {
private:
	static int numTecplotFileWritten;

public:
	FieldWriter() {};

	static void writeTecplotFile_GPU(
		double t_current,
		std::string filePath,
		std::string title,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::OutputNodeFieldSoA& output_node_field	);

	static void writeContinueFile_GPU(
		int i_step,
		double t_current,
		std::string filePath,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::FieldSoA& elementField);
	// �����Ux Uy
	static void writeContinueFile_2(
		int i_step,
		double t_current,
		std::string filePath,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::FieldSoA& elementField);

	// ��������ļ�����
	static int getNumTecplotFileWritten();
};

#endif