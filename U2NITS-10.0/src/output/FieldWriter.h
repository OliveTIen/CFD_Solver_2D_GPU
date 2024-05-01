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
#include "../gpu/datatype/EdgeSoA.h"
#include "../gpu/datatype/BoundaryV2.h"
/*
һ��û�����ݵ��࣬���ȫ���Ǿ�̬��Ա��������ô���������ռ��Ч����ͬ��
*/

class FieldWriter {
private:
	// �ο���������ȫ�����õ����ݣ�����Ҫÿ���ڵ㶼����һ��
	struct ReferenceData {
		myfloat p_inf = 0.0;// ����ѹ��
		myfloat p_dynamic_inf = 0.0;// ������ѹ

		std::vector<myfloat> CD_vector;// ��edgeSet������
		std::vector<myfloat> CL_vector;// ��edgeSet������

		// ���²ο����ݵ�����ѹ���Ͷ�ѹ��ÿ�����tecplotǰҪ����һ��
		void update_farfield();
		// �����edgeSet��������������m_referenceData
		void calculate_edgeSet_force(GPU::EdgeSoA& edge_host, GPU::OutputNodeFieldSoA& output_node_field, GPU::BoundaryV2& boundary_host_new);
	};
	struct HistData {
	public:
		int iteration = 0;// ��������
		myfloat physical_time = 0.0;// ����ʱ��
		std::vector<myfloat> residual;// �غ����в�
		myfloat CD_wall;// wall�߽������ϵ����x������
		myfloat CL_wall;// wall�߽������ϵ����y������
		// ������������
		void update(int _iteration, myfloat _physical_time, const std::vector<myfloat>& _res, ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new);
		void update(int _iteration, myfloat _physical_time, const double _res[4], ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new);
		// ����wall������ϵ��
		void set_CDCL_wall(ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new);
	};

	// �������Ԥ��
	enum OutputSchemePresetEnum {
		type_ruvp, // ���������ruvp
		type_ruvp_Cp,
		type_ruvp_Ma,
		type_full
	};
public:
	// �������
	class FieldOutputScheme {
	public:
		int preset = 0;// Ԥ��
		bool rho = true;
		bool u = true;
		bool v = true;
		bool p = true;
		bool Cp = true;//ѹ��ϵ��
		bool Ma = true;
		
		FieldOutputScheme(int _preset) {
			preset = _preset;
			updateBoolDataByPreset();
		}
		// ��Ԥ�������������ֵ
		void updateBoolDataByPreset();
		// ���������
		std::string get_variable_names();
		// �������ֵ�����ڵ����ֵת��Ϊstring���ÿո�ָ����������з�
		std::string get_variable_value_string(myint nodeID, GPU::NodeSoA& node_host, GPU::OutputNodeFieldSoA& output_node_field);
	};
	class HistOutputScheme {
	public:

		// ���������
		std::string get_variable_names();
		std::string get_variable_value_string(HistData& histData);
	};

private:
	static FieldWriter* p_instance;// ����ָ��

	int m_numTecplotFileWritten = 0;// ����������ļ���
	bool m_initialized = false;// �Ƿ��ѳ�ʼ��
	bool m_called_write_tecplot_hist_file = false;// �Ƿ��ѵ���write tecplot hist file��Ĭ��Ϊfalse
	FieldOutputScheme m_outputScheme_volume = FieldOutputScheme(type_ruvp_Ma);// �����������������
	FieldOutputScheme m_outputScheme_boundary = FieldOutputScheme(type_ruvp_Cp);// �������������������
	HistOutputScheme m_outputScheme_hist;// ��������������

	ReferenceData m_referenceData;// �ο����ݣ�������
	HistData m_histData;
	
private:
	FieldWriter() {};
	
public:
	// ��ȡ��ǰʵ��
	static FieldWriter* getInstance();

	// ��ʼ�����������ֻ��readConfigʱ����
	void initialize_outputScheme_usingConfig();

	void write_tecplot_volume_file(
		double t_current,
		std::string filePath,
		std::string title,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		GPU::OutputNodeFieldSoA& output_node_field
	);

	void write_tecplot_boundary_file(
		double t_current,
		std::string filePath,
		GPU::NodeSoA& node_host,
		GPU::EdgeSoA& edge_host,
		GPU::OutputNodeFieldSoA& output_node_field, 
		GPU::BoundaryV2& boundary_host_new
	);

	void write_tecplot_hist_file(
		std::string filePath,
		int iteration,
		double t_physics,
		double residual[4],
		GPU::EdgeSoA& edge_host, GPU::OutputNodeFieldSoA& output_node_field, GPU::BoundaryV2& boundary_host_new
	);

	void writeContinueFile_1(
		int i_step,
		double t_current,
		std::string filePath,
		GPU::NodeSoA& nodes,
		GPU::ElementSoA& elements,
		myfloat* elementField_U[4]);

	// ��������ļ�����
	int getNumTecplotFileWritten() {
		return m_numTecplotFileWritten;
	}

	// ����������������ڴ棬����ruvp������������ֻ��GPUSolver2��allocateHostMemory����
	void allocNodeFieldDataUsingOutputScheme(GPU::OutputNodeFieldSoA& nodeField, myint num_node);
	// �ͷ�outputNodeField���ֶ�����Ķ��ڴ�ruvp��ֻ��GPUSolver2����
	void freeNodeFieldData(GPU::OutputNodeFieldSoA& nodeField);
	// ����nodeField
	void update_nodeField();
	// ����nodeField.ruvp������������element_vruvpת��Ϊ�������nodeField��ruvp
private:
	// 
	void update_nodeField_ruvp(GPU::ElementSoA& element, GPU::OutputNodeFieldSoA& nodeField, myfloat* element_vruvp[4]);
	// ����nodeField�������������������μ�write_tecplot_boundary_file��Ҫ�����߽�edge������ȫ��node
	void update_nodeField_other_variables(GPU::OutputNodeFieldSoA& nodeField);
	
};

#endif