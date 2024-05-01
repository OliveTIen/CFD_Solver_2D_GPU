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
一个没有数据的类，如果全部是静态成员函数，那么它与命名空间的效果相同。
*/

class FieldWriter {
private:
	// 参考物理量。全场共用的数据，不需要每个节点都计算一遍
	struct ReferenceData {
		myfloat p_inf = 0.0;// 来流压力
		myfloat p_dynamic_inf = 0.0;// 来流动压

		std::vector<myfloat> CD_vector;// 各edgeSet的阻力
		std::vector<myfloat> CL_vector;// 各edgeSet的升力

		// 更新参考数据的来流压力和动压。每次输出tecplot前要调用一次
		void update_farfield();
		// 计算各edgeSet的升阻力，存入m_referenceData
		void calculate_edgeSet_force(GPU::EdgeSoA& edge_host, GPU::OutputNodeFieldSoA& output_node_field, GPU::BoundaryV2& boundary_host_new);
	};
	struct HistData {
	public:
		int iteration = 0;// 迭代步数
		myfloat physical_time = 0.0;// 物理时间
		std::vector<myfloat> residual;// 守恒量残差
		myfloat CD_wall;// wall边界的阻力系数。x正方向
		myfloat CL_wall;// wall边界的升力系数。y正方向
		// 设置所有数据
		void update(int _iteration, myfloat _physical_time, const std::vector<myfloat>& _res, ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new);
		void update(int _iteration, myfloat _physical_time, const double _res[4], ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new);
		// 设置wall升阻力系数
		void set_CDCL_wall(ReferenceData& referenceData, GPU::BoundaryV2& boundary_host_new);
	};

	// 输出方案预设
	enum OutputSchemePresetEnum {
		type_ruvp, // 最基本数据ruvp
		type_ruvp_Cp,
		type_ruvp_Ma,
		type_full
	};
public:
	// 输出方案
	class FieldOutputScheme {
	public:
		int preset = 0;// 预设
		bool rho = true;
		bool u = true;
		bool v = true;
		bool p = true;
		bool Cp = true;//压力系数
		bool Ma = true;
		
		FieldOutputScheme(int _preset) {
			preset = _preset;
			updateBoolDataByPreset();
		}
		// 用预设更新其他布尔值
		void updateBoolDataByPreset();
		// 输出变量名
		std::string get_variable_names();
		// 输出变量值，将节点变量值转换为string，用空格分隔，不带换行符
		std::string get_variable_value_string(myint nodeID, GPU::NodeSoA& node_host, GPU::OutputNodeFieldSoA& output_node_field);
	};
	class HistOutputScheme {
	public:

		// 输出变量名
		std::string get_variable_names();
		std::string get_variable_value_string(HistData& histData);
	};

private:
	static FieldWriter* p_instance;// 单例指针

	int m_numTecplotFileWritten = 0;// 已输出流场文件数
	bool m_initialized = false;// 是否已初始化
	bool m_called_write_tecplot_hist_file = false;// 是否已调用write tecplot hist file。默认为false
	FieldOutputScheme m_outputScheme_volume = FieldOutputScheme(type_ruvp_Ma);// 体流场数据输出方案
	FieldOutputScheme m_outputScheme_boundary = FieldOutputScheme(type_ruvp_Cp);// 表面流场数据输出方案
	HistOutputScheme m_outputScheme_hist;// 输出方案。待完成

	ReferenceData m_referenceData;// 参考数据，物理量
	HistData m_histData;
	
private:
	FieldWriter() {};
	
public:
	// 获取当前实例
	static FieldWriter* getInstance();

	// 初始化输出方案，只在readConfig时调用
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

	// 已输出的文件数量
	int getNumTecplotFileWritten() {
		return m_numTecplotFileWritten;
	}

	// 根据输出方案申请内存，包括ruvp和其他变量。只被GPUSolver2的allocateHostMemory调用
	void allocNodeFieldDataUsingOutputScheme(GPU::OutputNodeFieldSoA& nodeField, myint num_node);
	// 释放outputNodeField的手动申请的堆内存ruvp。只被GPUSolver2调用
	void freeNodeFieldData(GPU::OutputNodeFieldSoA& nodeField);
	// 更新nodeField
	void update_nodeField();
	// 更新nodeField.ruvp。将格心数据element_vruvp转化为格点数据nodeField的ruvp
private:
	// 
	void update_nodeField_ruvp(GPU::ElementSoA& element, GPU::OutputNodeFieldSoA& nodeField, myfloat* element_vruvp[4]);
	// 计算nodeField其他场变量。升阻力参见write_tecplot_boundary_file，要遍历边界edge而不是全场node
	void update_nodeField_other_variables(GPU::OutputNodeFieldSoA& nodeField);
	
};

#endif