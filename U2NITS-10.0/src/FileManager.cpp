#include "FileManager.h"

void FileManager::useDefaultInputPara_and_writeInputPara() {
	global::writeLogAndCout("use default file: input.json\n");
	inputpara.SetObject();
	rapidjson::Document::AllocatorType& allocator = inputpara.GetAllocator();

	rapidjson::Value basic(rapidjson::kObjectType);
	basic.AddMember("dimension", 2, allocator);
	basic.AddMember("continue", 1, allocator);
	inputpara.AddMember("basic", basic, allocator);//需要在basic.AddMember之后使用，否则报错

	rapidjson::Value space(rapidjson::kObjectType);
	rapidjson::Value space_1D(rapidjson::kObjectType);
	space_1D.AddMember("nElement", 300, allocator);
	space_1D.AddMember("x1", 0.0, allocator);
	space_1D.AddMember("x2", 1.0, allocator);
	rapidjson::Value space_2D(rapidjson::kObjectType);
	space.AddMember("space_1D", space_1D, allocator);
	space.AddMember("space_2D", space_2D, allocator);
	inputpara.AddMember("space", space, allocator);

	rapidjson::Value time(rapidjson::kObjectType);
	time.AddMember("CFL", 0.6, allocator);
	time.AddMember("T", 0.01, allocator);
	inputpara.AddMember("time", time, allocator);

	rapidjson::Value physicsModel(rapidjson::kObjectType);
	physicsModel.AddMember("gamma", 1.4, allocator);
	physicsModel.AddMember("equation:1-Eluer,2-NS", 1, allocator);
	inputpara.AddMember("physicsModel", physicsModel, allocator);

	rapidjson::Value boundaryCondition(rapidjson::kObjectType);
	rapidjson::Value bC_1D(rapidjson::kObjectType);
	rapidjson::Value bC_2D(rapidjson::kObjectType);
	rapidjson::Value bC_2D_inlet(rapidjson::kObjectType);
	bC_2D_inlet.AddMember("rho", 1.0, allocator);
	bC_2D_inlet.AddMember("u", 500.0, allocator);
	bC_2D_inlet.AddMember("v", 0.0, allocator);
	bC_2D_inlet.AddMember("p", 1.1e5, allocator);
	bC_2D_inlet.AddMember("Ma", 0.5, allocator);
	bC_2D_inlet.AddMember("AOA[degree]", 0.0, allocator);
	rapidjson::Value bC_2D_outlet(rapidjson::kObjectType);
	bC_2D_outlet.AddMember("rho", 1.0, allocator);
	bC_2D_outlet.AddMember("u", 300.0, allocator);
	bC_2D_outlet.AddMember("v", 0.0, allocator);
	bC_2D_outlet.AddMember("p", 1.0e5, allocator);
	bC_2D_outlet.AddMember("Ma", 1.5, allocator);
	bC_2D_outlet.AddMember("AOA[degree]", 0.0, allocator);
	rapidjson::Value bC_2D_inf(rapidjson::kObjectType);
	bC_2D_inf.AddMember("rho", 1.0, allocator);
	bC_2D_inf.AddMember("u", 300.0, allocator);
	bC_2D_inf.AddMember("v", 0.0, allocator);
	bC_2D_inf.AddMember("p", 1.0e5, allocator);
	bC_2D_inf.AddMember("Ma", 0.8, allocator);
	bC_2D_inf.AddMember("AOA[degree]", 1.25, allocator);
	bC_2D.AddMember("inlet", bC_2D_inlet, allocator);
	bC_2D.AddMember("outlet", bC_2D_outlet, allocator);
	bC_2D.AddMember("inf", bC_2D_inf, allocator);
	boundaryCondition.AddMember("1D", bC_1D, allocator);
	boundaryCondition.AddMember("2D", bC_2D, allocator);
	inputpara.AddMember("boundaryCondition", boundaryCondition, allocator);

	rapidjson::Value output(rapidjson::kObjectType);
	output.AddMember("step_per_output", 50, allocator);
	rapidjson::Value output_var(rapidjson::kObjectType);
	output_var.AddMember("rho", 1, allocator);
	output_var.AddMember("u", 1, allocator);
	output_var.AddMember("v", 1, allocator);
	output_var.AddMember("p", 1, allocator);
	output.AddMember("output_var:bool", output_var, allocator);
	inputpara.AddMember("output", output, allocator);

	createFolder("input");
	writeInputPara(global::exePath + "input\\input.json");

}

void FileManager::writeInputPara(std::string filepath_name) {
	rapidjson::StringBuffer buffer;
	rapidjson::PrettyWriter<rapidjson::StringBuffer> prettyWriter(buffer);  //格式化
	inputpara.Accept(prettyWriter);
	std::string content = buffer.GetString();
	std::ofstream outfile(filepath_name);
	if (outfile.is_open()) {
		outfile << content;
		outfile.close();
		//std::cout << "输出文件内容：\n";
		//std::cout << content;
	}
	else {
		std::cout << "(FileManager::writeInputPara)Error: fail to open " + filepath_name + "\n";
	}
}

void FileManager::readInputPara_initGlobalPara(std::string filepath_name) {
	//读取初始参数，初始化GlobalPara

	//读取输入参数，存入json树
	std::ifstream infile(filepath_name);
	if (infile.is_open()) {
		global::writeLogAndCout("read input parameter from " + filepath_name + "\n");
		std::string json_content((std::istreambuf_iterator<char>(infile)), std::istreambuf_iterator<char>()); //将文件的数据流转换位std::string类型
		infile.close();
		inputpara.Parse(json_content.c_str());
	}
	else {
		global::writeLogAndCout("(FileManager::readInputPara)Error: fail to open " + filepath_name + "\n");
		useDefaultInputPara_and_writeInputPara();
	}

	//用json树初始化GlobalPara
	useJsonTreeToUpdateGlobalPara();
}

void FileManager::checkIntegrity() {
}

void FileManager::useJsonTreeToUpdateGlobalPara() {
	using namespace GlobalPara;
	const char* err = "Error: Unexisted member, in function FileManager::updateData()\n";
	try {
		if (inputpara.HasMember("basic")) {
			rapidjson::Value& basic = inputpara["basic"];
			getMemberValue(basic, "dimension", basic::dimension);
			getMemberValue(basic, "continue", basic::_continue);
		}
		else throw err;

		if (inputpara.HasMember("space") && inputpara["space"].HasMember("space_1D")) {
			rapidjson::Value& space_1D = inputpara["space"]["space_1D"];
			getMemberValue<int>(space_1D, "nElement", space::_1D::nElement);
			getMemberValue<double>(space_1D, "x1", space::_1D::x1);
			getMemberValue<double>(space_1D, "x2", space::_1D::x2);
		}
		else throw err;

		if (inputpara.HasMember("time")) {
			rapidjson::Value& time = inputpara["time"];
			getMemberValue(time, "CFL", time::CFL);
			getMemberValue(time, "T", time::T);
		}
		else throw err;

		if (inputpara.HasMember("physicsModel")) {
			rapidjson::Value& physicsModel = inputpara["physicsModel"];
			getMemberValue(physicsModel, "gamma", Constant::gamma);
			int equation_;
			getMemberValue(physicsModel, "equation:1-Eluer,2-NS", equation_);
			if (equation_ == 1)physicsModel::equation = _EQ_euler;//变更为physicsModel
			else if (equation_ == 2)physicsModel::equation = _EQ_NS;
		}
		else throw err;

		if (inputpara.HasMember("boundaryCondition") && inputpara["boundaryCondition"].HasMember("2D")
			&& inputpara["boundaryCondition"]["2D"].HasMember("inlet") && inputpara["boundaryCondition"]["2D"].HasMember("outlet")
			&& inputpara["boundaryCondition"]["2D"].HasMember("inf")
			) {
			rapidjson::Value& boundaryCondition = inputpara["boundaryCondition"];
			getMemberValue(boundaryCondition["2D"]["inlet"], "use ruvp", boundaryCondition::_2D::inlet::use_ruvp);
			getMemberValue(boundaryCondition["2D"]["inlet"], "rho", boundaryCondition::_2D::inlet::ruvp[0]);
			getMemberValue(boundaryCondition["2D"]["inlet"], "u", boundaryCondition::_2D::inlet::ruvp[1]);
			getMemberValue(boundaryCondition["2D"]["inlet"], "v", boundaryCondition::_2D::inlet::ruvp[2]);
			getMemberValue(boundaryCondition["2D"]["inlet"], "p", boundaryCondition::_2D::inlet::ruvp[3]);
			getMemberValue(boundaryCondition["2D"]["inlet"], "Ma", boundaryCondition::_2D::inlet::Ma);
			getMemberValue(boundaryCondition["2D"]["inlet"], "AOA[degree]", boundaryCondition::_2D::inlet::AOA);

			getMemberValue(boundaryCondition["2D"]["outlet"], "use ruvp", boundaryCondition::_2D::outlet::use_ruvp);
			getMemberValue(boundaryCondition["2D"]["outlet"], "rho", boundaryCondition::_2D::outlet::ruvp[0]);
			getMemberValue(boundaryCondition["2D"]["outlet"], "u", boundaryCondition::_2D::outlet::ruvp[1]);
			getMemberValue(boundaryCondition["2D"]["outlet"], "v", boundaryCondition::_2D::outlet::ruvp[2]);
			getMemberValue(boundaryCondition["2D"]["outlet"], "p", boundaryCondition::_2D::outlet::ruvp[3]);
			getMemberValue(boundaryCondition["2D"]["outlet"], "Ma", boundaryCondition::_2D::outlet::Ma);
			getMemberValue(boundaryCondition["2D"]["outlet"], "AOA[degree]", boundaryCondition::_2D::outlet::AOA);

			getMemberValue(boundaryCondition["2D"]["inf"], "use ruvp", boundaryCondition::_2D::inf::use_ruvp);
			getMemberValue(boundaryCondition["2D"]["inf"], "rho", boundaryCondition::_2D::inf::ruvp[0]);
			getMemberValue(boundaryCondition["2D"]["inf"], "u", boundaryCondition::_2D::inf::ruvp[1]);
			getMemberValue(boundaryCondition["2D"]["inf"], "v", boundaryCondition::_2D::inf::ruvp[2]);
			getMemberValue(boundaryCondition["2D"]["inf"], "p", boundaryCondition::_2D::inf::ruvp[3]);
			getMemberValue(boundaryCondition["2D"]["inf"], "Ma", boundaryCondition::_2D::inf::Ma);
			getMemberValue(boundaryCondition["2D"]["inf"], "AOA[degree]", boundaryCondition::_2D::inf::AOA);
		}
		else throw err;

		if (inputpara.HasMember("initialCondition")) {
			rapidjson::Value& initialCondition = inputpara["initialCondition"];
			getMemberValue(initialCondition, "type", initialCondition::type);
		}
		else throw err;

		if (inputpara.HasMember("output") && inputpara["output"].HasMember("output_var:bool")) {
			rapidjson::Value& output = inputpara["output"];
			getMemberValue(output, "step_per_output", output::step_per_output);
			getMemberValue(output["output_var:bool"], "rho", output::output_var_ruvp[0]);
			getMemberValue(output["output_var:bool"], "u", output::output_var_ruvp[1]);
			getMemberValue(output["output_var:bool"], "v", output::output_var_ruvp[2]);
			getMemberValue(output["output_var:bool"], "p", output::output_var_ruvp[3]);
		}
		else throw err;


	}
	catch (const char* e) {
		global::writeLogAndCout(e);
	}

}

void FileManager::createFolder(std::string foldername) {
	std::string path = global::exePath;
	std::vector<std::string> files = global::ls(path);
	for (int i = 0; i < files.size(); i++) {
		if (files[i] == foldername)return;//文件夹已经存在
	}
	foldername = "mkdir " + path + foldername;
	system(foldername.c_str());
}
