#include "FVM_2D.h"
#include "output/ConsolePrinter.h"
#include "output/ResidualCalculator.h"
#include "output/LogWriter.h"
#include "global/StringProcessor.h"
#include "global/SystemInfo.h"
#include "global/FilePathManager.h"
#include "output/HistWriter.h"
#include "input/GmeshMeshReader.h"

FVM_2D* FVM_2D::pFVM2D = nullptr;

FVM_2D::FVM_2D() {
	pFVM2D = this;
}

void FVM_2D::run() {
	//尝试读取续算文件。若读取失败或者_continue==false则从头开始
	bool startFromZero = false;
	if (GlobalPara::basic::_continue) {
		int flag_readContinue = readContinueFile();
		if (flag_readContinue == -1) {
			std::string error_msg = "Warning: Fail to read previous mesh(pause_*.dat). ";
			error_msg += "Will try to start from zero again.\n";
			LogWriter::writeLogAndCout(error_msg);
			startFromZero = true;
		}
	}
	else {
		startFromZero = true;
	}
	//从头开始。读取网格、初始化
	if (startFromZero) {
		int flag_readMesh = GmeshMeshReader::readMeshFile(".inp");
		if (flag_readMesh == -1) {
			std::string error_msg = "Error: Fail to read mesh. Program will exit.\n";
			LogWriter::writeLogAndCout(error_msg);
			return;//退出
		}
		setInitialCondition();
	}

	
	//将边界edge打上标签，并检查周期边界完整性
	{
		int ret = boundaryManager.iniBoundaryEdgeSetID_and_iniBoundaryType(this);
		if (ret != 0)return;
	}

	//日志记录边界参数
	logBoundaryCondition();

	//求解器
	if (GlobalPara::physicsModel::equation == _EQ_euler) {
		solve_("_Euler.out", "_Euler.info");
	}
	else {
		std::cout << "Error: Invalid equation type";
	}

}

void FVM_2D::setInitialCondition() {
	switch (GlobalPara::initialCondition::type) {
	case 1:
	{
		//1.常数
		using namespace GlobalPara::boundaryCondition::_2D;
		for (int ie = 0; ie < elements.size(); ie++) {
			Math_2D::ruvp_2_U(inf::ruvp, elements[ie].U, Constant::gamma);
			//if(elements[ie].ID==173)system()
		}
		std::string str;
		str += "InitialCondition:\nUniform flow, ruvp\t" + StringProcessor::DoubleArray_2_String(inf::ruvp, 4) + "\n";
		LogWriter::writeLog(str);
	}
		break;

	case 2:
	{
		//2.均匀流和等熵涡的叠加

		using namespace GlobalPara::boundaryCondition::_2D;
		std::string str;
		str += "InitialCondition:\nUniform flow + isentropicVortex, ruvp of Uniform flow:" + StringProcessor::DoubleArray_2_String(inf::ruvp, 4) + "\n";
		LogWriter::writeLogAndCout(str);
		isentropicVortex_2(5, 5, 5, inf::ruvp);
		//for (int ie = 0; ie < elements.size(); ie++) {
		//	//均匀流+等熵涡 Isentropic vortex
		//	double deltau, deltav, deltaT;
		//	double x = elements[ie].x; double y = elements[ie].y;
		//	double xc = 5; double yc = 5;
		//	isentropicVortex(x, y, xc, yc, 1, deltau, deltav, deltaT);
		//	double ruvp[4]{};
		//	ruvp[1] = inf::ruvp[1] + deltau;
		//	ruvp[2] = inf::ruvp[2] + deltav;
		//	double rho0 = inf::ruvp[0];
		//	double p0 = inf::ruvp[3];
		//	double T0 = p0 / rho0 / Constant::gamma;//p=rho*R*T
		//	double deltap = deltaT * rho0 / (1 - Constant::R * rho0 * T0 / pow(p0, Constant::gamma));
		//	double deltarho = deltap * rho0 / pow(p0, Constant::gamma);
		//	ruvp[0] = rho0 + deltarho;
		//	ruvp[3] = p0 + deltap;
		//	Math_2D::ruvp_2_U(ruvp, elements[ie].U, Constant::gamma);
		//}

	}
		break;

	}

}

void FVM_2D::logBoundaryCondition() {
	//日志记录边界参数
	std::string str;
	str += "BoundaryCondition:\n";
	str += "inlet::ruvp\t" + StringProcessor::DoubleArray_2_String(GlobalPara::boundaryCondition::_2D::inlet::ruvp, 4)
		+ "\noutlet::ruvp\t" + StringProcessor::DoubleArray_2_String(GlobalPara::boundaryCondition::_2D::outlet::ruvp, 4)
		+ "\ninf::ruvp\t" + StringProcessor::DoubleArray_2_String(GlobalPara::boundaryCondition::_2D::inf::ruvp, 4)
		+ "\n";
	LogWriter::writeLog(str);
}

void FVM_2D::solve_(std::string suffix_out, std::string suffix_info) {
	FilePathManager* filePathManager = FilePathManager::getInstance();
	HistWriter histWriter(filePathManager->outputFolder_path + "\\" + GlobalPara::basic::filename + "_hist.dat");
	const int residual_vector_size = 4;
	const double big_value = 1e8;
	double residual_vector[residual_vector_size]{1,1,1,1};
	histWriter.writeHead();

	//输出格式(tecplot)，在文件夹中输出，每一帧存成一个文件
	std::vector<Element_T3> elements_old;

	//变量
	double t = t_previous;//当前物理时间。
	const double T = GlobalPara::time::T;//目标物理时间。用于控制是否终止计算
	int nFiles = 0;//输出文件个数，用于显示
	bool signal_pause = 0;//暂停信号，用于控制
	clock_t start_t;//计时变量，用于计算求解时间(基于CPU周期)
	start_t = clock();
	double calculateTime = 0.0;
	double calculateSpeed = 0.0;
	std::string solveInfo;
	//等熵涡误差计算文件头
	if (GlobalPara::initialCondition::type == 2)writeFileHeader_isentropicVortex();
	COORD p1 = ConsolePrinter::getCursorPosition();
	ConsolePrinter::drawProgressBar(int(t / T));	//绘制进度条
	for (int istep = istep_previous; istep < 1e6 && t < T; istep++) {
		elements_old = elements;
		//计算时间
		caldt();
		dt = (t + dt > T) ? (T - t) : dt;//if (t + dt > T)dt = T - t;
		t += dt;
		//计算通量、时间推进
		solver.evolve(dt);

		//定期输出进度
		if (istep % GlobalPara::output::step_per_print == 1) {
			ConsolePrinter::clearDisplay(p1, ConsolePrinter::getCursorPosition());
			ConsolePrinter::setCursorPosition(p1);
			ConsolePrinter::drawProgressBar(int(t / T * 100 + 2));//+x是为了防止停在99%

			calculateTime = (double)(clock() - start_t) / CLOCKS_PER_SEC;//
			calculateSpeed = ((istep - istep_previous)/calculateTime);
			std::cout;
			solveInfo = ConsolePrinter::printSolveInfo(calculateTime, istep, calculateSpeed, nFiles, t, T, residual_vector);
		}
		//检查非法值
		if (isNan()|| residual_vector[0] > big_value) {
			std::cout << "\nWarning: \"NaN\" detected. ";
			std::cout << "Possible Reason: \n"
				<< "  1. Last time, this program terminated abnormally, leading to broken autosave files.\n"
				<< "  2. Invalid boundary condition.\n"
				<< "  3. When you continue to compute, you use a different boundary condition.\n";
			std::cout << "Computation stopped.\n";
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
			writeContinueFile(filePathManager->getExePath_withSlash() + "output\\nan_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t, istep);
			break;
		}
		//定期输出流场
		if (istep % GlobalPara::output::step_per_output == 1) {
			//.dat 流场显示文件 tecplot格式
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);//若istep小于4位数，则补0
			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t);
			nFiles++;
			//.autosave 自动保存文件 续算文件
			updateAutoSaveFile(t, istep);

		}
		//定期输出残差
		if (istep % GlobalPara::output::step_per_output_hist == 1) {
			//计算残差，结果保存在residual_vector中
			ResidualCalculator::cal_residual(elements_old, elements, ResidualCalculator::NORM_INF, residual_vector);
			histWriter.writeData(istep, residual_vector, residual_vector_size);

			//等熵涡的误差文件 
			if (GlobalPara::initialCondition::type == 2) {
				//求解时间
				//end_t = clock();
				double time_used = (double)(clock() - start_t) / CLOCKS_PER_SEC;
				//输出误差文件
				ResidualCalculator::cal_error_isentropicVortex(0, 0, 10, 10, 5, t, istep, time_used, GlobalPara::boundaryCondition::_2D::inf::ruvp);
			}
		}
		//按esc终止
		if (_kbhit()) {
			char ch = _getch();
			if (ch == 27) {
				char szBuffer[20];
				sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
				writeContinueFile(filePathManager->getExePath_withSlash() + "output\\pause_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat",t,istep);
				std::cout << "[ESC]: Computation Interrupted\n";
				signal_pause = true;
			}
		}
		//达到规定迭代时间，终止
		else if (t >= T) {
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
			writeContinueFile(filePathManager->getExePath_withSlash() + "output\\pause_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t, istep);
			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t);
			std::cout << "Computation finished\n";
			signal_pause = true;
		}
		//残差足够小，终止
		else if (residual_vector[0] <= Constant::epsilon) {
			char szBuffer[20];
			sprintf_s(szBuffer, _countof(szBuffer), "%04d", istep);
			writeContinueFile(filePathManager->getExePath_withSlash() + "output\\pause_" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t, istep);
			writeTecplotFile(filePathManager->getExePath_withSlash() + "output\\" + GlobalPara::basic::filename + "[" + szBuffer + "].dat", t);
			std::cout << "Computation finished as the field is already stable\n";
#ifdef _WIN32
			MessageBox(NULL, "Computation finished", "U2NITS", MB_OK);
#endif // _WIN32
			signal_pause = true;
		}
		//终止
		if (signal_pause) {
			LogWriter::writeLog(solveInfo);
			break;
		}
	}

}

void FVM_2D::caldt() {
	//以下为新代码
	for (int ie = 0; ie < elements.size(); ie++) {
		double LambdaC_i = elements[ie].calLambdaFlux(this);
		double Omega_i = elements[ie].calArea(this);
		double dt_i = GlobalPara::time::CFL * Omega_i / LambdaC_i;
		dt = (std::min)(dt, dt_i);
	}
}

void FVM_2D::writeTecplotFile(std::string f_name, double t_current) {
	//std::cout << "wirtefile:" << f_name << std::endl;
	calculateNodeValue();
	std::ofstream f(f_name);
	if (!f.is_open())std::cout << "Error: Fail to open file " << f_name << std::endl;

	//header
	f << R"(TITLE = "Example: 2D Finite Element Data")" << "\n";
	f << R"(VARIABLES = "X", "Y")";
	if (GlobalPara::output::output_var_ruvp[0])		f << R"(, "rho")";
	if (GlobalPara::output::output_var_ruvp[1])		f << R"(, "u")";
	if (GlobalPara::output::output_var_ruvp[2])		f << R"(, "v")";
	if (GlobalPara::output::output_var_ruvp[3])		f << R"(, "p")";
	f << "\n";

	f << "ZONE NODES="<< nodes.size() 
		<<", ELEMENTS = "<< elements.size() 
		<<", DATAPACKING = POINT, ZONETYPE = FEQUADRILATERAL"
		<<", SOLUTIONTIME="<<std::scientific<< t_current<<std::fixed
		<<std::endl;
	
	//nodes
	for (int i = 0; i < nodes.size(); i++) {
		f << nodes[i].x << ' ' << nodes[i].y;
		if (GlobalPara::output::output_var_ruvp[0])		f << ' ' << rho_nodes[i];
		if (GlobalPara::output::output_var_ruvp[1])		f << ' ' << u_nodes[i];
		if (GlobalPara::output::output_var_ruvp[2])		f << ' ' << v_nodes[i];
		if (GlobalPara::output::output_var_ruvp[3])		f << ' ' << p_nodes[i];
		f << std::endl;
	}

	//elements
	for (int ie = 0; ie < elements.size(); ie++) {
		f << elements[ie].nodes[0] << " ";
		f << elements[ie].nodes[1] << " ";
		f << elements[ie].nodes[2] << " ";
		f << elements[ie].nodes[2] << std::endl;
	}

	f.close();
}

void FVM_2D::writeContinueFile(std::string f_name, double t, int istep) {
	//std::cout << "wirtefile:" << f_name << std::endl;
	int flag_OUT = -1;
	std::ofstream outfile(f_name);
	if (!outfile.is_open())std::cout << "Error: Fail to open file " << f_name << std::endl;

	outfile << "t, istep" << std::endl;
	outfile << t << " " << istep << std::endl;

	outfile << "nodes: ID, x, y" << std::endl;
	for (int in = 0; in < nodes.size(); in++) {
		std::vector<double> uu_i = nodes[in].calNodeValue(this);
		outfile << nodes[in].ID << " ";
		outfile << nodes[in].x << " ";
		outfile << nodes[in].y << std::endl;
	}

	outfile << "elements: ID, nodes, U" << std::endl;
	for (int ie = 0; ie < elements.size(); ie++) {
		outfile << elements[ie].ID << " ";
		outfile << elements[ie].nodes[0] << " ";
		outfile << elements[ie].nodes[1] << " ";
		outfile << elements[ie].nodes[2] << " ";
		outfile << elements[ie].U[0] << " ";
		outfile << elements[ie].U[1] << " ";
		outfile << elements[ie].U[2] << " ";
		outfile << elements[ie].U[3] << std::endl;
	}

	outfile << "boundaries: ID, name; edgeIDs" << std::endl;
	for (int ib = 0; ib < boundaryManager.vBoundarySets.size(); ib++) {
		outfile << boundaryManager.vBoundarySets[ib].ID << " ";
		outfile << boundaryManager.vBoundarySets[ib].name << std::endl;
		//一行太长会导致读取时buffer长度不够，从而导致难以预料的结果(死循环、读取失败等)
		int nEdge = (int)boundaryManager.vBoundarySets[ib].pEdges.size();//控制输出' '还是'\n'
		bool check = 0;//控制
		for (int iEdge = 0; iEdge < nEdge;iEdge++) {
			outfile << boundaryManager.vBoundarySets[ib].pEdges[iEdge]->ID;
			if (iEdge % 10 == 9){
				outfile << std::endl; 
				check = 1;
			}
			else {
				outfile << ' ';
				check = 0;
			}
		}
		if(!check)outfile << std::endl;//若上次输出的' '
	}
	outfile.close();

}

int FVM_2D::readContinueFile() {
	std::string m_path = FilePathManager::getInstance()->getExePath_withSlash() + "output\\";
	std::vector<std::string> files = FilePathManager::ls(m_path);//files为字符串向量，存储了路径下所有文件名
	int index_maxNstep = -1;
	//搜寻是否有pause_filename[xxx].dat文件，若有，则找xxx最大的
	int maxNstep = 0;
	int tmp_index = 0;//用于指示"[", "]"
	std::string tmp_str;
	for (int i = 0; i < files.size(); i++) {
		//搜寻以pause_filename开头的文件名
		std::string str_match = "pause_" + GlobalPara::basic::filename;
		if (files[i].substr(0, int(str_match.size())) == str_match) {
			//剔除"pause_filename["
			int pre_length = int(str_match.size() + 1);//"pause_filename["的长度
			std::string str = files[i].substr(pre_length);
			//剔除"].dat"
			int post_length = 5;//"].dat"的长度
			for (int i = 0; i < post_length; i++) {
				str.pop_back();
			}
			//将xxx存入num
			int num = std::stoi(str);
			if (num > maxNstep) {
				index_maxNstep = i;
				maxNstep = num;
			}
		}
	}
	//若无，则返回-1
	if (index_maxNstep == -1)return -1;
	std::ifstream infile(m_path + files[index_maxNstep]);
	if (!infile) {
		return -1;
	}

	std::cout << "Continue from: "<< files[index_maxNstep] << std::endl;

	int state = -1;//1-Node, 2-Element
	int maxnodeID = 1;
	int maxedgeID = 1;
	int maxelementID = 1;
	int iline_set = 0;//临时变量，用于记录读取set时经过了多少行
	std::vector<std::vector<int>> edges_of_all_sets;//临时变量，存储各set的edgeID
	const int bufferLength = 600;//!!!小心长度不够
	char buffer[bufferLength];
	std::string tLine;
	std::vector<std::string> tWords;
	std::vector<std::vector<std::string>> tWordsMatrix;
	std::vector<int> edge_ID_long;
	int nNullRow = 0;
	bool loop = 1;
	while (loop) {
		//get words and set state
		infile.getline(buffer, bufferLength);
		if (infile.eof())loop=0;
		tLine = buffer;
		tWords = StringProcessor::splitString(tLine);
		if (nNullRow >= 10)loop = 0;
		if (tWords.size() == 0) {
			nNullRow++;
			continue;
		}//空行，强迫开始下一次循环，防止tWords[0]内存错误
		else {
			if (tWords[0] == "t")state = 0;
			if (tWords[0] == "nodes:")state = 1;
			if (tWords[0] == "edges:")state = 2;
			if (tWords[0] == "elements:")state = 3;
			if (tWords[0] == "boundaries:")state = 4;
		}

		//translate
		if (state == 0 && tWords[0].substr(0, 1) != "t") {//t_solve, istep_solve
			t_previous = std::stod(tWords[0]);
			istep_previous = (int)std::stod(tWords[1]);
		}
		else if (state == 1 && tWords[0].substr(0, 1) != "n") {//nodes: ID, x, y
			Node_2D node;
			node.ID = (int)std::stod(tWords[0]);
			node.x = std::stod(tWords[1]);
			node.y = std::stod(tWords[2]);
			nodes.push_back(node);
			maxnodeID = (std::max)(maxnodeID, node.ID);
		}
		else if (state == 3 && tWords[0].substr(0, 1) != "e") {//elements: ID, nodes, U
			Element_T3 ele;
			ele.ID = (int)std::stod(tWords[0]);
			ele.nodes[0] = (int)std::stod(tWords[1]);
			ele.nodes[1] = (int)std::stod(tWords[2]);
			ele.nodes[2] = (int)std::stod(tWords[3]);
			ele.U[0] = std::stod(tWords[4]);
			ele.U[1] = std::stod(tWords[5]);
			ele.U[2] = std::stod(tWords[6]);
			ele.U[3] = std::stod(tWords[7]);
			elements.push_back(ele);
			maxelementID = (std::max)(maxelementID, ele.ID);
		}
		else if (state == 4 && tWords[0].substr(0, 1) != "b") {//boundaries: ID, name, edgeIDs
			iline_set++;
			if (!isdigit(tWords[1][0])) {//若不是数字 1.vBoundarySets新增条目、该条目的部分初始化 2.edges_of_all_sets新增条目
				//set的ID,name初始化
				VirtualBoundarySet_2D vb;
				vb.ID = (int)std::stod(tWords[0]);
				vb.name = tWords[1];
				boundaryManager.vBoundarySets.push_back(vb);
				//edges_of_all_sets新增条目
				edges_of_all_sets.push_back(std::vector<int>());
			}
			else {//若是数字 1.edges_of_all_sets的初始化
				//引用最后一个BoundarySet
				std::vector<int>& edges_of_current_set = edges_of_all_sets[edges_of_all_sets.size() - 1];
				std::vector<int> intVector = StringProcessor::Words2Ints(tWords);
				boundaryManager.attachToVector(edges_of_current_set, intVector);
			}
		}
	}

	infile.close();

	iniPNodeTable(maxnodeID);
	iniEdges();
	iniPEdgeTable();
	iniPElementTable(maxelementID);
	iniElement_xy_pEdges();
	iniNode_neighborElements();
	boundaryManager.iniBoundarySetPEdges_in_readContinueFile(this, edges_of_all_sets);

	return 0;
}

void FVM_2D::updateAutoSaveFile(double t, int istep) {
	//最多保存global::autosaveFileNum个autosave文件，再多则覆写之前的
	writeContinueFile(FilePathManager::getInstance()->getExePath_withSlash() + "output\\autosave" + std::to_string(iCurrentAutosaveFile) + "_" + GlobalPara::basic::filename + ".dat", t, istep);
	iCurrentAutosaveFile++;
	if (iCurrentAutosaveFile >= GlobalPara::output::autosaveFileNum)iCurrentAutosaveFile = 0;
}

void FVM_2D::calculateNodeValue() {

	if (GlobalPara::output::output_var_ruvp[0]) 		rho_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[1]) 		u_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[2]) 		v_nodes.resize(nodes.size());
	if (GlobalPara::output::output_var_ruvp[3]) 		p_nodes.resize(nodes.size());

	for (int i = 0; i < nodes.size(); i++) {
		std::vector<double> U_node = nodes[i].calNodeValue(this);
		double U[4]{ U_node[0],U_node[1],U_node[2],U_node[3] };
		double ruvp[4];
		math.U_2_ruvp(U, ruvp, Constant::gamma);
		if (GlobalPara::output::output_var_ruvp[0]) 		rho_nodes[i] = ruvp[0];
		if (GlobalPara::output::output_var_ruvp[1]) 		u_nodes[i] = ruvp[1];
		if (GlobalPara::output::output_var_ruvp[2]) 		v_nodes[i] = ruvp[2];
		if (GlobalPara::output::output_var_ruvp[3]) 		p_nodes[i] = ruvp[3];
	}

}

void FVM_2D::iniPNodeTable(int maxnodeID) {
	pNodeTable.resize(maxnodeID + 1);
	for (int in = 0; in < nodes.size(); in++) {
		pNodeTable[nodes[in].ID] = &(nodes[in]);
	}
}

void FVM_2D::iniEdges() {
	for (int ie = 0; ie < elements.size(); ie++) {
		int n0, n1, n2;
		n0 = elements[ie].nodes[0];
		n1 = elements[ie].nodes[1];
		n2 = elements[ie].nodes[2];
		iniEdges_registerSingle(n0, n1, &(elements[ie]));
		iniEdges_registerSingle(n1, n2, &(elements[ie]));
		iniEdges_registerSingle(n2, n0, &(elements[ie]));
	}
}

void FVM_2D::iniPEdgeTable() {
	pEdgeTable.resize(edges.size() + 1);
	for (int i = 0; i < edges.size(); i++) {
		pEdgeTable[edges[i].ID] = &(edges[i]);
	}
}

void FVM_2D::iniNode_neighborElements() {
	for (int ie = 0; ie < elements.size(); ie++) {
		elements[ie].nodes;
		for (int jn = 0; jn < 3; jn++) {
			Node_2D* pN = getNodeByID(elements[ie].nodes[jn]);
			pN->neighborElements.push_back(&(elements[ie]));
		}
	}
}

Edge_2D* FVM_2D::isEdgeExisted(int n0, int n1) {
	//n0, n1表示节点ID
	if(edges.size()==0)return nullptr;
	for (int i = 0; i < edges.size(); i++) {
		if (edges[i].nodes[0] == n0 && edges[i].nodes[1] == n1)return &(edges[i]);
		if (edges[i].nodes[0] == n1 && edges[i].nodes[1] == n0)return &(edges[i]);
	}
	return nullptr;
}

void FVM_2D::iniEdges_registerSingle(int n0, int n1, Element_T3* pE) {
	Edge_2D* pEdge = isEdgeExisted(n0, n1);
	if (pEdge == nullptr) {
		Edge_2D tmp_edge;
		tmp_edge.ID = (int)edges.size() + 1;
		tmp_edge.nodes[0] = n0;
		tmp_edge.nodes[1] = n1;
		tmp_edge.pElement_L = pE;
		edges.push_back(tmp_edge);
	}
	else {
		pEdge->pElement_R = pE;
	}
}

void FVM_2D::iniPElementTable(int maxelementID) {
	pElement_T3Table.resize(maxelementID + 1);
	for (int ie = 0; ie < elements.size(); ie++) {
		pElement_T3Table[elements[ie].ID] = &(elements[ie]);
	}
}

void FVM_2D::iniElement_xy_pEdges() {
	for (int ie = 0; ie < elements.size(); ie++) {
		elements[ie].inixy(this); //已经calxy，以后不必担心。但是必须要先读取node再读取element
		elements[ie].iniPEdges(this);
	}
}

Node_2D* FVM_2D::getNodeByID(int ID) {
	//if (pNodeTable.size() == 0) {
	//	std::cout << "Error: uninitialized pNodeTable. Initializing..." << std::endl;
	//	iniPNodeTable();
	//	//return nullptr;
	//}
	return pNodeTable[ID];
}

void FVM_2D::isentropicVortex(double x, double y, double xc, double yc, double chi, double& deltau, double& deltav, double& deltaT) {
	double xbar = x - xc;
	double ybar = y - yc;
	double r2 = xbar * xbar + ybar * ybar;
	const double gamma = Constant::gamma;
	const double PI = Constant::PI;
	deltau = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * (-ybar);
	deltav = chi / 2.0 / PI * exp(0.5 * (1 - r2)) * xbar;
	deltaT = -(gamma - 1) / chi / chi * 8 * gamma * PI * PI * exp(1 - r2);
}

void FVM_2D::isentropicVortex_2(double xc, double yc, double chi, const double* ruvp0) {
	//ruvp0：均匀流参数
	for (int ie = 0; ie < elements.size(); ie++) {
		Element_T3& e = elements[ie];
		double rho, u, v, p;
		double xbar, ybar, r2, du, dv, dT;
		const double PI = Constant::PI;
		const double gamma = Constant::gamma;
		xbar = e.x - xc;
		ybar = e.y - yc;
		r2 = xbar * xbar + ybar * ybar;
		du = chi / 2. / PI * exp(0.5 * (1. - r2)) * (-ybar);
		dv = chi / 2. / PI * exp(0.5 * (1. - r2)) * xbar;
		u = ruvp0[1] + du;
		v = ruvp0[2] + dv;
		dT = -(gamma - 1.) * chi * chi / (8. * gamma * PI * PI) * exp(1. - r2);
		rho = pow(ruvp0[3] + dT, 1. / (gamma - 1.));
		p = rho * (ruvp0[3] + dT);
		double ruvp[4]{ rho,u,v,p };
		Math_2D::ruvp_2_U(ruvp, e.U, gamma);
	}
	//lambda表达式 [ capture ] ( params ) opt -> ret { body; }; http://c.biancheng.net/view/3741.html
	//例如 auto f = [](int a) -> int { return a + 1; };
	//auto fun_WA = [](const double* xy) {
}


void FVM_2D::writeFileHeader_isentropicVortex() {
	std::ofstream outfile(FilePathManager::getInstance()->getExePath_withSlash() + "output\\" + "error_isentropicVortex_" + GlobalPara::basic::filename + ".txt", std::ios::app);//追加模式
	outfile << SystemInfo::getCurrentDateTime() << "\n";
	outfile
		<< "istep" << "\t"
		<< "cpu_time[s]" << "\t"
		<< "error_L1" << "\t"
		<< "error_L2" << "\t"
		<< "error_max"<< "\n";
	outfile.close();
}

bool FVM_2D::isStable(std::vector<Element_T3> old) {
	double sum = 0;
	double res;
	for (int ie = 0; ie < elements.size(); ie++) {
		//sum += abs(elements[ie].calculateu() - old[ie].calculateu());
		res = elements[ie].U[1] - old[ie].U[1];
		sum += abs(res);
	}
	//std::cout << sum << std::endl;
	if (sum <= 1e-12)return 1;
	else return 0;
}

bool FVM_2D::isNan() {
	bool is_nan = 0;
	for (int ie = 0; ie < elements.size(); ie++) {
		//只需要检查1项即可 rho
		if (isnan(elements[ie].U[0])) {
			is_nan = 1;
			const double& x = elements[ie].x;
			const double& y = elements[ie].y;
			std::string str = 
				"Warning: \"NaN\" detected in element (x=" + std::to_string(x) + ", y=" + std::to_string(y) 
				+ ", U[0,1,2,3]=" + std::to_string(elements[ie].U[0]) + ", " + std::to_string(elements[ie].U[1])
				+ ", " + std::to_string(elements[ie].U[2]) + ", " + std::to_string(elements[ie].U[3]) + "\n";
			LogWriter::writeLog(str, 1);
			//break;
		}
	}
	return is_nan;
}
