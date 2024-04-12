#include "Element_2D.h"
#include "FVM_2D.h"
#include "Edge_2D.h"
#include "output/LogWriter.h"

//用函数对静态成员初始化，需放进局部namespace
namespace {
    Eigen::MatrixXi iniSiti_Q4() {
        Eigen::MatrixXi si_ti(2, 4);//(维数, 顶点数)
        si_ti <<
            -1, +1, +1, -1,//si
            -1, -1, +1, +1;//ti
        return si_ti;
    }

    Eigen::MatrixXd iniGaussPoint_Q4() {
        //高斯积分点 1/sqrt(3)=0.577350269189625(7,8), double取0.577350269189626
        Eigen::MatrixXd GaussPointMatrix(8, 3);//(顶点数, 维数)
        double ops3 = 0.577350269189626;//onePerSqrt3
        GaussPointMatrix <<
            -ops3, -ops3,
            ops3, -ops3,
            ops3, ops3,
            -ops3, ops3;
        return GaussPointMatrix;
    }
}

//Eigen::MatrixXi Element_Q4::si_ti = iniSiti_Q4();
//Eigen::MatrixXd Element_Q4::GaussPointMatrix = iniGaussPoint_Q4();

double Element_2D::calArea(FVM_2D* f) {
    // 叉乘法计算三角形面积
    double xn[3]{}, yn[3]{};
    for (int i = 0; i < 3; i++) {
        xn[i] = f->getNodeByID(nodes[i])->x;
        yn[i] = f->getNodeByID(nodes[i])->y;
    }
    return 0.5 * abs(xn[0] * (yn[1] - yn[2]) + xn[1] * (yn[2] - yn[0]) + xn[2] * (yn[0] - yn[1]));

    //// 但事实上总有单元其节点顺序是逆时针，这是不可避免的
    //double area = 0.5 * (xn[0] * (yn[1] - yn[2]) + xn[1] * (yn[2] - yn[0]) + xn[2] * (yn[0] - yn[1]));
    //if (area <= 0) {
    //    // 要求节点顺序必须是逆时针，否则计算通量时会出错
    //    std::stringstream ss;
    //    ss << "Element area <= 0. Please check if the node order is counterclockwise.\n ";
    //    ss << "At element ID=" << ID << ", GPUID=" << GPUID << ", (x,y)=(" << x << "," << y << ")";
    //    LogWriter::logAndPrintError(ss.str());
    //    exit(-1);
    //}
    //return area;
}


std::vector<Element_2D*> Element_2D::findNeighbor() {
    //使用前需保证pEdges是最新的，且pEdges的pElement_L/R是最新的
    std::vector<Element_2D*> pNeighborElements(3);
    for (int iedge = 0; iedge < 3; iedge++) {
        if (pEdges[iedge] == nullptr) {
            std::cout << "Error: Element_2D::pEdges[" << iedge << "] is uninitialized(Element ID: " << this->ID << ")" << std::endl;
            pNeighborElements[iedge] = nullptr;
            //continue;//结束本次循环
        }
        else {
            if (pEdges[iedge]->pElement_L == this) {
                if (pEdges[iedge]->pElement_R == nullptr) {
                    pNeighborElements[iedge] = nullptr;
                }
                else {
                    pNeighborElements[iedge] = pEdges[iedge]->pElement_R;
                }
            }
            else if(pEdges[iedge]->pElement_R == this) {
                if (pEdges[iedge]->pElement_L == nullptr) {
                    pNeighborElements[iedge] = nullptr;
                }
                else {
                    pNeighborElements[iedge] = pEdges[iedge]->pElement_L;
                }
            }
        }
    }
    return pNeighborElements;
}

std::vector<Element_2D*> Element_2D::findNeighbor_withoutNullptr() {
    std::vector<Element_2D*>neighbors_origin = this->findNeighbor();
    std::vector<Element_2D*>neighbors;
    //清除nullptr
    for (int i = 0; i < 3; i++) {
        if (neighbors_origin[i] != nullptr) {
            //neighbors_origin[i]->calxy(f);//初始化坐标 calxy!!!
            neighbors.push_back(neighbors_origin[i]);
        }
    }
    return neighbors;
}

double Element_2D::calDistanceFromNearestNeighbor(FVM_2D* f) {
    std::vector<Element_2D*> n = findNeighbor();
    double dis = 1e10;
    double tmpx, tmpy, tmpdis;
    //this->calxy(f);读取文件时已经calxy了
    for (int in = 0; in < n.size(); in++) {
        if (n[in] != nullptr) {
            //n[in]->calxy(f);读取文件时已经calxy了
            tmpx = n[in]->x;
            tmpy = n[in]->y;
            tmpdis = sqrt((tmpx - x) * (tmpx - x) + (tmpy - y) * (tmpy - y));
            dis = (std::min)(dis, tmpdis);
        }
    }
    return dis;
}

void Element_2D::updateSlope_Barth(FVM_2D* f) {
    //计算UxUy并进行限制

    //3个邻居 超静定
    //   |x0-x y0-y|       |U0i-Ui|
    //   |x1-x y1-y| |Uxi|=|U1i-Ui|
    //   |x2-x y2-y| |Uyi| |U2i-Ui|
    //   
    
    //2个邻居
    //   
    //   |x0-x y0-y| |Uxi|=|U0i-Ui|
    //   |x1-x y1-y| |Uyi| |U1i-Ui|
    //   
    
    //1个邻居 不定 补充条件
    //   
    //   |x0-x y0-y| |Uxi|=|U0i-Ui|
    //               |Uyi|
    //   

    std::vector<Element_2D*>neighbors_origin = findNeighbor();
    std::vector<Element_2D*>neighbors;
    //清除nullptr
    for (int i = 0; i < 3; i++) {
        if (neighbors_origin[i] != nullptr) {
            //neighbors_origin[i]->calxy(f);//初始化坐标 calxy!!!
            neighbors.push_back(neighbors_origin[i]);
        }
    }
    //this->calxy(f);// 读取文件时已经calxy
    const int nNeighbor = (int)neighbors.size();

    for (int i_ = 0; i_ < 4; i_++) {
        Eigen::MatrixXd xy(nNeighbor, 2);
        Eigen::VectorXd Uxy(2);//待求量
        Eigen::VectorXd UU(nNeighbor);
        //initialize xy, UU
        for (int j = 0; j < nNeighbor; j++) {
            xy(j, 0) = neighbors[j]->x - this->x;//xj-x 已经calxy
            xy(j, 1) = neighbors[j]->y - this->y;//yj-y
            UU(j) = neighbors[j]->U[i_] - this->U[i_];//Uji-Ui
        }
        //solve
        if (nNeighbor == 3) {
            Uxy = xy.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(UU);//奇异值分解
            this->Ux[i_] = Uxy[0];
            this->Uy[i_] = Uxy[1];
        }
        else if (nNeighbor == 2) {
            Uxy = xy.householderQr().solve(UU);//解线性方程组
            this->Ux[i_] = Uxy[0];
            this->Uy[i_] = Uxy[1];
        }
        else {
            //该情况很罕见，除非是corner单元
            LogWriter::logAndPrint("未完成：Element_T3::calSlope_Barth(FVM_2D* f) nNeighbor==1 \n", LogWriter::Fatal);
            //指定梯度方向
        }

    }
    //Barth限制器
    restructor_in_updateSlope_Barth(f);
}

void Element_2D::restructor_in_updateSlope_Barth(FVM_2D* f) {
    //限制器，用来修正Ux,Uy。其理论依据是若三个顶点处的Ux,Uy满足有界性，则面内全满足有界性
    //计算偏差上下界
    std::vector<Element_2D*> neighbors = findNeighbor();
    double UU[3][4]{};//邻居函数值与自身函数值的差
    double UUup = 0;
    double UUdown = 0;
    for (int i = 0; i < neighbors.size(); i++) {
        if (neighbors[i] != nullptr) {
            for (int j = 0; j < 4; j++) {
                UU[i][j] = neighbors[i]->U[j] - this->U[j];
                UUup = (std::max)(UUup, UU[i][j]);
                UUdown = (std::min)(UUdown, UU[i][j]);
            }
        }
    }
    //计算顶点偏差
    double U_node[3][4]{};
    double UU_node[3][4]{};//顶点函数值与自身函数值的差
    double UUup_node = 0;
    double UUdown_node = 0;
    for (int i_node = 0; i_node < 3; i_node++) {
        Node_2D* pNode = f->pNodeTable[nodes[i_node]];
        get_U(pNode->x, pNode->y, U_node[i_node]);
        for (int j = 0; j < 4; j++) {
            UU_node[i_node][j] = U_node[i_node][j] - this->U[j];
            UUup_node = (std::max)(UUup_node, UU_node[i_node][j]);
            UUdown_node = (std::min)(UUdown_node, UU_node[i_node][j]);
        }
    }
    //修正Ux, Uy
    double ratio = 1;
    if (UUup_node > UUup) {
        ratio = (std::max)(ratio, UUup_node / UUup);
    }
    if (UUdown_node < UUdown) {
        ratio = (std::max)(ratio, UUdown_node / UUdown);
    }

    for (int j = 0; j < 4; j++) {
        Ux[j] /= ratio;
        Uy[j] /= ratio;
    }
}

void Element_2D::get_U(double xpoint, double ypoint, double* _U) {
    //常量重构
    if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_constant) {
        for (int j = 0; j < 4; j++) {
            _U[j] = U[j];//此处U为成员变量
        }
    }
    //线性重构
    //根据梯度Ux、Uy计算点(xpoint,ypoint)处_U值
    else if (GlobalPara::inviscid_flux_method::flag_reconstruct == _REC_linear) {
        for (int j = 0; j < 4; j++) {
            _U[j] = U[j] + Ux[j] * (xpoint - x) + Uy[j] * (ypoint - y);
        }
    }
}

Eigen::Vector4d Element_2D::get_U(double xpoint, double ypoint) {
    //初始化_U数组
    double _U[4];
    get_U(xpoint, ypoint, _U);//另一个重载函数
    //数组转向量
    Eigen::Vector4d _U_vector;
    for (int j = 0; j < 4; j++) {
        _U_vector[j] = _U[j];
    }
    return _U_vector;
}

std::vector<double> Element_2D::U2uv(const Eigen::Vector4d& Uc) {
    //守恒量ρ,ρu,ρv,ρE
    std::vector<double> uv;//非守恒量u v
    uv[0] = Uc[1]/Uc[0];
    uv[1] = Uc[2]/Uc[0];
    return uv;
}

void Element_2D::generateElementEdge(FVM_2D* f) {
    generateElementEdge_registerSingle(f, nodes[0], nodes[1], 0);
    generateElementEdge_registerSingle(f, nodes[1], nodes[2], 1);
    generateElementEdge_registerSingle(f, nodes[2], nodes[0], 2);
}

void Element_2D::generateElementEdge_registerSingle(FVM_2D* f, int ID_0, int ID_1, int iEdge) {
    //输入：待注册edge的两个node编号
    //操作：
    //1.添加edge或新增edge的Element_R信息，并修改edges和pEdgeTable
    //2.在Element.pEdges注册
    
    //遍历f_edges，若某元素的nodes与当前待确定nodes相同(不分先后)，说明已经注册了
    std::vector<Edge_2D>& f_edges = f->edges;//引用的本质是指针常量
    std::vector<Edge_2D*>& f_pEdgeTable = f->pEdgeTable;
    for (int i_edge = 0; i_edge < f_edges.size(); i_edge++) {
        Edge_2D& edge = f_edges[i_edge];//不必担心局部变量的反复创建释放的开销，因为编译器很聪明。
        //if ((ID_0 == edge.nodes[0] && ID_1 == edge.nodes[1]) || (ID_0 == edge.nodes[1] && ID_1 == edge.nodes[0])) {
        //已经注册过，因此只需补充element_R信息
        //按理说，若已经注册过，且都是逆时针旋转，则新注册的edge与已注册的edge必然是方向相反的
        if (ID_0 == edge.nodes[1] && ID_1 == edge.nodes[0]) {
            edge.pElement_R = this;
            //更新pEdges
            pEdges[iEdge] = &(f_edges[i_edge]);
            //直接终止
            return;
        }
    }
    //未搜寻到，说明是新注册的 L R H ID nodes
    Edge_2D new_edge;
    new_edge.pElement_L = this;
    new_edge.ID = (int)f_edges.size() + 1;//不严格按照inp文件中的ID来
    int edgeID = new_edge.ID;
    new_edge.nodes[0] = ID_0;
    new_edge.nodes[1] = ID_1;
    f_edges.push_back(new_edge);
    Edge_2D* pNewEdge = &(f_edges[f_edges.size() - 1]);
    //更新table
    if (f_pEdgeTable.size() < edgeID + 1) {
        f_pEdgeTable.resize(edgeID + 1);//第ID号，需要ID+1的空间
        f_pEdgeTable[edgeID] = pNewEdge;
    }
    //更新pEdges
    pEdges[iEdge] = pNewEdge;//取自&(f_edges[f_edges.size() - 1])，因此不会因局部变量而销毁
}


double Element_2D::calLambda(const double gamma) {
    double rho = U[0];
    double rho_u = U[1];//rho*u
    double rho_v = U[2];
    double rho_E = U[3];
    double E = rho_E / rho;
    double u = rho_u / rho;
    double v = rho_v / rho;
    double V2 = u * u + v * v;//|U|^2
    double p = 0.5 * rho * (gamma - 1) * (2 * E - V2);
    return sqrt(V2) + sqrt(gamma * p / rho);
}

double Element_2D::calLambdaFlux(FVM_2D* f) {
    double LambdaC = 0;
    // 对于每条边，计算其坐标，然后获取边上的LambdaC
    // 目前这段代码有问题
    for (int ie = 0; ie < 3; ie++) {
        // 边坐标、U
        double ex, ey;
        double eU[4]{};//边中点处ρ,ρu,ρv,ρE
        pEdges[ie]->getxy(f, ex, ey);
        get_U(ex, ey, eU);
        // 边法向量
        std::vector<double> en = pEdges[ie]->getDirectionN();
        //eabs
        double u = eU[1] / eU[0];
        double v = eU[2] / eU[0];
        double E = eU[3] / eU[0];
        double eabs = abs(u * en[0] + v * en[1]);//|euv·en|
        //ec
        double V2 = u * u + v * v;
        double& rho = eU[0];
        double p = 0.5 * rho * (GlobalPara::constant::gamma - 1) * (2 * E - V2);
        double ec = sqrt(GlobalPara::constant::gamma * p / rho);
        double dl = pEdges[ie]->getLength();
        LambdaC += (eabs + ec) * dl;
    }
    return LambdaC;
}
