#include "MaterialManager.h"
#include "../output/LogWriter.h"
#include "TomlFileManager.h"
#include "../global/GlobalPara.h"
#include "../global/CExit.h"

MaterialManager* MaterialManager::p_instance = nullptr;

MaterialManager* MaterialManager::getInstance() {
    if (p_instance == nullptr) {
        p_instance = new MaterialManager();
    }
    return p_instance;
}

void MaterialManager::addMaterial(FluidMaterial material) {
    if (m_materials.empty()) {
        m_materials.push_back(material);
    }
    else {
        // 目前只允许添加一种材料，不允许重复添加
        LogWriter::logAndPrintError("can only add one material\n");
        exit(-1);
    }
    
}


const MaterialManager::FluidMaterial& MaterialManager::getMaterial(size_t index) {
    if (m_materials.empty()) {
        LogWriter::logAndPrintError("MaterialManager.m_materials.empty()\n");
        exit(-1);
    }
    size_t size = m_materials.size();
    if (index >= 0 && index < size) {
        return m_materials[index];
    }
    else {
        LogWriter::logAndPrintError("MaterialManager.m_materials.out_of_range\n");
        exit(-1);
    }
}

void MaterialManager::initialize_using_config(double rho_ref, double U_ref, double L_ref) {
    /*
    放在FieldInitializer中，因为要等到确定初始化方式后才知道用哪个u作为参考速度
    功能：
    读取toml文件，先判断是否用Re修正mu0，然后添加Material
    
    读取策略
    若模式为：用Re修正粘度系数
    */
    
    TomlFileManager* t = TomlFileManager::getInstance();
    int calculate_mu0_by_Re = 0;
    double Re = 1e6;
    
    t->getValueIfExists("constant.calculate_mu0_by_Re", calculate_mu0_by_Re);
    if (calculate_mu0_by_Re) {
        t->getValue("constant.Re", Re);
        addMaterial(MaterialManager::FluidMaterial("air", rho_ref, U_ref, L_ref, Re));
    }
    else {
        double mu0 = 17.9e-6;
        t->getValueIfExists("constant.mu0", mu0);
        addMaterial({ "air", mu0 });
    }
}

