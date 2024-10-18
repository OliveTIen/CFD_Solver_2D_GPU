#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H
#include <string>
#include <vector>
/*
目前不想用类，只用静态函数
先只考虑一种材料。后续会用vector
[注意]后面修改为GPU程序时，GPU传参时，SDevicePara要传最新的mu
*/

class MaterialManager {
private:
	struct FluidMaterial {
	public:
		std::string name = "air";
		double mu0 = 17.9e-6; // 动力粘性系数参考值
		FluidMaterial(std::string _name, double _mu0_dynamic_viscosity) {
			name = _name;
			mu0 = _mu0_dynamic_viscosity;
		}
		// 给定参考密度rho、参考速度u、参考长度L、雷诺数Re，计算mu0并存储
		FluidMaterial(std::string _name, double rho, double u, double L, double Re) {
			name = _name;
			mu0 = rho * u * L / Re;
		}
		FluidMaterial() {};
		//FluidMaterial(FluidMaterial&) = delete;
		//FluidMaterial(FluidMaterial&&) = delete;
	};

	static MaterialManager* p_instance;
	std::vector<FluidMaterial> m_materials;

	MaterialManager() {};

public:
	static MaterialManager* getInstance();
	void addMaterial(FluidMaterial material);
	// 返回常引用，不可作为左值，不可被修改
	const FluidMaterial& getMaterial(size_t index = 0);
	void initialize_using_config(double rho_ref, double U_ref, double L_ref);
};

#endif // !MATERIAL_MANAGER_H
