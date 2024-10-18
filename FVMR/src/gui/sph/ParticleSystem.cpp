#include <glad/glad.h> // 和glad.c对应
#include <GLFW/glfw3.h>// 和glfw3.lib对应
#include <glm/gtc/type_ptr.hpp>
#include <unordered_set>
#include "ParticleSystem.h"


void ParticleSystem::setParameters() {
	auto& para = physicalParameters;
	para.m_standard_distance = 0.05;
	para.smooth_length = 0.08;


	para.bounce1_damping_coefficient = 0.224;
	para.m_standard_distance = 0.295;
	para.smooth_length = 0.110;
	para.gravity_direction = glm::normalize(glm::vec3(-0.1, -1, 0.1));
	
}

void ParticleSystem::resetParticles() {

	auto& para = physicalParameters;


	float space = para.m_standard_distance/2;
	//float space = 0.1;
	glm::vec3 initialCubeSize = initialCube[1] - initialCube[0];
	initialNum = initialCubeSize / space;
	vertices.clear();
	for (float z = initialCube[1].z; z > initialCube[0].z; z -= space)
		for (float x = initialCube[0].x; x < initialCube[1].x; x += space)
			for (float y = initialCube[0].y; y < initialCube[1].y; y += space) {
				vertices.push_back({ x,y,z });
			}

	velocities.resize(vertices.size(), { 0,0,0 });
	accelerations.resize(vertices.size(), { 0,0,0 });
	densities.resize(vertices.size());
	pressures.resize(vertices.size());
	particleNeighbors.resize(vertices.size());


	return;

	// resize vectors
	size_t nx = initialNum.x;
	size_t ny = initialNum.y;
	size_t nz = initialNum.z;
	vertices.resize(nx * ny * nz);

	velocities.resize(vertices.size());
	accelerations.resize(vertices.size());
	densities.resize(vertices.size());
	pressures.resize(vertices.size());
	particleNeighbors.resize(vertices.size());

	// set particle positions and velocities
	auto getDx = [](size_t nx, float x)->float {
		// if nx<=1, then return 0, else return dx
		if (x < 0)x = -x;
		return (nx <= 1) ? (0.0f) : (x / (nx - 1));
	};
	float dx = getDx(nx, initialCube[1].x - initialCube[0].x);
	float dy = getDx(ny, initialCube[1].y - initialCube[0].y);
	float dz = getDx(nz, initialCube[1].z - initialCube[0].z);
	for (int i = 0; i < nx; i++) {
		for (int j = 0; j < ny; j++) {
			for (int k = 0; k < nz; k++) {
				float rand_f = (float(rand() % 1000) / 1000.0f) * 0.1 * physicalParameters.m_standard_distance;
				int particle = i + j * nx + k * nx * ny;
				vertices[particle] = glm::vec3(i * dx + rand_f, j * dy, k * dz + rand_f) + initialCube[0];
				velocities[particle] = glm::vec3(0, 0, 0);
			}
		}
	}
	// 
	//para.normalization_density = getDensity(para.particle_mass, para.smooth_length);
	//para.normalization_pressure_force = getPressureForce(para.particle_mass, para.smooth_length);
	//para.normalization_viscous_force = getViscousForce(para.particle_mass, para.dynamic_viscosity, para.smooth_length);

	para.m_rest_density = para.particle_mass / pow(para.m_standard_distance, 3);
}

int ParticleSystem::find_biggest_pow_of_2(int n) {
	// find the biggest power of 2 no more than n
	/*
	algorithm: https://blog.csdn.net/dreamispossible/article/details/91162847
	1. for example, n=19, that is 0b00010011
	2. the result must have the format of 0b.10...0
	3. so 0b00010011 -> 0b00011111 -> 0b00010000
	*/
	n |= n >> 1;// 0b00011011
	n |= n >> 2;// 0b00011111
	n |= n >> 4;
	n |= n >> 8;
	n |= n >> 16; 
	return (n + 1) >> 1;// 0b00010000
}

void ParticleSystem::initializeBgGrid() {
	auto& elemNum = bgGridElemNum;
	// sort cube points
	auto Sort = [](float& x, float& y)->void {
		// swap x,y if x>y
		if (x > y) {
			float tmp = y;
			y = x;
			x = tmp;
		}
	};
	for (int i = 0; i < 3; i++) {
		Sort(bgGridCube[0][i], bgGridCube[1][i]);
	}

	// get background grid element number
	glm::vec3 cubeSize = bgGridCube[1] - bgGridCube[0];
	float refSize = physicalParameters.smooth_length * 2;// diameter of smooth sphere
	size_t particleNum = vertices.size();
	int root3_particleNum = (int)powf(float(particleNum), 1.0 / 3.0);
	for (int i = 0; i < 3; i++) {
		int refNum = find_biggest_pow_of_2(int(cubeSize[i] / refSize));
		// O(n*n) -> O(m*n), so m should <= n
		if (refNum > root3_particleNum)refNum = root3_particleNum;
		if (refNum < 1)refNum = 1;
		elemNum[i] = refNum;
	}
	
	// get background grid element size
	for (size_t i = 0; i < 3; i++) {
		bgGridElemSize[i] = cubeSize[i] / elemNum[i];
	}

	particleGridIds.resize(vertices.size());
	gridParticleIds.resize(elemNum.x * elemNum.y * elemNum.z);
	updateParticleGridIds();
}

void ParticleSystem::updateParticleGridIds() {
	// clear
	for (auto& grid : gridParticleIds) {
		grid.clear();// set size to zero
	}

	// update grid ids according to position
	for (size_t i = 0; i < vertices.size(); i++) {
		const auto& position = vertices[i];
		auto relative_position = position - bgGridCube[0];
		auto& gridId = particleGridIds[i];
		const auto& elemNum = bgGridElemNum;
		// update particleGridIds
		for (int j = 0; j < 3; j++) {
			int id = int(relative_position[j] / bgGridElemSize[j]);
			if (id < 0)id = 0;
			if (id >= elemNum[j])id = elemNum[j];
			gridId[j] = id;
		}
		// update gridParticleIds
		int gridId_sequence = gridId.x + gridId.y * elemNum.x
			+ gridId.z * elemNum.x * elemNum.y;
		gridParticleIds[gridId_sequence].push_back(i);
	}
}

std::vector<std::pair<int, float>> ParticleSystem::searchNeighborsInRadius(int particleId, float radius) {
	if (particleId < 0 || particleId >= vertices.size())throw "particleId out of range";

	const auto& elemNum = bgGridElemNum;
	const auto& currentGrid = particleGridIds[particleId];
	// expand searching range to 9 grids
	std::unordered_set<int> particleToSearch;
	for (int i = -1; i <= 1; i++) {
		for (int j = -1; j <= 1; j++) {
			for (int k = -1; k <= 1; k++) {
				glm::ivec3 grid = glm::ivec3(i, j, k);
				int gridId_sequence = grid.x + grid.y * elemNum.x
					+ grid.z * elemNum.x * elemNum.y;
				if (gridId_sequence < 0 || gridId_sequence >= gridParticleIds.size())continue;// filter grids to avoid subscript out of range
				const auto& particles_in_grid = gridParticleIds[gridId_sequence];
				for (const auto& p : particles_in_grid) {
					particleToSearch.insert(p);
				}
			}
		}
	}

	// filter neighbors in candidates using radius
	std::vector<std::pair<int, float>> neighbors;
	for (const auto& particle : particleToSearch) {
		if (particle < 0 || particle >= vertices.size())throw "particle out of range";
		float distance = glm::distance(vertices[particle], vertices[particleId]);
		if (distance < radius)neighbors.push_back(std::pair<int, float>(particle, distance));
	}
	// return value optimization: https://blog.csdn.net/wangkai6666/article/details/135183961
	return neighbors;// 编译器会自动优化，避免拷贝复制，此处不用管
}

bool ParticleSystem::in_bgGridCube(const glm::vec3& point, float ball_radius) {
	const auto& p = point;
	glm::vec3 p0 = bgGridCube[0];
	glm::vec3 p1 = bgGridCube[1];
	p0 += ball_radius;
	p1 -= ball_radius;
	bool condition = p.x > p0.x && p.x < p1.x
		&& p.y > p0.y && p.y < p1.y
		&& p.z > p0.z && p.z < p1.z;
	if (condition)return true;
	return false;
}

void ParticleSystem::in_bgGridCube(const glm::vec3& point, float ball_radius, int& flag_x, int& flag_y, int& flag_z) {
	const auto& p = point;
	glm::vec3 p0 = bgGridCube[0];
	glm::vec3 p1 = bgGridCube[1];
	p0 += ball_radius;
	p1 -= ball_radius;
	auto getFlag = [](float px, float p0x, float p1x)->int {
		if (px > p1x)return 1;
		else if (px < p0x)return -1;
		else return 0;
	};
	flag_x = getFlag(p.x, p0.x, p1.x);
	flag_y = getFlag(p.y, p0.y, p1.y);
	flag_z = getFlag(p.z, p0.z, p1.z);
}

void ParticleSystem::setBoundaryCondition_Wall(int particle, float ball_radius, float dt) {




	auto& p = vertices[particle];
	auto& v = velocities[particle];
	auto& a = accelerations[particle];
	glm::vec3 p0 = bgGridCube[0];
	glm::vec3 p1 = bgGridCube[1];
	p0 += ball_radius;
	p1 -= ball_radius;
	
	int bounceType = physicalParameters.bounceType;
	float bounce0_spring = physicalParameters.bounce0_spring;
	float bounce0_norm_drag = physicalParameters.bounce0_norm_drag;
	float bounce0_wall_friction = physicalParameters.bounce0_wall_friction;
	float epsilon = 0.5 * ball_radius;

	if (bounceType == 0) {
		// apply bouncy acceleration, using linear spring
		glm::vec3 acc{0, 0, 0};
		auto addBouncyForce = [bounce0_spring,epsilon,&v, bounce0_norm_drag, bounce0_wall_friction](glm::vec3 p_p0, glm::vec3 N, glm::vec3& acc)->void {
			float a_max = 1e5f;
			float phi = glm::dot(p_p0, N);// projection of p_p0 onto N. For example, for plane Z+, phi = p.z - p0.z
			if (phi < epsilon) {
				acc += bounce0_spring * (epsilon - phi) * N;
				// apply drag in norm
				glm::vec3 vn = glm::dot(v, N) * N;
				glm::vec3 vt = v - vn;
				glm::vec3 drag_acc = -vn * bounce0_norm_drag;
				acc += drag_acc;
				// apply friction in shear
				if (glm::length(vt) > 0.001) {// glm::length: sqrt(dot(v, v))
					acc += -vt * bounce0_wall_friction;
				}
			}
		};
		addBouncyForce(p - p0, { 0, 0, 1 }, acc);
		addBouncyForce(p - p0, { 0, 1, 0 }, acc);
		addBouncyForce(p - p0, { 1, 0, 0 }, acc);
		addBouncyForce(p - p1, { 0, 0,-1 }, acc);
		addBouncyForce(p - p1, { 0,-1, 0 }, acc);
		addBouncyForce(p - p1, {-1, 0, 0 }, acc);
		a += acc;
	}
	else {// bounce 1 
		// limit bounce_damping range
		float bounce_damping = -physicalParameters.bounce1_damping_coefficient;
		if (bounce_damping < -1)bounce_damping = -1;
		if (bounce_damping > 0)bounce_damping = 0;
		glm::vec3 v_old = v;
		auto setBCWall = [bounce_damping](float& px, float& vx, float& ax, float p0x, float p1x)->int {
			int flag = 0;
			if (px == p1x || px == p0x) {
				vx = 0;
			}
			if (px > p1x) {
				px = p1x;
				vx *= bounce_damping;
				//ax = 0;
				flag = 1;
			}
			else if (px < p0x) {
				px = p0x;
				vx *= bounce_damping;
				//ax = 0;
				flag = -1;
			}
			return flag;
		};
		for (int i = 0; i < 3; i++) {
			//float v2_old = glm::dot(v[i], v[i]);
			int flag = setBCWall(p[i], v[i], a[i], p0[i], p1[i]);
			if (flag != 0) {

			}
		}

		// convert normal speed to shear speed
		//未完待续;
	}

}



float ParticleSystem::getDensity(float particleMass, float smoothingLength) {
	return ((315 * particleMass) / (64 * PI * powf(smoothingLength, 9.0f)));
}

float ParticleSystem::getKernelPoly6(float h) {
	return 315 / (64 * PI * powf(h, 9.0f));
}

float ParticleSystem::getWPoly6(float h, float r) {
	// equation ref: https://www.thecodeway.com/?p=161
	if (r > h)return 0;
	if (r < 0)return 0;
	float K = getKernelPoly6(h);
	float sd = h * h - r * r;// square difference
	return K * sd * sd * sd;
}

float ParticleSystem::getPressureForce(float particleMass, float smoothingLength) {
	return (-(45 * particleMass) / (PI * powf(smoothingLength, 6.0f)));
}

float ParticleSystem::getViscousForce(float particleMass, float dynamicViscosity, float smoothingLength) {
	return ((45 * dynamicViscosity * particleMass) / (PI * powf(smoothingLength, 6.0f)));
}


void ParticleSystem::init() {
	if (buffer_allocated) {
		printf("buffer already allocated.\n");
		return;
	}
	
	setParameters();
	resetParticles();
	initializeBgGrid();

	// rendering
	shader = new ShaderParticle("./shader/particle.vert", "./shader/particle.frag");
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	glGenBuffers(1, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), glm::value_ptr(vertices[0]), GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	buffer_allocated = true;
	current_cursor_position = GUI::ConsoleProxy::getCursorPosition();
}

void ParticleSystem::update(float dt) {
	if (vertices.empty())return;
	//if (!enable_rendering)return;
	
	{
		// update dt queue (push new dt to dt_queue[0])
		for (size_t i = dt_queue.size() - 1; i >= 1; i--) {
			dt_queue[i] = dt_queue[i - 1];
		}
		dt_queue[0] = dt;
		// remove max dt, then get average
		int max_index = 0;
		int size = dt_queue.size();
		float max_dt = dt_queue[0];
		for (int i = 1; i < size; i++) {
			if (dt_queue[i] > max_dt) {
				max_index = i;
				max_dt = dt_queue[i];
			}
		}
		float dt_sum = 0;
		int num = 0;
		for (int i = 0; i < size; i++) {
			if (i == max_index)continue;
			if (dt_queue[i] <= 0)continue;
			dt_sum += dt_queue[i];
			num += 1;
		}
		if (num == 0)num = 1;
		dt = dt_sum /= float(num);
	}


	const auto& para = physicalParameters;
	if (para.printOnConsole) {
		GUI::ConsoleProxy::setCursorPosition(current_cursor_position);
		GUI::ConsoleProxy::showCursor(false);
		std::cout << "vertices.size:" << vertices.size() << "\n";
		int printLineNum = 10;
		for (int i = 0; i < vertices.size() && i < printLineNum; i++) {
			printf("vertex %2d: ", i);
			printf("x=%-+5.2f, y=%-+10.2f, z=%-+5.2f", vertices[i].x, vertices[i].y, vertices[i].z);
			printf(", ");
			printf("u=%-+7.2f, v=%-+10.2f, w=%-+7.2f", velocities[i].x, velocities[i].y, velocities[i].z);
			printf("\n");
			printf("p=%-+7.2f", densities[i]);
			printf("\n");
		}
	}

	for (int particle = 0; particle < vertices.size(); particle++) {
		const int& i = particle;
		if (para.useBgGrid) {
			particleNeighbors[i] = searchNeighborsInRadius(i, para.smooth_length);
		}
		else {
			std::vector<std::pair<int, float>> pN;
			for (int j = 0; j < vertices.size(); j++) {
				float d = glm::distance(vertices[i], vertices[j]);
				pN.push_back(std::pair<int, float>(j, d));
			}
			particleNeighbors[i] = pN;
		}
	}
	
	const float& h = para.smooth_length;// smooth length
	const float& m = para.particle_mass;// particle mass
	float rigidCollision_radius = para.rigidCollision_radius;
	float h2 = h * h;
	float h3 = h2 * h;
	float h6 = h3 * h3;
	float h9 = h6 * h3;
	float W_poly = m * 315.0 / (64.0 * PI * h9);// coefficient;  = m * 315.0 / (64.0 * PI * h9)
	//Spiky Kernel
	float W_kernelSpiky = -45.0f / (PI * h6);
	//Viscosity Kernel
	float W_kernelViscosity = 45.0f / (PI * h6);

	// pressure
	for (int i = 0; i < vertices.size(); i++) {
		pressures[i] = densities[i] * 0.01;// rho*R*T
		float sum = 0.0f;
		for (const auto& neighbor : particleNeighbors[i]) {
			int j = neighbor.first;
			float r = neighbor.second;// distance
			float r2 = r * r;
			float h2_r2 = h2 - r2;
			if (j == i) {
				// do nothing
			}
			if (h2_r2 < 0 && para.m_clip_smooth_function)h2_r2 = 0;
			sum += pow(h2_r2, 3.f);
		}
		densities[i] = W_poly * m * sum;
		pressures[i] = (densities[i] - para.m_rest_density) * para.m_gasConstantK;
	}
	// clear force
	for (int i = 0; i < vertices.size(); i++) {
		accelerations[i] = glm::vec3(0, 0, 0);
	}
	// rigid collision
	if (para.applyRigidCollision) {
		for (int i = 0; i < vertices.size(); i++) {
			glm::vec3 aci{0, 0, 0};
			for (const auto& neighbor : particleNeighbors[i]) {
				int j = neighbor.first;
				if (j == i)continue;
				glm::vec3 rj_ri = vertices[j] - vertices[i];
				float r = neighbor.second;// distance
				glm::vec3 nij{0,0,0};// norm, from i to j
				if (r > 0.01 * rigidCollision_radius) {
					nij = glm::normalize(rj_ri);
				}
				float phi = r - rigidCollision_radius;
				if (phi < 0) {
					aci += -para.rigidCollision_spring * (0 - phi) * nij;
				}
			}
			accelerations[i] += aci;
		}
	}

	// pressure force
	const float C_pa = m * 45 / (PI * h6) * 1e1;// m * 45 / (PI * h6) * ...
	for (int i = 0; i < vertices.size(); i++) {
		if (!para.applyPressureForce)break;
		glm::vec3 api = glm::vec3(0, 0, 0);
		for (const auto& neighbor : particleNeighbors[i]) {
			int j = neighbor.first;
			if (j == i)continue;
			float r = neighbor.second;// distance
			float h_r_2 = (h - r) * (h - r);
			if (h < r && para.m_clip_smooth_function)h_r_2 = 0;

			//const float K_Spiky = 15 / (PI * h6);
			//float W_Spiky = K_Spiky * powf(h - r, 3);
			//if (r > h)W_Spiky = 0;

			const auto& pi = pressures[i];
			const auto& pj = pressures[j];
			const auto& rhoi = densities[i];
			const auto& rhoj = densities[j];
			const auto& ri = vertices[i];
			const auto& rj = vertices[j];
			float A = (pi + pj) / (2 * rhoi * rhoj);
			if (r < para.m_epsilon)r = para.m_epsilon;
			//api += C_pa * A * B;

			float pterm = -m * W_kernelSpiky * h_r_2 * A;
			api += (ri - rj) * pterm / r;


			float limitValue = 1.0e3;
			for (int k = 0; k < 3; k++) {
				if (api[k] > limitValue)api[k] = limitValue;
				if (api[k] < -limitValue)api[k] = -limitValue;
			}

			if (isnan(api[0])) {
				std::cout << "ap_i isnan\n";
			}


		}
		accelerations[i] += api;
	}
	// viscous force
	for (int i = 0; i < vertices.size(); i++) {
		if (!para.applyViscousForce)break;

		glm::vec3 avi = glm::vec3(0, 0, 0);
		for (const auto& neighbor : particleNeighbors[i]) {
			int j = neighbor.first;
			if (j == i)continue;
			float r = neighbor.second;// distance

			const auto& ui = velocities[i];
			const auto& uj = velocities[j];
			const auto& rhoi = densities[i];
			const auto& rhoj = densities[j];

			float rhoirhoj = rhoi * rhoj;
			if (rhoirhoj < para.m_epsilon)rhoirhoj = para.m_epsilon;
			glm::vec3 A = (uj - ui) / rhoirhoj;
			float h_r = h - r;
			if (h_r < 0)h_r = 0;
			//avi += mu * C_pa * A * h_r;
			avi += (uj - ui) * W_kernelViscosity * para.m_viscosity * h_r * m / rhoirhoj;
		}

		if (isnan(avi[0])) {
			std::cout << "isnan\n";
		}
		accelerations[i] += avi;
	}
	// other force
	for (int i = 0; i < vertices.size(); i++){
		// gravity
		auto gravity_force = para.gravity * para.gravity_direction;
		accelerations[i] += gravity_force;

		// drag
		accelerations[i] += -velocities[i] * physicalParameters.air_drag;

		// boundary conditions
		setBoundaryCondition_Wall(i, para.sphere_radius, dt);


		//float frac = glm::length(ap) / glm::length(av);
		//printf("%2d: ", i);
		//printf("ap.x=%f, ap.y=%f, ap.z=%f, ap/av=%f, p=%f, rho=%f,\n", ap.x, ap.y, ap.z, frac, pressures[i], densities[i]);

	}
	// time advance
	for (int i = 0; i < vertices.size(); i++) {
		velocities[i] += dt * accelerations[i];
		vertices[i] += dt * velocities[i];
	}


	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), glm::value_ptr(vertices[0]), GL_DYNAMIC_DRAW);

	return;
	/*
	old 
	*/

	//// calculate
	//for (int particle = 0; particle < vertices.size(); particle++) {
	//	// search neighbors
	//	particleNeighbors[particle] = searchNeighborsInRadius(particle, para.smooth_length);

	//	const auto& neighbors = particleNeighbors[particle];
	//	for (const auto& neighbor : neighbors) {
	//		int neighborId = neighbor.first;
	//		float distance = neighbor.second;
	//		float dL2 = smooth_length_square - distance * distance;
	//		densities[particle] += para.normalization_density * dL2 * dL2 * dL2;
	//		// pressures = ISOTROPIC_EXPONENT * (densities - BASE_DENSITY)
	//		pressures[particle] = para.isotropic_exponent * (densities[particle] - para.base_density);
	//		forces[particle] = glm::vec3(0.0f, 0.0f, 0.0f);
	//	}
	//}

	//// forces
	//for (int particle = 0; particle < vertices.size(); particle++) {
	//	const auto& i = particle;
	//	auto& v = velocities;
	//	const auto& rho = densities;
	//	const auto& neighbors = particleNeighbors[particle];
	//	// force from neighbors
	//	for (const auto& neighbor : neighbors) {
	//		int neighborId = neighbor.first;
	//		float distance = neighbor.second;
	//		if (neighborId == particle)continue;// exclude the particle itself

	//		const auto& p = pressures;
	//		const auto& j = neighborId;
	//		const float pn = para.normalization_pressure_force;

	//		glm::vec3 dX = vertices[j] - vertices[i];
	//		float dL = para.smooth_length - distance;
	//		// pressure force
	//		forces[i] += pn * (-dX) / distance * (p[j] + p[i]) / (2.0f * rho[j]) * dL * dL;
	//		// viscous force
	//		forces[i] += para.normalization_viscous_force * ((v[j] - v[i]) / rho[j] * (para.smooth_length - distance));
	//	}
	//	// gravity force
	//	forces[i] += para.gravity_force * rho[i];

	//	// boundary conditions
	//	setBoundaryCondition_Wall(particle, para.sphere_radius, dt, forces[i]);

	//	// Euler Step
	//	float dt = para.time_step_length;
	//	v[i] += dt * forces[i] / rho[i];
	//	vertices[i] += dt * v[i];

	//	
	//}
	//

	//glBindBuffer(GL_ARRAY_BUFFER, vbo);
	//glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), glm::value_ptr(vertices[0]), GL_DYNAMIC_DRAW);

	
}

void ParticleSystem::render(const glm::mat4& projection, const glm::mat4& view) const {
	if (vertices.empty())return;
	if (!enable_rendering)return;

	glEnable(GL_PROGRAM_POINT_SIZE);
	glBindVertexArray(vao);
	shader->use();
	shader->setMat4("projection", projection);
	shader->setMat4("view", view);
	shader->setMat4("model", matrix_model);
	shader->setVec4("userColor", userColor);
	shader->setFloat("pointRadius", pointRadius_render);
	shader->setFloat("refDistance", refDistance);
	//glPointSize(10);

	glDrawArrays(GL_POINTS, 0, vertices.size());
	glBindVertexArray(0);
	glDisable(GL_PROGRAM_POINT_SIZE);
}

void ParticleSystem::cleanup() {
	if (!buffer_allocated) {
		printf("buffer not allocated, cannot delete buffer\n");
		return;
	}

	glDeleteVertexArrays(1, &vao);
	glDeleteBuffers(1, &vbo);
	delete shader;
}




