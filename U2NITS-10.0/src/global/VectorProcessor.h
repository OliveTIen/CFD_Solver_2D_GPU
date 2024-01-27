#ifndef VECTOR_PROCESSOR_H
#define VECTOR_PROCESSOR_H

#include <vector>

class VectorProcessor {
public:
	// 拼接数组 会修改accepter
	static void appendToVector(std::vector<int>& accepter, const std::vector<int>& giver);

};


#endif