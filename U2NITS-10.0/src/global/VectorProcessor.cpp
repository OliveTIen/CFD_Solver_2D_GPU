#include "VectorProcessor.h"

void VectorProcessor::appendToVector(std::vector<int>& accepter, const std::vector<int>& giver) {
	// ƴ������
	for (int i = 0; i < giver.size(); i++) {
		accepter.push_back(giver[i]);
	}
}
