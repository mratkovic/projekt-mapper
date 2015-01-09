/*
 * lis.h
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 */

#ifndef SRC_LIS_H_
#define SRC_LIS_H_

#include <vector>
#include <stdint.h>

namespace bioutil {

class Lis {
public:
	static void calcLIS(std::vector<int>* result, const std::vector<std::pair<uint32_t, uint32_t> >& elements);
	static int estimateBeginingPosFromLIS(std::vector<std::pair<uint32_t, uint32_t> >& positions,
			std::vector<int>& lis);
	static void reconstructLIS(std::vector<int>* result, int lastPos, int32_t* dpPath, int n);
};

} /* namespace bioutil */

#endif /* SRC_LIS_H_ */