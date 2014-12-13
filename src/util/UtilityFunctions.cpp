#include "Fenwick.h"
#include <cassert>
#include <cstdio>
#include <algorithm>
#include <utility>
#include <vector>

#include "UtilityFunctions.h"

namespace {

}

template<typename type>
void freeVectorOfPtrs(std::vector<type>& vec) {
	typename std::vector<type>::iterator myIter = vec.begin();
	while (myIter != vec.end()) {
		if (*myIter) {
			delete (*myIter);
		}
		++myIter;
	}
	vec.clear();
}


