/*
 * BioUtil.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef BIOUTIL_H_
#define BIOUTIL_H_

#include <cstring>
#include <cassert>

static inline bool isValidBaseACGT(char base) {
	if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
		return true;
	} else {
		return false;
	}
}
inline int baseToInt(char base) {
	base = toupper(base);
	assert(isValidBaseACGT(base));
	if (base == 'A')
		return 1;
	if (base == 'C')
		return 2;
	if (base == 'G')
		return 3;
	if (base == 'T')
		return 4;
	// ???
	assert(false);
	return -1;
}

inline void arrayBaseToInt(char *part, size_t len) {
	for(size_t i = 0; i < len; ++i) {
		part[i] = baseToInt(part[i]);
	}
}

#endif /* UTIL_H_ */
