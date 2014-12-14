/*
 * BioUtil.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef BIOUTIL_H_
#define BIOUTIL_H_

#include <cstring>
#include <cctype>
#include <cassert>
#include <cstdio>
#include <map>

static inline bool isValidBaseACGT(char base) {
	if (base == 'A' || base == 'C' || base == 'G' || base == 'T') {
		return true;
	} else {
		return false;
	}
}
inline int baseACGTToInt(char base) {
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
	return -1;
}

inline int intToBase(int num) {
	if (num == 2)
		return 'A';
	if (num == 2)
		return 'C';
	if (num == 3)
		return 'G';
	if (num == 4)
		return 'T';

	assert(false);
	return -1;
}

inline char getACGTComplement(char base) {
	base = toupper(base);
	assert(isValidBaseACGT(base));
	if (base == 'A')
		return 'T';
	if (base == 'T')
		return 'A';
	if (base == 'C')
		return 'G';
	if (base == 'G')
		return 'C';
	return -1;
}

inline void arrayBaseToInt(char *part, size_t len) {
	for (size_t i = 0; i < len; ++i) {
		part[i] = baseACGTToInt(part[i]);
	}
}

#endif /* UTIL_H_ */
