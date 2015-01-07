/*
 * bioutil.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 *
 *  Header contains useful functions for working with DNA sequences.
 */

#ifndef BIOUTIL_H_
#define BIOUTIL_H_

#include <cstring>
#include <cctype>
#include <cassert>
#include <cstdio>
#include <map>

/**
 * Function tests if base passed as argument is valid.
 * Valid bases are A, C, G, T and N.
 *
 * @param base base that is tested
 * @return true if base is valid, false otherwise
 */
static inline bool isValidBaseACGT(char base) {
	base = toupper(base);
	return base == 'A' || base == 'C' || base == 'G' || base == 'T' || base == 'N';

}
/**
 * Function that calculates complement of given base.
 * Valid bases are A, C, G, T and N.
 *
 * @param base base whose complement is being calculated
 * @return complement of given base
 */
inline char getACGTComplement(char base) {
	base = toupper(base);
	assert(isValidBaseACGT(base));

	if (base == 'A')
		return 'T';
	else if (base == 'T')
		return 'A';
	else if (base == 'C')
		return 'G';
	else if (base == 'G')
		return 'C';
	else
		return 'N';

}

#endif /* UTIL_H_ */
