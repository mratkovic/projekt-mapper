/*
 * lcsk.h
 *
 *  Created on: Jan 6, 2015
 *      Author: marko
 *  Static class that defines method for calculation of LCSk between two strings.
 *  LCSk is number of common substrings of length k between two given strings.
 */

#ifndef SRC_LCSK_H_
#define SRC_LCSK_H_

#include <vector>
#include <stdint.h>

namespace bioutil {

class LCSk {

public:
	/**
	 * Method that is used for calculation of LCSk between two given strings with already
	 * filled vector of positions of every common substring of length k in those two strings.
	 *
	 * @param k parameter k of LCSk
	 * @result vector for keeping reconstruction of LCSk. It contains pairs of integers that
	 * correspond with positions in strings.
	 * @elements vector filled with pairs (i, j) that represent starting positions of common substrings
	 * of length k. i is starting position in first string, j starting position in second string
	 * @return value of LCSk between two strings
	 */
	static uint32_t calcLCSk(uint32_t k, std::vector<std::pair<uint32_t, uint32_t> >* result,
			std::vector<std::pair<uint32_t, uint32_t> >& elements);

};

} /* namespace bioutil */

#endif /* SRC_LCSK_H_ */
