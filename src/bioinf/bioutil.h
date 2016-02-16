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
#include <cstdint>
#include <climits>
#include <map>

namespace bioinf {

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
 * Function transforms base to small int number.
 *	'0' - 'A'
 *	'1' - 'C'
 *	'2' - 'G'
 *	'3' - 'T'
 *	'4' - 'N'
 *
 * @param base base
 * @return number assigned to base
 */
static inline uint8_t baseToInt(char base) {
  base = toupper(base);
  switch (base) {
    case 'A':
      return 0;
    case 'C':
      return 1;
    case 'G':
      return 2;
    case 'T':
      return 3;
    case 'N':
      return 4;
    default:
      return UCHAR_MAX;
  }

}
/**
 * Function transforms base from small int number to ACTG base.
 *	'0' - 'A'
 *	'1' - 'C'
 *	'2' - 'G'
 *	'3' - 'T'
 *	'4' - 'N'
 *
 * @param base base
 * @return number assigned to base
 */
inline char intToBase(const uint8_t num) {
  switch (num) {
    case 0:
      return 'A';
    case 1:
      return 'C';
    case 2:
      return 'G';
    case 3:
      return 'T';
    case 4:
      return 'N';
    default:
      return 0;
  }
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

  switch (base) {
    case 'A':
      return 'T';
    case 'C':
      return 'G';
    case 'G':
      return 'C';
    case 'T':
      return 'A';
    default:
      return 'N';
  }
}

/**
 * Function that calculates complement of base as small int number.
 * Valid bases are A, C, G, T and N given as small int numbers 0, 1, 2, 3
 *
 * @param base base whose complement is being calculated given as small int number
 * @return complement of given base also as small int
 */
inline uint8_t getACGTComplementAsSmallInt(uint8_t base) {

  switch (base) {
    case 0:
      // A to T
      return 3;
    case 1:
      // C to G
      return 2;
    case 2:
      // G to C
      return 1;
    case 3:
      // T to A
      return 0;
    default:
      return 4;
  }

}

}  // end namespace

#endif /* UTIL_H_ */
