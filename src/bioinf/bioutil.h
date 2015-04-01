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

  if (base == 'A')
    return 0;
  if (base == 'C')
    return 1;
  if (base == 'G')
    return 2;
  if (base == 'T')
    return 3;
  if (base == 'N')
    return 4;

  return UCHAR_MAX;

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
  if (num == 0)
    return 'A';
  if (num == 1)
    return 'C';
  if (num == 2)
    return 'G';
  if (num == 3)
    return 'T';

  if (num == 4)
    return 'N';

  // should never happen
  return 0;
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

/**
 * Function that calculates complement of base as small int number.
 * Valid bases are A, C, G, T and N given as small int numbers 0, 1, 2, 3
 *
 * @param base base whose complement is being calculated given as small int number
 * @return complement of given base also as small int
 */
inline uint8_t getACGTComplementAsSmallInt(uint8_t base) {
  // A to T
  if (base == 0)
    return 3;
  else if (base == 3)
    return 0;

  // C to G
  else if (base == 1)
    return 2;
  else if (base == 2)
    return 1;

  // N
  else
    return 4;

}

}  // end namespace

#endif /* UTIL_H_ */
