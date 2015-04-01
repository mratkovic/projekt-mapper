/*
 * position.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_POSITION_H_
#define SRC_UTIL_POSITION_H_

#include <stdint.h>

class Position {

 public:
  Position();
  Position(uint32_t score, uint32_t start, uint32_t end, bool complemented,
           const char* cigar, uint32_t cigarLen);

  virtual ~Position();
  bool isComplement();
  uint32_t score() const;
  uint32_t end();
  uint32_t start();

  void end(uint32_t end);
  void start(uint32_t start);
  void cigar(const char* cigar, uint32_t len);
  void setComplement(bool complemented);

  const char* cigar();
  bool operator <(const Position &other) const {
    return score_ <= other.score();

  }

 private:
  uint32_t score_;
  uint32_t start_;
  uint32_t end_;

  bool complement_;
  char* cigar_;

};

#endif /* SRC_UTIL_POSITION_H_ */
