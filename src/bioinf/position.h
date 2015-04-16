/*
 * position.h
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#ifndef SRC_UTIL_POSITION_H_
#define SRC_UTIL_POSITION_H_

#include <stdint.h>
#include <cstdlib>

class Position {

 public:
  Position();
  Position(uint32_t score,  uint32_t secondaryScore,uint32_t start, uint32_t end, bool complemented,
           const char* cigar, uint32_t cigarLen);

  virtual ~Position();
  bool isComplement();
  uint32_t score() const;
  uint32_t secondaryScore() const;
  uint32_t end();
  uint32_t start() const;

  void end(uint32_t end);
  void start(uint32_t start);
  void cigar(const char* cigar, uint32_t len);
  void setComplement(bool complemented);
  void score(uint32_t score);
  void secondaryScore(uint32_t secondary_score);

  const char* cigar();
  bool operator <(const Position &other) const {
    if (score_ != other.score()) {
      return score_ <= other.score();
    }

    if (secondary_score_ != other.secondaryScore()) {
      return secondary_score_ <= other.secondaryScore();
    }

    if (abs(start_ - other.start()) > 20) {
      return start_ - other.start();

    }

    return 0;

  }

 private:
  uint32_t score_;
  uint32_t secondary_score_;
  uint32_t start_;
  uint32_t end_;

  bool complement_;
  char* cigar_;

};

#endif /* SRC_UTIL_POSITION_H_ */
