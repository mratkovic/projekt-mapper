/*
 * read.h
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#ifndef SRC_CORE_READ_H_
#define SRC_CORE_READ_H_

#include <zlib.h>
#include <set>

#include "sequence.h"
#include "bioutil.h"
#include "../external/kseq.h"
#include "../core/suffix_array.h"
#include "../core/mapping.h"


KSEQ_INIT(int, read)

template<typename T, typename Pred = std::less<T> >
struct ptr_compare : Pred {
  ptr_compare(Pred const & p = Pred())
      : Pred(p) {
  }

  bool operator()(T const * p1, T const * p2) const {
    return Pred::operator()(*p1, *p2);
  }
};

class Read {
 public:
  Read();
  ~Read();

  void clear();
  bool readNextFromFASTQ(kseq_t *seq);
  void allBasesToSmallInt();
  void allBasesToLetters();

  Read* getReverseComplement();
  void printReadSAM(FILE* outFile, bioutil::Sequence* seq);

  const char* data();
  const char* id();
  uint32_t dataLen();

  void addMapping(double score, uint32_t start, uint32_t end, bool isComplement,
                  const char* cigar, uint32_t cigarLen);

  std::multiset<Mapping*, ptr_compare<Mapping> >& mappings();
  Mapping* bestMapping(uint32_t index);
  uint32_t mappingsSize();
  bool basesInt();

 private:
  char* id_;
  char* data_;
  uint32_t dataLen_;

  char* optional_identifier_;
  char* quality_;

  bool basesAsInt_;
  std::multiset<Mapping*, ptr_compare<Mapping> > mappings_;

};

#endif /* SRC_CORE_READ_H_ */
