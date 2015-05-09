/*
 * incremental_lcsk_solver.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#include <core/incremental_lcsk_solver.h>
#include <external/edlib.h>
#include <metrics_algorithm/lcskppV2.h>
#include <stdio.h>
#include <algorithm>

IncrementalLCSkSolver::IncrementalLCSkSolver(bioinf::Sequence* seq)
    : LCSkSolver(seq) {

  maxMatchNum_ = MAX_MATCH_NUM;
  minMatchNum_ = MIN_MATCH_NUM;
}
IncrementalLCSkSolver::~IncrementalLCSkSolver() {
}

void IncrementalLCSkSolver::fillPositions(bioinf::Read* read) {
  std::vector<triplet_t<uint32_t>> pos;
  assert(read->dataLen() >= kmerK_);

  uint32_t len = kmerK_;
  for (uint32_t i = kmerK_; i < read->dataLen(); ++i) {
    len = getKmerPositions(read, pos, i - kmerK_, len);
  }

  std::sort(pos.begin(), pos.end());
  if (pos.size() == 0) {
    return;
  }

  uint32_t startIndex = 0;
  uint32_t endIndex = 0;
  while (true) {

    // prozor dug max WINDOW SIZE * duzina reada
    uint32_t windowEnd = pos[startIndex].first + windowSize_ * read->dataLen();
    for (; endIndex < pos.size() && pos[endIndex].first < windowEnd;
        ++endIndex) {
      ;
    }

    runLCSkpp(startIndex, endIndex, pos, read);

    uint32_t lastPosition = pos[startIndex].first;
    // pomakni pocetak prozora za duljinu reada od proslog pocetka
    for (;
        startIndex < pos.size()
            && pos[startIndex].first < lastPosition + read->dataLen();
        ++startIndex) {
      ;
    }
    if (startIndex >= pos.size()) {
      break;
    }

  }
}

void IncrementalLCSkSolver::runLCSkpp(int startIndex, int endIndex,
                                      std::vector<triplet_t<uint32_t>> &pos,
                                      bioinf::Read* read) {
  std::vector<triplet_t<uint32_t>> lcsKData;

  // TODO Provjera koliko kmera se nalazi izmedu start i end
  // te ukoliko je manje od mog minimuma preskoci
  // ovisno o kvaliteti readova

  // MOZE BEZ OVOG?????, dam iterator i velicinu
  // std::vector<int>::const_iterator a i to e to
  // jeli vridnoo

  uint32_t aproxMaxLen = 0;
  for (int i = startIndex; i < endIndex; ++i) {
    lcsKData.push_back(pos[i]);
    aproxMaxLen += pos[i].third;
  }

  if (read->bestPosition(0) != NULL) {
    uint32_t currentLen = read->bestPosition(0)->score() / KEEP_FACTOR;
    if (currentLen > aproxMaxLen) {
      // no way it could produce better position
      return;
    }
  }

  std::vector<std::pair<uint32_t, uint32_t> > result;
  uint32_t score = LCSkppV2::calcLCSkpp(result, lcsKData);

  //result contains pairs (refPos, readPos)
  int len = result.size();
  int32_t beginPos = result[0].first;
  int32_t endPos = result[len - 1].first;

  if (len > 2) {
    endPos = std::min<int32_t>(
        endPos,
        result[len - 1].first + read->dataLen() - result[len - 1].second);
  }

  beginPos = std::max(beginPos, 0);
  read->addPosition(score, beginPos, endPos, false, NULL, 0, lcsKData.size());
}

uint32_t IncrementalLCSkSolver::getKmerPositions(
    bioinf::Read* read, std::vector<triplet_t<uint32_t>> &positions,
    int kmerStart, uint32_t initialLen) {

  int numOfSolutions;
  int len;
  const int* matches = sa()->iterativeSearch(read->data() + kmerStart,
                                             read->dataLen() - kmerStart - 1,
                                             initialLen, &numOfSolutions,
                                             maxMatchNum_, minMatchNum_, &len);

  if (*matches == -1) {
    return kmerK_;
  }
  for (int i = 0; i < numOfSolutions; ++i) {
    positions.push_back(triplet_t<uint32_t>(matches[i], kmerStart, len));
  }

  return std::max<int>(kmerK_, len - 1);
}

void IncrementalLCSkSolver::printInfo() {
  fprintf(
      stderr,
      "IncrementalLCSkSolver: k:%d; window:%d; max_match:%d; min_match:%d\n",
      kmerK_, windowSize_, maxMatchNum_, minMatchNum_);
}
