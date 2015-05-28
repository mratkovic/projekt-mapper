/*
 * lcsk_solver.cpp
 *
 *  Created on: Mar 31, 2015
 *      Author: marko
 */

#include <algorithm>
#include <core/lcsk_solver.h>
#include "metrics_algorithm/lcskpp.h"
#include "external/edlib.h"
#include "external/ssw_cpp.h"

LCSkSolver::LCSkSolver(bioinf::Sequence* seq) {

  kmerK_ = KMER_K;
  windowSize_ = WINDOW_SIZE;
//  swStartOffset_ = SW_START_OFFSET;
//  swEndOffset_ = SW_END_OFFSET;
  seq_ = seq;
  sa_ = NULL;

}

LCSkSolver::~LCSkSolver() {
  delete sa_;
}

bioinf::Sequence* LCSkSolver::seq() {
  return seq_;

}
SuffixArray* LCSkSolver::sa() {
  return sa_;
}

void LCSkSolver::readSuffixArrayFromFile(const char* saInPath) {
  FILE* saIn = fopen(saInPath, "rb");

  fprintf(stderr, "Reading suffix array from file\n");
  sa_ = new SuffixArray(saIn, seq_->data(), seq_->dataLen());

  fprintf(stderr, "SuffixArray read\n");
  fclose(saIn);
}

void LCSkSolver::fillPositions(bioinf::Read* read) {
  std::vector<std::pair<uint32_t, uint32_t> > pos;
  assert(read->dataLen() >= kmerK_);

  for (uint32_t i = kmerK_; i < read->dataLen(); ++i) {
    getKmerPositions(read, pos, i - kmerK_);
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

void LCSkSolver::runLCSkpp(int startIndex, int endIndex,
                           std::vector<std::pair<uint32_t, uint32_t>> &pos,
                           bioinf::Read* read) {

  std::vector<std::pair<uint32_t, uint32_t> > lcsKData;

  // TODO Provjera koliko kmera se nalazi izmedu start i end
  // te ukoliko je manje od mog minimuma preskoci
  // ovisno o kvaliteti readova

  //TODO dodat reserve

  for (int i = startIndex; i < endIndex; ++i) {
    lcsKData.push_back(pos[i]);
  }

  std::vector<std::pair<uint32_t, uint32_t> > result;
  uint32_t score = LCSkpp::calcLCSkpp(kmerK_, result, lcsKData);

  //result contains pairs (refPos, readPos)
  int32_t beginPos = result[0].first - result[0].second;
  int32_t endPos = beginPos + read->dataLen() * 2l;

  int len = result.size();
  if (len > 2) {
    endPos = std::min<int32_t>(endPos, result[len - 1].first);
  }

  beginPos = std::max(beginPos, 0);
  read->addPosition(result.size(), beginPos, beginPos + read->dataLen());

}

void LCSkSolver::getKmerPositions(
    bioinf::Read* read, std::vector<std::pair<uint32_t, uint32_t>> &positions,
    int kmerStart) {

  int numOfSolutions;
  const int* matches = sa_->search(read->data() + kmerStart, kmerK_,
                                   &numOfSolutions);

  if (*matches == -1) {
    return;
  }
  for (int i = 0; i < numOfSolutions; ++i) {
    positions.push_back(std::make_pair(matches[i], kmerStart));
  }

}

void LCSkSolver::findReadPosition(bioinf::Read* read) {
  if (seq_->basesInt()) {
    read->allBasesToSmallInt();
  }
  bioinf::Read* reverse_complement = read->getReverseComplement();
  fillPositions(read);

  // add to keep best score and skip lcsk calculations
  // that will never give best score
  Position* p = read->bestPosition(0);
  if (p) {
    reverse_complement->addPosition(p->score(), p->start(), p->end(),
                                    p->isComplement());
  }
  fillPositions(reverse_complement);

  std::multiset<Position*, ptr_compare<Position> >::reverse_iterator it =
      reverse_complement->positions().rbegin();

  for (auto& it : reverse_complement->positions()) {
    if (p && it->score() == p->score() && it->start() == p->start()
        && it->end() == p->end()) {
      // skup position that was previously inserted to track best score
      continue;
    }

    it->setComplement(true);
    read->addPosition(it->score(), it->start(), it->end(), true);
  }

  std::set<Position*, ptr_compare<Position> > tmp_set = read->positions();
  read->positions().clear();
  read->keepRatio(1.00f);

  for (it = tmp_set.rbegin(); it != tmp_set.rend(); ++it) {
    int32_t start = std::max<int32_t>(0,
                                      (*it)->start() - 0.5 * read->dataLen());
    uint32_t end = std::min<uint32_t>(sa()->size() - 1,
                                      (*it)->end() + 0.5 * read->dataLen());

    int score, numLocations, alignmentLength;
    int* startLocations;
    int* endLocations;
    unsigned char* alignment;

    if ((*it)->isComplement()) {
      edlibCalcEditDistance((const unsigned char *) reverse_complement->data(),
                            reverse_complement->dataLen(),
                            (const unsigned char *) (sa()->text() + start),
                            end - start + 1, 5, -1, EDLIB_MODE_HW, true, true,
                            &score, &endLocations, &startLocations,
                            &numLocations, &alignment, &alignmentLength);

    } else {
      edlibCalcEditDistance((const unsigned char *) read->data(),
                            read->dataLen(),
                            (const unsigned char *) (sa()->text() + start),
                            end - start + 1, 5, -1, EDLIB_MODE_HW, true, true,
                            &score, &endLocations, &startLocations,
                            &numLocations, &alignment, &alignmentLength);
    }

    char* cigar;
    edlibAlignmentToCigar(alignment, alignmentLength, EDLIB_CIGAR_STANDARD,
                          &cigar);

    start = startLocations[0] + start;
    end = endLocations[0] + start;

//    int secondaryScore = read->dataLen()
//        - abs(read->dataLen() - (end - start + 1));
//
//    secondaryScore = (*it)->secondaryScore();
    int secondaryScore = (*it)->score();
    int newScore = 3 * read->dataLen() - score;

    read->addPosition(newScore, start, end, (*it)->isComplement(), cigar,
                      strlen(cigar), tmp_set.size());
    free(cigar);
    if (endLocations) {
      free(endLocations);
    }
    if (startLocations) {
      free(startLocations);
    }
    if (alignment) {
      free(alignment);
    }
    delete (*it);
  }

  delete reverse_complement;

}

// Pokusaj sa SSW i scored matrix

//void LCSkSolver::findReadPosition(bioinf::Read* read) {
//  if (seq_->basesInt()) {
//    read->allBasesToSmallInt();
//  }
//
//
//  bioinf::Read* reverse_complement = read->getReverseComplement();
//  fillPositions(read);
//  fillPositions(reverse_complement);
//
//  std::multiset<Position*, ptr_compare<Position> >::reverse_iterator it =
//      reverse_complement->positions().rbegin();
//
//  for (; it != reverse_complement->positions().rend(); ++it) {
//    (*it)->setComplement(true);
//    read->addPosition((*it)->score(), (*it)->start(), (*it)->end(), true);
//  }
//
//  std::set<Position*, ptr_compare<Position> > tmp_set = read->positions();
//  read->positions().clear();
//  read->keepRatio(1.00f);
//
//
//
//
//  StripedSmithWaterman::Aligner aligner(SSW_MATCH, SSW_MISMATCH, SSW_GAP_OPEN,
//  SSW_GAP_EXTEND);
//  StripedSmithWaterman::Filter filter;
//  uint32_t missmatchMaxNum = read->dataLen() * 0.5;
//
//  filter.score_filter = (read->dataLen() - missmatchMaxNum) * SSW_MATCH
//      - missmatchMaxNum * SSW_MISMATCH;
//
//  for (it = tmp_set.rbegin(); it != tmp_set.rend(); ++it) {
//    int32_t start = std::max<int32_t>(0,
//                                      (*it)->start() - 0.5 * read->dataLen());
//    uint32_t end = std::min<uint32_t>(sa()->size() - 1,
//                                      (*it)->end() + 0.5 * read->dataLen());
//
//    StripedSmithWaterman::Alignment alignment;
//
//    if ((*it)->isComplement()) {
//      aligner.Align(reverse_complement->data(), sa_->text() + start,
//                    end - start + 1, filter, &alignment);
//    } else {
//      aligner.Align(read->data(), sa_->text() + start, end - start + 1, filter,
//                    &alignment);
//    }
//
//    end = start + alignment.ref_end;
//    start += alignment.ref_begin;
//
////    int score, numLocations, alignmentLength;
////    int* startLocations;
////    int* endLocations;
////    unsigned char* alignment;
////
////    if ((*it)->isComplement()) {
////      edlibCalcEditDistance((const unsigned char *) reverse_complement->data(),
////                            reverse_complement->dataLen(),
////                            (const unsigned char *) (sa()->text() + start),
////                            end - start + 1, 5, -1, EDLIB_MODE_HW, true, true,
////                            &score, &endLocations, &startLocations,
////                            &numLocations, &alignment, &alignmentLength);
////
////    } else {
////      edlibCalcEditDistance((const unsigned char *) read->data(),
////                            read->dataLen(),
////                            (const unsigned char *) (sa()->text() + start),
////                            end - start + 1, 5, -1, EDLIB_MODE_HW, true, true,
////                            &score, &endLocations, &startLocations,
////                            &numLocations, &alignment, &alignmentLength);
////    }
////
////    char* cigar;
////    edlibAlignmentToCigar(alignment, alignmentLength, EDLIB_CIGAR_EXTENDED,
////                          &cigar);
////
////
////
////    start = startLocations[0] + start;
////    end = endLocations[0] + start;
//
//    int secondaryScore = read->dataLen()
//        - abs(read->dataLen() - (end - start + 1));
//
//    secondaryScore = (*it)->secondaryScore();
//    secondaryScore = 0;
//    int newScore = alignment.sw_score;
//
//    read->addPosition(newScore, start, end, (*it)->isComplement(),
//                      alignment.cigar_string.c_str(),
//                      alignment.cigar_string.size(), secondaryScore);
//
//    delete (*it);
//  }
//
//  delete reverse_complement;
//
//}

void LCSkSolver::printInfo() {
  fprintf(stderr, "LCSkSolver: k%d; window:%d;\n", kmerK_, windowSize_);

}

