/*
 * mapper.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <algorithm>
#include "../external/kseq.h"
#include "../external/edlib.h"
#include "../bioutil/bioutil.h"
#include "mapper.h"
#include "mapping.h"
#include "lis.h"
#include "../util/utility_functions.h"

using namespace bioutil;

void Mapper::getKmerPositions(
    Read* read, SuffixArray* sa,
    std::vector<std::pair<uint32_t, uint32_t> > &positions, int kmerStart) {

  int numOfSolutions;
  const int* matches = sa->search(read->data() + kmerStart, KMER_K,
                                  &numOfSolutions);

  if (*matches == -1) {
    return;
  }
  for (int i = 0; i < numOfSolutions; ++i) {
    positions.push_back(std::make_pair(matches[i], kmerStart));
  }
}

void Mapper::getKmerPositions3(Read* read, SuffixArray* sa,
                               std::vector<Triplet> &positions, int kmerStart) {

  int numOfSolutions;
  int len;
  const int* matches = sa->iterativeSearch(read->data() + kmerStart,
                                           read->dataLen() - kmerStart - 1,
                                           KMER_K,
                                           &numOfSolutions, MAX_MATCH_NUM,
                                           MIN_MATCH_NUM,
                                           &len);

  if (*matches == -1) {
    return;
  }
  for (int i = 0; i < numOfSolutions; ++i) {
    positions.push_back(Triplet(matches[i], kmerStart, len));
  }
}

void Mapper::mapReadToSuffixArray(Read* read, SuffixArray* sa,
                                  bool generateCIGAR) {
  read->allBasesToSmallInt();
  Read* reverse_complement = read->getReverseComplement();
  fillMappings3(read, sa);
  fillMappings3(reverse_complement, sa);

  std::multiset<Mapping*, ptr_compare<Mapping> >::reverse_iterator it =
      reverse_complement->mappings().rbegin();
  for (; it != reverse_complement->mappings().rend(); ++it) {
    (*it)->setComplement(true);
    read->addMapping((*it)->score(), (*it)->start(), (*it)->end(), true, NULL,
                     0);
  }

  if (!generateCIGAR) {
    delete reverse_complement;
    return;
  }

  std::multiset<Mapping*, ptr_compare<Mapping> > tmp_set = read->mappings();
  read->mappings().clear();

  for (it = tmp_set.rbegin(); it != tmp_set.rend(); ++it) {
    int32_t start = std::max<int32_t>(0, (*it)->start() - SW_START_OFFSET);
    uint32_t end = std::min<uint32_t>(sa->size() - 1,
                                      (*it)->end() + SW_END_OFFSET);

    int score, numLocations, alignmentLength;
    int* startLocations;
    int* endLocations;
    unsigned char* alignment;

    if ((*it)->isComplement()) {
      edlibCalcEditDistance((const unsigned char *) reverse_complement->data(),
                            reverse_complement->dataLen(),
                            (const unsigned char *) (sa->text() + start),
                            end - start + 1, 5, -1, EDLIB_MODE_HW, true, true,
                            &score, &endLocations, &startLocations,
                            &numLocations, &alignment, &alignmentLength);

    } else {
      edlibCalcEditDistance((const unsigned char *) read->data(),
                            read->dataLen(),
                            (const unsigned char *) (sa->text() + start),
                            end - start + 1, 5, -1, EDLIB_MODE_HW, true, true,
                            &score, &endLocations, &startLocations,
                            &numLocations, &alignment, &alignmentLength);
    }

    char* cigar;
    edlibAlignmentToCigar(alignment, alignmentLength, EDLIB_CIGAR_EXTENDED,
                          &cigar);

    int newScore = read->dataLen() - score;
    start = startLocations[0] + start;
    end = endLocations[0] + start;

    read->addMapping(newScore, start, end, (*it)->isComplement(), cigar,
                     strlen(cigar));
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

void Mapper::runLCSk3(int startIndex, int endIndex, std::vector<Triplet> &pos,
                      Read* read) {
  std::vector<Triplet> lcsKData;

  // TODO Provjera koliko kmera se nalazi izmedu start i end
  // te ukoliko je manje od mog minimuma preskoci
  // ovisno o kvaliteti readova

  for (int i = startIndex; i < endIndex; ++i) {
    lcsKData.push_back(pos[i]);
  }

  std::vector<std::pair<uint32_t, uint32_t> > result;
  uint32_t score = LCSk::calcLCSpp(result, lcsKData);

  //result contains pairs (refPos, readPos)
  int len = result.size();
  int32_t beginPos = result[0].first - result[0].second;
  int32_t endPos = beginPos + read->dataLen() * 2l;

  if (len > 2) {
    endPos = std::min<int32_t>(endPos, result[len - 1].first);
  }

  beginPos = std::max(beginPos, 0);
  read->addMapping(score, beginPos, endPos, false, NULL, 0);
}

void Mapper::fillMappings3(Read* read, SuffixArray* sa) {
  std::vector<Triplet> pos;
  assert(read->dataLen() >= KMER_K);

  for (uint32_t i = KMER_K; i < read->dataLen(); ++i) {
    getKmerPositions3(read, sa, pos, i - KMER_K);
  }

  std::sort(pos.begin(), pos.end());
  if (pos.size() == 0) {
    return;
  }

  uint32_t startIndex = 0;
  uint32_t endIndex = 0;
  while (true) {

    // prozor dug max WINDOW SIZE * duzina reada
    uint32_t windowEnd = pos[startIndex].first + WINDOW_SIZE * read->dataLen();
    for (; endIndex < pos.size() && pos[endIndex].first < windowEnd;
        ++endIndex) {
      ;
    }

    runLCSk3(startIndex, endIndex, pos, read);

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

void Mapper::fillMappings(Read* read, SuffixArray* sa) {
  std::vector<std::pair<uint32_t, uint32_t> > pos;
  assert(read->dataLen() >= KMER_K);

  for (uint32_t i = KMER_K; i < read->dataLen(); ++i) {
    getKmerPositions(read, sa, pos, i - KMER_K);
  }

  std::sort(pos.begin(), pos.end());
  if (pos.size() == 0) {
    return;
  }

  uint32_t startIndex = 0;
  uint32_t endIndex = 0;
  while (true) {

    // prozor dug max WINDOW SIZE * duzina reada
    uint32_t windowEnd = pos[startIndex].first + WINDOW_SIZE * read->dataLen();
    for (; endIndex < pos.size() && pos[endIndex].first < windowEnd;
        ++endIndex) {
      ;
    }

    runLCSk(startIndex, endIndex, pos, read);

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

void Mapper::runLIS(int startIndex, int endIndex,
                    std::vector<std::pair<uint32_t, uint32_t> > &pos,
                    Read* read) {
  std::vector<std::pair<uint32_t, uint32_t> > lisData;

  for (int i = startIndex; i < endIndex; ++i) {
    lisData.push_back(pos[i]);
  }

  std::vector<int> lisResult;
  Lis::calcLIS(&lisResult, lisData);

  int beginPos = Lis::estimateBeginingPosFromLIS(pos, lisResult);
  read->addMapping(lisResult.size(), beginPos, beginPos + read->dataLen(),
                   false, NULL, 0);
}

void Mapper::runLCSk(int startIndex, int endIndex,
                     std::vector<std::pair<uint32_t, uint32_t> > &pos,
                     Read* read) {
  std::vector<std::pair<uint32_t, uint32_t> > lcsKData;

  // TODO Provjera koliko kmera se nalazi izmedu start i end
  // te ukoliko je manje od mog minimuma preskoci
  // ovisno o kvaliteti readova

  for (int i = startIndex; i < endIndex; ++i) {
    lcsKData.push_back(pos[i]);
  }

  std::vector<std::pair<uint32_t, uint32_t> > result;
  uint32_t score = LCSk::calcLCSkpp(KMER_K, result, lcsKData);

  //result contains pairs (refPos, readPos)
  int32_t beginPos = result[0].first - result[0].second;
  int32_t endPos = beginPos + read->dataLen() * 2l;

  int len = result.size();
  if (len > 2) {
    endPos = std::min<int32_t>(endPos, result[len - 1].first);
  }

  beginPos = std::max(beginPos, 0);
  read->addMapping(result.size(), beginPos, beginPos + read->dataLen(), false,
  NULL,
                   0);
}

void copyFromTmpToFile(char* tmpFileName, FILE* src, FILE *dest) {
  char buffer[4096];

  while (fgets(buffer, sizeof buffer, src) != NULL) {
    fputs(buffer, dest);
  }

  fclose(src);
  remove(tmpFileName);
}

void fillSAMHeader(FILE* out, Sequence* seq) {
  fprintf(out, "@HD\tVN:1.4\tSQ:unsorted\n");

  for (uint32_t i = 0; i < seq->numOfSequences(); ++i) {
    fprintf(out, "@SQ\tSN:%s\tLN:%u\n", seq->info(i), seq->seqLen(i));
  }
  fprintf(out, "@PG\tID:mapper\tPN:mapper\n");
}

void Mapper::mapAllReads(char* readsInPath, char* solutionOutPath,
                         SuffixArray* sa, Sequence* seq, bool generateCIGAR) {
  FILE* readsIn = fopen(readsInPath, "r");
  FILE* solutionOut = fopen(solutionOutPath, "w");

  fillSAMHeader(solutionOut, seq);
  kseq_t *kseq = kseq_init(fileno(readsIn));
  Read* singleRead = new Read;

  int threadNum = omp_get_num_procs();
  //threadNum = 3;

  omp_set_dynamic(0);
  omp_set_num_threads(threadNum);

  FILE* tmpOutput[threadNum];
  char tmpFilesNames[threadNum][50];

  for (int i = 0; i < threadNum; ++i) {
    sprintf(tmpFilesNames[i], "tmp%d.out", i);
    assert(validateOutputFile(tmpFilesNames[i]));
    tmpOutput[i] = fopen(tmpFilesNames[i], "w");
  }
  printf("Tmp files created\n");

#pragma omp parallel
  {
#pragma omp single
    {

      uint32_t cntr = 0;

      while (singleRead->readNextFromFASTQ(kseq)) {
        Read* read = singleRead;

        // create task solve single read
#pragma omp task firstprivate(read) shared(tmpOutput) shared(sa) shared(seq)
        {
          mapReadToSuffixArray(read, sa, generateCIGAR);
          read->printReadSAM(tmpOutput[omp_get_thread_num()], seq);
          delete read;
        }

        if (cntr % 250 == 0) {
          printf("-%d-\n", cntr);
        }
        ++cntr;
        singleRead = new Read();
      }
    }
#pragma omp barrier
  }

  delete singleRead;

  for (int i = 0; i < threadNum; ++i) {
    fclose(tmpOutput[i]);
    tmpOutput[i] = fopen(tmpFilesNames[i], "r");
  }
  printf("Merging tmp files\n");
  for (int i = 0; i < threadNum; ++i) {
    copyFromTmpToFile(tmpFilesNames[i], tmpOutput[i], solutionOut);
  }

  fclose(solutionOut);
  fclose(readsIn);
  printf("Completed\n");
}

