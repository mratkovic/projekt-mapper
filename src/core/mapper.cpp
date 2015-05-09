/*
 * mapper.cpp
 *
 *  Created on: Dec 13, 2014
 *      Author: marko
 */

#include <algorithm>
#include "external/kseq.h"
#include "external/edlib.h"
#include "bioinf/bioutil.h"
#include "mapper.h"
#include <omp.h>

using namespace bioinf;


Mapper::Mapper(Sequence* seq, Solver* solver, uint32_t threadNum,
               float minKeepScoreRatio, uint32_t maxPositionsPerRead)
    : seq_(seq),
      solver_(solver),
      threadNum_(threadNum),
      minKeepScoreRatio_(minKeepScoreRatio),
      maxPositionsPerRead_(maxPositionsPerRead) {

}

void Mapper::copyFromTmpFileAndDelete(char* tmpFileName, FILE* src,
                                      FILE *dest) {
  char buffer[4096];

  while (fgets(buffer, sizeof buffer, src) != NULL) {
    fputs(buffer, dest);
  }

  fclose(src);
  remove(tmpFileName);
}

void Mapper::createTmpFiles(FILE* tempFiles[],
                            char tmpFileNames[][MAX_TMP_NAME_LEN],
                            uint32_t numOfFiles) {

  for (uint32_t i = 0; i < numOfFiles; ++i) {
    sprintf(tmpFileNames[i], "tmp%d.out", i);
    assert(validateOutputFile(tmpFileNames[i]));
    tempFiles[i] = fopen(tmpFileNames[i], "w");

  }
  printf("Tmp files created\n");

}

void Mapper::fillSAMHeader(FILE* out) {
  fprintf(out, "@HD\tVN:1.4\tSQ:unsorted\n");

  for (uint32_t i = 0; i < seq_->numOfSequences(); ++i) {
    fprintf(out, "@SQ\tSN:%s\tLN:%u\n", seq_->info(i), seq_->seqLen(i));
  }
  fprintf(out, "@PG\tID:mapper\tPN:mapper\n");
}

void Mapper::mergeTmpFiles(char fileNames[][MAX_TMP_NAME_LEN], FILE* tmpFiles[],
                           FILE* solutionFile, int numberOfFiles) {

  for (int i = 0; i < numberOfFiles; ++i) {
    fclose(tmpFiles[i]);
    tmpFiles[i] = fopen(fileNames[i], "r");
  }
  printf("Merging tmp files\n");
  for (int i = 0; i < numberOfFiles; ++i) {
    copyFromTmpFileAndDelete(fileNames[i], tmpFiles[i], solutionFile);
  }
}

void Mapper::mapAllReads(char* readsInPath, char* solutionOutPath) {
  fprintf(stderr, "Stats: t:%d\n", threadNum_);
  solver_->printInfo();

  FILE* readsIn = fopen(readsInPath, "r");
  FILE* solutionOut = fopen(solutionOutPath, "w");

  fillSAMHeader(solutionOut);
  kseq_t* kseq = kseq_init(fileno(readsIn));
  Read* singleRead = new Read(minKeepScoreRatio_, maxPositionsPerRead_);

  omp_set_dynamic(0);
  omp_set_num_threads(threadNum_);

  FILE* tmpOutput[threadNum_];
  char tmpFileNames[threadNum_][MAX_TMP_NAME_LEN];
  createTmpFiles(tmpOutput, tmpFileNames, threadNum_);

  // OPENMP workaround
  Sequence* seq = seq_;
  Solver* solver = solver_;

  #pragma omp parallel
  {
     #pragma omp single
    {

      uint32_t cntr = 0;

      while (singleRead->readNextFromFASTQ(kseq)) {
        Read* read = singleRead;

        #pragma omp task firstprivate(read) shared(tmpOutput) shared(seq) shared(solver)
        {
          // create task for  solving single read
          solver->findReadPosition(read);
          read->printReadSAM(tmpOutput[omp_get_thread_num()], seq);
          delete read;
        }

        if (cntr % 25 == 0) {
          fprintf(stderr, "-%d-\n", cntr);
        }
        ++cntr;
        singleRead = new Read(minKeepScoreRatio_, maxPositionsPerRead_);
      }
    }
  #pragma omp barrier
  }

  mergeTmpFiles(tmpFileNames, tmpOutput, solutionOut, threadNum_);

  delete singleRead;
  fclose(solutionOut);
  fclose(readsIn);
  fprintf(stderr, "Completed\n");
}

