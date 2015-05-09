/*
 * aligner.cpp
 *
 *  Created on: May 9, 2015
 *      Author: marko
 */

#include <core/aligner.h>
#include <bioinf/sequence.h>
#include <core/mapper.h>

#define MODE_LEN 25

Aligner::Aligner(uint32_t kmerK, uint32_t lowerMatchesLimit,
                 uint32_t upperMatchesLimit, uint32_t threadNum,
                 float keepFactor, uint32_t maxPos)
    : kmerK_(kmerK),
      lowerMatchesLimit_(lowerMatchesLimit),
      upperMatchesLimit_(upperMatchesLimit),
      threadNum_(threadNum),
      keepFactor_(keepFactor),
      maxPos_(maxPos) {
}

Aligner::~Aligner() {
}

void Aligner::run(int argc, char **argv) {
  int option;
  char mode[MODE_LEN];
  while ((option = getopt(argc, argv, "m:t:k:l:h:kf:pos:")) >= 0) {
    switch (option) {
      case 'm':
        strncpy(mode, optarg, MODE_LEN);
        break;
    }
  }

  if (!strcmp(mode, "index") && optind + 2 == argc) {
    constructSA(argv[optind], argv[optind + 1]);

  } else if (!strcmp(mode, "map") && optind + 4 == argc) {
    optind = 0;
    mapReads(argc, argv);

  } else {
    verboseUsageAndExit();
  }

  printf("\nEXIT 0\n");
}

void Aligner::verboseUsageAndExit() {

  fprintf(stderr, "\n");
  fprintf(
      stderr,
      "Usage - create index: mapper -m index <fastaFile> <suffixArrayOutputFile>\n");
  fprintf(
      stderr,
      "Usage - map reads: mapper -m map [options...] <fastaFile> <suffixArrayOutputFile> <readsFASTQ> <resultOutputFile>\n");

  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "  -t N  N is thread number. [default: number of cores]\n");
  fprintf(stderr, "  -k N  N is seed length. [default: 15]\n");
  fprintf(
      stderr,
      "  -l lowerLimit  lowerLimit is minimum required number of hits. [default: / ]\n");
  fprintf(
      stderr,
      "  -h upperLimit  upperLimit is maximum allowed number of hits. [default: / ]\n");

  fprintf(
      stderr,
      "  -kf factor  factor float value that represents minimum ratio bestScore/score of reads that are being kept. [default: 1.2 ]\n");

  fprintf(
      stderr,
      "  -pos N  N is maximum allowed number of positions that are kept per read [default: 80 ]\n");

  exit(-1);
}

void Aligner::constructSA(char* fastaInPath, char* saOutputPath) {
  assert(validateInputFile(fastaInPath));
  assert(validateOutputFile(saOutputPath));

  FILE* fastaIn = fopen(fastaInPath, "r");
  bioinf::Sequence *seq = new bioinf::Sequence;

  fprintf(stderr, "Reading sequence file %s\n", fastaInPath);
  seq->readSequencesFromFASTA(fastaIn);
  fprintf(stderr, "Number of sequences: %lu\n",
          (unsigned long) seq->numOfSequences());

  for (uint32_t i = 0; i < seq->numOfSequences(); ++i) {
    fprintf(stderr, "\t%d: %s: [%d]\n", i, seq->info(i), seq->seqLen(i));
  }
  fclose(fastaIn);

  fprintf(stderr, "Constructing suffix array started\n");
  SuffixArray *sa = new SuffixArray(seq->data(), seq->dataLen());
  fprintf(stderr, "Constructing suffix array completed\n");
  fprintf(stderr, "Size of array: %d\n", sa->size());

  FILE* saOut = fopen(saOutputPath, "wb");
  fprintf(stderr, "Saving suffix array to %s\n", saOutputPath);
  sa->saveSuffixArray(saOut);
  fprintf(stderr, "Completed\n");

  delete sa;
  delete seq;
  fclose(saOut);
}

void Aligner::mapReads(int argc, char** argv) {

  int option;
  while ((option = getopt(argc, argv, "m:t:k:l:h:kf:pos:")) >= 0) {
    switch (option) {
      case 't':
        threadNum_ = atoi(optarg);
        break;
      case 'k':
        kmerK_ = atoi(optarg);
        break;
      case 'l':
        lowerMatchesLimit_ = atoi(optarg);
        break;
      case 'h':
        upperMatchesLimit_ = atoi(optarg);
        break;
      case 'kf':
        keepFactor_ = atof(optarg);
        break;
      case 'pos':
        maxPos_ = atoi(optarg);
        break;
    }
  }

  char* fastaInPath = argv[optind];
  char* saFile = argv[optind + 1];
  char* readsInPath = argv[optind + 2];
  char* outputFilePath = argv[optind + 3];
  mapReads(fastaInPath, saFile, readsInPath, outputFilePath);

}

void Aligner::mapReads(char* fastaInPath, char* saInputPath, char* readsInPath,
                       char* SAMOutputPath) {

  assert(validateInputFile(fastaInPath));
  assert(validateInputFile(saInputPath));
  assert(validateInputFile(readsInPath));
  assert(validateOutputFile(SAMOutputPath));

  FILE* fastaIn = fopen(fastaInPath, "r");
  bioinf::Sequence *seq = new bioinf::Sequence;
  fprintf(stderr, "Reading sequence file %s\n", fastaInPath);
  seq->readSequencesFromFASTA(fastaIn);
  seq->allBasesToSmallInt();
  fclose(fastaIn);

  IncrementalLCSkSolver* solver = new IncrementalLCSkSolver(seq);

  solver->kmerK_ = kmerK_;
  solver->maxMatchNum_ = upperMatchesLimit_;
  solver->minMatchNum_ = lowerMatchesLimit_;

  solver->readSuffixArrayFromFile(saInputPath);
  Mapper* mapper = new Mapper(seq, solver, threadNum_, keepFactor_, maxPos_);
  fprintf(stderr, "Mapping reads to sequence\n");
  mapper->mapAllReads(readsInPath, SAMOutputPath);

  fprintf(stderr, "\nComplete\n");

  delete mapper;
  delete solver;
  delete seq;

}
