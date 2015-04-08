/*
 * main.cpp
 *
 *  Created on: Mar 23, 2015
 *      Author: marko
 */

/*
 * main.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 * File contains main function of program that maps short reads to reference gene.
 *
 * Usage:
 * mapper construct <fastaFile> <suffixArrayOutputFile>
 * Construction of suffix array for gene given in FASTA format file
 *
 *  mapper map <fastaFile> <suffixArrayOutputFile> <reads> <resultOutputFile>
 *  Using previously constructed suffix array reads are mapped to reference genome in FASTA format file
 *  and that information is stored in SAM format to output file .
 */

#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <vector>
#include <utility>
#include <algorithm>
#include <cstring>
#include <omp.h>

#include "core/mapper.h"
#include "core/suffix_array.h"
#include "bioinf/sequence.h"
#include "util/utility_functions.h"
#include "core/incremental_lcsk_solver.h"
#include "core/lcsk_solver.h"

#define MODE_LEN 25

using namespace bioinf;

/**
 * Function that displays valid command line arguments in case of wrong
 * program start and exits program with status -1.
 */
void verboseUsageAndExit();
/**
 * Method that constructs suffix array from FASTA file and stores
 * suffix array to binary file on given position.
 *
 * @param fastaInPath path to fasta file containing DNA sequence.
 * @param saOutputPath path to file where suffix array is going to be saved.
 */
void constructSA(char* fastaInPath, char* saOutputPath);

/**
 * Method that maps short reads to reference gene.
 *
 * @param fastaInPath path to fasta file containing DNA sequence.
 * @param saFile path to file containing previously constructed suffix array.
 * @param readsInPath path to fastq file containing short reads.
 * @param outputFilePath path to output file for storing mapping information.
 *
 */
void mapReads(int argc, char** argv);

/**
 * Main method. Entry point of this project.
 *
 * @param parameters are described in header of this file
 * @return exit status
 */
int main(int argc, char **argv) {

  int option;
  char mode[MODE_LEN];
  while ((option = getopt(argc, argv, "m:t:k:l:h:")) >= 0) {
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
  return 0;
}

void verboseUsageAndExit() {

  fprintf(stderr, "\n");
  fprintf(stderr,
          "Usage - create index: mapper -m index <fastaFile> <suffixArrayOutputFile>\n");
  fprintf(
      stderr,
      "Usage - map reads: mapper -m map [options...] <fastaFile> <suffixArrayOutputFile> <readsFASTQ> <resultOutputFile>\n");

  fprintf(stderr, "\n");
  fprintf(stderr, "Options:\n");
  fprintf(stderr,
          "  -t N  N is thread number. [default: number of processors]\n");
  fprintf(stderr, "  -k N  N is seed length. [default: 15]\n");
  fprintf(
      stderr,
      "  -l lowerLimit  lowerLimit is minimum required number of hits. [default: / ]\n");
  fprintf(
      stderr,
      "  -h upperLimit  upperLimit is maximum allowed number of hits. [default: / ]\n");

  exit(-1);
}

void constructSA(char* fastaInPath, char* saOutputPath) {
  assert(validateInputFile(fastaInPath));
  assert(validateOutputFile(saOutputPath));

  FILE* fastaIn = fopen(fastaInPath, "r");
  Sequence *seq = new Sequence;

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

void mapReads(int argc, char** argv) {

  uint32_t k = KMER_K;
  uint32_t l = MIN_MATCH_NUM;
  uint32_t h = MAX_MATCH_NUM;
  uint32_t t = omp_get_num_procs();

  int option;
  while ((option = getopt(argc, argv, "m:t:k:l:h:")) >= 0) {
    switch (option) {
      case 't':
        t = atoi(optarg);
        break;
      case 'k':
        k = atoi(optarg);
        break;
      case 'l':
        l = atoi(optarg);
        break;
      case 'h':
        h = atoi(optarg);
        break;
    }
  }

  char* fastaInPath = argv[optind];
  assert(validateInputFile(fastaInPath));

  char* saFile = argv[optind + 1];
  assert(validateInputFile(saFile));

  char* readsInPath = argv[optind + 2];
  assert(validateInputFile(readsInPath));

  char* outputFilePath = argv[optind + 3];
  assert(validateOutputFile(outputFilePath));

  FILE* fastaIn = fopen(fastaInPath, "r");
  Sequence *seq = new Sequence;
  fprintf(stderr, "Reading sequence file %s\n", fastaInPath);
  seq->readSequencesFromFASTA(fastaIn);
  seq->allBasesToSmallInt();
  fclose(fastaIn);

  IncrementalLCSkSolver* solver = new IncrementalLCSkSolver(seq);

  solver->kmerK_ = k;
  solver->maxMatchNum_ = h;
  solver->minMatchNum_ = l;

  solver->readSuffixArrayFromFile(saFile);
  Mapper* mapper = new Mapper(seq, solver, t);
  fprintf(stderr, "Mapping reads to sequence\n");
  mapper->mapAllReads(readsInPath, outputFilePath);

  fprintf(stderr, "\nComplete\n");

  delete mapper;
  delete solver;
  delete seq;

}

