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
#include <omp.h>

#include "core/mapper.h"
#include "core/suffix_array.h"
#include "bioinf/sequence.h"
#include "util/utility_functions.h"
#include "core/incremental_lcsk_solver.h"

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
void mapReads(char* fastaInPath, char* saFile, char* readsInPath,
              char* outputFilePath, uint32_t threadNum);

/**
 * Main method. Entry point of this project.
 *
 * @param parameters are described in header of this file
 * @return exit status
 */
int main(int argc, char **argv) {
  if (argc != 4 && argc != 6) {
    verboseUsageAndExit();
  }

  if (!strcmp(argv[1], "index") && argc == 4) {
    constructSA(argv[2], argv[3]);

  } else if (!strcmp(argv[1], "map") && argc == 6) {
    mapReads(argv[2], argv[3], argv[4], argv[5], omp_get_num_procs() / 2);  //TODO

  } else {
    verboseUsageAndExit();
  }

  printf("\nEXIT 0\n");
  return 0;
}

void verboseUsageAndExit() {
  printf("Invalid call\n");
  printf("mapper index <fastaFile> <suffixArrayOutputFile>\n");
  printf(
      "mapper map <fastaFile> <suffixArrayOutputFile> <reads> <resultOutputFile>\n");
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

void mapReads(char* fastaInPath, char* saFile, char* readsInPath,
              char* outputFilePath, uint32_t threadNum) {

  assert(validateInputFile(fastaInPath));
  assert(validateInputFile(saFile));
  assert(validateInputFile(readsInPath));
  assert(validateOutputFile(outputFilePath));

  FILE* fastaIn = fopen(fastaInPath, "r");
  Sequence *seq = new Sequence;
  fprintf(stderr, "Reading sequence file %s\n", fastaInPath);
  seq->readSequencesFromFASTA(fastaIn);
  seq->allBasesToSmallInt();
  fclose(fastaIn);

  IncrementalLCSkSolver* solver = new IncrementalLCSkSolver(seq);
  solver->kmerK_ = 8;
  solver->maxMatchNum_ = 18;
  solver->minMatchNum_ = 10;

  solver->readSuffixArrayFromFile(saFile);
  Mapper* mapper = new Mapper(seq, solver, 4);
  printf("Mapping reads to sequence\n");
  mapper->mapAllReads(readsInPath, outputFilePath);

  printf("\nComplete\n");

  delete mapper;
  delete solver;
  delete seq;

}

