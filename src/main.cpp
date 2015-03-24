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

#include "bioutil/bioutil.h"
#include "core/mapper.h"
#include "core/suffix_array.h"
#include "bioutil/read.h"
#include "bioutil/sequence.h"
#include "util/utility_functions.h"
#include "validator.h"

using namespace bioutil;

/**
 * Function that displays valid command line arguments in case of wrong
 * program start and exits program with status -1.
 */
void displayInvalidCallMsg();
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
void mapReads(char* fastaInPath, char* saFile, char* readsInPath, char* outputFilePath);

/**
 * Method that validates SAM output file. Validation is done by comparison between
 * test SAM file and reference SAM file.
 *
 * @param referenceFilePath path to SAM file containing valid information.
 * @param testFilePath path to SAM file containing data that needs to be tested
 */
void validate(char* referenceFilePath, char* testFilePath);

/**
 * Main method. Entry point of this project.
 *
 * @param parameters are described in header of this file
 * @return exit status
 */
int main(int argc, char **argv) {
	if (argc != 4 && argc != 6) {
		displayInvalidCallMsg();
	}

	if (!strcmp(argv[1], "construct") && argc == 4) {
		constructSA(argv[2], argv[3]);

	} else if (!strcmp(argv[1], "map") && argc == 6) {
		mapReads(argv[2], argv[3], argv[4], argv[5]);

	} else if (!strcmp(argv[1], "validate") && argc == 4) {
		validate(argv[2], argv[3]);
	} else {
		displayInvalidCallMsg();
	}
	return 0;
}

void displayInvalidCallMsg() {
	printf("Invalid call\n");
	printf("mapper construct <fastaFile> <suffixArrayOutputFile>\n");
	printf("mapper map <fastaFile> <suffixArrayOutputFile> <reads> <resultOutputFile>\n");
	printf("mapper validate <referenceSAMfile> <testSAMfile>\n");
	exit(-1);
}

void constructSA(char* fastaInPath, char* saOutputPath) {
	assert(validateInputFile(fastaInPath));
	assert(validateOutputFile(saOutputPath));

	FILE* fastaIn = fopen(fastaInPath, "r");
	Sequence *seq = new Sequence;
	printf("Reading sequence file %s\n", fastaInPath);
	seq->readSequencesFromFASTA(fastaIn);
	printf("Number of sequences: %lu\n", (unsigned long) seq->numOfSequences());

	for(uint32_t i = 0; i < seq->numOfSequences(); ++i) {
	  printf("\t%d: %s: [%d]\n", i, seq->info(i));
	}
	fclose(fastaIn);

	printf("Constructing suffix array started\n");
	SuffixArray *sa = new SuffixArray(seq->data(), seq->dataLen());

	printf("Constructing suffix array completed\n");
	printf("Size of array: %d\n", sa->size());

	FILE* saOut = fopen(saOutputPath, "wb");
	printf("Saving suffix array to %s\n", saOutputPath);
	sa->saveSuffixArray(saOut);
	printf("Completed\n");

	delete sa;
	delete seq;
	fclose(saOut);
}

void mapReads(char* fastaInPath, char* saFile, char* readsInPath, char* outputFilePath) {
	assert(validateInputFile(fastaInPath));
	assert(validateInputFile(saFile));
	assert(validateInputFile(readsInPath));
	assert(validateOutputFile(outputFilePath));

	FILE* fastaIn = fopen(fastaInPath, "r");
	Sequence *seq = new Sequence;
	printf("Reading sequence file %s\n", fastaInPath);
	seq->readSequencesFromFASTA(fastaIn);
	seq->allBasesToSmallInt();
	fclose(fastaIn);

	FILE* saIn = fopen(saFile, "rb");
	printf("Reading suffix array from file\n");
	SuffixArray *sa = new SuffixArray(saIn, seq->data(), seq->dataLen());
	printf("SuffixArray read\n");
	fclose(saIn);

	printf("Mapping reads to sequence\n");
	Mapper::mapAllReads(readsInPath, outputFilePath, sa, seq, true);

	delete sa;
	delete seq;

}

void validate(char* referenceFilePath, char* testFilePath) {
	assert(validateInputFile(referenceFilePath));
	assert(validateInputFile(testFilePath));

	FILE* refFile = fopen(referenceFilePath, "r");
	FILE* testFile = fopen(testFilePath, "r");

	Validator::validateSAM(refFile, testFile);
}

