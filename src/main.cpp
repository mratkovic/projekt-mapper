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
 * File contains main function of program that maps long reads to reference gene.
 *
 */

#include <core/aligner.h>

int main(int argc, char **argv) {
  Aligner aligner;
  aligner.run(argc, argv);
}

