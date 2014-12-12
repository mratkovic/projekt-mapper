/*
 * SuffixArray.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include "SuffixArray.h"
#include <cstdio>
#include <cstring>

#define NUM_OF_SYMBOLS 5
typedef unsigned long ulint;

inline bool leq(int a1, int a2, int b1, int b2) {
	return (a1 < b1 || (a1 == b1 && a2 <= b2));
}                                                   // and triples
inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) {
	return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
}

int SuffixArray::compare(char *pattern, size_t patternLen, ulint startPos) {
	for (ulint i = startPos; i < std::min(this->_size, startPos + patternLen); ++i) {
		int diff = this->_sequence[i] - pattern[i - startPos];
		if (diff != 0) {
			return diff;
		}
	}
	return 0;
}
ulint SuffixArray::findStartingPositions(char *pattern, size_t patternLen, ulint id,
		std::vector<std::pair<int, ulint> > dest) {
	ulint lo, hi, mid;
	lo = 0;
	hi = this->_size;

	while (lo + 1 < hi) {
		mid = (lo + hi) / 2;
		int cmp = this->compare(pattern, patternLen, _array[mid]);
		if (cmp == 0) {

			for (lo = mid; lo > 0 && this->compare(pattern, patternLen, _array[lo - 1]) == 0; lo--) {
				  dest.push_back(std::make_pair (id, this->_array[lo - 1]));
			}
			for (hi = mid; hi < this->_size && this->compare(pattern, patternLen, _array[hi + 1]) == 0; hi++) {
				 dest.push_back(std::make_pair (id, this->_array[hi + 1]));
			}

			return hi - lo + 1;
		} else if (cmp < 0) {
			lo = mid;
		} else {
			hi = mid;
		}
	}
	return 0;
}
void SuffixArray::constructFromGene(bioutil::Gene* gene) {
	if (gene->dataSize() % 3 == 1) {
		this->_size = gene->dataSize() + 2;
	} else if (gene->dataSize() % 3 == 2) {
		this->_size = gene->dataSize() + 1;
	} else {
		this->_size = gene->dataSize();
	}

	this->_array = (ulint *) malloc(this->_size * sizeof(ulint));
	this->_sequence = (char *) malloc(this->_size * sizeof(char));
	gene->turnBaseToInt(this->_sequence);
	construct();
}

void SuffixArray::printSuffixArray(FILE *file) {

	for (ulint i = 0; i < _size; ++i) {
		fprintf(file, "%5lu-- %5lu  ", i, _array[i]);

		for(ulint j = _array[i]; j < std::min(5+_array[i], _size); ++j) {
			fprintf(file, "%d", _sequence[j]);
		}

		fprintf(file, "\n");
	}
}

bool SuffixArray::construct() {
	ulint *tmp = (ulint *) malloc(this->_size * sizeof(ulint));
	for (ulint i = 0; i < this->_size; ++i) {
		tmp[i] = this->_sequence[i];
	}
	SuffixArray::constructArray(tmp, this->_array, this->_size, NUM_OF_SYMBOLS);
	free(tmp);
}

void SuffixArray::radixSortPass(ulint* in, ulint* out, ulint* sequence, ulint n, int numOfSymbols) {
	int* cntr = new int[numOfSymbols];
	for (int i = 0; i < numOfSymbols; i++) {
		cntr[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		cntr[sequence[in[i]]]++;
	}

	for (ulint i = 0, sum = 0; i < numOfSymbols; i++) {
		int t = cntr[i];
		cntr[i] = sum;
		sum += t;
	}
	for (ulint i = 0; i < n; i++) {
		out[cntr[sequence[in[i]]]++] = in[i];
	}
	for (ulint i = 0; i < n; i++) {
	}
	delete[] cntr;
}

void SuffixArray::constructArray(ulint* str, ulint* SA, ulint n, int K) {
	ulint n0 = (n + 2) / 3;
	ulint n1 = (n + 1) / 3;
	ulint n2 = n / 3;

	ulint n02 = n0 + n2;

	ulint* R12 = new ulint[n02 + 3];
	R12[n02] = R12[n02 + 1] = R12[n02 + 2] = 0;

	ulint* SA12 = new ulint[n02 + 3];
	SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;

// R12 indeksi pocetaka mod 1 i mod 2 sufiksa
// +(n0-n1) dodaje jedan viska mod 1 sufiks ako je n%3 == 1
	for (ulint i = 0, j = 0; i < n + (n0 - n1); i++) {
		if (i % 3 != 0) {
			R12[j++] = i;
		}
	}

// radix za sortirat triplete koji pocinju na R12 indeksima
	radixSortPass(R12, SA12, str + 2, n02, K);
	radixSortPass(SA12, R12, str + 1, n02, K);
	radixSortPass(R12, SA12, str, n02, K);

//S12[i] opocetak i-tog po redu tripleta

// dodjeliti tripletima lexikografska imena
	ulint name = 0, c0 = -1, c1 = -1, c2 = -1;
	for (ulint i = 0; i < n02; i++) {
		if (str[SA12[i]] != c0 || str[SA12[i] + 1] != c1 || str[SA12[i] + 2] != c2) {
			name++;
			c0 = str[SA12[i]];
			c1 = str[SA12[i] + 1];
			c2 = str[SA12[i] + 2];
		}
		if (SA12[i] % 3 == 1) {
			R12[SA12[i] / 3] = name;
		} else {
			// desna polovica
			R12[SA12[i] / 3 + n0] = name;
		}
	}

// rekurzija ako imena nisu jedinstvena
	if (name < n02) {
		constructArray(R12, SA12, n02, name + 1);
		// jedinstvena imena u R12 pomocu suffiksnog polja
		for (int i = 0; i < n02; i++)
			R12[SA12[i]] = i + 1;
	} else
		// sufiksno polje izravno iz imena
		for (int i = 0; i < n02; i++) {
			SA12[R12[i] - 1] = i;
		}

	ulint* R0 = new ulint[n0];
	ulint* SA0 = new ulint[n0];

	for (int i = 0, j = 0; i < n02; i++)
		if (SA12[i] < n0)
			R0[j++] = 3 * SA12[i];
	radixSortPass(R0, SA0, str, n0, K);

// merge SA0 i SA12
	for (ulint p = 0, t = n0 - n1, k = 0; k < n; k++) {
		//pozicija trenutnog 12 sufiksa
		int i = SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2;
		// pozicija trenutnog 0 sufiksa
		int j = SA0[p];
		if (SA12[t] < n0 ?
				leq(str[i], R12[SA12[t] + n0], str[j], R12[j / 3]) :
				leq(str[i], str[i + 1], R12[SA12[t] - n0 + 1], str[j], str[j + 1], R12[j / 3 + n0])) { // suffix from SA12 is smaller
			SA[k] = i;
			t++;
			if (t == n02) {
				// isprazni preostale 0 sufikse
				for (k++; p < n0; p++, k++) {
					SA[k] = SA0[p];
				}
			}
		} else {
			SA[k] = j;
			p++;
			if (p == n0) {
				// isprazni preostale 12 sufikse
				for (k++; t < n02; t++, k++) {
					SA[k] = SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2;
				}
			}
		}
	}
	delete[] R12;
	delete[] SA12;
	delete[] SA0;
	delete[] R0;
}

