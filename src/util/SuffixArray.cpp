/*
 * SuffixArray.cpp
 *
 *  Created on: Dec 12, 2014
 *      Author: marko
 */

#include "SuffixArray.h"
#include "BioUtil.h"
#include <cstdio>
#include <cstring>
#include <algorithm>

#define NUM_OF_SYMBOLS 5

inline bool leq(int a1, int a2, int b1, int b2) {
	return (a1 < b1 || (a1 == b1 && a2 <= b2));
}
inline bool leq(int a1, int a2, int a3, int b1, int b2, int b3) {
	return (a1 < b1 || (a1 == b1 && leq(a2, a3, b2, b3)));
}

int SuffixArray::compare(char *pattern, size_t patternLen, int startPos) {
	for (int i = startPos; i < std::min<int>(this->_size, (startPos + patternLen)); ++i) {
		int diff;

//		if(complement) {
//			diff = this->_sequence[i] - getACGTComplement(pattern[i - startPos]);
//		} else {
		diff = this->_sequence[i] - pattern[i - startPos];
		if (diff != 0) {
			return diff;
		}
	}
	return 0;
}

int SuffixArray::findStartingPositions(char *pattern, size_t patternLen, int id,
		std::vector<std::pair<int, int> > &dest) {

	int lo, hi, mid;
	lo = 0;
	hi = this->_size;

	while (lo + 1 < hi) {
		mid = (lo + hi) / 2;
		int cmp = this->compare(pattern, patternLen, _array[mid]);
		if (cmp == 0) {
			dest.push_back(std::make_pair(this->_array[mid], id));
			for (lo = mid; lo > 0 && this->compare(pattern, patternLen, _array[lo - 1]) == 0; lo--) {
				dest.push_back(std::make_pair(this->_array[lo - 1], id));
			}
			for (hi = mid; hi < this->_size && this->compare(pattern, patternLen, _array[hi + 1]) == 0; hi++) {
				dest.push_back(std::make_pair(this->_array[hi + 1], id));
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
//
//int SuffixArray::findStartingPositions(char *pattern, size_t patternLen, int id,
//		std::vector<std::pair<int, int> > &dest) {
//	int lo, hi, mid;
//	lo = 0;
//	hi = this->_size;
//
//	while (lo + 1 < hi) {
//		mid = (lo + hi) / 2;
//		int cmp = this->compare(pattern, patternLen, _array[mid]);
//		if (cmp == 0) {
//			dest.push_back(std::make_pair(this->_array[mid], id));
//			for (lo = mid; lo > 0 && this->compare(pattern, patternLen, _array[lo - 1]) == 0; lo--) {
//				dest.push_back(std::make_pair(this->_array[lo - 1], id));
//			}
//			for (hi = mid; hi < this->_size && this->compare(pattern, patternLen, _array[hi + 1]) == 0; hi++) {
//				dest.push_back(std::make_pair(this->_array[hi + 1], id));
//			}
//
//			return hi - lo + 1;
//		} else if (cmp < 0) {
//			lo = mid;
//		} else {
//			hi = mid;
//		}
//	}
//	return 0;
//}
void SuffixArray::constructFromSequence(bioutil::Sequence* seq) {
	if (seq->dataSize() % 3 == 1) {
		this->_size = seq->dataSize() + 2;
	} else if (seq->dataSize() % 3 == 2) {
		this->_size = seq->dataSize() + 1;
	} else {
		this->_size = seq->dataSize();
	}

	this->_array = (int *) malloc(this->_size * sizeof(int));
	seq->turnBaseToInt();
	this->_sequence = seq->data();
	construct();
}

bool SuffixArray::saveSuffixArray(FILE *file) {
	int cntr = 0;

	cntr += fwrite(&_size, sizeof(int), 1, file);
	cntr += fwrite(&_symbolNum, sizeof(int), 1, file);
	cntr += fwrite(_array, sizeof(int), this->_size, file);

	return cntr == this->_size + 2;

}
SuffixArray::SuffixArray(FILE *f) {
	int size;
	fread(&size, sizeof(int), 1, f);
	this->_size = size;

	int numOfSym;
	fread(&numOfSym, sizeof(int), 1, f);
	this->_symbolNum = numOfSym;

	this->_array = (int *) malloc(_size * sizeof(int));
	fread(_array, sizeof(int), _size, f);
	this->_sequence = 0;
}

bool SuffixArray::construct() {
	int *tmp = (int *) malloc(this->_size * sizeof(int));
	for (int i = 0; i < this->_size; ++i) {
		tmp[i] = this->_sequence[i];
	}
	SuffixArray::constructArray(tmp, this->_array, this->_size, NUM_OF_SYMBOLS);
	free(tmp);
}

void SuffixArray::radixSortPass(int* in, int* out, int* sequence, int n, int numOfSymbols) {
	int* cntr = new int[numOfSymbols];
	for (int i = 0; i < numOfSymbols; i++) {
		cntr[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		cntr[sequence[in[i]]]++;
	}

	for (int i = 0, sum = 0; i < numOfSymbols; i++) {
		int t = cntr[i];
		cntr[i] = sum;
		sum += t;
	}
	for (int i = 0; i < n; i++) {
		out[cntr[sequence[in[i]]]++] = in[i];
	}
	for (int i = 0; i < n; i++) {
	}
	delete[] cntr;
}

void SuffixArray::constructArray(int* str, int* SA, int n, int K) {
	int n0 = (n + 2) / 3;
	int n1 = (n + 1) / 3;
	int n2 = n / 3;

	int n02 = n0 + n2;

	int* R12 = new int[n02 + 3];
	R12[n02] = R12[n02 + 1] = R12[n02 + 2] = 0;

	int* SA12 = new int[n02 + 3];
	SA12[n02] = SA12[n02 + 1] = SA12[n02 + 2] = 0;

	// R12 indeksi pocetaka mod 1 i mod 2 sufiksa
	// +(n0-n1) dodaje jedan viska mod 1 sufiks ako je n%3 == 1
	for (int i = 0, j = 0; i < n + (n0 - n1); i++) {
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
	int name = 0, c0 = -1, c1 = -1, c2 = -1;
	for (int i = 0; i < n02; i++) {
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

	int* R0 = new int[n0 + 3];
	int* SA0 = new int[n0 + 3];

	for (int i = 0, j = 0; i < n02; i++)
		if (SA12[i] < n0)
			R0[j++] = 3 * SA12[i];
	radixSortPass(R0, SA0, str, n0, K);

	// merge SA0 i SA12
	for (int p = 0, t = n0 - n1, k = 0; k < n; k++) {
		//pozicija trenutnog 12 sufiksa
		int i = SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2;
		// pozicija trenutnog 0 sufiksa
		int j = SA0[p];
		if (SA12[t] < n0 ?
				leq(str[i], R12[SA12[t] + n0], str[j], R12[j / 3]) :
				leq(str[i], str[i + 1], R12[SA12[t] - n0 + 1], str[j], str[j + 1], R12[j / 3 + n0])) {
			// suffix od SA12 manji
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

