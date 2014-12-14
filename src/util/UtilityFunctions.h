#ifndef SRC_UTIL_UTILITY_H_
#define SRC_UTIL_UTILITY_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <dirent.h>
#include <unistd.h>
#include <utility>
#include <vector>

template<typename type>
void freeVectorOfPtrs(std::vector<type>& vec);

inline size_t trimEnd(char *data) {
	int pos = strlen(data) - 1;

	while (pos >= 0 && isspace(data[pos])) {
		data[pos] = 0;
		--pos;
	}
	return pos + 1;
}

inline int lobit(const int& a) {
	return a & -a;
}
inline size_t trimBegin(char *data) {
	int cntr = 0, len = strlen(data);
	while (isspace(data[cntr])) {
		cntr++;
	}
	for (int i = 0; i < len - cntr; ++i) {
		data[i] = data[i + cntr];
	}
	data[len - cntr] = 0;
	return len - cntr;
}

inline size_t trim(char *data) {
	trimBegin(data);
	return trimEnd(data);
}

inline bool canWriteToFIle(char*filePath) {
	printf("W %s\n", filePath);
	FILE* f = fopen(filePath, "w");
	if (f != NULL) {
		fclose(f);
		return true;
	}
	return false;
}

inline bool canReadFromFile(char*filePath) {
	printf("R %s\n", filePath);
	FILE* f = fopen(filePath, "r");

	if (f != NULL) {
		fclose(f);
		return true;
	}
	return false;
}

inline bool isFolder(char*folderPath) {
	DIR* dir = opendir(folderPath);
	if (dir != NULL) {
		closedir(dir);
		return true;
	}
	return false;
}

#endif
