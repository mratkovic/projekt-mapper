#ifndef SRC_UTIL_UTILITY_H_
#define SRC_UTIL_UTILITY_H_

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <dirent.h>
#include <unistd.h>

inline bool validateOutputFile(char*filePath) {
  printf("W %s\n", filePath);
  FILE* f = fopen(filePath, "w");
  if (f != NULL) {
    fclose(f);
    return true;
  }
  return false;
}

inline bool validateInputFile(char*filePath) {
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
