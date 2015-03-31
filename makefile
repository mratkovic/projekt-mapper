PROJECT := mapper

CURDIR := .

INCLUDES := -I$(CURDIR)/src/

SRC_DIR := $(CURDIR)/src/
EXTERNAL_DIR := $(CURDIR)/src/external/
CORE_DIR := $(CURDIR)/src/core/
UTIL_DIR := $(CURDIR)/src/util/
BIOUTIL_DIR := $(CURDIR)/src/bioutil/


SRC_OBJ_DIR := $(CURDIR)/obj/src/
EXTERNAL_OBJ_DIR := $(CURDIR)/obj/external/
CORE_OBJ_DIR := $(CURDIR)/obj/core/
UTIL_OBJ_DIR := $(CURDIR)/obj/util/
BIOUTIL_OBJ_DIR := $(CURDIR)/obj/bioutil/


EXTERNAL_H_FILES := $(wildcard $(EXTERNAL_DIR)/*.h)
EXTERNAL_CPP_FILES := $(wildcard $(EXTERNAL_DIR)/*.cpp)
EXTERNAL_OBJ_FILES := $(addprefix $(EXTERNAL_OBJ_DIR),$(notdir $(EXTERNAL_CPP_FILES:.cpp=.o)))

SRC_H_FILES := $(wildcard $(SRC_DIR)/*.h)
SRC_CPP_FILES := $(wildcard $(SRC_DIR)/*.cpp)
SRC_OBJ_FILES := $(addprefix $(SRC_OBJ_DIR),$(notdir $(SRC_CPP_FILES:.cpp=.o)))


CORE_H_FILES := $(wildcard $(CORE_DIR)/*.h)
CORE_CPP_FILES := $(wildcard $(CORE_DIR)/*.cpp)
CORE_OBJ_FILES := $(addprefix $(CORE_OBJ_DIR),$(notdir $(CORE_CPP_FILES:.cpp=.o)))

UTIL_H_FILES := $(wildcard $(UTIL_DIR)/*.h)
UTIL_CPP_FILES := $(wildcard $(UTIL_DIR)/*.cpp)
UTIL_OBJ_FILES := $(addprefix $(UTIL_OBJ_DIR),$(notdir $(UTIL_CPP_FILES:.cpp=.o)))

BIOUTIL_H_FILES := $(wildcard $(BIOUTIL_DIR)/*.h)
BIOUTIL_CPP_FILES := $(wildcard $(BIOUTIL_DIR)/*.cpp)
BIOUTIL_OBJ_FILES := $(addprefix $(BIOUTIL_OBJ_DIR),$(notdir $(BIOUTIL_CPP_FILES:.cpp=.o)))

CXX := g++
CXXFLAGS := -g -Wall -O2 -std=c++11 -Wno-unused-function -fopenmp

CC := $(CXX)
CC_FLAGS := $(CXXFLAGS)

LD_FLAGS := -lz -fopenmp -ldivsufsort
LD_LIBS := $(INCLUDES)



all: $(PROJECT)

$(PROJECT): $(BIOUTIL_OBJ_FILES) $(CORE_OBJ_FILES) $(UTIL_OBJ_FILES)  $(SRC_OBJ_FILES) $(EXTERNAL_OBJ_FILES)
	$(CC) -o $@ $^  $(LD_FLAGS) $(LD_LIBS)

$(SRC_OBJ_FILES): $(SRC_CPP_FILES) $(SRC_H_FILES)
	@mkdir -p $(SRC_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(SRC_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(CORE_OBJ_FILES): $(CORE_CPP_FILES) $(CORE_H_FILES)
	@mkdir -p $(CORE_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(CORE_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(UTIL_OBJ_FILES): $(UTIL_CPP_FILES) $(UTIL_H_FILES)
	@mkdir -p $(UTIL_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(UTIL_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(BIOUTIL_OBJ_FILES): $(BIOUTIL_CPP_FILES) $(BIOUTIL_H_FILES)
	@mkdir -p $(BIOUTIL_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(BIOUTIL_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(EXTERNAL_OBJ_FILES): $(EXTERNAL_CPP_FILES) $(EXTERNAL_H_FILES)
	@mkdir -p $(EXTERNAL_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(EXTERNAL_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

run: $(PROJECT)
	./$(PROJECT) $(COMMANDLINE_OPTIONS)

.PHONY: clean

clean:
	rm -rf obj/
	rm -f $(PROJECT)
