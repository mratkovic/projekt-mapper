#Custom make file

PROJECT := mapper

CURDIR := .

INCLUDES := -I$(CURDIR)/src/

SRC_DIR := $(CURDIR)/src/
EXTERNAL_DIR := $(CURDIR)/src/external/
CORE_DIR := $(CURDIR)/src/core/
UTIL_DIR := $(CURDIR)/src/util/
BIOINF_DIR := $(CURDIR)/src/bioinf/
METRICS_DIR := $(CURDIR)/src/metrics_algorithm/


SRC_OBJ_DIR := $(CURDIR)/obj/src/
EXTERNAL_OBJ_DIR := $(CURDIR)/obj/external/
CORE_OBJ_DIR := $(CURDIR)/obj/core/
UTIL_OBJ_DIR := $(CURDIR)/obj/util/
BIOINF_OBJ_DIR := $(CURDIR)/obj/bioutil/
METRICS_OBJ_DIR := $(CURDIR)/obj/metrics_algorithm/


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

BIOINF_H_FILES := $(wildcard $(BIOINF_DIR)/*.h)
BIOINF_CPP_FILES := $(wildcard $(BIOINF_DIR)/*.cpp)
BIOINF_OBJ_FILES := $(addprefix $(BIOINF_OBJ_DIR),$(notdir $(BIOINF_CPP_FILES:.cpp=.o)))


METRICS_H_FILES := $(wildcard $(METRICS_DIR)/*.h)
METRICS_CPP_FILES := $(wildcard $(METRICS_DIR)/*.cpp)
METRICS_OBJ_FILES := $(addprefix $(METRICS_OBJ_DIR),$(notdir $(METRICS_CPP_FILES:.cpp=.o)))

CXX := g++
CXXFLAGS := -g -Wall -O2 -std=c++11 -Wno-unused-function -fopenmp $(INCLUDES) #-w

CC := $(CXX)
CC_FLAGS := $(CXXFLAGS)

LD_FLAGS := -lz -fopenmp -ldivsufsort
LD_LIBS := $(INCLUDES)



all: $(PROJECT)

$(PROJECT): $(BIOINF_OBJ_FILES) $(METRICS_OBJ_FILES) $(CORE_OBJ_FILES) $(UTIL_OBJ_FILES)  $(SRC_OBJ_FILES) $(EXTERNAL_OBJ_FILES)
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

$(BIOINF_OBJ_FILES): $(BIOINF_CPP_FILES) $(BIOINF_H_FILES)
	@mkdir -p $(BIOINF_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(BIOINF_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))


$(METRICS_OBJ_FILES): $(METRICS_CPP_FILES) $(METRICS_H_FILES)
	@mkdir -p $(METRICS_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(METRICS_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(EXTERNAL_OBJ_FILES): $(EXTERNAL_CPP_FILES) $(EXTERNAL_H_FILES)
	@mkdir -p $(EXTERNAL_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(EXTERNAL_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

run: $(PROJECT)
	./$(PROJECT) $(COMMANDLINE_OPTIONS)

.PHONY: clean

clean:
	rm -rf obj/
	rm -f $(PROJECT)
