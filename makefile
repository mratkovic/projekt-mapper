PROJECT := mapper

INCLUDES := -I$(CURDIR)/src/

SRC_DIR := $(CURDIR)/src/
EXTERNAL_DIR := $(CURDIR)/src/external/

SRC_OBJ_DIR := $(CURDIR)/obj/src/
EXTERNAL_OBJ_DIR := $(CURDIR)/obj/external/

EXTERNAL_H_FILES := $(wildcard $(EXTERNAL_DIR)/*.h)
EXTERNAL_CPP_FILES := $(wildcard $(EXTERNAL_DIR)/*.cpp)
EXTERNAL_OBJ_FILES := $(addprefix $(EXTERNAL_OBJ_DIR),$(notdir $(EXTERNAL_CPP_FILES:.cpp=.o)))

SRC_H_FILES := $(wildcard $(SRC_DIR)/*.h)
SRC_CPP_FILES := $(wildcard $(SRC_DIR)/*.cpp)
SRC_OBJ_FILES := $(addprefix $(SRC_OBJ_DIR),$(notdir $(SRC_CPP_FILES:.cpp=.o)))

CXX := g++
CC := $(CXX)
CXXFLAGS := -g -Wall -Wno-unused-function  $(INCLUDES) -fopenmp
CC_FLAGS := $(CXXFLAGS)
LD_FLAGS := -lz -fopenmp -ldivsufsort
LD_LIBS := $(INCLUDES)

all: $(PROJECT) 

$(PROJECT): $(SRC_OBJ_FILES) $(EXTERNAL_OBJ_FILES)
	$(CC) -o $@ $^  $(LD_FLAGS) $(LD_LIBS)

$(SRC_OBJ_FILES): $(SRC_CPP_FILES) $(SRC_H_FILES)
	mkdir -p $(SRC_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(SRC_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

$(EXTERNAL_OBJ_FILES): $(EXTERNAL_CPP_FILES) $(EXTERNAL_H_FILES)
	mkdir -p $(EXTERNAL_OBJ_DIR)
	$(CC) $(CC_FLAGS) -c -o $@ $(EXTERNAL_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))

run: $(PROJECT)
	./$(PROJECT) $(COMMANDLINE_OPTIONS)

.PHONY: clean

clean:
	rm -rf obj/
	rm -f $(PROJECT)
