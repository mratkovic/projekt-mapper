export INCLUDES := -I$(CURDIR)/src/
export MAIN_DIR := $(CURDIR)/src/main/
export UTIL_DIR := $(CURDIR)/src/util/

UTIL_H_FILES := $(wildcard $(UTIL_DIR)/*.h)
UTIL_CPP_FILES := $(wildcard $(UTIL_DIR)/*.cpp)
UTIL_OBJ_FILES := $(addprefix obj/util/,$(notdir $(UTIL_CPP_FILES:.cpp=.o)))


MAIN_H_FILES := $(wildcard $(MAIN_DIR)/*.h)
MAIN_CPP_FILES := $(wildcard $(MAIN_DIR)/*.cpp)
MAIN_OBJ_FILES := $(addprefix obj/main/,$(notdir $(MAIN_CPP_FILES:.cpp=.o)))

CC := g++
LD_FLAGS := -Wall
CC_FLAGS :=  $(INCLUDES)

all: mapper 

forceall: clean all

mapper: $(MAIN_OBJ_FILES)  $(UTIL_OBJ_FILES)
	$(CC) $(CC_FLAGS) -o $@ obj/util/*.o  obj/main/*.o  $(LD_FLAGS)


$(UTIL_OBJ_FILES): $(UTIL_CPP_FILES) $(UTIL_H_FILES)
	mkdir -p obj/util
	$(CC) $(CC_FLAGS) -c -o $@ $(UTIL_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))


$(MAIN_OBJ_FILES): $(MAIN_CPP_FILES) $(MAIN_H_FILES)
	mkdir -p obj/main
	$(CC) $(CC_FLAGS) -c -o $@ $(MAIN_DIR)/$(notdir $(patsubst %.o, %.cpp, $@))


clean:
	rm -rf obj/
	rm -f mapper
