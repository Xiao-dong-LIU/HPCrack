CXX      := -mpicxx
CXXFLAGS := -O3 -fopenmp -std=c++11 -std=c++17 
#CXXFLAGS := -pedantic-errors -Wall -Wextra -Werror -O3 -fopenmp -std=c++11
LDFLAGS  := -L/usr/lib -lstdc++ -llapack -lblas 
BUILD    := ./build
OBJ_DIR  := $(BUILD)/objects
TARGET   := program
INCLUDE  := -I ./include
SRC      := $(wildcard src/*.cpp)        

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)
DEPENDS := $(SRC:%.cpp=$(OBJ_DIR)/%.d)

all: $(TARGET)

-include $(DEPENDS)

$(OBJ_DIR)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS)  $(INCLUDE)  -MD -MP -o $@ -c $<

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS)  $(INCLUDE) -o $(TARGET) $(OBJECTS)

.PHONY: all build clean debug release

build:
	@mkdir -p $(OBJ_DIR)


debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(TARGET)

