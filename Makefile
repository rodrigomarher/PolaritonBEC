detected_OS := $(shell uname)
use_MPI := $(USE_MPI)

ifeq ($(detected_OS),Linux)
	CC = g++
endif
ifeq ($(detected_OS), Darwin)
	CC = g++-12
endif

ifeq ($(use_MPI), Yes)
	CC =  mpic++
endif

CFLAGS = -lm -O3 -fopenmp -ffast-math -g -Wall

SRC = ./src/

INCLUDE = ./include/
BUILD = ./build/
LIB = ./lib/

SOURCES := $(wildcard $(SRC)*.cpp)
OBJECTS := $(patsubst $(SRC)%.cpp,$(BUILD)%.o,$(SOURCES))

all: PolaritonBEC

test: $(SRC)/test.cpp
	$(CC) -I$(INCLUDE) $(CFLAGS) $(SRC)/test.cpp -o test

PolaritonBEC: $(LIB)libpolaritonbec.so
	mkdir -p results
	$(CC) -L$(LIB) -lpolaritonbec -I$(INCLUDE) $(CFLAGS) $(SRC)/main.cpp -o PolaritonBEC

$(LIB)libpolaritonbec.so: $(OBJECTS)
	mkdir -p $(LIB)
	$(CC) -shared $^ -fopenmp -o $(LIB)libpolaritonbec.so

$(BUILD)%.o: $(SRC)/%.cpp
	mkdir -p ./build/
	$(CC) -c -fPIC -I$(INCLUDE) $(CFLAGS) $< -o $@

clean:
	rm -rf results/
	rm -rf PolaritonBEC
	rm -rf output_*
	rm -rf slurm*
	rm -rf $(LIB)
	rm -rf $(BUILD)

