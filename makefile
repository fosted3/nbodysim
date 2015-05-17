CC=g++
NVCC=nvcc
CFLAGS=-c -Wall -std=c++11 -Wextra -march=native -mtune=native -DFLOAT
CUDA_CFLAGS=-c -DFLOAT -m64 -arch=compute_30 -code=sm_30
LDFLAGS=-lpthread -lfreeimage
CUDA_LDFLAGS=-lcuda -lcudart
EXECUTABLE=nbodysim
DIRS=bin/debug/ bin/release/  build/debug/ build/release/ data/ img/
CUDA_INCLUDES=-I/opt/cuda/include

all:
	make debug
	make release

debug: CFLAGS += -g -Og
debug: CUDA_CFLAGS += -g -G -O0
debug: BIN_DIR=bin/debug
debug: BUILD_DIR=build/debug
debug: setup
debug: $(EXECUTABLE)

release: CFLAGS += -O2
release: CUDA_CFLAGS += -O2
release: BIN_DIR=bin/release
release: BUILD_DIR=build/release
release: setup
release: $(EXECUTABLE)

$(EXECUTABLE): $(BUILD_DIR)/main.o $(BUILD_DIR)/particle.o $(BUILD_DIR)/octree.o $(BUILD_DIR)/vector.o $(BUILD_DIR)/thread_functions.o $(BUILD_DIR)/cuda_code.o $(BUILD_DIR)/cuda_helper.o
	$(NVCC) $(BUILD_DIR)/main.o $(BUILD_DIR)/particle.o $(BUILD_DIR)/octree.o $(BUILD_DIR)/vector.o $(BUILD_DIR)/thread_functions.o $(BUILD_DIR)/cuda_code.o $(BUILD_DIR)/cuda_helper.o -o $(BIN_DIR)/$(EXECUTABLE) $(CUDA_LDFLAGS) $(LDFLAGS)

$(BUILD_DIR)/main.o: src/main.cpp
	$(CC) $(CFLAGS) src/main.cpp -o $(BUILD_DIR)/main.o

$(BUILD_DIR)/particle.o: src/particle.cpp
	$(CC) $(CFLAGS) src/particle.cpp -o $(BUILD_DIR)/particle.o

$(BUILD_DIR)/octree.o: src/octree.cpp
	$(CC) $(CFLAGS) src/octree.cpp -o $(BUILD_DIR)/octree.o

$(BUILD_DIR)/vector.o: src/vector.cpp
	$(CC) $(CFLAGS) src/vector.cpp -o $(BUILD_DIR)/vector.o

$(BUILD_DIR)/thread_functions.o: src/thread_functions.cpp
	$(CC) $(CFLAGS) src/thread_functions.cpp -o $(BUILD_DIR)/thread_functions.o

$(BUILD_DIR)/cuda_code.o: src/cuda_code.cu
	$(NVCC) $(CUDA_CFLAGS) src/cuda_code.cu -o $(BUILD_DIR)/cuda_code.o $(CUDA_INCLUDES)

$(BUILD_DIR)/cuda_helper.o: src/cuda_helper.cpp
	$(CC) $(CFLAGS) src/cuda_helper.cpp -o $(BUILD_DIR)/cuda_helper.o $(CUDA_INCLUDES)

setup:
	mkdir -p $(DIRS)

clean:
	rm -f build/debug/*.o build/release/*.o bin/debug/* bin/release/*

clean_data:
	rm -f data/* img/*

test:
	make debug
	valgrind --leak-check=full --show-reachable=yes --read-var-info=yes --track-origins=yes ./bin/debug/nbodysim
