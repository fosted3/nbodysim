CC=g++
NVCC=nvcc
CFLAGS=-c -Wall -std=c++11 -Wextra -march=native -mtune=native -DFLOAT
CUDA_CFLAGS=-c -DFLOAT -m64 -arch=compute_30 -code=sm_30
LDFLAGS=-lpthread -lfreeimage
CUDA_LDFLAGS=-lcuda -lcudart
EXECUTABLE=bin/nbodysim
DIRS=bin/ build/
CUDA_INCLUDES=-I/opt/cuda/include

all: $(EXECUTABLE)

debug: CFLAGS += -g -Og
debug: CUDA_CFLAGS += -g -G -O0
debug: $(EXECUTABLE)

release: CFLAGS += -O2
release: CUDA_CFLAGS += -O2
release: $(EXECUTABLE)

$(EXECUTABLE): build/main.o build/particle.o build/octree.o build/vector.o build/thread_functions.o build/cuda_code.o build/cuda_helper.o
	@mkdir -p $(DIRS)
	$(NVCC) build/main.o build/particle.o build/octree.o build/vector.o build/thread_functions.o build/cuda_code.o build/cuda_helper.o -o $(EXECUTABLE) $(CUDA_LDFLAGS) $(LDFLAGS)

build/main.o: src/main.cpp
	@mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/main.cpp -o build/main.o

build/particle.o: src/particle.cpp
	@mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/particle.cpp -o build/particle.o

build/octree.o: src/octree.cpp
	@mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/octree.cpp -o build/octree.o

build/vector.o: src/vector.cpp
	@mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/vector.cpp -o build/vector.o

build/thread_functions.o: src/thread_functions.cpp
	@mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/thread_functions.cpp -o build/thread_functions.o

build/cuda_code.o: src/cuda_code.cu
	@mkdir -p $(DIRS)
	$(NVCC) $(CUDA_CFLAGS) src/cuda_code.cu -o build/cuda_code.o $(CUDA_INCLUDES)

build/cuda_helper.o: src/cuda_helper.cpp
	@mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/cuda_helper.cpp -o build/cuda_helper.o $(CUDA_INCLUDES)

clean:
	rm -f build/*.o bin/*
