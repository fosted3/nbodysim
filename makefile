CC=g++
CFLAGS=-c -Wall -g -O3 -std=c++11 -Wextra -march=native -mtune=native
LDFLAGS=-lpthread -lfreeimage
EXECUTABLE=bin/nbodysim
DIRS=bin/ build/

all: $(EXECUTABLE)

$(EXECUTABLE): build/main.o build/particle.o build/octree.o build/vector.o build/thread_functions.o
	mkdir -p $(DIRS)
	$(CC) build/main.o build/particle.o build/octree.o build/vector.o build/thread_functions.o -o $(EXECUTABLE) $(LDFLAGS)

build/main.o: src/main.cpp
	mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/main.cpp -o build/main.o

build/particle.o: src/particle.cpp
	mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/particle.cpp -o build/particle.o

build/octree.o: src/octree.cpp
	mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/octree.cpp -o build/octree.o

build/vector.o: src/vector.cpp
	mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/vector.cpp -o build/vector.o

build/thread_functions.o: src/thread_functions.cpp
	mkdir -p $(DIRS)
	$(CC) $(CFLAGS) src/thread_functions.cpp -o build/thread_functions.o

clean:
	rm -rf build/*.o bin/nbodysim
