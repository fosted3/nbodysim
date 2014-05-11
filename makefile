CC=g++
CFLAGS=-c -Wall -g -O3 -std=c++0x -Wextra
LDFLAGS=-lpthread
EXECUTABLE=bin/nbodysim

all: $(EXECUTABLE)

$(EXECUTABLE): build/main.o build/particle.o build/octree.o build/vector.o
	$(CC) build/main.o build/particle.o build/octree.o build/vector.o -o $(EXECUTABLE) $(LDFLAGS)

build/main.o: src/main.cpp
	$(CC) $(CFLAGS) src/main.cpp -o build/main.o

build/particle.o: src/particle.cpp
	$(CC) $(CFLAGS) src/particle.cpp -o build/particle.o

build/octree.o: src/octree.cpp
	$(CC) $(CFLAGS) src/octree.cpp -o build/octree.o

build/vector.o: src/vector.cpp
	$(CC) $(CFLAGS) src/vector.cpp -o build/vector.o

clean:
	rm -rf build/*.o bin/nbodysim
