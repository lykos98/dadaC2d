LIBRARIES=-lm -fopenmp 
OPTIM=-O4 -mavx2 -march=native  -Wall -Wextra -DUSE_NORM
DEBUG=
SRC="src"
VERBOSE=-DVERBOSE

CC=gcc

all: lib 

lib: bin/libdadac.so

bin/libdadac.so: bin/dadac.o bin/kdtree.o bin/heap.o bin/vptree.o
	${CC} -shared bin/dadac.o bin/kdtree.o bin/heap.o bin/vptree.o ${DEBUG} ${OPTIM} ${LIBRARIES} -o bin/libdadac.so 

bin/dadac.o: src/dadac.c
	${CC} -c src/dadac.c -o bin/dadac.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

bin/kdtree.o: src/kdtree.c
	${CC} -c src/kdtree.c -o bin/kdtree.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

bin/vptree.o: src/vptree.c
	${CC} -c src/vptree.c -o bin/vptree.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}

bin/heap.o: src/heap.c
	${CC} -c src/heap.c -o bin/heap.o ${OPTIM} -fopenmp ${DEBUG} -fpic ${VERBOSE}


clean:
	rm bin/*.so bin/*.o *.o driver test
