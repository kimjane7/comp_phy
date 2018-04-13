CC = g++
CFLAGS = -std=c++11


all: run 



run: main.o SIRS.o
	${CC} ${CFLAGS} main.o SIRS.o -o run ${LIB}



main.o: main.cpp 
	${CC} ${CFLAGS} -c main.cpp
SIRS.o: SIRS.cpp
	${CC} ${CFLAGS} -c SIRS.cpp ${LIB}



clean: 
	rm run
	rm *o *.dat

