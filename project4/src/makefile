CC = g++
CFLAGS = -std=c++11

all: run 



run: main.o SIRS.o
	${CC} ${CFLAGS} main.o SIRS.o -o run



main.o: main.cpp 
	${CC} ${CFLAGS} -c main.cpp
SIRS.o: SIRS.cpp
	${CC} ${CFLAGS} -c SIRS.cpp



clean: 
	rm run
	rm *o *.dat

