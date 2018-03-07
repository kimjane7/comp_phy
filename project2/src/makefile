CC = g++
CFLAGS = -std=c++11
LIB = -larmadillo
PROG = testcode

all: toeplitz harmonic quantum ${PROG}

toeplitz: toeplitz.cpp jacobi.cpp
	${CC} ${CFLAGS} toeplitz.cpp jacobi.cpp -o toeplitz $(LIB)

harmonic: harmonic_oscillator.cpp jacobi.cpp
	${CC} ${CFLAGS} harmonic_oscillator.cpp jacobi.cpp -o harmonic $(LIB)

quantum: quantum_dots.cpp jacobi.cpp
	${CC} ${CFLAGS} quantum_dots.cpp jacobi.cpp -o quantum $(LIB)

${PROG}: test_main.o test_jacobi.o jacobi.o
	${CC} ${CFLAGS} test_main.o test_jacobi.o jacobi.o -o ${PROG} ${LIB}

test_main.o: test_main.cpp
	${CC} ${CFLAGS} -c test_main.cpp

test_jacobi.o: test_jacobi.cpp
	${CC} ${CFLAGS} -c test_jacobi.cpp

jacobi.o: jacobi.cpp
	${CC} ${CFLAGS} -c jacobi.cpp

clean: 
	rm toeplitz harmonic quantum testcode