CC = g++
CFLAGS = -std=c++11

all: main

main: main.cpp solar_system.cpp planet.cpp
	${CC} ${CFLAGS} main.cpp solar_system.cpp planet.cpp -o run $(LIB)

clean: 
	rm run *.dat

