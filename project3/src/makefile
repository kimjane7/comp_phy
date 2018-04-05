CC = g++
CFLAGS = -std=c++11

all: run testcode



run: main.o solar_system.o planet.o
	${CC} ${CFLAGS} main.o solar_system.o planet.o -o run
testcode: test_main.o solar_system.o planet.o
	${CC} ${CFLAGS} test_main.o solar_system.o planet.o -o testcode



main.o: main.cpp 
	${CC} ${CFLAGS} -c main.cpp
test_main.o: test_main.cpp
	${CC} ${CFLAGS} -c test_main.cpp
solar_system.o: solar_system.cpp
	${CC} ${CFLAGS} -c solar_system.cpp
planet.o: planet.cpp
	${CC} ${CFLAGS} -c planet.cpp



clean: 
	rm run testcode
	rm *o *.dat
	rm ../report/report.aux ../report/report.log 
	rm ../report/reportNotes.bib ../report/report.out 
	rm ../report/report.synctex.gz

