matrixIntArray: Makefile main.o matrixIntArray.o
	g++-9 -std=c++17 -Wall -O2 -o MAIN_CPP main.o matrixIntArray.o

main.o: Makefile main.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g main.cpp

matrixIntArray.o: Makefile matrixIntArray.cpp matrixIntArray.h
	g++-9 -std=c++17 -Wall -O2 -c -g matrixIntArray.cpp

run:
	./MAIN_CPP

clean:
	rm *.o MAIN_CPP
