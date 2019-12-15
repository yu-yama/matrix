test: Makefile matrix.o test.o
	g++-9 -std=c++17 -Wall -O2 -o TEST_CPP test.o matrix.o

matrix.o: Makefile matrix.h matrix.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g matrix.cpp

test.o: Makefile test.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g test.cpp

run:
	./TEST_CPP

clean:
	rm *.o TEST_CPP
