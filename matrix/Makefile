matrix: Makefile matrix.o main.o
	g++-9 -std=c++17 -Wall -O2 -o MAIN_CPP main.o matrix.o

matrix.o: Makefile matrix.h matrix.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g matrix.cpp

main.o: Makefile main.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g main.cpp

run:
	./MAIN_CPP

clean:
	rm *.o MAIN_CPP
