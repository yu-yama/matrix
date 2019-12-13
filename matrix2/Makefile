matrix2: Makefile main.o matrix2.o
	g++-9 -std=c++17 -Wall -O2 -o MAIN_CPP main.o matrix2.o

main.o: Makefile main.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g main.cpp

matrix2.o: Makefile matrix2.h matrix2.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g matrix2.cpp

run:
	./MAIN_CPP

clean:
	rm *.o MAIN_CPP
