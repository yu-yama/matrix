matrixInt: main.o matrixInt.o
	g++-9 -std=c++17 -Wall -O2 -o MAIN_CPP main.o matrixInt.o

main.o: main.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g main.cpp

matrixInt.o: matrixInt.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g matrixInt.cpp

run:
	./MAIN_CPP

clean:
	rm *.o MAIN_CPP
