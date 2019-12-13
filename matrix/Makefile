matrix: main.o matrix.o
	g++-9 -std=c++17 -Wall -O2 -o MAIN_CPP main.o matrix.o

main.o: main.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g main.cpp

matrix.o: matrix.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g matrix.cpp

run:
	./MAIN_CPP

clean:
	rm *.o MAIN_CPP
