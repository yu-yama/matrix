matrix2dv: main.o matrix2dv.o
	g++-9 -std=c++17 -Wall -O2 -o MAIN_CPP main.o matrix2dv.o

main.o: main.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g main.cpp

matrix2dv.o: matrix2dv.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g matrix2dv.cpp

run:
	./MAIN_CPP

clean:
	rm *.o MAIN_CPP
