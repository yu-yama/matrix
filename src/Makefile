test: Makefile matrix.o test.o
	g++-9 -std=c++17 -Wall -O2 -o TEST_CPP test.o matrix.o || (echo "Compilation Failed!"; exit 1)

matrix.o: Makefile matrix.h matrix.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g matrix.cpp || (echo "Failed to build matrix.o!"; exit 1)

test.o: Makefile matrix.h matrix.cpp test.cpp
	g++-9 -std=c++17 -Wall -O2 -c -g test.cpp || (echo "Failed to build test.o!"; exit 1)

test-github: matrix.o-github test.o-github
	g++ -std=c++17 -Wall -O2 -o TEST_CPP test.o matrix.o || (echo "Compilation Failed!"; exit 1)

matrix.o-github:
	g++ -std=c++17 -Wall -O2 -c -g matrix.cpp || (echo "Failed to build matrix.o!"; exit 1)

test.o-github:
	g++ -std=c++17 -Wall -O2 -c -g test.cpp || (echo "Failed to build test.o!"; exit 1)

run:
	./TEST_CPP > output.txt || (echo "Execution Failed!"; exit 1)
	diff output.txt expected_output.txt || echo "Unexpected Output!"
	rm output.txt

clean:
	rm *.o TEST_CPP
