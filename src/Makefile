# MIT License
#
# Copyright (c) 2019 yu-yama
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
