.PHONY: all clean check tests gcov_report dvi install uninstall dist
XX=g++

FILES= s21_graph.a \
	   s21_graph_algorithms.a \
	   console_interface.cpp \
	   main.cpp \

GCOVFILES= s21_graph.cpp \
	   s21_graph_algorithms.cpp \

LDFLAGS= -Wall -Wextra -Werror -std=c++17 -g
OS := $(shell uname -s)
TEST_FLAGS= -lgtest -lgtest_main -lpthread
CFLAGS= -fprofile-arcs -ftest-coverage
BD= build/

ifeq ($(OS), Darwin)
	CC += -D OS_MAC
else
	CC += -D OS_LINUX
endif

all: build

build: clean s21_graph.a s21_graph_algorithms.a
	$(XX) $(LDFLAGS) $(FILES) -o navigator.out

s21_graph.a:
	$(XX) -c s21_graph.cpp
	ar rcs s21_graph.a *.o
	ranlib $@
	cp $@ lib$@

s21_graph_algorithms.a:
	$(XX) -c s21_graph_algorithms.cpp
	ar rcs s21_graph_algorithms.a *.o
	ranlib $@
	cp $@ lib$@

clean:
	rm -rf build test *.a *.o *.cfg *.gcda *.gcno *.html *.out *.info *.dSYM report

check: 
	cp ../materials/linters/CPPLINT.cfg CPPLINT.cfg
	cp ../materials/linters/suppressions.cfg suppressions.cfg
	python3 ../materials/linters/cpplint.py *.cpp *.h
	cppcheck -q --std=c++17 --enable=all --language=c++ --suppressions-list=suppressions.cfg *.cpp *.h

gcov_report: clean
	$(XX) $(FLAGS) $(CFLAGS) $(GCOVFILES) tests.cpp -lstdc++ $(TEST_FLAGS) -o test
	./test
	
	lcov --no-external -c -d . -o tests.info
	genhtml tests.info -o report
ifeq ($(OS), Darwin)
	open -a "Google Chrome" report/index.html
else
	xdg-open report/index.html
endif

tests: test

test: clean s21_graph.a s21_graph_algorithms.a
	$(XX) $(LDFLAGS) $(FLAGS) tests.cpp s21_graph.a s21_graph_algorithms.a -lstdc++ $(TEST_FLAGS) -o test.out
	./test.out

leaks: clean tests
	leaks -atExit -- ./test.out

valgrind:
	valgrind --leak-check=full  --show-leak-kinds=all --track-origins=yes --verbose --log-file=valgrind.txt ./test.out
