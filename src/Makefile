CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c++17 -I h -I /usr/local/include/gtest/ -c
LFLAGS = -std=c++17 -I h -pthread -lstdc++


all: test

clean:
	rm -rf s21_matrix_oop.o tests.o s21_matrix_oop.a s21_matrix_oop.out

test: clean objects s21_matrix_oop.a
	$(CC) $(LFLAGS) -o s21_matrix_oop.out s21_matrix_oop.a tests.o /usr/local/lib/libgtest.a
	./s21_matrix_oop.out

objects:
	${CC} ${CFLAGS} tests.cpp -o tests.o

s21_matrix_oop.a:
	$(CC) ${CFLAGS} s21_matrix_oop.cpp -o s21_matrix_oop.o
	ar rcs s21_matrix_oop.a s21_matrix_oop.o
	ranlib s21_matrix_oop.a

git:
	git status 
	git add s21_matrix_oop.cpp s21_matrix_oop.h test.cpp Makefile tests.cpp
	git commit -m "$m"
	git push origin develop
	git status 

check:
	cp ../materials/linters/.clang-format .
	clang-format -n *.cpp *.h



