CC = gcc
CFLAGS = -O0 -g

all: test1 test2 quadratic


test1: test1.c
	$(CC) $(CFLAGS) -o test1 test1.c

test2: test2.c
	$(CC) $(CFLAGS) -o test2 test2.c

quadratic: quadratic.c
	$(CC) $(CFLAGS) -o quadratic quadratic.c

clean :
	rm -f test1 test2 quadratic
