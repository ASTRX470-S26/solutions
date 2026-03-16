CC = gcc
CFLAGS = -O3
LDLIBS = -lm
HEAD = methods.h

all: rootfind

# this says look for any .c file and compile it to a .o file
# of the same name.  This will rerun if any of the .c files
# have been edited
%.o: %.c $(HEAD)
	$(CC) -c -o $@ $< $(CFLAGS)

# this takes the .o files and links them to generate the executable
# with the name rootfind
rootfind: main.o bisection.o secant.o newton.o
	$(CC) -o rootfind main.o bisection.o secant.o newton.o $(LDLIBS)

clean :
	rm -f *.o
	rm -f rootfind
