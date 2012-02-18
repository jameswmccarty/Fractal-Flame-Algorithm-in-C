#Makefile for fractal generator

CC = gcc
CFLAGS = -ansi -pedantic -Werror -Wall -O3 -D_GNU_SOURCE
LFLAGS = -lm -ltiff -lpthread

all: fractal

fractal: fractal.c fractal.h
	$(CC) $(CFLAGS) fractal.c -o fractal $(LFLAGS)

profile: fractal.c fractal.h
	$(CC) $(CFLAGS) -pg fractal.c -o fractal $(LFLAGS)

debug: fractal.c fractal.h
	$(CC) -ansi -pedantic -Werror -Wall -g -O0 -D_GNU_SOURCE fractal.c -o fractal $(LFLAGS)

clean:
	\rm fractal
