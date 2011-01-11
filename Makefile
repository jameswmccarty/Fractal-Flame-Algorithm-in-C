#Makefile for fractal generator
# Author: James McCarty

CC = gcc
CFLAGS = -ansi -pedantic -D_GNU_SOURCE -Werror -Wall -O2 
LFLAGS = -lm -ljpeg -ltiff

fractal: fractal.c fractal.h
	$(CC) $(CFLAGS) fractal.c -o fractal $(LFLAGS)

profile: fractal.c fractal.h
	$(CC) $(CFLAGS) -pg fractal.c -o fractal $(LFLAGS)

clean:
	\rm fractal
