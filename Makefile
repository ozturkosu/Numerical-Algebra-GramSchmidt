
CC = cc

#Note: Change -O to -g when using debugger.
#CFLAGS = -g -Wall
CFLAGS = -O -Wall

all:	hw2_6

hw2_6.o: hw2_6.c
	$(CC) $(CFLAGS) -c $<

hw2_6: hw2_6.o
	$(CC) $< -o $@ -lm

plot: plotresults.m results.m
	octave -q -f --eval plotresults
	epstopdf errors.eps
	epstopdf times.eps

clean:
	rm -f hw2_6.o hw2_6
