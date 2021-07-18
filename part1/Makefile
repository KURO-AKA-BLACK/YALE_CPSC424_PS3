#
# This Makefile assumes that the following module files are loaded:
#
#   iomkl
#
# Set the compiler command
CC = mpicc

# Set the compiler options
CFLAGS = -g -O3 -xHost -fno-alias -std=c99

all: step1 

task1: task1.o timing.o rwork.o
	$(CC) -o $@ $(CFLAGS) $^
	
task2: task2.o timing.o rwork.o
	$(CC) -o $@ $(CFLAGS) $^
	
task3: task3.o timing.o rwork.o
	$(CC) -o $@ $(CFLAGS) $^

.c.o:
	$(CC) $(CFLAGS) -c $<

# Be careful not to remove rwork.o and timing.o
clean:
	rm -f task1 task1.o task2 task2.o task3 task3.o

