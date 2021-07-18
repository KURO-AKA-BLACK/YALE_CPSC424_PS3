# NOTE: This Makefile presumes that the following command has been executed:

#                module load iomkl

TIMINGDIR = /home/cpsc424_ahs3/shared/utils/timing
CC = mpicc
CFLAGS = -g -O3 -xHost -fno-alias -std=c99 -I$(TIMINGDIR)
EXECUTABLES = serial task5 task6 task7 task8

serial:	serial.o matmul.o $(TIMINGDIR)/timing.o
	icc -o $@ $(CFLAGS) $^
	
task5:	task5.o $(TIMINGDIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^
	
task6:	task6.o $(TIMINGDIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^

task7:	task7.o $(TIMINGDIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^

task8:	task8.o $(TIMINGDIR)/timing.o
	$(CC) -o $@ $(CFLAGS) $^
	
.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(EXECUTABLES) *.o
	