CC = gcc
MPICC = mpicc

FLAGS = -Wall -Werror -pedantic -std=gnu99
LIBS = -lm

all: run-serial run-parallel

run-serial:	main-serial.o
	$(CC) main-serial.o -lplsa-serial $(LIBS) -o run-serial

run-parallel:	main-parallel.o
	$(MPICC) main-parallel.o -lplsa-parallel $(LIBS) -o run-parallel

main-serial.o:	main.c
	$(CC) $(FLAGS) -c main.c -o main-serial.o

main-parallel.o: main.c
	$(MPICC) $(FLAGS) -c main.c -DMPI -o main-parallel.o

clean:
	rm -fr  *.o run-serial run-parallel final_score plsa.log *.state logs/
