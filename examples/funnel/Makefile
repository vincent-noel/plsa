all:	run-serial run-parallel

run-serial:	main-serial.o
	gcc main-serial.o -lm -lplsa-serial -o run-serial

run-parallel:	main-parallel.o
	mpicc main-parallel.o -lm -lplsa-parallel -o run-parallel

main-serial.o:	main.c
	gcc -c main.c -o main-serial.o

main-parallel.o: main.c
	mpicc -c main.c -DMPI -o main-parallel.o

clean: 
	rm -fr  *.o run-serial run-parallel pid final_score plsa.log *.state logs/

