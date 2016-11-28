# Copyright (C) 2016 Vincent Noel (vincent.noel@butantan.gov.br)
#
# plsa is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# plsa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with plsa. If not, see <http://www.gnu.org/licenses/>.

CC = gcc
MPICC = mpicc

FLAGS = -Wall -Werror -pedantic -std=gnu11
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
