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

SRCDIR = src

CC = gcc
MPICC = mpicc

FLAGS = -Wall -Werror -pedantic -std=gnu99
LIBS = -lm

OBJ = config.o error.o distributions.o random.o
SOBJ = plsa.o lsa.o moves.o savestate.o score.o
POBJ = plsa_p.o lsa_p.o moves_p.o savestate_p.o score_p.o


all: libplsa-serial.so libplsa-parallel.so

install: libplsa-serial.so libplsa-parallel.so $(SRCDIR)/config.h $(SRCDIR)/global.h $(SRCDIR)/sa.h
	cp libplsa-serial.so libplsa-parallel.so /usr/lib/
	mkdir -p /usr/include/plsa
	cp $(SRCDIR)/config.h $(SRCDIR)/global.h $(SRCDIR)/sa.h	/usr/include/plsa
	rm -f *.so *.o

uninstall:
	rm -f /usr/lib/libplsa-serial.so /usr/lib/libplsa-parallel.so
	rm -fr /usr/include/plsa

clean:
	rm -f *.o *.so

examples: run-funnel-serial run-funnel-parallel
	rm *.o

clean_examples:
	rm -f *.o *.so run-funnel-serial run-funnel-parallel
	rm -fr logs/ final_score plsa.log *.state

run-funnel-serial: main-funnel-serial.o
	$(CC) main-funnel-serial.o -lplsa-serial -lm -o run-funnel-serial

main-funnel-serial.o: examples/funnel/main.c
	$(CC) $(FLAGS) -c examples/funnel/main.c -o main-funnel-serial.o

run-funnel-parallel: main-funnel-parallel.o
	$(MPICC) main-funnel-parallel.o -lplsa-parallel -lm -o run-funnel-parallel

main-funnel-parallel.o:	examples/funnel/main.c
	$(MPICC) $(FLAGS) -c examples/funnel/main.c -DMPI -o main-funnel-parallel.o


# Shared objects
libplsa-serial.so:	$(OBJ) $(SOBJ)
	$(CC) -shared -o libplsa-serial.so $(OBJ) $(SOBJ) $(LIBS)

libplsa-parallel.so: $(OBJ) $(POBJ)
	$(MPICC) -shared -o libplsa-parallel.so $(OBJ) $(POBJ) $(LIBS)


# Objects
# parallel ones
plsa_p.o: $(SRCDIR)/plsa.c
	$(MPICC) $(FLAGS) -fpic -DMPI -o plsa_p.o -c $(SRCDIR)/plsa.c

lsa_p.o: $(SRCDIR)/lsa.c
	$(MPICC) $(FLAGS) -fpic -DMPI -o lsa_p.o -c $(SRCDIR)/lsa.c

moves_p.o: $(SRCDIR)/moves.c $(SRCDIR)/moves.h
	$(MPICC) $(FLAGS) -fpic -DMPI -o moves_p.o -c $(SRCDIR)/moves.c

savestate_p.o: $(SRCDIR)/savestate.c
	$(MPICC) $(FLAGS) -fpic -DMPI -o savestate_p.o -c $(SRCDIR)/savestate.c

score_p.o: $(SRCDIR)/score.c $(SRCDIR)/score.h
	$(MPICC) $(FLAGS) -fpic -DMPI -o score_p.o -c $(SRCDIR)/score.c



# serial ones
plsa.o: $(SRCDIR)/plsa.c
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/plsa.c

lsa.o: $(SRCDIR)/lsa.c
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/lsa.c

moves.o: $(SRCDIR)/moves.c $(SRCDIR)/moves.h
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/moves.c

savestate.o: $(SRCDIR)/savestate.c
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/savestate.c

score.o: $(SRCDIR)/score.c $(SRCDIR)/score.h
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/score.c


#Other ones
config.o: $(SRCDIR)/config.c $(SRCDIR)/config.h
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/config.c

error.o: $(SRCDIR)/error.c $(SRCDIR)/error.h
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/error.c

distributions.o: $(SRCDIR)/distributions.c $(SRCDIR)/distributions.h
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/distributions.c

random.o: $(SRCDIR)/random.c $(SRCDIR)/random.h
	$(CC) $(FLAGS) -fpic -c $(SRCDIR)/random.c
