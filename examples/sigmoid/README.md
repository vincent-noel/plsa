##A sigmoid function

This is a model with one variable, and four parameters :

	f(x) = k.x^n/(x^n + theta^n) + basal

This examples fits this model to some data.


##Compiling
	make all

##Running the serial version
	./run-serial

##Running the parallel version
	mpirun -np 2 ./run-parallel

##License
	Copyright (C) 2016 Vincent Noel (vincent.noel@butantan.gov.br)

	plsa is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	plsa is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with plsa. If not, see <http://www.gnu.org/licenses/>.
