# plsa
C library for non-linear optimization, using a parallelized Lam-Delosme simulated annealing.



##Dependencies
You will need MPI libraries to execute C code in parallel

	libopenmpi-dev openmpi-bin



##Installation
To compile the shared library:

	make all

To compile and install it in your /usr/lib and /usr/include

	sudo make install



##References
This library is based on the following publication :

	King-Wai Chu, Yuefan Deng, and John Reinitz. "Parallel simulated annealing by mixing of states.",
	Journal of Computational Physics 148.2 (1999): 646-662.

Thanks to John Reinitz for sending me the original code.



##License
	Copyright (C) 2016 Vincent Noel
	
	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.
	
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.
	
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
