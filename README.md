# plsa [![Build Status](https://travis-ci.org/vincent-noel/plsa.svg?branch=master)](https://travis-ci.org/vincent-noel/plsa)
C library for non-linear optimization, using a parallelized Lam-Delosme simulated annealing.



## Dependencies
You will need MPI libraries to execute C code in parallel

	libopenmpi-dev openmpi-bin



## Installation
To compile and install it in your /usr/lib and /usr/include

	sudo make install



## Examples
The library is provided with two simples examples:

	make examples


- A funnel : the score function is the euclidian distance to a number.

		./run-funnel-serial
		mpirun -np 2 ./run-funnel-parallel



- A sigmoid function : the score function is the distance to data.

		./run-sigmoid-serial
		mpirun -np 2 ./run-sigmoid-parallel



## References
This library is based on the following publication :

	King-Wai Chu, Yuefan Deng, and John Reinitz. "Parallel simulated annealing by mixing of states.",
	Journal of Computational Physics 148.2 (1999): 646-662.

Thanks to John Reinitz for sending me the original code.
Originally written by Jimmy Lam, Dan Greening, John Reinitz, King-Wai Chu, Johannes Jaeger.



## License
	Copyright (C) 2016 Vincent Noel

	plsa is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
	License as published by the Free Software Foundation, either version 3 of the License, or (at your option)
	any later version.

	plsa is distributed in the hope that it will be useful,	but WITHOUT ANY WARRANTY; without even the implied
	warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
	more details.

	You should have received a copy of the GNU General Public License along with plsa.
	If not, see <http://www.gnu.org/licenses/>.

## Financial support

	This library was developed within the CeTICS project, at the Butantan Institute.

<p align="center">
	<a href="http://cetics.butantan.gov.br"><img src="docs/logos/cetics.png" align="middle" hspace="50"></a>
	<a href="http://www.butantan.gov.br"><img src="docs/logos/butantan.png" width="300" align="middle" hspace="50"></a>
</p>

	The work was supported by grants #13/07467-1, and #13/24212-7
	of the São Paulo Research Foundation (FAPESP)


<p align="center">
	<a href="http://www.fapesp.br"><img src="docs/logos/FAPESP.jpg" width="300" align="middle"></a>
</p>
