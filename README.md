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
