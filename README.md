# Fox-Algorithm
Fox‘s algorithm is a parallel matrix multiplication function, which distributes the matrix using a checkerboard scheme.

Most parallel matrix multiplication functions use a checkerboard distribution of the matrices. This means that the processes are viewed as a grid, and, rather than assigning entire rows or entire columns to each process, we assign small sub-matrices. For example, if we have four processes, we might assign the element of a 4x4 matrix as shown below, checkerboard mapping of a 4x4 matrix to four processes.

Fox‘s algorithm is a one that distributes the matrix using a checkerboard scheme like the above. Here is a example:

To simplify the disscussion, we assume have 3x3(nxn) matrix and p = 9 processes(or cpus). aij , bij and cij are assigned to process(i, j), or  let's say number i*n+j.
