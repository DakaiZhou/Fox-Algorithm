# Fox-Algorithm
Fox‘s algorithm is a parallel matrix multiplication function, which distributes the matrix using a checkerboard scheme.

Most parallel matrix multiplication functions use a checkerboard distribution of the matrices. This means that the processes are viewed as a grid, and, rather than assigning entire rows or entire columns to each process, we assign small sub-matrices. For example, if we have four processes, we might assign the element of a 4x4 matrix as shown below, checkerboard mapping of a 4x4 matrix to four processes.

<img src="https://github.com/DakaiZhou/Fox-Algorithm/blob/master/checkerboard.png" height="150" width="400">

Fox‘s algorithm is a one that distributes the matrix using a checkerboard scheme like the above. Here is a example:

To simplify the disscussion, we assume have 3x3(nxn) matrix and p = 9 processes(or cpus). aij , bij and cij are assigned to process(i, j), or  let's say number i*n+j.

<img src="https://github.com/DakaiZhou/Fox-Algorithm/blob/master/fox.png" height="200" width="500">

• stage 0 on process (i, j) : cij = aii*bij

– broadcast aii across the i-th row of processes

– do the local multiplication with bij

– ”shift” the elements of B up one process row, the elements in the top row are shifted to bottom row (circular shift)

• stage 1 on process (i, j) : cij = cij + ai,i+1*bi+1,j in the last row: i + 1 -> (i + 1)mod n

– broadcast ai,i+1 across the i-th row of processes

– do the local multiplication, after the circular shift in stage 0, now the local element of B is bi+1,j

– circular shift of the elements of B up one process row

• ...

• stage k on process (i, j) : cij = cij + ai,k*bk,j with k = (i + k) mod n

• ...

• after stage k = n  1, full multiplication on process (i, j):

=> cij = aii*bij + ai,i+1*bi+1,j + ... + ai,n-1*bn-1,j + ai0*b0j + ... + ai,i-1*bi-1,j

Reference: http://csis.uni-wuppertal.de/courses/scripts/lab2_chapters/PP4.pdf

Test data sets:
