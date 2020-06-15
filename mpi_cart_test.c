/*
 * Example program to show the use of the MPI functions
 * to create and use a cartesian grid of processes.
 *
 * Use:
 *    * if the chosen number of processes is 6, then it will
 *      demonstrate the properties of a 2-dimensional (2x3) grid
 *
 *    * if the chosen number of processes is 24, then it will
 *      demonstrate the properties of a 3-dimensional (2x3x4) grid
 *
 *    * else: it prints an error message
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"





/* FUNCTION TO SETUP THE GRID */
void test_2D_grid()
{
  int old_rank;
  int dimensions[2] = {2,3};    /* assuming p = 6  */
  int wrap_around[2];
  int coordinates[2];
  int free_coords[2];

  int nrows, ncols;

  MPI_Comm cart_comm, row_comm, col_comm;
  int my_cart_rank;

  /* set up global grid information */
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

  /* circular shift in second dimension, also in first just because */
  wrap_around[0] = 1;
  wrap_around[1] = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &cart_comm);
  MPI_Comm_rank(cart_comm, &my_cart_rank);

  /* get process coordinates in grid communicator */
  MPI_Cart_coords(cart_comm, my_cart_rank, 2, coordinates);


  /* set up row communicator */
  free_coords[0] = 0;
  free_coords[1] = 1;
  MPI_Cart_sub(cart_comm, free_coords, &row_comm);


  /* set up column communicator */
  free_coords[0] = 1;
  free_coords[1] = 0;
  MPI_Cart_sub(cart_comm, free_coords, &col_comm);

  MPI_Comm_size(row_comm, &ncols);
  MPI_Comm_size(col_comm, &nrows);


  if( old_rank == 0 )
  {
    printf("\n2-dimensional cartesian grid, with dimensions[2] = {2,3}\n");
    printf("nr of processes in a row: %d\n", ncols);
    printf("nr of processes in a column: %d\n\n\n(see source code for further details)\n", nrows);
  }

  /* print grid info */
  //printf("old rank = %2d,\tCart. rank = %2d,\tcoords = (%2d, %2d)\n", old_rank, my_cart_rank, coordinates[0], coordinates[1]);
}



/* FUNCTION TO SETUP THE GRID */
void test_3D_grid()
{
  int old_rank;
  int  dimensions[3] = {2,3,4};   /* assuming p = 24  */
  int wrap_around[3];
  int coordinates[3];
  int free_coords[3];

  int x_size, y_size, z_size, xy_size, xz_size, yz_size;

  MPI_Comm cart_comm, x_comm, y_comm, z_comm, xy_comm, xz_comm, yz_comm;
  int my_cart_rank;

  /* set up global grid information */
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

  /* circular shift in second dimension, also in first just because */
  wrap_around[0] = 1;
  wrap_around[1] = 1;
  wrap_around[2] = 1;

  MPI_Cart_create(MPI_COMM_WORLD, 3, dimensions, wrap_around, 1, &cart_comm);
  MPI_Comm_rank(cart_comm, &my_cart_rank);

  /* get process coordinates in grid communicator */
  MPI_Cart_coords(cart_comm, my_cart_rank, 3, coordinates);


  /* set up communicator at fixed X coordinate, i.e. on the YZ-plane */
  free_coords[0] = 0;
  free_coords[1] = 1;
  free_coords[2] = 1;
  MPI_Cart_sub(cart_comm, free_coords, &yz_comm);

  /* set up communicator at fixed Y coordinate, i.e. on the XZ-plane */
  free_coords[0] = 1;
  free_coords[1] = 0;
  free_coords[2] = 1;
  MPI_Cart_sub(cart_comm, free_coords, &xz_comm);

  /* set up communicator at fixed Z coordinate, i.e. on the XY-plane */
  free_coords[0] = 1;
  free_coords[1] = 1;
  free_coords[2] = 0;
  MPI_Cart_sub(cart_comm, free_coords, &xy_comm);



  /* set up communicator over the X coordinate, i.e. keeping fixed Y and Z coords */
  free_coords[0] = 1;
  free_coords[1] = 0;
  free_coords[2] = 0;
  MPI_Cart_sub(cart_comm, free_coords, &x_comm);

  /* set up communicator over the Y coordinate, i.e. keeping fixed X and Z coords */
  free_coords[0] = 0;
  free_coords[1] = 1;
  free_coords[2] = 0;
  MPI_Cart_sub(cart_comm, free_coords, &y_comm);

  /* set up communicator over the Z coordinate, i.e. keeping fixed X and Y coords */
  free_coords[0] = 0;
  free_coords[1] = 0;
  free_coords[2] = 1;
  MPI_Cart_sub(cart_comm, free_coords, &z_comm);


  /* get sizes of all communicators to print them */
  MPI_Comm_size(x_comm, &x_size);
  MPI_Comm_size(y_comm, &y_size);
  MPI_Comm_size(z_comm, &z_size);
  MPI_Comm_size(xy_comm, &xy_size);
  MPI_Comm_size(xz_comm, &xz_size);
  MPI_Comm_size(yz_comm, &yz_size);


  if( old_rank == 0 )
  {
    printf("\n3-dimensional cartesian grid, with dimensions[3] = {2,3,4}\n");
    printf("nr of processes in x_comm: %d\n",      x_size);
    printf("nr of processes in y_comm: %d\n",      y_size);
    printf("nr of processes in z_comm: %d\n",      z_size);
    printf("nr of processes in xy_comm: %d\n",     xy_size);
    printf("nr of processes in xz_comm: %d\n",     xz_size);
    printf("nr of processes in yz_comm: %d\n\n\n(see source code for further details)\n", yz_size);
  }

  /* print grid info */
  //printf("old rank = %2d,\tCart. rank = %2d,\tcoords = (%2d, %2d, %2d)\n", old_rank, my_cart_rank, coordinates[0], coordinates[1], coordinates[2]);
}









 /****************   MAIN   *****************/
int main(int argc, char* argv[])
{
  int p, r;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &r);

  if(p == 6) { test_2D_grid(); }
  else if(p == 24) { test_3D_grid(); }
  else { if(r == 0) {printf("this program only works with 6 or 24 processes!!!\n"); } }

  MPI_Finalize();

  return 0;
}
