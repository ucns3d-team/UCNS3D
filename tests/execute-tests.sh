#!/bin/bash

chmod +X ucns3d_p
mpi_tasks=2
echo "ucns3d_p will run on $mpi_tasks MPI tasks."

mpirun -np "$mpi_tasks" ./ucns3d_p