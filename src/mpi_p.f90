MODULE MPIINFO
USE MPI
IMPLICIT NONE

 !INCLUDE "mpif.h"
!$omp declare target(N)
INTEGER:: N !THE NUMBER OF RANK THAT I HAVE(FOR EACH PROCESSOR!
INTEGER:: ISIZE !THE TOTAL NUMBER OF RANKS(SIZE OF)
! INTEGER:: ICOMMUNICATOR !THE COMMUNICATOR OF COMM_wORLD
INTEGER::IERROR,provided
INTEGER::STATUS(MPI_STATUS_SIZE)
 END MODULE
