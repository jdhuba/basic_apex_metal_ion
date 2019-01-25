!-----------------------------------------------------------------------
      subroutine shutdown(msg)
!
! An fatal error has occurred -- shut down the model, including MPI.
!
      implicit none
#ifdef MPI
#include <mpif.h>
#endif
!
! Args:
      character(len=*) :: msg
!
! Local:
      integer :: ier
      character(len=80) :: errorcode
!
      write(6,"(/,28('>'),' MODEL SHUTDOWN ',28('<'))")
      write(6,"('Shutdown: stop message: ',a)") trim(msg)
#ifdef MPI
      write(6,"('Shutdown calling mpi_abort..')")
      write(6,"(72('>'))")
      call mpi_abort(MPI_COMM_WORLD,errorcode,ier) 
#endif
      stop
      end subroutine shutdown
!-----------------------------------------------------------------------
