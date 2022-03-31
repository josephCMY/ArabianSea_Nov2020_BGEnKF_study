












module da_par_util1

   use da_control, only : rootproc, ierr, comm, root




   
   
   
   
   

   implicit none

   include 'mpif.h'
   integer, parameter :: true_mpi_real    = mpi_real8
   integer, parameter :: true_mpi_complex = mpi_double_complex

   contains

subroutine da_proc_sum_int (value)

   !---------------------------------------------------------------------------
   !  Purpose: Do MPI sum operation across processors to get the global sum of
   !           an integer value. The sum is returned only on the root processor,
   !           i.e., processor 0. (In this way, we do not have to do all-to-all 
   !           communication, unlike wrf_dm_sum_int, which does)
   !
   ! The routine generates a MPI barrier
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(inout) :: value     ! Value on processor.

   integer :: sum                   ! Sum across processors.

   ! Don't trace, as called within trace routines
   ! if (trace_use_frequent) call da_trace_entry("da_proc_sum_int")

   sum=0
   call mpi_reduce(value, sum, 1, mpi_integer, mpi_sum, root, &
      comm, ierr)

   if (rootproc) value = sum

   ! if (trace_use_frequent) call da_trace_exit("da_proc_sum_int")

end subroutine da_proc_sum_int


subroutine da_proc_sum_ints (values)

   !---------------------------------------------------------------------------
   !  Purpose: Do MPI sum operation across processors to get the global sum of
   !           an integer array. The sum is returned only on the root processor,
   !           i.e., processor 0. (In this way, we do not have to do all-to-all 
   !           communication, unlike wrf_dm_sum_ints, which does)
   !
   ! The routine generates a MPI barrier
   !---------------------------------------------------------------------------

   implicit none

   integer, intent(inout) :: values(:) ! Values

   integer, allocatable :: sums(:) ! Sum across processors.


   ! Don't trace, as called within trace routines
   ! if (trace_use_frequent) call da_trace_entry("da_proc_sum_ints")

   allocate (sums(size(values)))
   sums(:)=0
   call mpi_reduce(values(:), sums(:), size(values), mpi_integer, mpi_sum, &
      root, comm, ierr)

   if (rootproc) values(:) = sums(:)
   deallocate(sums)

   ! if (trace_use_frequent) call da_trace_exit("da_proc_sum_ints")

end subroutine da_proc_sum_ints


subroutine da_proc_sum_real (value)

   !---------------------------------------------------------------------------
   ! Purpose: Do MPI reduction operation across processors to sum a real 
   ! vector.  
   ! On return, each of the N components of the vector "value" represents
   ! the sum of parts from each processor. The result is stored only on 
   ! the root processor, i.e., processor 0. (In this way, we do not have 
   ! to do all-to-all communication, unlike wrf_dm_sum_real, which does)
   !
   ! The routine generates a MPI barrier
   !---------------------------------------------------------------------------

   implicit none

   real, intent(inout) :: value(:) ! Array to be summed componentwise.

   real              :: apsum(size(value))             ! Sum across processors.

   ! Don't trace, as called within trace routines
   !if (trace_use_frequent) call da_trace_entry("da_proc_sum_real")

   apsum(:)=0
   call mpi_reduce(value, apsum, SIZE(value), true_mpi_real, mpi_sum, root, &
      comm, ierr)

   if (rootproc) value = apsum

   !if (trace_use_frequent) call da_trace_exit("da_proc_sum_real")

end subroutine da_proc_sum_real



end module da_par_util1
