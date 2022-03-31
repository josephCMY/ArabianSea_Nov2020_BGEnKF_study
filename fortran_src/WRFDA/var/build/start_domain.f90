















SUBROUTINE start_domain ( grid , allowed_to_read )

   USE module_domain
   USE module_configure

   IMPLICIT NONE

   
   TYPE (domain)          :: grid
   LOGICAL, INTENT(IN)    :: allowed_to_read
   
   INTEGER :: idum1, idum2

   INTERFACE
   END INTERFACE

   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )



END SUBROUTINE start_domain

