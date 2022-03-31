












MODULE da_wavelet



 USE da_control,    ONLY: ime,ims,ite,its,jme,jms,jts,kme,kms,trace_use_dull
 USE module_domain_type, ONLY: domain, vp_type
 IMPLICIT NONE
 CHARACTER            :: namw		
 INTEGER              :: lf		
 INTEGER              :: nb		
 INTEGER, ALLOCATABLE :: nij(:,:,:)	
 REAL, ALLOCATABLE    :: ws(:)		
 REAL, ALLOCATABLE    :: wsd(:,:)	
CONTAINS
 SUBROUTINE da_transform_through_wavelet(grid,mz,wsd,v,vv)
!
! Purpose: Horizontal control-variable transform v -> vv through idwtai_w.c.
! Author: Aime' Fournier
!
! vv(:ni*nj,m) = MATMUL(uh(:,:,m),v(:nq,m)), m the vertical mode,
! uh(ij,q,m) = w(q,ij)*wsd(q,m)/p_x(ij) the m'th horiz. CVT matrix,
! w(q,:ni*nj) the q'th of nq horizontal 2D wavelets,
! wsd(q,m) the standard deviation of the q'th wavelet coefficient,
! p_x(ij)**2 the horizontal area in ij. 
!
! Method:
!  1) Diagonal multiply by wsd(:nq,m);
!  2) MatMul by Transpose(w(:nq,:ni*nj)) using idwtai_w();
!  3) Diagonal divide by p_x(:ni*nj), cf. gen_be_stage4_regional.
!
 IMPLICIT NONE
 TYPE(domain),  INTENT( IN)::grid
 REAL,          INTENT(OUT)::vv(ims:ime,jms:jme,kms:kme)	! Grid point/EOF equivalent.
 INTEGER,       INTENT( IN)::mz					! Vertical truncation.
 REAL,          INTENT( IN)::wsd(nij(0,1,2),nij(0,0,2),mz)	! Wavelet standard deviation:
 REAL,          INTENT( IN)::v(nij(0,0,2)*nij(0,1,2)*mz)	! Field to be transformed.
 INTEGER                   ::dj(0:1),dk(0:1),dv,i,j,m,n
 LOGICAL, SAVE             ::call1=.TRUE.
 REAL                      ::u(nij(0,1,2),nij(0,0,2))		! Since INTENT(IN)::v in caller.

    PRINT'(a,": must compile after setenv WAVELET 1")',"da_transform_through_wavelet.inc"
    CALL wrf_abort

ENDSUBROUTINE da_transform_through_wavelet
SUBROUTINE da_transform_through_wavelet_adj(grid,mz,wsd,v,vv)
!
! Purpose: Horizontal control-var. transform-adjoint vv -> v through dwtai_w.c.
! Author: Aime' Fournier
!
! v(:nq,m) = MATMUL(TRANSPOSE(uh(:,:,m)),vv(:ni*nj,m)), m the vert. mode,
! uh(ij,q,m) = w(q,ij)*wsd(q,m)/p_x(ij) the m'th horiz. CVT matrix,
! w(q,:ni*nj) the q'th of nq horizontal 2D wavelets,
! wsd(q,m) the standard deviation of the q'th wavelet coefficient,
! p_x(ij)**2 the horizontal area in ij. 
!
! Method:
!  3) Diagonal divide by p_x(:ni*nj), cf. gen_be_stage4_regional.
!  2) MatMul by w(:nq,:ni*nj) using dwtai_w();
!  1) Diagonal multiply by wsd(:nq,m);
!
 IMPLICIT NONE
 TYPE(domain), INTENT( IN)::grid
 REAL,         INTENT( IN)::vv(ims:ime,jms:jme,kms:kme)		! Grid point/EOF equivalent.
 INTEGER,      INTENT( IN)::mz					! Vertical truncation.
 REAL,         INTENT( IN)::wsd(nij(0,1,2),nij(0,0,2),mz)	! Wavelet standard deviation.
 REAL,         INTENT(OUT)::v(nij(0,0,2)*nij(0,1,2)*mz)		! Field to be transformed.
 INTEGER                  ::dj(0:1),dk(0:1),dv,i,j,m,n
 LOGICAL, SAVE            ::call1=.TRUE.

    PRINT'(a,": must compile after setenv WAVELET 1")',"da_transform_through_wavelet_adj.inc"
    CALL wrf_abort

ENDSUBROUTINE da_transform_through_wavelet_adj
END MODULE da_wavelet
