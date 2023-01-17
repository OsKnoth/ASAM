PROGRAM TestReshape
  IMPLICIT NONE
  INTEGER ,PARAMETER :: nx = 4
  INTEGER ,PARAMETER :: ny = 3
  INTEGER ,PARAMETER :: nz = 2
  REAL(8) :: xx(nx,ny,nz)
  REAL(8) :: x(nz*ny*nx)

  INTEGER :: ix,iy,iz
  INTEGER :: ind

  ind=1
  DO ix=1,nx
    DO iy=1,ny
      DO iz=1,nz
        xx(ix,iy,iz)=ind
        ind=ind+1
      END DO  
    END DO  
  END DO  
  x=PACK(RESHAPE(xx,(/nz,ny,nx/),ORDER=(/3,2,1/)),MASK=.TRUE.)
  DO ix=1,SIZE(x)
    WRITE(*,*) ix,x(ix)
  END DO  

ENDPROGRAM TestReshape
