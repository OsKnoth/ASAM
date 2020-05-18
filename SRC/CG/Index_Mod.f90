MODULE Index_Mod

  IMPLICIT NONE
  INTERFACE IndexC
    MODULE PROCEDURE IndexC3,IndexC4
  END INTERFACE


CONTAINS

FUNCTION IndexC3(ix,iy,iz,nx,ny,nz)

  INTEGER :: IndexC3
  INTEGER :: ix,iy,iz,nx,ny,nz

! SELECT CASE(GridType)
  IndexC3=ix+((iy-1)+(iz-1)*ny)*nx 
! IndexC3=iz+((iy-1)+(ix-1)*ny)*nz 
! IndexC3=iz+((ix-1)+(iy-1)*nx)*nz 
! IndexC3=iy+((iz-1)+(ix-1)*nz)*ny 

END FUNCTION IndexC3
FUNCTION IndexC4(ip,ix,iy,iz,np,nx,ny,nz)

  INTEGER :: IndexC4
  INTEGER :: ip,ix,iy,iz,np,nx,ny,nz

  IndexC4=ip+((ix-1)+((iy-1)+(iz-1)*ny)*nx)*np
! IndexC4=iz+((iy-1)+(ix-1)*ny)*nz 
! IndexC4=ip+((iz-1)+((ix-1)+(iy-1)*nx)*nz)*np 
! IndexC4=iy+((iz-1)+(ix-1)*nz)*ny 

END FUNCTION IndexC4
END MODULE Index_Mod
