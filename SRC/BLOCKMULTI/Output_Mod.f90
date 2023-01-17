MODULE Output_Mod
USE netcdf         
USE Kind_Mod
USE DomDecomp_Mod, ONLY: read_blocks_from_netcdf, SubDomain, block
USE List_Mod, ONLY: RArray

INTEGER :: nc_err

CONTAINS

SUBROUTINE output_serial(subdom, outfields, outtypes, outnames, outpath)

  IMPLICIT NONE

  TYPE(SubDomain), INTENT(IN) :: subdom
  TYPE(RArray), INTENT(in) :: outfields(:)
  INTEGER, INTENT(in) :: outtypes(:)
  CHARACTER(LEN=*), INTENT(in) :: outnames(:)
  CHARACTER(LEN=*), INTENT(IN) :: outpath

  INTEGER :: iblock, jblock, nblocks, nfields, ifield, icomp_fld, icomp, icell, iface, iface_tmp
  INTEGER :: ncid, blockdimid, blockdimid2, blockvarid, blockvarid2
  INTEGER, DIMENSION(COUNT(subdom%blockiscomp)) :: xdimid, ydimid, zdimid, &
                                                   x2dimid, y2dimid, z2dimid, &
                                                    xvarid, yvarid, zvarid, &
                                                    x2varid, y2varid, z2varid, blockids
  INTEGER :: fldvarid(SUM(outtypes, 1), COUNT(subdom%blockiscomp))
  INTEGER, DIMENSION(3) :: dimids
  INTEGER :: ii_st, i, j, k
  CHARACTER(LEN=8) :: numb2str

  !REAL, ALLOCATABLE :: fld_sub(:,:,:)

  INTEGER :: nx, ny, nz

  nfields = SIZE(outfields, 1)

  nblocks = COUNT(subdom%blockiscomp)

  WRITE(*,*) "Write file ", outpath

  nc_err = nf90_create(trim(outpath), nf90_clobber, ncid)
  nc_err = nf90_def_dim(ncid, 'nblocks', 1, blockdimid)
  nc_err = nf90_def_dim(ncid, 'block', COUNT(subdom%blockiscomp), blockdimid2)
  nc_err = nf90_def_var(ncid, 'nblocks', NF90_INT, blockdimid, blockvarid)
  nc_err = nf90_def_var(ncid, 'blockids', NF90_INT, blockdimid2, blockvarid2)

  jblock = 1
  DO iblock = 1, subdom%nblocks
    IF (subdom%blockiscomp(iblock) .EQV. .TRUE.) THEN
      blockids(jblock) = subdom%blockids_compghst(iblock)
      WRITE(numb2str, '(I8)') subdom%blockids_compghst(iblock)
      nx = subdom%blocks(iblock)%fld_shape(3)
      ny = subdom%blocks(iblock)%fld_shape(2)
      nz = subdom%blocks(iblock)%fld_shape(1)
      nc_err = nf90_def_dim(ncid, 'xb'//ADJUSTL(numb2str), nx, xdimid(jblock))
      nc_err = nf90_def_dim(ncid, 'yb'//ADJUSTL(numb2str), ny, ydimid(jblock))
      nc_err = nf90_def_dim(ncid, 'zb'//ADJUSTL(numb2str), nz, zdimid(jblock))
      nc_err = nf90_def_dim(ncid, 'x2b'//ADJUSTL(numb2str), nx + 1, x2dimid(jblock))
      nc_err = nf90_def_dim(ncid, 'y2b'//ADJUSTL(numb2str), ny + 1, y2dimid(jblock))
      nc_err = nf90_def_dim(ncid, 'z2b'//ADJUSTL(numb2str), nz + 1, z2dimid(jblock))

      nc_err = nf90_def_var(ncid, 'xb'//ADJUSTL(numb2str), NF90_REAL, xdimid(jblock), xvarid(jblock))
      nc_err = nf90_def_var(ncid, 'yb'//ADJUSTL(numb2str), NF90_REAL, ydimid(jblock), yvarid(jblock))
      nc_err = nf90_def_var(ncid, 'zb'//ADJUSTL(numb2str), NF90_REAL, zdimid(jblock), zvarid(jblock))
      nc_err = nf90_def_var(ncid, 'x2b'//ADJUSTL(numb2str), NF90_REAL, x2dimid(jblock), x2varid(jblock))
      nc_err = nf90_def_var(ncid, 'y2b'//ADJUSTL(numb2str), NF90_REAL, y2dimid(jblock), y2varid(jblock))
      nc_err = nf90_def_var(ncid, 'z2b'//ADJUSTL(numb2str), NF90_REAL, z2dimid(jblock), z2varid(jblock))

      dimids = (/xdimid(jblock), ydimid(jblock), zdimid(jblock)/)

      icomp = 1
      DO ifield = 1, nfields
        IF (outtypes(ifield) .EQ. 1) THEN
          nc_err = nf90_def_var(ncid, trim(outnames(ifield))//'_b'//ADJUSTL(numb2str), NF90_REAL, dimids, fldvarid(icomp, jblock))  
          icomp = icomp + 1
        ELSE
          dimids = (/x2dimid(jblock), ydimid(jblock), zdimid(jblock)/)
          nc_err = nf90_def_var(ncid, trim(outnames(ifield))//'_u_b'//ADJUSTL(numb2str), &
                                NF90_REAL, dimids, fldvarid(icomp, jblock))
          icomp = icomp + 1
          dimids = (/xdimid(jblock), y2dimid(jblock), zdimid(jblock)/)
          nc_err = nf90_def_var(ncid, trim(outnames(ifield))//'_v_b'//ADJUSTL(numb2str), &
                                NF90_REAL, dimids, fldvarid(icomp, jblock))
          icomp = icomp + 1
          dimids = (/xdimid(jblock), ydimid(jblock), z2dimid(jblock)/)
          nc_err = nf90_def_var(ncid, trim(outnames(ifield))//'_w_b'//ADJUSTL(numb2str), &
                                NF90_REAL, dimids, fldvarid(icomp, jblock))
          icomp = icomp + 1
        END IF
      END DO
    
      jblock = jblock + 1
    END IF
  END DO

  nc_err = nf90_enddef(ncid)
  nc_err = nf90_put_var(ncid, blockvarid, nblocks)
  nc_err = nf90_put_var(ncid, blockvarid2, blockids)

  jblock = 1
  icell = 1
  iface = 1
  DO iblock = 1, subdom%nblocks
    IF (subdom%blockiscomp(iblock) .EQV. .TRUE.) THEN
      nx = subdom%blocks(iblock)%fld_shape(3)
      ny = subdom%blocks(iblock)%fld_shape(2)
      nz = subdom%blocks(iblock)%fld_shape(1)
      nc_err = nf90_put_var(ncid, xvarid(jblock), subdom%blocks(iblock)%x)
      nc_err = nf90_put_var(ncid, yvarid(jblock), subdom%blocks(iblock)%y)
      nc_err = nf90_put_var(ncid, zvarid(jblock), subdom%blocks(iblock)%z)
      nc_err = nf90_put_var(ncid, x2varid(jblock), subdom%blocks(iblock)%x2)
      nc_err = nf90_put_var(ncid, y2varid(jblock), subdom%blocks(iblock)%y2)
      nc_err = nf90_put_var(ncid, z2varid(jblock), subdom%blocks(iblock)%z2)

      icomp = 1
      DO ifield = 1, nfields
        IF (outtypes(ifield) .EQ. 1) THEN

          nc_err = nf90_put_var(ncid, fldvarid(icomp, jblock), RESHAPE(outfields(ifield)%data(icell:icell + &
                                subdom%blocks(iblock)%ncells - 1), (/nx, ny, nz/), order=(/3, 2, 1/)))
          icomp = icomp + 1
        ELSE
          iface_tmp = iface

          nc_err = nf90_put_var(ncid, fldvarid(icomp, jblock), &
                                RESHAPE(outfields(ifield)%data(iface_tmp:iface_tmp + (nx + 1) * ny * nz - 1), &
                                (/nx + 1, ny, nz/), order=(/3, 2, 1/)))

          icomp = icomp + 1
          iface_tmp = iface_tmp + (nx + 1) * ny * nz

          nc_err = nf90_put_var(ncid, fldvarid(icomp, jblock), &
                                RESHAPE(outfields(ifield)%data(iface_tmp:iface_tmp + nx * (ny + 1) * nz - 1), &
                                (/nx, ny + 1, nz/), order=(/3, 2, 1/)))

          icomp = icomp + 1
          iface_tmp = iface_tmp + (ny + 1) * nx * nz

          nc_err = nf90_put_var(ncid, fldvarid(icomp, jblock), &
                                RESHAPE(outfields(ifield)%data(iface_tmp:iface_tmp + nx * ny * (nz + 1) - 1), &
                                (/nx, ny, nz + 1/), order=(/3, 2, 1/)))
          icomp = icomp + 1
        END IF
      END DO
      icell = icell + subdom%blocks(iblock)%ncells
      iface = iface + subdom%blocks(iblock)%nfaces
      jblock = jblock + 1
    END IF
  END DO

  nc_err = nf90_close(ncid)

END SUBROUTINE output_serial


END MODULE Output_Mod



