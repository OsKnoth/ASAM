MODULE Gradient_block_interpolation_Mod

USE Kind_Mod
USE Block_Mod, ONLY: Block


IMPLICIT NONE

CONTAINS

SUBROUTINE boundint_weights_yz(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                               x, y, z, y_bnd, z_bnd, y2_bnd, z2_bnd, volseff_bnd, arseffx_bnd, &
                               j, k, j_bnds, k_bnds, i, point_dest)

  IMPLICIT NONE

  REAL(Realkind), INTENT(inout) :: weights(3), points_volsarea(3, 2)
  INTEGER, INTENT(inout) :: points_loc(3, 4)
  INTEGER, INTENT(inout) :: ipoint
  TYPE(Block), POINTER, INTENT(in) :: blocks(:)
  INTEGER, INTENT(in) :: bndblock_id
  REAL(Realkind), POINTER, INTENT(in) :: x(:), y(:), z(:), y_bnd(:), z_bnd(:), y2_bnd(:), z2_bnd(:)
  REAL(Realkind), POINTER, INTENT(in) :: volseff_bnd(:,:,:), arseffx_bnd(:,:,:)
  INTEGER, INTENT(in) :: j, k, i
  INTEGER, INTENT(in) :: j_bnds(2), k_bnds(2)
  REAL(Realkind), INTENT(inout) :: point_dest(2)
  
  REAL(Realkind) :: points(3, 2)
  REAL(Realkind) :: tmp_arr(16)
  INTEGER :: nz_bnd, ny_bnd, nx_bnd
  INTEGER :: nz_tmp, ny_tmp, nx_tmp
  REAL(Realkind), POINTER :: volseff_tmp(:,:,:), arseffx_tmp(:,:,:), arseffy_tmp(:,:,:), arseffz_tmp(:,:,:)
  REAL(Realkind), POINTER :: z_tmp(:), y_tmp(:), x_tmp(:)
  INTEGER :: bndblock_id2
  INTEGER :: k_bnds2(2), j_bnds2(2), i_bnds2(2)
  INTEGER :: l, kk, jj, ii, k_bnd_2, j_bnd_2, i_bnd_2
  
  REAL(Realkind), PARAMETER :: eps=1e-20
  

  nz_bnd = blocks(bndblock_id)%fld_shape(1)
  ny_bnd = blocks(bndblock_id)%fld_shape(2)
  nx_bnd = blocks(bndblock_id)%fld_shape(3)

  point_dest(1) = y(j)
  point_dest(2) = z(k)


  points_loc(1, 1) = bndblock_id
  points_loc(1, 2) = lubnd(volseff_bnd, 3, i)
  points_loc(1, 3) = j_bnds(1)
  points_loc(1, 4) = k_bnds(1)
  points_volsarea(1, 1) = volseff_bnd(k_bnds(1), j_bnds(1), lubnd(volseff_bnd, 3, i))
  points_volsarea(1, 2) = arseffx_bnd(k_bnds(1), j_bnds(1), lubnd(arseffx_bnd, 3, i))

  points(1, 1) = y_bnd(j_bnds(1))
  points(1, 2) = z_bnd(k_bnds(1))

  ipoint = 1


  IF (ABS(y(j) - y2_bnd(j_bnds(1))) .LT. ABS(y(j) - y2_bnd(j_bnds(1) + 1))) THEN

    IF (j_bnds(1) == 1) THEN
      IF (ANY(blocks(bndblock_id)%cface_s)) THEN
        CALL find_indices_ik(blocks, blocks(bndblock_id)%cface_s, lubnd(volseff_bnd, 3, i), &
                             i_bnds2, k_bnds(1), k_bnds2, bndblock_id, bndblock_id2)

        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(k) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0

        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(lubnd(volseff_bnd, 3, i)) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = y2_bnd(j_bnds(1)) - 0.5 * volseff_tmp(k_bnd_2, ny_tmp, i_bnd_2) / &
                            (arseffy_tmp(k_bnd_2, ny_tmp + 1, i_bnd_2) + eps)
        points(ipoint, 2) = z_tmp(k_bnd_2)
        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, ny_tmp, i_bnd_2)
        IF (i .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffx_tmp(k_bnd_2, ny_tmp, i_bnd_2 + 1)
        ELSE
          points_volsarea(ipoint, 2) = arseffx_tmp(k_bnd_2, ny_tmp, i_bnd_2)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = ny_tmp
        points_loc(ipoint, 4) = k_bnd_2

      END IF

    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = y_bnd(j_bnds(1) - 1)
      points(ipoint, 2) = z_bnd(k_bnds(1))
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1), j_bnds(1) - 1, lubnd(volseff_bnd, 3, i))
      points_volsarea(ipoint, 2) = arseffx_bnd(k_bnds(1), j_bnds(1) - 1, lubnd(arseffx_bnd, 3, i))

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = lubnd(volseff_bnd, 3, i)
      points_loc(ipoint, 3) = j_bnds(1) - 1
      points_loc(ipoint, 4) = k_bnds(1)

    END IF

  ELSE
    IF (j_bnds(1) .EQ. ny_bnd) THEN

      IF (ANY(blocks(bndblock_id)%cface_n)) THEN
        CALL find_indices_ik(blocks, blocks(bndblock_id)%cface_n, lubnd(volseff_bnd, 3, i), &
                             i_bnds2, k_bnds(1), k_bnds2, bndblock_id, bndblock_id2)

        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(k) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(lubnd(volseff_bnd, 3, i)) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1


        ipoint = ipoint + 1
        points(ipoint, 1) = y2_bnd(j_bnds(1) + 1) + 0.5 * volseff_tmp(k_bnd_2, 1, i_bnd_2) / &
                            (arseffy_tmp(k_bnd_2, 1, i_bnd_2) + eps)
        points(ipoint, 2) = z_tmp(k_bnd_2)

        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, 1, i_bnd_2)
        IF (i .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffx_tmp(k_bnd_2, 1, i_bnd_2)
        ELSE
          points_volsarea(ipoint, 2) = arseffx_tmp(k_bnd_2, 1, i_bnd_2 + 1)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = 1
        points_loc(ipoint, 4) = k_bnd_2

      END IF
    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = y_bnd(j_bnds(1) + 1)
      points(ipoint, 2) = z_bnd(k_bnds(1))
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1), j_bnds(1) + 1, lubnd(volseff_bnd, 3, 1))
      points_volsarea(ipoint, 2) = arseffx_bnd(k_bnds(1), j_bnds(1) + 1, lubnd(arseffx_bnd, 3, 1))

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = lubnd(volseff_bnd, 3, i)
      points_loc(ipoint, 3) = j_bnds(1) + 1
      points_loc(ipoint, 4) = k_bnds(1)

    END IF
  END IF

  IF (ABS(z(k) - z2_bnd(k_bnds(1))) .LT. ABS(z(k) - z2_bnd(k_bnds(1) + 1))) THEN
    IF (k_bnds(1) .EQ. 1) THEN
      IF (ANY(blocks(bndblock_id)%cface_b)) THEN
        CALL find_indices_ij(blocks, blocks(bndblock_id)%cface_b, lubnd(volseff_bnd, 3, i), &
                             i_bnds2, j_bnds(1), j_bnds2, bndblock_id, bndblock_id2)     

        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(j) - y_tmp(jj)))
          l = l + 1
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(lubnd(volseff_bnd, 3, i)) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1
 
        ipoint = ipoint + 1                                
        points(ipoint, 1) = y_tmp(j_bnd_2)
        points(ipoint, 2) = z2_bnd(k_bnds(1)) - 0.5 * volseff_tmp(nz_tmp, j_bnd_2, i_bnd_2) / &
                            (arseffz_tmp(nz_tmp + 1, j_bnd_2, i_bnd_2) + eps)
        points_volsarea(ipoint, 1) = volseff_tmp(nz_tmp, j_bnd_2, i_bnd_2)
        IF (i .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffx_tmp(nz_tmp, j_bnd_2, i_bnd_2 + 1)
        ELSE
          points_volsarea(ipoint, 2) = arseffx_tmp(nz_tmp, j_bnd_2, i_bnd_2)
        END IF

               
        points_loc(ipoint, 1) = bndblock_id2 
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = nz_tmp

      END IF

    ELSE

      ipoint = ipoint + 1            
      points(ipoint, 1) = y_bnd(j_bnds(1))
      points(ipoint, 2) = z_bnd(k_bnds(1) - 1)
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1) - 1, j_bnds(1), lubnd(volseff_bnd, 3, 1))
      points_volsarea(ipoint, 2) = arseffx_bnd(k_bnds(1) - 1, j_bnds(1), lubnd(arseffx_bnd, 3, 1))
            
      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = lubnd(volseff_bnd, 3, i)
      points_loc(ipoint, 3) = j_bnds(1)
      points_loc(ipoint, 4) = k_bnds(1) - 1

    END IF
  ELSE   
    IF (k_bnds(1) .EQ. nz_bnd) THEN
      IF (ANY(blocks(bndblock_id)%cface_t)) THEN
        CALL find_indices_ij(blocks, blocks(bndblock_id)%cface_t, lubnd(volseff_bnd, 3, i), &
                             i_bnds2, j_bnds(1), j_bnds2, bndblock_id, bndblock_id2)

        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(j) - y_tmp(jj)))
          l = l + 1
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(lubnd(volseff_bnd, 3, i)) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = y_tmp(j_bnd_2)
        points(ipoint, 2) = z2_bnd(k_bnds(1) + 1) + 0.5 * volseff_tmp(1, j_bnd_2, i_bnd_2) / &
                            (arseffz_tmp(1, j_bnd_2, i_bnd_2) + eps)
        points_volsarea(ipoint, 1) = volseff_tmp(1, j_bnd_2, i_bnd_2)
        IF (i .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffx_tmp(1, j_bnd_2, i_bnd_2 + 1)
        ELSE
          points_volsarea(ipoint, 2) = arseffx_tmp(1, j_bnd_2, i_bnd_2)
        END IF  
              
        points_loc(ipoint, 1) = bndblock_id2 
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = 1

      END IF

    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = y_bnd(j_bnds(1))
      points(ipoint, 2) = z_bnd(k_bnds(1) + 1)
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1) + 1, j_bnds(1), lubnd(volseff_bnd, 3, 1))
      points_volsarea(ipoint, 2) = arseffx_bnd(k_bnds(1) + 1, j_bnds(1), lubnd(arseffx_bnd, 3, 1))
            
      points_loc(ipoint, 1) = bndblock_id 
      points_loc(ipoint, 2) = lubnd(volseff_bnd, 3, i)
      points_loc(ipoint, 3) = j_bnds(1)
      points_loc(ipoint, 4) = k_bnds(1) + 1

    END IF
  END IF
  
  CALL linear_weights(weights, points, ipoint, point_dest)


END SUBROUTINE boundint_weights_yz


SUBROUTINE boundint_weights_xz(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                               x, y, z, x_bnd, z_bnd, x2_bnd, z2_bnd, volseff_bnd, arseffy_bnd, &
                               i, k, i_bnds, k_bnds, j, point_dest)

  IMPLICIT NONE

  REAL(Realkind), INTENT(inout) :: weights(3), points_volsarea(3, 2)
  INTEGER, INTENT(inout) :: points_loc(3, 4)
  INTEGER, INTENT(inout) :: ipoint
  TYPE(Block), POINTER, INTENT(in) :: blocks(:)
  INTEGER, INTENT(in) :: bndblock_id
  REAL(Realkind), POINTER, INTENT(in) :: x(:), y(:), z(:), x_bnd(:), z_bnd(:), x2_bnd(:), z2_bnd(:)
  REAL(Realkind), POINTER, INTENT(in) :: volseff_bnd(:,:,:), arseffy_bnd(:,:,:)
  INTEGER, INTENT(in) :: j, k, i
  INTEGER, INTENT(in) :: i_bnds(2), k_bnds(2)
  REAL(Realkind), INTENT(inout) :: point_dest(2)

  REAL(Realkind) :: points(3, 2)
  REAL(Realkind) :: tmp_arr(16)
  INTEGER :: nz_bnd, ny_bnd, nx_bnd
  INTEGER :: nz_tmp, ny_tmp, nx_tmp
  REAL(Realkind), POINTER :: volseff_tmp(:,:,:), arseffx_tmp(:,:,:), arseffy_tmp(:,:,:), arseffz_tmp(:,:,:)
  REAL(Realkind), POINTER :: z_tmp(:), y_tmp(:), x_tmp(:)
  INTEGER :: bndblock_id2
  INTEGER :: k_bnds2(2), j_bnds2(2), i_bnds2(2)
  INTEGER :: l, kk, jj, ii, k_bnd_2, j_bnd_2, i_bnd_2

  REAL(Realkind), PARAMETER :: eps=1e-20

  nz_bnd = blocks(bndblock_id)%fld_shape(1)
  ny_bnd = blocks(bndblock_id)%fld_shape(2)
  nx_bnd = blocks(bndblock_id)%fld_shape(3)

  point_dest(1) = x(i)
  point_dest(2) = z(k)
  
  points_loc(1, 1) = bndblock_id
  points_loc(1, 2) = i_bnds(1)
  points_loc(1, 3) = lubnd(volseff_bnd, 2, j)
  points_loc(1, 4) = k_bnds(1)
  points_volsarea(1, 1) = volseff_bnd(k_bnds(1), lubnd(volseff_bnd, 2, j), i_bnds(1))
  points_volsarea(1, 2) = arseffy_bnd(k_bnds(1), lubnd(arseffy_bnd, 2, j), i_bnds(1))
  
  points(1, 1) = x_bnd(i_bnds(1))
  points(1, 2) = z_bnd(k_bnds(1))
      
  ipoint = 1
        
  IF (ABS(x(i) - x2_bnd(i_bnds(1))) .LT. ABS(x(i) - x2_bnd(i_bnds(1) + 1))) THEN
      
    IF (i_bnds(1) == 1) THEN
      IF (ANY(blocks(bndblock_id)%cface_w)) THEN 
        CALL find_indices_jk(blocks, blocks(bndblock_id)%cface_w, lubnd(volseff_bnd, 2, j), &
                             j_bnds2, k_bnds(1), k_bnds2, bndblock_id, bndblock_id2)
        
        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)
        
        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(k) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(lubnd(volseff_bnd, 2, j)) - y_tmp(jj)))
          l = l + 1
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = x2_bnd(i_bnds(1)) - 0.5 * volseff_tmp(k_bnd_2, j_bnd_2, nx_tmp) / &
                            (arseffx_tmp(k_bnd_2, j_bnd_2, nx_tmp + 1) + eps)
        points(ipoint, 2) = z_tmp(k_bnd_2)
        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, j_bnd_2, nx_tmp)
        IF (j .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffy_tmp(k_bnd_2, j_bnd_2 + 1, nx_tmp)
        ELSE
          points_volsarea(ipoint, 2) = arseffy_tmp(k_bnd_2, j_bnd_2, nx_tmp)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = nx_tmp
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = k_bnd_2

      END IF

    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = x_bnd(i_bnds(1) - 1)
      points(ipoint, 2) = z_bnd(k_bnds(1))
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1), lubnd(volseff_bnd, 2, j), i_bnds(1) - 1)
      points_volsarea(ipoint, 2) = arseffy_bnd(k_bnds(1), lubnd(arseffy_bnd, 2, j), i_bnds(1) - 1)

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = i_bnds(1) - 1
      points_loc(ipoint, 3) = lubnd(volseff_bnd, 2, j)
      points_loc(ipoint, 4) = k_bnds(1)

    END IF

  ELSE
    IF (i_bnds(1) .EQ. nx_bnd) THEN

      IF (ANY(blocks(bndblock_id)%cface_e)) THEN
        CALL find_indices_jk(blocks, blocks(bndblock_id)%cface_e, lubnd(volseff_bnd, 2, j), &
                             j_bnds2, k_bnds(1), k_bnds2, bndblock_id, bndblock_id2)

        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(k) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(lubnd(volseff_bnd, 2, j)) - y_tmp(jj)))
          l = l + 1         
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = x2_bnd(i_bnds(1) + 1) + 0.5 * volseff_tmp(k_bnd_2, j_bnd_2, 1) / &
                            (arseffx_tmp(k_bnd_2, j_bnd_2, 1) + eps)
        points(ipoint, 2) = z_tmp(k_bnd_2)
        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, j_bnd_2, 1)
        IF (j .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffy_tmp(k_bnd_2, j_bnd_2 + 1, 1)
        ELSE
          points_volsarea(ipoint, 2) = arseffy_tmp(k_bnd_2, j_bnd_2, 1)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = 1
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = k_bnd_2

      END IF
    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = x_bnd(i_bnds(1) + 1)
      points(ipoint, 2) = z_bnd(k_bnds(1))
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1), lubnd(volseff_bnd, 2, 1), i_bnds(1) + 1)
      points_volsarea(ipoint, 2) = arseffy_bnd(k_bnds(1), lubnd(arseffy_bnd, 2, 1), i_bnds(1) + 1)

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = i_bnds(1) + 1
      points_loc(ipoint, 3) = lubnd(volseff_bnd, 2, j)
      points_loc(ipoint, 4) = k_bnds(1)

    END IF
  END IF

  IF (ABS(z(k) - z2_bnd(k_bnds(1))) .LT. ABS(z(k) - z2_bnd(k_bnds(1) + 1))) THEN
    IF (k_bnds(1) .EQ. 1) THEN
      IF (ANY(blocks(bndblock_id)%cface_b)) THEN 
        CALL find_indices_ij(blocks, blocks(bndblock_id)%cface_b, i_bnds(1), i_bnds2, &
                             lubnd(volseff_bnd, 2, j), j_bnds2, bndblock_id, bndblock_id2)
        
        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)
        
        tmp_arr(:) = 0.0
        l = 0
        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(i) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1
        
        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(lubnd(volseff_bnd, 2, j)) - y_tmp(jj)))
          l = l + 1
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        ipoint = ipoint + 1                                
        points(ipoint, 1) = x_tmp(i_bnd_2)
        points(ipoint, 2) = z2_bnd(k_bnds(1)) - 0.5 * volseff_tmp(nz_tmp, j_bnd_2, i_bnd_2) / &
                            (arseffz_tmp(nz_tmp + 1, j_bnd_2, i_bnd_2) + eps)
        points_volsarea(ipoint, 1) = volseff_tmp(nz_tmp, j_bnd_2, i_bnd_2)
        IF (j .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffy_tmp(nz_tmp, j_bnd_2 + 1, i_bnd_2)
        ELSE
          points_volsarea(ipoint, 2) = arseffy_tmp(nz_tmp, j_bnd_2, i_bnd_2)
        END IF             
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = nz_tmp

      END IF

    ELSE
      
      ipoint = ipoint + 1            
      points(ipoint, 1) = x_bnd(i_bnds(1))
      points(ipoint, 2) = z_bnd(k_bnds(1) - 1)
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1) - 1, lubnd(volseff_bnd, 2, 1), i_bnds(1))
      points_volsarea(ipoint, 2) = arseffy_bnd(k_bnds(1) - 1, lubnd(arseffy_bnd, 2, 1), i_bnds(1))
            
      points_loc(ipoint, 1) = bndblock_id 
      points_loc(ipoint, 2) = i_bnds(1)
      points_loc(ipoint, 3) = lubnd(volseff_bnd, 2, j)
      points_loc(ipoint, 4) = k_bnds(1) - 1

    END IF

  ELSE  

    IF (k_bnds(1) .EQ. nz_bnd) THEN
      IF (ANY(blocks(bndblock_id)%cface_t)) THEN
        CALL find_indices_ij(blocks, blocks(bndblock_id)%cface_t, i_bnds(1), i_bnds2, &
                             lubnd(volseff_bnd, 2, j), j_bnds2, bndblock_id, bndblock_id2)

        z_tmp => blocks(bndblock_id2)%z
        y_tmp => blocks(bndblock_id2)%y
        x_tmp => blocks(bndblock_id2)%x
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(i) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(lubnd(volseff_bnd, 2, j)) - y_tmp(jj)))
          l = l + 1
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = x_tmp(i_bnd_2)
        points(ipoint, 2) = z2_bnd(k_bnds(1) + 1) + 0.5 * volseff_tmp(1, j_bnd_2, i_bnd_2) / &
                            (arseffz_tmp(1, j_bnd_2, i_bnd_2) + eps)
        points_volsarea(ipoint, 1) = volseff_tmp(1, j_bnd_2, i_bnd_2)
        IF (j .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffy_tmp(1, j_bnd_2 + 1, i_bnd_2)
        ELSE
          points_volsarea(ipoint, 2) = arseffy_tmp(1, j_bnd_2, i_bnd_2)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = 1

      END IF

    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = x_bnd(i_bnds(1))
      points(ipoint, 2) = z_bnd(k_bnds(1) + 1)
      points_volsarea(ipoint, 1) = volseff_bnd(k_bnds(1) + 1, lubnd(volseff_bnd, 2, 1), i_bnds(1))
      points_volsarea(ipoint, 2) = arseffy_bnd(k_bnds(1) + 1, lubnd(arseffy_bnd, 2, 1), i_bnds(1))

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = i_bnds(1)
      points_loc(ipoint, 3) = lubnd(volseff_bnd, 2, j)
      points_loc(ipoint, 4) = k_bnds(1) + 1

    END IF

  END IF

  CALL linear_weights(weights, points, ipoint, point_dest)

END SUBROUTINE boundint_weights_xz


SUBROUTINE boundint_weights_xy(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                               x, y, z, x_bnd, y_bnd, x2_bnd, y2_bnd, volseff_bnd, arseffz_bnd, &
                               i, j, i_bnds, j_bnds, k, point_dest)

  IMPLICIT NONE

  REAL(Realkind), INTENT(inout) :: weights(3), points_volsarea(3, 2)
  INTEGER, INTENT(inout) :: points_loc(3, 4)
  INTEGER, INTENT(inout) :: ipoint
  TYPE(Block), POINTER, INTENT(in) :: blocks(:)
  INTEGER, INTENT(in) :: bndblock_id
  REAL(Realkind), POINTER, INTENT(in) :: x(:), y(:), z(:), x_bnd(:), y_bnd(:), x2_bnd(:), y2_bnd(:)
  REAL(Realkind), POINTER, INTENT(in) :: volseff_bnd(:,:,:), arseffz_bnd(:,:,:)
  INTEGER, INTENT(in) :: j, k, i
  INTEGER, INTENT(in) :: i_bnds(2), j_bnds(2)
  REAL(Realkind), INTENT(inout) :: point_dest(2)

  REAL(Realkind) :: points(3, 2)
  REAL(Realkind) :: tmp_arr(16)
  INTEGER :: nz_bnd, ny_bnd, nx_bnd
  INTEGER :: nz_tmp, ny_tmp, nx_tmp
  REAL(Realkind), POINTER :: volseff_tmp(:,:,:), arseffx_tmp(:,:,:), arseffy_tmp(:,:,:), arseffz_tmp(:,:,:)
  REAL(Realkind), POINTER :: z_tmp(:), y_tmp(:), x_tmp(:)
  INTEGER :: bndblock_id2
  INTEGER :: k_bnds2(2), j_bnds2(2), i_bnds2(2)
  INTEGER :: l, kk, jj, ii, k_bnd_2, j_bnd_2, i_bnd_2

  REAL(Realkind), PARAMETER :: eps=1e-20

  nz_bnd = blocks(bndblock_id)%fld_shape(1)
  ny_bnd = blocks(bndblock_id)%fld_shape(2)
  nx_bnd = blocks(bndblock_id)%fld_shape(3)

  point_dest(1) = x(i)
  point_dest(2) = y(j)

  points_loc(1, 1) = bndblock_id
  points_loc(1, 2) = i_bnds(1)
  points_loc(1, 3) = j_bnds(1)
  points_loc(1, 4) = lubnd(volseff_bnd, 1, k)
  points_volsarea(1, 1) = volseff_bnd(lubnd(volseff_bnd, 1, k), j_bnds(1), i_bnds(1))
  points_volsarea(1, 2) = arseffz_bnd(lubnd(arseffz_bnd, 1, k), j_bnds(1), i_bnds(1))

  points(1, 1) = x_bnd(i_bnds(1))
  points(1, 2) = y_bnd(j_bnds(1))

  ipoint = 1

  IF (ABS(x(i) - x2_bnd(i_bnds(1))) .LT. ABS(x(i) - x2_bnd(i_bnds(1) + 1))) THEN

    IF (i_bnds(1) == 1) THEN
      IF (ANY(blocks(bndblock_id)%cface_w)) THEN
        CALL find_indices_jk(blocks, blocks(bndblock_id)%cface_w, i_bnds(1), i_bnds2, &
                             lubnd(volseff_bnd, 1, k), k_bnds2, bndblock_id, bndblock_id2)

        x_tmp => blocks(bndblock_id2)%x
        y_tmp => blocks(bndblock_id2)%y
        z_tmp => blocks(bndblock_id2)%z

        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(j) - y_tmp(jj)))
          l = l + 1
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(lubnd(volseff_bnd, 1, k)) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = x2_bnd(i_bnds(1)) - 0.5 * volseff_tmp(k_bnd_2, j_bnd_2, nx_tmp) / &
                            (arseffy_tmp(k_bnd_2, j_bnd_2, nx_tmp + 1) + eps)
        points(ipoint, 2) = y_tmp(j_bnd_2)
        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, j_bnd_2, nx_tmp)
        IF (k .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffx_tmp(k_bnd_2 + 1, j_bnd_2, nx_tmp)
        ELSE
          points_volsarea(ipoint, 2) = arseffx_tmp(k_bnd_2, j_bnd_2, nx_tmp)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = nx_tmp
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = k_bnd_2

      END IF

    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = x_bnd(i_bnds(1) - 1)
      points(ipoint, 2) = y_bnd(j_bnds(1))
      points_volsarea(ipoint, 1) = volseff_bnd(lubnd(volseff_bnd, 1, k), j_bnds(1), i_bnds(1) - 1)
      points_volsarea(ipoint, 2) = arseffz_bnd(lubnd(arseffz_bnd, 1, k), j_bnds(1), i_bnds(1) - 1)

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = i_bnds(1) - 1
      points_loc(ipoint, 3) = j_bnds(1)
      points_loc(ipoint, 4) = lubnd(volseff_bnd, 1, k)

    END IF

  ELSE
    IF (i_bnds(1) .EQ. nx_bnd) THEN

      IF (ANY(blocks(bndblock_id)%cface_e)) THEN
        CALL find_indices_jk(blocks, blocks(bndblock_id)%cface_e, i_bnds(1), i_bnds2, &
                             lubnd(volseff_bnd, 1, k), k_bnds2, bndblock_id, bndblock_id2)

        x_tmp => blocks(bndblock_id2)%x
        y_tmp => blocks(bndblock_id2)%y
        z_tmp => blocks(bndblock_id2)%z
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffx_tmp => blocks(bndblock_id2)%arseffx
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO jj = j_bnds2(1), j_bnds2(2)
          tmp_arr(l + 1) = ABS((y(j) - y_tmp(jj)))
          l = l + 1
        END DO
        j_bnd_2 = MINLOC(tmp_arr(1:l), 1) + j_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(lubnd(volseff_bnd, 1, k)) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = x2_bnd(i_bnds(1) + 1) + 0.5 * volseff_tmp(k_bnd_2, j_bnd_2, 1) / &
                            (arseffx_tmp(k_bnd_2, j_bnd_2, 1) + eps)
        points(ipoint, 2) = y_tmp(j_bnd_2)
        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, j_bnd_2, 1)
        IF (k .EQ. -1) THEN
           points_volsarea(ipoint, 2) = arseffz_tmp(k_bnd_2 + 1, j_bnd_2, 1)
        ELSE
           points_volsarea(ipoint, 2) = arseffz_tmp(k_bnd_2, j_bnd_2, 1)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = 1
        points_loc(ipoint, 3) = j_bnd_2
        points_loc(ipoint, 4) = k_bnd_2

      END IF
    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = x_bnd(i_bnds(1) + 1)
      points(ipoint, 2) = y_bnd(j_bnds(1))
      points_volsarea(ipoint, 1) = volseff_bnd(lubnd(volseff_bnd, 1, 1), j_bnds(1), i_bnds(1) + 1)
      points_volsarea(ipoint, 2) = arseffz_bnd(lubnd(arseffz_bnd, 1, 1), j_bnds(1), i_bnds(1) + 1)

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = i_bnds(1) + 1
      points_loc(ipoint, 3) = j_bnds(1)
      points_loc(ipoint, 4) = lubnd(volseff_bnd, 1, k)

    END IF
  END IF

  IF (ABS(y(j) - y2_bnd(j_bnds(1))) .LT. ABS(y(j) - y2_bnd(j_bnds(1) + 1))) THEN
    IF (j_bnds(1) .EQ. 1) THEN
      IF (ANY(blocks(bndblock_id)%cface_s)) THEN
        CALL find_indices_ik(blocks, blocks(bndblock_id)%cface_s, i_bnds(1), i_bnds2, &
                             lubnd(volseff_bnd, 1, k), k_bnds2, bndblock_id, bndblock_id2)

        x_tmp => blocks(bndblock_id2)%x
        y_tmp => blocks(bndblock_id2)%y
        z_tmp => blocks(bndblock_id2)%z
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(i) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(lubnd(volseff_bnd, 1, k)) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = x_tmp(i_bnd_2)
        points(ipoint, 2) = y2_bnd(j_bnds(1)) - 0.5 * volseff_tmp(k_bnd_2, ny_tmp, i_bnd_2) / &
                            (arseffy_tmp(k_bnd_2, ny_tmp + 1, i_bnd_2) + eps)
        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, ny_tmp, i_bnd_2)
        IF (k .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffz_tmp(k_bnd_2 + 1, ny_tmp, i_bnd_2)
        ELSE
          points_volsarea(ipoint, 2) = arseffz_tmp(k_bnd_2, ny_tmp, i_bnd_2)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = ny_tmp
        points_loc(ipoint, 4) = k_bnd_2

      END IF

    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = x_bnd(i_bnds(1))
      points(ipoint, 2) = y_bnd(j_bnds(1) - 1)
      points_volsarea(ipoint, 1) = volseff_bnd(lubnd(volseff_bnd, 1, 1), j_bnds(1) - 1, i_bnds(1))
      points_volsarea(ipoint, 2) = arseffz_bnd(lubnd(arseffz_bnd, 1, 1), j_bnds(1) - 1, i_bnds(1))

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = i_bnds(1)
      points_loc(ipoint, 3) = j_bnds(1) - 1
      points_loc(ipoint, 4) = lubnd(volseff_bnd, 1, k)

    END IF

  ELSE

    IF (j_bnds(1) .EQ. ny_bnd) THEN
      IF (ANY(blocks(bndblock_id)%cface_n)) THEN
        CALL find_indices_ik(blocks, blocks(bndblock_id)%cface_n, i_bnds(1), i_bnds2, &
                             lubnd(volseff_bnd, 1, k), k_bnds2, bndblock_id, bndblock_id2)

        x_tmp => blocks(bndblock_id2)%x
        y_tmp => blocks(bndblock_id2)%y
        z_tmp => blocks(bndblock_id2)%z
        volseff_tmp => blocks(bndblock_id2)%volseff
        arseffy_tmp => blocks(bndblock_id2)%arseffy
        arseffz_tmp => blocks(bndblock_id2)%arseffz
        nz_tmp = blocks(bndblock_id2)%fld_shape(1)
        ny_tmp = blocks(bndblock_id2)%fld_shape(2)
        nx_tmp = blocks(bndblock_id2)%fld_shape(3)

        tmp_arr(:) = 0.0
        l = 0
        DO ii = i_bnds2(1), i_bnds2(2)
          tmp_arr(l + 1) = ABS((x(i) - x_tmp(ii)))
          l = l + 1
        END DO
        i_bnd_2 = MINLOC(tmp_arr(1:l), 1) + i_bnds2(1) - 1

        tmp_arr(:) = 0.0
        l = 0
        DO kk = k_bnds2(1), k_bnds2(2)
          tmp_arr(l + 1) = ABS((z(lubnd(volseff_bnd, 1, k)) - z_tmp(kk)))
          l = l + 1
        END DO
        k_bnd_2 = MINLOC(tmp_arr(1:l), 1) + k_bnds2(1) - 1

        ipoint = ipoint + 1
        points(ipoint, 1) = x_tmp(i_bnd_2)
        points(ipoint, 2) = y2_bnd(j_bnds(1) + 1) + 0.5 * volseff_tmp(k_bnd_2, 1, i_bnd_2) / &
                            (arseffy_tmp(k_bnd_2, 1, i_bnd_2) + eps)
        points_volsarea(ipoint, 1) = volseff_tmp(k_bnd_2, 1, i_bnd_2)
        IF (k .EQ. -1) THEN
          points_volsarea(ipoint, 2) = arseffz_tmp(k_bnd_2 + 1, 1, i_bnd_2)
        ELSE
          points_volsarea(ipoint, 2) = arseffz_tmp(k_bnd_2, 1, i_bnd_2)
        END IF
        points_loc(ipoint, 1) = bndblock_id2
        points_loc(ipoint, 2) = i_bnd_2
        points_loc(ipoint, 3) = ny_tmp
        points_loc(ipoint, 4) = k_bnd_2

      END IF

    ELSE

      ipoint = ipoint + 1
      points(ipoint, 1) = x_bnd(i_bnds(1))
      points(ipoint, 2) = y_bnd(j_bnds(1) + 1)
      points_volsarea(ipoint, 1) = volseff_bnd(lubnd(volseff_bnd, 1, 1), j_bnds(1) + 1, i_bnds(1))
      points_volsarea(ipoint, 2) = arseffz_bnd(lubnd(arseffz_bnd, 1, 1), j_bnds(1) + 1, i_bnds(1))

      points_loc(ipoint, 1) = bndblock_id
      points_loc(ipoint, 2) = i_bnds(1)
      points_loc(ipoint, 3) = j_bnds(1) + 1
      points_loc(ipoint, 4) = lubnd(volseff_bnd, 1, k)

    END IF

  END IF

  CALL linear_weights(weights, points, ipoint, point_dest)

END SUBROUTINE boundint_weights_xy



SUBROUTINE find_indices_jk(blocks, bounding, j, j_bnds, k, k_bnds,  myblock_id, bndblock_id)

  IMPLICIT NONE

  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  INTEGER, INTENT(in) :: j, k
  INTEGER, INTENT(in) :: myblock_id
  INTEGER, INTENT(inout) :: j_bnds(2), k_bnds(2)
  INTEGER, INTENT(inout) :: bndblock_id

  INTEGER :: iblock, nblocks
  INTEGER :: jj, kk

  nblocks = SIZE(blocks)

  bndblock_id = 100000000

  DO iblock = 1, nblocks
    IF (bounding(iblock) .EQV. .TRUE.) THEN
      IF (ALL((/blocks(myblock_id)%y2(j + 1) .gt. blocks(iblock)%y2(1), &
                blocks(myblock_id)%y2(j) .lt. blocks(iblock)%y2(SIZE(blocks(iblock)%y2)), &
                blocks(myblock_id)%z2(k + 1) .gt. blocks(iblock)%z2(1), &
                blocks(myblock_id)%z2(k) .lt. blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) &
               /))) THEN
        j_bnds(1) = 1
        DO jj = 1, SIZE(blocks(iblock)%y2) - 1
          IF (blocks(iblock)%y2(jj) .gt. blocks(myblock_id)%y2(j)) EXIT
          j_bnds(1) = jj
        END DO
        DO jj = j_bnds(1), SIZE(blocks(iblock)%y2) - 1
          j_bnds(2) = jj
          IF (blocks(iblock)%y2(jj + 1) .ge. blocks(myblock_id)%y2(j + 1)) EXIT
        END DO

        k_bnds(1) = 1
        DO kk = 1, SIZE(blocks(iblock)%z2) - 1
          IF (blocks(iblock)%z2(kk) .gt. blocks(myblock_id)%z2(k)) EXIT
          k_bnds(1) = kk
        END DO
        DO kk = k_bnds(1), SIZE(blocks(iblock)%z2) - 1
          k_bnds(2) = kk
          IF (blocks(iblock)%z2(kk + 1) .ge. blocks(myblock_id)%z2(k + 1)) EXIT
        END DO

        bndblock_id = iblock
        EXIT
      END IF
    END IF
  END DO

END SUBROUTINE find_indices_jk


SUBROUTINE find_indices_ik(blocks, bounding, i, i_bnds, k, k_bnds,  myblock_id, bndblock_id)

  IMPLICIT NONE

  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  INTEGER, INTENT(in) :: i, k
  INTEGER, INTENT(in) :: myblock_id
  INTEGER, INTENT(inout) :: i_bnds(2), k_bnds(2)
  INTEGER, INTENT(inout) :: bndblock_id

  INTEGER :: iblock, nblocks
  INTEGER :: ii, kk

  nblocks = SIZE(blocks)

  bndblock_id = 100000000

  DO iblock = 1, nblocks
    IF (bounding(iblock) .EQV. .TRUE.) THEN
      IF (ALL((/blocks(myblock_id)%x2(i + 1) .gt. blocks(iblock)%x2(1), &
                blocks(myblock_id)%x2(i) .lt. blocks(iblock)%x2(SIZE(blocks(iblock)%x2)), &
                blocks(myblock_id)%z2(k + 1) .gt. blocks(iblock)%z2(1), &
                blocks(myblock_id)%z2(k) .lt. blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) &
               /))) THEN
        i_bnds(1) = 1
        DO ii = 1, SIZE(blocks(iblock)%x2) - 1
          IF (blocks(iblock)%x2(ii) .gt. blocks(myblock_id)%x2(i)) EXIT
          i_bnds(1) = ii
        END DO
        DO ii = i_bnds(1), SIZE(blocks(iblock)%x2) - 1
          i_bnds(2) = ii
          IF (blocks(iblock)%x2(ii + 1) .ge. blocks(myblock_id)%x2(i + 1)) EXIT
        END DO

        k_bnds(1) = 1
        DO kk = 1, SIZE(blocks(iblock)%z2) - 1
          IF (blocks(iblock)%z2(kk) .gt. blocks(myblock_id)%z2(k)) EXIT
          k_bnds(1) = kk
        END DO
        DO kk = k_bnds(1), SIZE(blocks(iblock)%z2) - 1
          k_bnds(2) = kk
          IF (blocks(iblock)%z2(kk + 1) .ge. blocks(myblock_id)%z2(k + 1)) EXIT
        END DO

        bndblock_id = iblock
        EXIT
      END IF
    END IF
  END DO

END SUBROUTINE find_indices_ik


SUBROUTINE find_indices_ij(blocks, bounding, i, i_bnds, j, j_bnds,  myblock_id, bndblock_id)

  IMPLICIT NONE

  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  INTEGER, INTENT(in) :: i, j
  INTEGER, INTENT(in) :: myblock_id
  INTEGER, INTENT(inout) :: i_bnds(2), j_bnds(2)
  INTEGER, INTENT(inout) :: bndblock_id

  INTEGER :: iblock, nblocks
  INTEGER :: ii, jj

  nblocks = SIZE(blocks)

  bndblock_id = 100000000

  DO iblock = 1, nblocks
    IF (bounding(iblock) .EQV. .TRUE.) THEN
      IF (ALL((/blocks(myblock_id)%x2(i + 1) .gt. blocks(iblock)%x2(1), &
                blocks(myblock_id)%x2(i) .lt. blocks(iblock)%x2(SIZE(blocks(iblock)%x2)), &
                blocks(myblock_id)%y2(j + 1) .gt. blocks(iblock)%y2(1), &
                blocks(myblock_id)%y2(j) .lt. blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) &
               /))) THEN
        i_bnds(1) = 1
        DO ii = 1, SIZE(blocks(iblock)%x2) - 1
          IF (blocks(iblock)%x2(ii) .gt. blocks(myblock_id)%x2(i)) EXIT
          i_bnds(1) = ii
        END DO
        DO ii = i_bnds(1), SIZE(blocks(iblock)%x2) - 1
          i_bnds(2) = ii
          IF (blocks(iblock)%x2(ii + 1) .ge. blocks(myblock_id)%x2(i + 1)) EXIT
        END DO

        j_bnds(1) = 1
        DO jj = 1, SIZE(blocks(iblock)%y2) - 1
          IF (blocks(iblock)%y2(jj) .gt. blocks(myblock_id)%y2(j)) EXIT
          j_bnds(1) = jj
        END DO
        DO jj = j_bnds(1), SIZE(blocks(iblock)%y2) - 1
          j_bnds(2) = jj
          IF (blocks(iblock)%y2(jj + 1) .ge. blocks(myblock_id)%y2(j + 1)) EXIT
        END DO

        bndblock_id = iblock
        EXIT
      END IF
    END IF
  END DO

END SUBROUTINE find_indices_ij



SUBROUTINE linear_weights(weights, points, ipoint, point_dest)

  IMPLICIT NONE

  REAL(Realkind), INTENT(inout) :: weights(3)
  REAL(Realkind), INTENT(in) :: points(3, 2)
  INTEGER, INTENT(in) :: ipoint
  REAL(Realkind), INTENT(in) :: point_dest(2)
  REAL(Realkind) :: mat(3, 3), mat_inv(3, 3), b_vecs(3, ipoint), b(3), w_matrix(ipoint, 3)

  INTEGER :: i, j

  REAL(Realkind), PARAMETER :: eps=1e-20
  
  weights(:) = 0.0
  b_vecs(:,:) = 0.0

  IF (ipoint .EQ. 1) THEN
    weights(1) = 1.0
  ELSE
    mat(:, 1) = 1.0
    DO i = 1, ipoint
      b_vecs(i, i) = 1.0
      DO j = 2, 3
        mat(i, j) = points(i, j - 1)
      END DO
    END DO
    IF (ipoint .EQ. 2) THEN
      mat(3, 2) = points(1, 1) + 1.0
      mat(3, 3) = points(1, 2) - (points(1, 1) - points(2, 1)) / (points(1, 2) - points(2, 2) + eps)
      b_vecs(3, 1) = 1.0
    END IF

    CALL inv(mat, mat_inv)

    w_matrix = TRANSPOSE(MATMUL(mat_inv, b_vecs))

    b(1) = 1.0
    b(2) = point_dest(1)
    b(3) = point_dest(2)
   
    weights(1:ipoint) = MATMUL(w_matrix, b)
!    weights(1) = 1.0
!    weights(2:3) = 0.0
!  DO i = 1, 3
!    IF (weights(i) .GT. 1) WRITE(*,*) weights, mat, id, ipoint, points, b
!  END DO
  END IF

END SUBROUTINE linear_weights


SUBROUTINE inv(A, Ainv)
 
  IMPLICIT NONE
 
  REAL(Realkind), INTENT(in) :: A(:,:)
  REAL(Realkind), INTENT(inout) :: Ainv(:,:)

  REAL(Realkind) :: work(size(A,1)) ! work array for LAPACK
  INTEGER :: n, info, ipiv(size(A,1))     ! pivot indices

  ! Store A in Ainv to prevent it from being overwritten by LAPACK
  Ainv = A
  n = size(A,1)
  ! (S)DGETRF computes an LU factorization of a general M-by-N matrix A
  ! using partial pivoting with row interchanges.
  SELECT CASE (Realkind)

  CASE(4) !Single precision
     CALL SGETRF(n, n, Ainv, n, ipiv, info)
     IF (info .ne. 0) stop 'Matrix is numerically singular!'
     CALL SGETRI(n, Ainv, n, ipiv, work, n, info)
     IF (info .ne. 0) stop 'Matrix inversion failed!'
  CASE(8) !Double precision
     CALL DGETRF(n, n, Ainv, n, ipiv, info)
     IF (info .ne. 0) stop 'Matrix is numerically singular!'
     CALL DGETRI(n, Ainv, n, ipiv, work, n, info)
     IF (info .ne. 0) stop 'Matrix inversion failed!'
  END SELECT

END SUBROUTINE inv


FUNCTION lubnd(arr, dim, l) RESULT(i)

  IMPLICIT NONE

  REAL(Realkind), INTENT(in) :: arr(:,:,:)
  INTEGER, INTENT(in) :: dim, l
  INTEGER i

  IF (l .EQ. 1) THEN
    i = LBOUND(arr, dim)
  ELSE IF (l .EQ. -1) THEN
    i = UBOUND(arr, dim)
  ELSE
    i = -1 !Error
  END IF

END FUNCTION lubnd

END MODULE Gradient_block_interpolation_Mod
