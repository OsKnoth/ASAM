  program PreGrid

!---------------------------------------
!--  Preprozessor for grid generation
!---------------------------------------

     implicit none

     TYPE block
        INTEGER ib 
        INTEGER ixa, ixe
        INTEGER iya, iye
        INTEGER ref
        TYPE (block), POINTER :: next
     END TYPE block
 
     TYPE (block), POINTER :: fein
     TYPE (block), POINTER :: neu
     TYPE (block), POINTER :: start_fein

     character*50 :: InpFile, OutFile
     integer, parameter :: inr = 72, onr = 73
     integer, parameter :: nSqrtBand = 2

     integer, allocatable :: grid(:,:)
     integer :: nx, ny, nc, nxa, nxe, nya, nye
     integer :: nxa1, nxe1, nya1, nye1
     integer :: MetRef, iRef, MRef, RefLoc
     integer :: ixa, ixe, iya, iye, ix, iy  
     integer :: nhmin, nh, nxy, nhstp
     integer :: ios, narg, pos, pos1, iargc
     integer :: lx, ly, ibx, iby, irx, iry
     integer :: NHeap, nBlocks, iBlock, nb, ib

     real*8 :: rHeap

!---------------------------------
!--  Open Input and Output Files
!---------------------------------

     narg = iargc()

     if (narg >= 2)  then                ! Set file names
        CALL getarg(1,InpFile)
        CALL getarg(2,OutFile)
     else if (narg == 1)  then
        CALL getarg(1,InpFile)

        InpFile = adjustl(InpFile)
        pos     = len(InpFile) 

        pos1    = index(InpFile,' ')-1
        if (pos1 >= 1)  pos = pos1
        pos1    = index(InpFile,'.')-1
        if (pos1 >= 1)  pos = pos1

        OutFile = InpFile(1:pos)//'.grid'
     else
        write(*,*) ' Input File (*.preg):'
        read (*,*)  InpFile

        InpFile = adjustl(InpFile)
        pos     = len(InpFile) 

        pos1    = index(InpFile,' ')-1
        if (pos1 >= 1)  pos = pos1
        pos1    = index(InpFile,'.')-1
        if (pos1 >= 1)  pos = pos1

        OutFile = InpFile(1:pos)

        write(*,*) ' Output File (*.grid):',InpFile(1:pos)//'.grid',' ?'
        read (*,*)  OutFile
        if (OutFile == ' ')  OutFile = InpFile(1:pos)//'.grid'
     end if

     InpFile = adjustl(InpFile)          ! Input File
     open (inr, file=trim(InpFile), status='old',iostat=ios)

     if (ios /= 0) then 
        open (inr, file=trim(InpFile)//'.preg', status='old',iostat=ios)
        if (ios /= 0)  then
           write(*,*) ' Input File >',trim(InpFile),'< not found!'
           stop ' PreGrid: Input File not found!'
        else
           write(*,*) ' Input File :',trim(InpFile)//'.preg'
        end if
     else
        write(*,*) ' Input File :',trim(InpFile)
     end if

     OutFile = adjustl(OutFile)          ! Output File
     open (onr, file=trim(OutFile), status='unknown')
     write(*,*) ' Output File :',trim(OutFile)
           
     NHeap = 25
     write(*,*) ' NHeap (maximum number of columns per block): [25] '
     read (*,*)  nHeap
  
!------------------------
!--  Set whole domain
!------------------------

     read(inr,*)  nxa1, nxe1
     read(inr,*)  nya1, nye1
     read(inr,*)  MetRef
     
     if (MetRef >= 0)  then
        MRef = 2**MetRef
        nxa  = nxa1 * MRef
        nxe  = nxe1 * MRef
        nya  = nya1 * MRef
        nye  = nye1 * MRef
     else
        MRef = 2**(-MetRef)
        nxa  = nxa1 / MRef
        nxe  = nxe1 / MRef
        nya  = nya1 / MRef
        nye  = nye1 / MRef

        if ((MRef*nxa /= nxa1) .OR. (MRef*nxe /= nxe1))  then
           write(*,8018) MetRef, nxa, nxa1, nxe, nxe1
           stop 'Inconsistency in domain refinement level (x-direction) !'
        end if

        if ((MRef*nya /= nya1) .OR. (MRef*nye /= nye1))  then
           write(*,8019) MetRef, nya, nya1, nye, nye1
           stop 'Inconsistency in domain refinement level (y-direction) !'
        end if

     end if

     write(onr,8001)  nxa, nxe
     write(onr,8002)  nya, nye
     write(onr,*) 
  
     read(inr,*)  
     read(inr,*)  nBlocks

!------------------------
!--  Loop for subdomains
!------------------------

     ALLOCATE (fein)
     NULLIFY  (fein%next)
     NULLIFY  (start_fein)

     ib = 0

     do iBlock=1,nBlocks

        read(inr,*)  
        read(inr,*)  ixa, ixe
        read(inr,*)  iya, iye
        read(inr,*)  iRef
  
        if (MetRef >= 0)  then
           ixa  = ixa * MRef
           ixe  = ixe * MRef
           iya  = iya * MRef
           iye  = iye * MRef
        else
           ixa  = ixa / MRef
           ixe  = ixe / MRef
           iya  = iya / MRef
           iye  = iye / MRef
        end if
     
        iRef = iRef - MetRef
        if (iRef >= 0)  then
           RefLoc = 2**iRef
           nx   = (ixe-ixa) * RefLoc
           ny   = (iye-iya) * RefLoc
        else
           RefLoc = 2**(-iRef)
           nx   = (ixe-ixa) / RefLoc
           ny   = (iye-iya) / RefLoc
        end if

        nc   = nx * ny

!--------------------------------------
!--  optimize decomposition of subblock
        
!-- take original block
        lx = 0
        if (nx*ny <= NHeap) then
           lx  = nx
           ibx = 1
           irx = nx

           ly  = ny
           iby = 1
           iry = ny
        
!-- whole fine block greater than NHeap
        else if ((iRef>=1) .AND. (RefLoc*RefLoc>=NHeap)) then
           lx  = nx
           ibx = 1
           irx = nx

           ly  = ny
           iby = 1
           iry = ny
           write(*,8017) nHeap, iBlock, iRef, RefLoc*RefLoc
        
!-- small dimension
        else if (nx == 1) then
           lx  = 1
           ibx = 1
           irx = 1

           ly  = min(NHeap/lx,ny)
           iby = (ny-1)/ly + 1
           iry = ny - ly*(iby-1)

           if (nxe-nxa >= 2)  then
              write(*,8015) iBlock, nx
           end if
        else if (ny == 1) then
           ly  = 1
           iby = 1
           iry = 1

           lx  = min(NHeap/ly,nx)
           ibx = (nx-1)/lx + 1
           irx = nx - lx*(ibx-1)

           if (nye-nya >= 2)  then
              write(*,8016) iBlock, ny
           end if
        else if (nx == 2) then
           lx  = 2
           ibx = 1
           irx = 2

           ly  = min(NHeap/lx,ny)
           iby = (ny-1)/ly + 1
           iry = ny - ly*(iby-1)
        else if (ny == 2) then
           ly  = 2
           iby = 1
           iry = 2

           lx  = min(NHeap/ly,nx)
           ibx = (nx-1)/lx + 1
           irx = nx - lx*(ibx-1)
        else if (nx == 3) then
           lx  = 3
           ibx = 1
           irx = 3

           ly  = min(NHeap/lx,ny)
           iby = (ny-1)/ly + 1
           iry = ny - ly*(iby-1)
        else if (ny == 3) then
           ly  = 3
           iby = 1
           iry = 3

           lx  = min(NHeap/ly,nx)
           ibx = (nx-1)/lx + 1
           irx = nx - lx*(ibx-1)
        
!-- large dimension
        else
           rHeap = nHeap + 1.e-2
           nh    = sqrt(rHeap)
           nhmin = max(nh-nSqrtBand,4)
           nb    = nx * ny

           nhstp = 1
           if (iRef >= 1)  then            ! save: nhstp = RefLoc * xxxx
              nhstp = RefLoc
              nhmin = max(nhmin/RefLoc,1)
              nhmin = nhmin * RefLoc
           end if
          
           do ix=nhmin,nh+nhstp,nhstp 
              ibx = (nx-1)/ix + 1
              iy   = min(NHeap/ix,ny)
              if (iRef >= 1)  then         ! save: iy = RefLoc * xxxx
                 iy = iy / RefLoc
                 iy = iy * RefLoc
              end if
              if (iy <= 0)  exit
              iby = (ny-1)/iy + 1
            
              nxy = ibx * iby
              if (nb > nxy)  then
                 nb = nxy
                 lx = ix
                 ly = iy
              end if
           end do
          
!--------------------------------------
           do iy=nhmin,nh+nhstp,nhstp 
              iby = (ny-1)/iy + 1
              ix  = min(NHeap/iy,nx)
              if (iRef >= 1)  then         ! save: ix = RefLoc * xxxx
                 ix = ix / RefLoc
                 ix = ix * RefLoc
              end if
              if (ix <= 0)  exit
              ibx = (nx-1)/ix  + 1
            
              nxy = ibx * iby
              if (nb > nxy)  then
                 nb = nxy
                 lx = ix 
                 ly = iy
              end if
           end do

           ibx = (nx-1)/lx + 1
           irx = nx - lx*(ibx-1)

           iby = (ny-1)/ly + 1
           iry = ny - ly*(iby-1)
        end if

!--------------------------------------
!--- modification of blocks with irx=1 or iry=1
!--- (set length of boundary blocks to lx+1 or ly+1)
!
        if (irx == 1)  then
           irx = lx  + 1
           ibx = ibx - 1
        end if

        if (iry == 1)  then
           iry = ly  + 1
           iby = iby - 1
        end if

!--------------------------------------
!--- check consistency
        if (iRef >= 0)  then
           ix = ixa + (lx*(ibx-1)+irx) / RefLoc
           iy = iya + (ly*(iby-1)+iry) / RefLoc
        else
           ix = ixa + (lx*(ibx-1)+irx) * RefLoc
           iy = iya + (ly*(iby-1)+iry) * RefLoc
        end if

        if (ix /= ixe)  then
           write(*,8011) iBlock, ib, ix, ixe
           stop 'Inconsistency in block definition (in x-direction)!'
        end if

        if (iy /= iye)  then
           write(*,8012) iBlock, ib, iy, iye
           stop 'Inconsistency in block definition (in y-direction)!'
        end if

!--------------------------------------
!--------------------------------------
!--  define new partitions             
        do ix=1,ibx-1                   ! full blocks
           do iy=1,iby-1
              ib = ib + 1

              fein%ib  = iBlock
              if (iRef >= 0)  then
                 fein%ixa = ixa + lx*(ix-1) / RefLoc
                 fein%ixe = ixa + lx*ix     / RefLoc
                 fein%iya = iya + ly*(iy-1) / RefLoc
                 fein%iye = iya + ly*iy     / RefLoc
              else
                 fein%ixa = ixa + lx*(ix-1) * RefLoc
                 fein%ixe = ixa + lx*ix     * RefLoc
                 fein%iya = iya + ly*(iy-1) * RefLoc
                 fein%iye = iya + ly*iy     * RefLoc
              end if
              fein%ref = iRef

!--- allocate next block 
              ALLOCATE (neu)
              NULLIFY  (neu%next)

              IF (.not.associated(start_fein)) THEN
                 start_fein => fein
              END if

              fein%next => neu
              fein      => neu
           end do

           ib = ib  + 1                  ! y-boundary block
           iy = iby 

           fein%ib  = iBlock
           if (iRef >= 0)  then
              fein%ixa = ixa + lx*(ix-1) / RefLoc
              fein%ixe = ixa + lx*ix     / RefLoc
              fein%iya = iya + ly*(iy-1) / RefLoc
              fein%iye = iye 
           else
              fein%ixa = ixa + lx*(ix-1) * RefLoc
              fein%ixe = ixa + lx*ix     * RefLoc
              fein%iya = iya + ly*(iy-1) * RefLoc
              fein%iye = iye 
           end if
           fein%ref = iRef

!--- allocate next block 
           ALLOCATE (neu)
           NULLIFY  (neu%next)

           IF (.not.associated(start_fein)) THEN
              start_fein => fein
           END IF

           fein%next => neu
           fein      => neu
        end do

!--------------------------------------
        ix = ibx
        do iy=1,iby-1                  ! x-boundary blocks
           ib = ib + 1

           fein%ib  = iBlock
           if (iRef >= 0)  then
              fein%ixa = ixa + lx*(ix-1) / RefLoc
              fein%ixe = ixe 
              fein%iya = iya + ly*(iy-1) / RefLoc
              fein%iye = iya + ly*iy     / RefLoc
           else
              fein%ixa = ixa + lx*(ix-1) * RefLoc
              fein%ixe = ixe 
              fein%iya = iya + ly*(iy-1) * RefLoc
              fein%iye = iya + ly*iy     * RefLoc
           end if
           fein%ixe = ixe 
           fein%ref = iRef

!--- allocate next block 
           ALLOCATE (neu)
           NULLIFY  (neu%next)

           IF (.not.associated(start_fein)) THEN
              start_fein => fein
           END IF

           fein%next => neu
           fein      => neu
        end do
        ib = ib + 1                     ! y-boundary block
        iy = iby

        fein%ib  = iBlock
        if (iRef >= 0)  then
           fein%ixa = ixa + lx*(ix-1) / RefLoc
           fein%ixe = ixe 
           fein%iya = iya + ly*(iy-1) / RefLoc
           fein%iye = iye
        else
           fein%ixa = ixa + lx*(ix-1) * RefLoc
           fein%ixe = ixe 
           fein%iya = iya + ly*(iy-1) * RefLoc
           fein%iye = iye
        end if
        fein%ref = iRef

!--- allocate next block 
        ALLOCATE (neu)
        NULLIFY  (neu%next)

        IF (.not.associated(start_fein)) THEN
           start_fein => fein
        END IF

        fein%next => neu
        fein      => neu
                 
     end do

!-----------------------------------
!--  Write new grid decomposition
!-----------------------------------

     nb = ib                            ! number of blocks
     write(onr,8003)  nb

     allocate (grid(nya+1:nye, nxa+1:nxe))
     grid = 0

     fein      => start_fein
     do ib = 1,nb
        iBlock = fein%ib  
        ixa    = fein%ixa 
        ixe    = fein%ixe 
        iya    = fein%iya 
        iye    = fein%iye 
        iRef   = fein%ref 
        
        write(onr,8004)  ib, iBlock
        write(onr,8005)  ixa, ixe
        write(onr,8006)  iya, iye
        write(onr,8007)  iRef

        grid(iya+1:iye, ixa+1:ixe) = grid(iya+1:iye, ixa+1:ixe) + 1

        fein      => fein%next

     end do

     close (inr)
     close (onr)

!-----------------------------------
!--  Check consistency        
!-----------------------------------
     
     do ix=nxa+1,nxe
        do iy=nya+1,nye
           if (grid(iy,ix) <= 0)  then
              write(*,8013)  iy, ix, grid(iy,ix)
           end if
           if (grid(iy,ix) >= 2)  then
              write(*,8014)  iy, ix, grid(iy,ix)
           end if
        end do
     end do

!-----------------------------------
!--  Format statements
!-----------------------------------
!
!--  output file
8001  format( i5,i5,'    ! Domain: X-direction ... nxa, nxe')
8002  format( i5,i5,'    !         Y-direction ... nya, nye')

8003  format( i5,5x,'    ! Number of blocks    ... nb      ')
8004  format( 5x,5x,'    ! Block ',i4,' (',i3,')')        

8005  format( i5,i5,'    !        ixa, ixe')
8006  format( i5,i5,'    !        iya, iye')
8007  format( i5,5x,'    !        iRef')

!-----------------------------------
!--  error messages

8011  format( 1x//' Inconsistency in block ',i4,' (',i3,') '/ &
&                 ' X-direction:  IX=',i4,'not equal to  IXE=',i4/1x)
8012  format( 1x//' Inconsistency in block ',i4,' (',i3,') '/ &
&                 ' Y-direction:  IY=',i4,'not equal to  IYE=',i4/1x)

8013  format( 1x//' Inconsistency in Block Structure:'/ &
&                 '    (',i3,',',i3,') ... not included. GRID=',i3/1x)
8014  format( 1x//' Inconsistency in Block Structure:'/ &
&                 '    (',i3,',',i3,') ...  2x included. GRID=',i3/1x)

8015  format( 1x//' Warning: X-dimension of Block ',i4,' = ',i2 /  &
&                 '          AdvOrd = 2 is not realized !'      /1x)
8016  format( 1x//' Warning: Y-dimension of Block ',i4,' = ',i2 /  &
&                 '          AdvOrd = 2 is not realized !'      /1x)

8017  format( 1x//' Warning: nHeap =',i4,'to small for Block',i4 /  &
&                 '          iRef  =',i3,' ==> ',i5,                &
&                 ' columns have to be in block',i4,' !'      /1x)

8018  format( 1x//' Inconsistency in domain refinement level ',      &
&                 ' ... MetRef =',i2                               / &
&                 ' X-direction:  nxa =',i4,', nxa1 (given) =',i4, / &
&                 ',    nxe =',i4,', nxe1 (given) =',i4/1x)

8019  format( 1x//' Inconsistency in domain refinement level ',      &
&                 ' ... MetRef =',i2                               / &
&                 ' Y-direction:  nya =',i4,', nya1 (given) =',i4, / &
&                 ',    nye =',i4,', nye1 (given) =',i4/1x)

  end program PreGrid
