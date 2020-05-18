MODULE Grid_Mod

  USE Kind_Mod
  USE FLoor_Mod
  USE Parallel_Mod
  USE BoundaryCond_Mod
  USE Physics_Mod
  USE Control_Mod
  USE Output_Mod

  IMPLICIT NONE

! Metis Parameter
  INTEGER, PRIVATE, PARAMETER :: METIS_NOPTIONS=40

  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_PTYPE=     1
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_OBJTYPE=   2
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_CTYPE=     3
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_IPTYPE=    4
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_RTYPE=     5
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_DBGLVL=    6
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_NITER=     7
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_NCUTS=     8
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_SEED=      9 
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_NO2HOP=   10
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_MINCONN=  11
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_CONTIG=   12
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_COMPRESS= 13
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_CCORDER=  14
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_PFACTOR=  15
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_NSEPS=    16
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_UFACTOR=  17
  INTEGER, PRIVATE, PARAMETER :: METIS_OPTION_NUMBERING=18
! Partitioning Schemes 
  INTEGER, PRIVATE, PARAMETER :: METIS_PTYPE_RB=0       ! Multilevel recursive sectioning
  INTEGER, PRIVATE, PARAMETER :: METIS_PTYPE_KWAY=1     ! Multilevel k-way partitioning            
! Objective  
  INTEGER, PRIVATE, PARAMETER :: METIS_OBJTYPE_CUT=0    ! Edge-cut minimization
  INTEGER, PRIVATE, PARAMETER :: METIS_OBJTYPE_VOL=1    ! Total communicatio volume minimization
! Matching scheme  
  INTEGER, PRIVATE, PARAMETER :: METIS_CTYPE_RM=0       ! Edge-cut minimization
  INTEGER, PRIVATE, PARAMETER :: METIS_CTYPE_SHEM=1     ! Total communicatio volume minimization
! Initial partitioning
  INTEGER, PRIVATE, PARAMETER :: METIS_IPTYPE_GROW=0    ! Edge-cut minimization
  INTEGER, PRIVATE, PARAMETER :: METIS_IPTYPE_RANDOM=1  ! Edge-cut minimization
  INTEGER, PRIVATE, PARAMETER :: METIS_IPTYPE_EDGE=2    ! Edge-cut minimization
  INTEGER, PRIVATE, PARAMETER :: METIS_IPTYPE_NODE=3    ! Edge-cut minimization
! Refinement  
  INTEGER, PRIVATE, PARAMETER :: METIS_RTYPE_FM=0    ! 
  INTEGER, PRIVATE, PARAMETER :: METIS_RTYPE_GREEDY=1    ! 
  INTEGER, PRIVATE, PARAMETER :: METIS_RTYPE_SEP2SIDED=2    ! 
  INTEGER, PRIVATE, PARAMETER :: METIS_RTYPE_SEP1SIDED=3    ! 

  

CONTAINS

SUBROUTINE refine_grid

!---------------------------------------------------------------
!   Refine blocks

END SUBROUTINE refine_grid

SUBROUTINE ngbr_bounds
!
!---------------------------------------------------------------
!---  Determine Neighbor boundary data
!---------------------------------------------------------------

    INTEGER :: ib1, in
 
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO in=1,AnzahlNachbar
        ib1=Nachbars(in)%ib
        Nachbars(in)%ix0=MAX(igx0,Floor(ib1)%igx0)*2.e0**RefineX
        Nachbars(in)%ix1=MIN(igx1,Floor(ib1)%igx1)*2.e0**RefineX
        Nachbars(in)%iy0=MAX(igy0,Floor(ib1)%igy0)*2.e0**RefineY
        Nachbars(in)%iy1=MIN(igy1,Floor(ib1)%igy1)*2.e0**RefineY
        Nachbars(in)%iz0=MAX(igz0,Floor(ib1)%igz0)*2.e0**RefineZ
        Nachbars(in)%iz1=MIN(igz1,Floor(ib1)%igz1)*2.e0**RefineZ
        Nachbars(in)%iNx0=MAX(igx0,Floor(ib1)%igx0)*2.e0**Floor(ib1)%RefineX
        Nachbars(in)%iNx1=MIN(igx1,Floor(ib1)%igx1)*2.e0**Floor(ib1)%RefineX
        Nachbars(in)%iNy0=MAX(igy0,Floor(ib1)%igy0)*2.e0**Floor(ib1)%RefineY
        Nachbars(in)%iNy1=MIN(igy1,Floor(ib1)%igy1)*2.e0**Floor(ib1)%RefineY
        Nachbars(in)%iNz0=MAX(igz0,Floor(ib1)%igz0)*2.e0**Floor(ib1)%RefineZ
        Nachbars(in)%iNz1=MIN(igz1,Floor(ib1)%igz1)*2.e0**Floor(ib1)%RefineZ
        IF (Nachbars(in)%nType=='iw') THEN
          Nachbars(in)%ix0=Floor(ib1)%igx1*2.e0**RefineX
          Nachbars(in)%ix1=Floor(ib1)%igx1*2.e0**RefineX
          Nachbars(in)%iNx0=Floor(ib1)%igx1*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%iNx1=Floor(ib1)%igx1*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%ixO=ix0
          Nachbars(in)%ixI=ix0+1
          Nachbars(in)%iNxO=Nachbars(in)%iNx1+1
          Nachbars(in)%iNxI=Nachbars(in)%iNx1
        ELSE IF (Nachbars(in)%nType=='ie') THEN
          Nachbars(in)%ix0=Floor(ib1)%igx0*2.e0**RefineX
          Nachbars(in)%ix1=Floor(ib1)%igx0*2.e0**RefineX
          Nachbars(in)%iNx0=Floor(ib1)%igx0*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%iNx1=Floor(ib1)%igx0*2.e0**Floor(ib1)%RefineX
          Nachbars(in)%ixO=ix1+1
          Nachbars(in)%ixI=ix1
          Nachbars(in)%iNxO=Nachbars(in)%iNx0
          Nachbars(in)%iNxI=Nachbars(in)%iNx0+1
        ELSE IF (Nachbars(in)%nType=='is') THEN
          Nachbars(in)%iy0=Floor(ib1)%igy1*2.e0**RefineY
          Nachbars(in)%iy1=Floor(ib1)%igy1*2.e0**RefineY
          Nachbars(in)%iNy0=Floor(ib1)%igy1*2.e0**Floor(ib1)%RefineY
          Nachbars(in)%iNy1=Floor(ib1)%igy1*2.e0**Floor(ib1)%RefineY
          Nachbars(in)%iyO=iy0
          Nachbars(in)%iyI=iy0+1
          Nachbars(in)%iNyO=Nachbars(in)%iNy1+1
          Nachbars(in)%iNyI=Nachbars(in)%iNy1
        ELSE IF (Nachbars(in)%nType=='in') THEN
          Nachbars(in)%iy0=Floor(ib1)%igy0*2.e0**RefineY
          Nachbars(in)%iy1=Floor(ib1)%igy0*2.e0**RefineY
          Nachbars(in)%iNy0=Floor(ib1)%igy0*2.e0**Floor(ib1)%RefineY
          Nachbars(in)%iNy1=Floor(ib1)%igy0*2.e0**Floor(ib1)%RefineY
          Nachbars(in)%iyO=iy1+1
          Nachbars(in)%iyI=iy1
          Nachbars(in)%iNyO=Nachbars(in)%iNy0
          Nachbars(in)%iNyI=Nachbars(in)%iNy0+1
        ELSE IF (Nachbars(in)%nType=='ib') THEN
          Nachbars(in)%iz0=Floor(ib1)%igz1*2.e0**RefineZ
          Nachbars(in)%iz1=Floor(ib1)%igz1*2.e0**RefineZ
          Nachbars(in)%iNz0=Floor(ib1)%igz1*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%iNz1=Floor(ib1)%igz1*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%izO=iz0
          Nachbars(in)%izI=iz0+1
          Nachbars(in)%iNzO=Nachbars(in)%iNz1+1
          Nachbars(in)%iNzI=Nachbars(in)%iNz1
        ELSE IF (Nachbars(in)%nType=='it') THEN
          Nachbars(in)%iz0=Floor(ib1)%igz0*2.e0**RefineZ
          Nachbars(in)%iz1=Floor(ib1)%igz0*2.e0**RefineZ
          Nachbars(in)%iNz0=Floor(ib1)%igz0*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%iNz1=Floor(ib1)%igz0*2.e0**Floor(ib1)%RefineZ
          Nachbars(in)%izO=iz1+1
          Nachbars(in)%izI=iz1
          Nachbars(in)%iNzO=Nachbars(in)%iNz0
          Nachbars(in)%iNzI=Nachbars(in)%iNz0+1
        END IF
      END DO      ! in
    END DO         ! ib

 END SUBROUTINE ngbr_bounds

SUBROUTINE getbufsize
      
      INTEGER :: ib1, ib_n, il, in, ip, ip_n
      INTEGER :: area
      TYPE(gliedT), POINTER :: SingleG

      
!     Determine nlink

      ALLOCATE(chain(0:NumProcs-1))
!     
      chain(:)%nlink=0 
      DO ib1=1,nb                                    
         ip=blMPI(ib1)%proc                        
         IF (MyId == ip) THEN                     
            DO in=1,Floor(ib1)%AnzahlNachbar    
               ib_n=Floor(ib1)%Nachbars(in)%ib   
               ip_n=blMPI(ib_n)%proc
               IF (MyId /= ip_n) THEN  
                  chain(ip_n)%nlink=chain(ip_n)%nlink+1
                  Floor(ib1)%Nachbars(in)%ntype(1:1)='p'
               END IF
            END DO                                  
         END IF
      END DO                                        

!     Allocate space
      NumberNeiProc=0
      DO ip_n=0,NumProcs-1
         ALLOCATE(chain(ip_n)%glied(chain(ip_n)%nlink))
         IF (chain(ip_n)%nlink > 0) THEN
            NumberNeiProc=NumberNeiProc+1
         END IF 
      END DO

!     Determine neighbours

      chain(:)%nlink=0 
      DO ib1=1,nb
         ip=blMPI(ib1)%proc
         IF (MyId == ip) THEN                     
            DO in=1,Floor(ib1)%AnzahlNachbar
               ib_n=Floor(ib1)%Nachbars(in)%ib
               ip_n=blMPI(ib_n)%proc
               IF (MyId /= ip_n) THEN  
                  chain(ip_n)%nlink=chain(ip_n)%nlink+1
                  il=chain(ip_n)%nlink
                  chain(ip_n)%glied(il)%OwnBlock=ib1
                  chain(ip_n)%glied(il)%NeiBlock(1)=ib_n
                  chain(ip_n)%glied(il)%NeiBlock(2)=in
                  chain(ip_n)%glied(il)%side=Floor(ib1)%Nachbars(in)%ntype
                  chain(ip_n)%glied(il)%RefOwn=Floor(ib1)%Refine
                  chain(ip_n)%glied(il)%RefNei=Floor(ib_n)%Refine
                  chain(ip_n)%glied(il)%IncrX=Floor(ib1)%Nachbars(in)%IncrX
                  chain(ip_n)%glied(il)%IncrY=Floor(ib1)%Nachbars(in)%IncrY
                  chain(ip_n)%glied(il)%IncrZ=Floor(ib1)%Nachbars(in)%IncrZ
                  chain(ip_n)%glied(il)%Nachbar=>Floor(ib1)%Nachbars(in)
               END IF
            END DO                                  
         END IF
      END DO                                        

!     Sort chain

      DO ip=0,NumProcs-1
         IF (ip < MyId) THEN
            CALL sortOwn(chain(ip)%glied)
         ELSE IF (ip > MyId) THEN
            CALL sortNei(chain(ip)%glied)
         END IF
      END DO

!     Fill chain

      DO ip=0,NumProcs-1                             
         chain(ip)%LenPack=0
         chain(ip)%LenUnPack=0
         DO il=1,chain(ip)%nlink                    
            SingleG=>chain(ip)%glied(il)
            ib1=SingleG%OwnBlock
            ib_n=SingleG%NeiBlock(1)
            in=SingleG%NeiBlock(2)

!           -- Default coordinates --
            SingleG%ix0U=Floor(ib1)%Nachbars(in)%ix0 ! coord. of bound. piece
            SingleG%ix1U=Floor(ib1)%Nachbars(in)%ix1
            SingleG%iy0U=Floor(ib1)%Nachbars(in)%iy0
            SingleG%iy1U=Floor(ib1)%Nachbars(in)%iy1
            SingleG%iz0U=Floor(ib1)%Nachbars(in)%iz0
            SingleG%iz1U=Floor(ib1)%Nachbars(in)%iz1

!           -- Adapt coordinates to sides and compute buffer sizes --

            IF (SingleG%side == 'pw' .OR. SingleG%side == 'pe') THEN
               IF (SingleG%RefNei < SingleG%RefOwn) THEN
                  area=     (SingleG%iy1U-SingleG%iy0U)/SingleG%IncrY
                  area=area*((SingleG%iz1U-SingleG%iz0U)/SingleG%IncrZ)
               ELSE
                  area=(SingleG%iy1U-SingleG%iy0U)*(SingleG%iz1U-SingleG%iz0U)
               END IF

            ELSE IF (SingleG%side == 'ps' .OR. SingleG%side == 'pn') THEN
               IF (SingleG%RefNei < SingleG%RefOwn) THEN
                  area=     (SingleG%ix1U-SingleG%ix0U)/SingleG%IncrX
                  area=area*((SingleG%iz1U-SingleG%iz0U)/SingleG%IncrZ)
               ELSE
                  area=(SingleG%ix1U-SingleG%ix0U)*(SingleG%iz1U-SingleG%iz0U)
               END IF

            ELSE IF (SingleG%side == 'pb' .OR. SingleG%side == 'pt') THEN
               IF (SingleG%RefNei < SingleG%RefOwn) THEN
                  area=     (SingleG%ix1U-SingleG%ix0U)/SingleG%IncrX
                  area=area*((SingleG%iy1U-SingleG%iy0U)/SingleG%IncrY)
               ELSE
                  area=(SingleG%ix1U-SingleG%ix0U)*(SingleG%iy1U-SingleG%iy0U)
               END IF

            END IF
            
            chain(ip)%LenUnPack = chain(ip)%LenUnPack + area            

         END DO                                      ! chainlinks
         
         chain(ip)%LenPack = chain(ip)%LenUnPack
          
      END DO                                         ! processors           
      MaxPutBuf = 0
      MaxGetBuf = 0
      DO ip=0,NumProcs-1
         MaxPutBuf = MAX(MaxPutBuf,chain(ip)%LenPack)
         MaxGetBuf = MAX(MaxGetBuf,chain(ip)%LenUnPack)
      END DO                                         ! processors

      RETURN

 CONTAINS
 SUBROUTINE sortOwn(glied)

   IMPLICIT NONE

   TYPE (gliedT), POINTER :: glied(:)
   TYPE (gliedT) :: Temp

   INTEGER :: i, j, len

   len=size(glied)

   DO i=1,len
     DO j=1,len-i
       IF (glied(j)%OwnBlock > glied(j+1)%OwnBlock) THEN
         Temp=glied(j)
         glied(j)=glied(j+1)
         glied(j+1)=Temp
       ELSE IF (glied(j)%OwnBlock == glied(j+1)%OwnBlock) THEN
         IF (glied(j)%NeiBlock(1) > glied(j+1)%NeiBlock(1)) THEN
           Temp=glied(j)
           glied(j)=glied(j+1)
           glied(j+1)=Temp
         ELSE IF (glied(j)%NeiBlock(1) == glied(j+1)%NeiBlock(1)) THEN
           IF (glied(j)%NeiBlock(2) < glied(j+1)%NeiBlock(2)) THEN  
             Temp=glied(j)
             glied(j)=glied(j+1)
             glied(j+1)=Temp
           END IF
         END IF
       END IF
     END DO
   END DO
 END SUBROUTINE sortOwn

 SUBROUTINE sortNei(glied)

   IMPLICIT NONE

   TYPE (gliedT), POINTER :: glied(:)
   TYPE (gliedT) :: Temp

   INTEGER :: i, j, len

   len=size(glied)

   DO i=1,len
     DO j=1,len-i
       IF (glied(j)%NeiBlock(1) > glied(j+1)%NeiBlock(1)) THEN
         Temp=glied(j)
         glied(j)=glied(j+1)
         glied(j+1)=Temp
       ELSE IF (glied(j)%NeiBlock(1) == glied(j+1)%NeiBlock(1)) THEN
         IF (glied(j)%OwnBlock > glied(j+1)%OwnBlock) THEN
           Temp=glied(j)
           glied(j)=glied(j+1)
           glied(j+1)=Temp
         ELSE IF (glied(j)%OwnBlock == glied(j+1)%OwnBlock) THEN
           IF (glied(j)%NeiBlock(2) > glied(j+1)%NeiBlock(2)) THEN  
             Temp=glied(j)
             glied(j)=glied(j+1)
             glied(j+1)=Temp
           END IF
         END IF
       END IF
     END DO
   END DO
 END SUBROUTINE sortNei

END SUBROUTINE getbufsize              

SUBROUTINE assignment
!
!=======================================================
!---   Initialization of Chemistry-Transport Model
!=======================================================
!

!-----------------------------------------------------------

       INTEGER :: VertOrder(nb+1), Edge(nb*MaxBlkNgbr)
       INTEGER :: VertWgt(nb), EdgeWgt(nb*MaxBlkNgbr)
       INTEGER :: part(nb)
       INTEGER :: ObjVal
       INTEGER, POINTER :: vsize(:)=>NULL()
       REAL(RealKind), POINTER :: tpwgts(:,:)=>NULL()
       REAL(RealKind), POINTER :: ubvec(:)=>NULL()
       INTEGER, POINTER :: options(:)=>NULL()

!      INTEGER :: options(5)
       INTEGER :: WgtFlag
       INTEGER :: NumFlag
       INTEGER :: nparts, edgecut
       INTEGER :: ib1
       INTEGER :: max_ib, ip

       INTEGER, PARAMETER :: puni = 76

       INTEGER :: i !Temporär

!-----------------------------------------------------------
!---  Determination of Partitions
!-----------------------------------------------------------
       IF (ivar > 5) ivar = 0

       IF (NumProcs == 1) THEN
          part(:) = 1
          
!--- METIS-partition
       ELSE IF (ivar == 0) THEN
       
!---  Set graph for block structure
          CALL graph_metis(MaxBlkNgbr,VertOrder,Edge,VertWgt,EdgeWgt)

!   IF (MyId==0) THEN
!     WRITE(400,*) 'nb ',nb,'NumEdge',VertOrder(nb+1)
!     WRITE(400,*) 'NumProc',NumProc
!     DO i=1,nb
!       WRITE(400,*) 'Len ',VertOrder(i+1)-VertOrder(i)
!       WRITE(400,*) 'i',i,'j',Edge(VertOrder(i)+1:VertOrder(i+1))+1
!       WRITE(400,*) 'v',VertWgt(i),'e',EdgeWgt(VertOrder(i)+1:VertOrder(i+1))
!     END DO  
!   END IF  
!---  call to METIS
          ALLOCATE(options(METIS_NOPTIONS))
          options=0
          options(METIS_OPTION_OBJTYPE)=METIS_OBJTYPE_VOL
          options(METIS_OPTION_CTYPE)=METIS_CTYPE_RM
          options(METIS_OPTION_IPTYPE)=METIS_IPTYPE_EDGE
          options(METIS_OPTION_RTYPE)=METIS_RTYPE_GREEDY
          options(METIS_OPTION_NO2HOP)=1
          options(METIS_OPTION_NCUTS)=1
          options(METIS_OPTION_NITER)=10
          options(METIS_OPTION_UFACTOR)=1
          options(METIS_OPTION_MINCONN)=1
          options(METIS_OPTION_CONTIG)=1
          options(METIS_OPTION_SEED)=12345
          options(METIS_OPTION_NUMBERING)=0
          options(METIS_OPTION_DBGLVL)=32
          WgtFlag=1
          NumFlag=1


          IF (NumProcs >= 2)  THEN
!            CALL METIS_PartGraphKway(nb,VertOrder,Edge,VertWgt,EdgeWgt, &
!                                     WgtFlag,NumFlag,NumProcs,options,edgecut,part)

    CALL METIS_PartGraphKway(nb,1,VertOrder,Edge &
                            ,VertWgt,vsize,EdgeWgt &
                            ,NumProcs &
                            ,tpwgts,ubvec,options &
                            ,ObjVal &
                            ,part)
    part=part+1
!         ELSE IF (NumProcs >= 2)  THEN
!            CALL METIS_PartGraphRecursive(nb,VertOrder,Edge,VertWgt,EdgeWgt, &
!                                          3,1,NumProcs,options,edgecut,part)
          END IF

!---  manual partition
       ELSE IF (ivar == 1) THEN    ! first nb/np blocks to proc 0, next: 1, ...
          max_ib=nb/NumProcs
          ip=1
          DO ib1=1,nb
             part(ib1)=ip
             IF (ib1 >= max_ib) THEN
                ip=ip+1                
                IF (ip == NumProcs) THEN
                   max_ib=nb
                ELSE
                   max_ib=max_ib+nb/NumProcs
                END IF
             END IF
          END DO

       ELSE IF (ivar == 2) THEN   ! 1st block to proc 0, 2nd: 1, ... , np'th: 0
          ip=1
          DO ib1=1,nb
             part(ib1)=ip
             ip=ip+1
             IF (ip == NumProcs+1) ip=1
          END DO

       ELSE IF (ivar == 3) THEN   ! 1st two blocks to proc 0, next two: 1, ...
          ip=1 
          DO ib1=1,nb
             part(ib1)=ip
             IF (MOD(ib1,2) == 0)  ip=ip+1
             IF (ip == NumProcs+1) ip=1
          END DO          
          
       ELSE IF (ivar == 4) THEN   ! 1st block: 0, 2nd: 2, ((nb+1)/2): 1, ...
          ip=1
          DO ib1=1,nb
             part(ib1)=ip
             ip=ip+2
             IF (ip >= NumProcs+1) ip=1
             IF (ib1 == nb/2) ip=2
          END DO
                 
       ELSE IF (ivar == 5) THEN   ! 1st block: np, 2nd: np-2, (nb/2): np-1, ...
          ip=NumProcs
          DO ib1=1,nb
             part(ib1)=ip
             ip=ip-2
             IF (ip <= 0) ip=NumProcs
             IF (ib1 >= nb/2) ip=ip-1
             IF (ip<1) ip=1
          END DO
       
       END IF                     ! NumProcs=1 ; ivar

!-----------------------------------------------------------
!---  Distribution of blocks
!
       ALLOCATE(blMPI(nb))
       DO ib1=1,nb
          blMPI(ib1)%proc = part(ib1)-1
       END DO

!---  protocol
       IF (MyId == 0)  THEN
          open (puni,FILE='Part.dat',STATUS='unknown')
          WRITE(puni,443)  NumProcs, nb
          DO ib1=1,nb
             WRITE(puni,444) ib1, blMPI(ib1)%proc
          END DO
          close(puni)
       END IF
!---  Test ---
       DO ip=0,NumProcs-1
          DO ib1=1,nb
            IF (ip==blMPI(ib1)%proc) THEN
              EXIT
            END IF
            IF (ib1==nb) THEN
              IF (MyId==0) THEN
                WRITE(*,*) 'Processor ',ip,'has no block'
              END IF  
!             CALL MPI_Finalize(MPIErr)
!             STOP
            END IF  
          END DO  
        END DO  

!-----------------------------------------------------------
!---  Set local block structure
!
       nbLoc=0
       DO ib1=1,nb
         IF (MyId==blMPI(ib1)%proc) THEN
           nbLoc=nbLoc+1
         END IF
       END DO
       ALLOCATE(LocGlob(nbLoc))
       LocGlob=0
       nbLoc=0

       DO ib1=1,nb
         IF (MyId == blMPI(ib1)%proc) THEN
           nbLoc=nbLoc+1
           blMPI(ib1)%ibLoc=nbLoc
           LocGlob(nbLoc)=ib1
         END IF
       END DO

!-----------------------------------------------------------
!---  FORMAT Statement
!-----------------------------------------------------------

443    FORMAT(' ================  PARTITION - PROTOCOL  =========='// &
&             ' NumProcs  : ',i4  /  ' NumBlocks : ',i4  //           &
&             '      Block     Processor '/                           &
&             ' ------------------------------')
444    FORMAT(8i10)

!===========================================================

CONTAINS

SUBROUTINE graph_metis(maxngbr,VertOrder,Edge,VertWgt,EdgeWgt)

  IMPLICIT NONE

  INTEGER :: maxngbr
  INTEGER :: VertOrder(:), Edge(:), VertWgt(:), EdgeWgt(:)
  INTEGER :: EdgeLoc(maxngbr),EdgeWgtLoc(maxngbr)

  INTEGER :: i,iLoc,ib,in,ix0,ix1,iy0,iy1,iz0,iz1,j,k

!----------------------------------------------------------------

    i=1
    DO ib=1,nb
      IF (Floor(ib)%AnzahlNachbar .GT. maxngbr) THEN
         write(*,*) 'Error GRAPH_METIS: Number of neighbors of block',  &
                    ib,'is',Floor(ib)%AnzahlNachbar,  &
                    '- so it exceeds assumed max. number = ',maxngbr,'!'
         stop 'GRAPH_METIS: Wrong maximum number of neighbors of block!'
      END IF

      VertWgt(ib)   = Floor(ib)%nc
      VertOrder(ib) = i

      EdgeLoc=0
      EdgeWgtLoc=0
      iLoc=1
      DO in=1,Floor(ib)%AnzahlNachbar

         ix0=Floor(ib)%Nachbars(in)%ix0               ! global indices
         iy0=Floor(ib)%Nachbars(in)%iy0
         iz0=Floor(ib)%Nachbars(in)%iz0
         ix1=Floor(ib)%Nachbars(in)%ix1
         iy1=Floor(ib)%Nachbars(in)%iy1
         iz1=Floor(ib)%Nachbars(in)%iz1
         IF (Floor(ib)%RefineX<Floor(ib)%Nachbars(in)%RefineX) THEN
           ix0=ix0*2
           ix1=ix1*2
         END IF
         IF (Floor(ib)%RefineY<Floor(ib)%Nachbars(in)%RefineY) THEN
           iy0=iy0*2
           iy1=iy1*2
         END IF
         IF (Floor(ib)%RefineZ<Floor(ib)%Nachbars(in)%RefineZ) THEN
           iz0=iz0*2
           iz1=iz1*2
         END IF

         IF (Floor(ib)%Nachbars(in)%ntype.eq.'iw'.or. &
             Floor(ib)%Nachbars(in)%ntype.eq.'ie') THEN
           EdgeLoc(iLoc)=Floor(ib)%Nachbars(in)%ib
           EdgeWgtLoc(iLoc)=(iy1-iy0)*(iz1-iz0)
           iLoc=iLoc+1
         ELSE IF (Floor(ib)%Nachbars(in)%ntype.eq.'is'.or. &
                  Floor(ib)%Nachbars(in)%ntype.eq.'in') THEN
           EdgeLoc(iLoc)=Floor(ib)%Nachbars(in)%ib
           EdgeWgtLoc(iLoc)=(ix1-ix0)*(iz1-iz0)
           iLoc=iLoc+1
         ELSE IF (Floor(ib)%Nachbars(in)%ntype.eq.'ib'.or. &
                  Floor(ib)%Nachbars(in)%ntype.eq.'it') THEN
           EdgeLoc(iLoc)=Floor(ib)%Nachbars(in)%ib
           EdgeWgtLoc(iLoc)=(ix1-ix0)*(iy1-iy0)
           iLoc=iLoc+1
         END IF
      END DO                                                       ! in
      iLoc=iLoc-1
      DO j=1,iLoc
        DO k=j+1,iLoc
          IF (EdgeLoc(j)==EdgeLoc(k)) THEN
            EdgeWgtLoc(j)=EdgeWgtLoc(j)+EdgeWgtLoc(k)
            EdgeLoc(k)=-1
          END IF
        END DO 
      END DO  
      DO j=1,iLoc
        IF (EdgeLoc(j)>0) THEN
          Edge(i)=EdgeLoc(j)
          EdgeWgt(i)=EdgeWgtLoc(j)
          i=i+1
        END IF
      END DO 
   END DO                                                          ! ib
   VertOrder(nb+1) = i
   VertOrder=VertOrder-1
   Edge=Edge-1

 END SUBROUTINE graph_metis

!-----------------------------------------------------------
  END SUBROUTINE assignment
SUBROUTINE read_grid2(FileName)


!---------------------------------------------------------------
!---  Determine Multiblock Structure
!---------------------------------------------------------------

    CHARACTER(*) :: FileName
    INTEGER :: i
    CHARACTER(300) :: Line
    CHARACTER(4)   :: indata_type
    REAL(RealKind) :: deg2rad
    LOGICAL :: SoilDepthManually=.FALSE.
    LOGICAL :: SelfMultiBlock=.TRUE.
    LOGICAL :: MultiBlock=.TRUE.
    INTEGER :: xRes,yRes,zRes,BX,BY,BZ,ix,iy,iz
    INTEGER :: xIncr,yIncr,zIncr


!---------------------------------------------------------------
!   Read data from file

  Domain%nrsoillayers=7
  ALLOCATE(Domain%zSDepth(0:Domain%nrsoillayers))
  Domain%zSDepth(0)=0.0d0
  DO i=1,Domain%nrsoillayers
    Domain%zSDepth(i)=0.01d0*3.0d0**(i-1.d0)
  END DO
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'#Gitter')>0) THEN
      READ(InputUnit,*) indata_type
      READ(InputUnit,*) Domain%x0,Domain%x1,Domain%nx
      READ(InputUnit,*) Domain%y0,Domain%y1,Domain%ny
      READ(InputUnit,*) Domain%z0,Domain%z1,Domain%nz
      IF (indata_type=='geo') THEN
        deg2rad=Pi/180.0d0
        Domain%x0=Domain%x0*deg2rad
        Domain%x1=Domain%x1*deg2rad
        Domain%y0=Domain%y0*deg2rad
        Domain%y1=Domain%y1*deg2rad
      END IF
    ELSE IF (INDEX(Line,'#NrSoilLayers')>0) THEN
      DEALLOCATE(Domain%zSDepth)
      READ(InputUnit,*) Domain%nrsoillayers,SoilDepthManually
      ALLOCATE(Domain%zSDepth(0:Domain%nrsoillayers))
      IF (SoilDepthManually) THEN
        READ(InputUnit,*) (Domain%zSDepth(i), i=0,Domain%nrsoillayers)
      ELSE
        Domain%zSDepth(0)=0.0d0
        DO i=1,Domain%nrsoillayers
          Domain%zSDepth(i)=0.01d0*3.0d0**(i-1.d0)
        END DO
      END IF
    ELSE IF (INDEX(Line,'#Multiblock')>0.AND.MultiBlock) THEN
      SelfMultiBlock=.FALSE.
      Domain%ix0=0
      Domain%iy0=0
      Domain%iz0=0
      READ(InputUnit,*) Domain%ix0,Domain%ix1
      READ(InputUnit,*) Domain%iy0,Domain%iy1
      READ(InputUnit,*) Domain%iz0,Domain%iz1
      Domain%igx0=Domain%ix0
      Domain%igx1=Domain%ix1
      Domain%igy0=Domain%iy0
      Domain%igy1=Domain%iy1
      Domain%igz0=Domain%iz0
      Domain%igz1=Domain%iz1

      Domain%nx = Domain%ix1 - Domain%ix0
      Domain%ny = Domain%iy1 - Domain%iy0
      Domain%nz = Domain%iz1 - Domain%iz0
      Domain%nc = Domain%nx * Domain%ny * Domain%nz

      READ(InputUnit,*)
      READ(InputUnit,*) nb
      ALLOCATE(Floor(nb))

!---------------------------------------------------------------
!   Block structure

      DO ib=1,nb
        READ(InputUnit,*)
        READ(InputUnit,*) Floor(ib)%igx0,Floor(ib)%igx1
        READ(InputUnit,*) Floor(ib)%igy0,Floor(ib)%igy1
        READ(InputUnit,*) Floor(ib)%igz0,Floor(ib)%igz1

        READ(InputUnit,*) Floor(ib)%Refine,Floor(ib)%RefineX,Floor(ib)%RefineY,Floor(ib)%RefineZ,Floor(ib)%RefLevel &
                        ,Floor(ib)%JacAdvection,Floor(ib)%JacDiffusion
        IF (MultiEx) THEN
          Floor(ib)%RefLevel=1
        END IF  
        IF (Floor(ib)%Refine>=0) THEN
          Floor(ib)%ix0=Floor(ib)%igx0*2**Floor(ib)%RefineX
          Floor(ib)%ix1=Floor(ib)%igx1*2**Floor(ib)%RefineX
          Floor(ib)%iy0=Floor(ib)%igy0*2**Floor(ib)%RefineY
          Floor(ib)%iy1=Floor(ib)%igy1*2**Floor(ib)%RefineY
          Floor(ib)%iz0=Floor(ib)%igz0*2**Floor(ib)%RefineZ
          Floor(ib)%iz1=Floor(ib)%igz1*2**Floor(ib)%RefineZ
          Floor(ib)%xShift=0
          Floor(ib)%yShift=0
          Floor(ib)%zShift=0
        ELSE
          Floor(ib)%ix0=Floor(ib)%igx0/2**(-Floor(ib)%RefineX)
          Floor(ib)%ix1=Floor(ib)%igx1/2**(-Floor(ib)%RefineX)
          Floor(ib)%iy0=Floor(ib)%igy0/2**(-Floor(ib)%RefineY)
          Floor(ib)%iy1=Floor(ib)%igy1/2**(-Floor(ib)%RefineY)
          Floor(ib)%iz0=Floor(ib)%igz0/2**(-Floor(ib)%RefineZ)
          Floor(ib)%iz1=Floor(ib)%igz1/2**(-Floor(ib)%RefineZ)
          Floor(ib)%xShift=mod(Floor(ib)%igx0,2**(-Floor(ib)%RefineX))
          Floor(ib)%yShift=mod(Floor(ib)%igy0,2**(-Floor(ib)%RefineY)) ! Hinneburg (vorher RefineX)
          Floor(ib)%zShift=mod(Floor(ib)%igz0,2**(-Floor(ib)%RefineZ)) ! Hinneburg (vorher RefineX)
        END IF
        Floor(ib)%ib=ib

        Floor(ib)%nx = Floor(ib)%ix1 - Floor(ib)%ix0
        Floor(ib)%ny = Floor(ib)%iy1 - Floor(ib)%iy0
        Floor(ib)%nz = Floor(ib)%iz1 - Floor(ib)%iz0
        Floor(ib)%nc = Floor(ib)%nx * Floor(ib)%ny * Floor(ib)%nz

        Floor(ib)%TypeW='iw'
        IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West/='Period') THEN
          Floor(ib)%TypeW='ow'
        END IF 
        Floor(ib)%TypeE='ie'
        IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East/='Period') THEN
          Floor(ib)%TypeE='oe'
        END IF 
        Floor(ib)%TypeS='is'
        IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South/='Period') THEN
          Floor(ib)%TypeS='os'
        END IF 
        Floor(ib)%TypeN='in'
        IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North/='Period') THEN
          Floor(ib)%TypeN='on'
        END IF 
        Floor(ib)%TypeB='ib'
        IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom/='Period') THEN
          Floor(ib)%TypeB='ob'
        END IF 
        Floor(ib)%TypeT='it'
        IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top/='Period') THEN
          Floor(ib)%TypeT='ot'
        END IF 
      END DO  ! ib
    ELSE IF (INDEX(Line,'#SelfMultiblock1')>0.AND.SelfMultiBlock) THEN
      MultiBlock=.FALSE.
      Domain%ix0=0
      Domain%iy0=0
      Domain%iz0=0
      READ(InputUnit,*) Domain%ix0,Domain%ix1,BX
      READ(InputUnit,*) Domain%iy0,Domain%iy1,BY
      READ(InputUnit,*) Domain%iz0,Domain%iz1,BZ
      Domain%igx0=Domain%ix0
      Domain%igx1=Domain%ix1
      Domain%igy0=Domain%iy0
      Domain%igy1=Domain%iy1
      Domain%igz0=Domain%iz0
      Domain%igz1=Domain%iz1

      Domain%nx = Domain%ix1 - Domain%ix0
      Domain%ny = Domain%iy1 - Domain%iy0
      Domain%nz = Domain%iz1 - Domain%iz0
      Domain%nc = Domain%nx * Domain%ny * Domain%nz
      nb=BX*BY*BZ
      xRes=MOD(Domain%nx,BX)
      xIncr=Domain%nx/BX
      yRes=MOD(Domain%ny,BY)
      yIncr=Domain%ny/BY
      zRes=MOD(Domain%nz,BZ)
      zIncr=Domain%nz/BZ
      ALLOCATE(Floor(nb))
      iz0=0
      ib=1
      DO iz=1,BZ
        IF (iz<=zRes) THEN
          iz1=iz0+zIncr+1
        ELSE  
          iz1=iz0+zIncr
        END IF  
        iy0=0
        DO iy=1,BY
          IF (iy<=yRes) THEN
            iy1=iy0+yIncr+1
          ELSE  
            iy1=iy0+yIncr
          END IF  
          ix0=0 
          DO ix=1,BX
            IF (ix<=xRes) THEN
              ix1=ix0+xIncr+1
            ELSE  
              ix1=ix0+xIncr
            END IF  
            Floor(ib)%igx0=ix0
            Floor(ib)%igy0=iy0
            Floor(ib)%igz0=iz0
            Floor(ib)%igx1=ix1
            Floor(ib)%igy1=iy1
            Floor(ib)%igz1=iz1

            Floor(ib)%Refine=0
            Floor(ib)%RefineX=0
            Floor(ib)%RefineY=0
            Floor(ib)%RefineZ=0
            Floor(ib)%RefLevel=4
            Floor(ib)%JacAdvection=.TRUE.
            Floor(ib)%JacDiffusion=.TRUE.

            IF (Floor(ib)%Refine>=0) THEN
              Floor(ib)%ix0=Floor(ib)%igx0*2**Floor(ib)%RefineX
              Floor(ib)%ix1=Floor(ib)%igx1*2**Floor(ib)%RefineX
              Floor(ib)%iy0=Floor(ib)%igy0*2**Floor(ib)%RefineY
              Floor(ib)%iy1=Floor(ib)%igy1*2**Floor(ib)%RefineY
              Floor(ib)%iz0=Floor(ib)%igz0*2**Floor(ib)%RefineZ
              Floor(ib)%iz1=Floor(ib)%igz1*2**Floor(ib)%RefineZ
              Floor(ib)%xShift=0
              Floor(ib)%yShift=0
              Floor(ib)%zShift=0
            ELSE
              Floor(ib)%ix0=Floor(ib)%igx0/2**(-Floor(ib)%RefineX)
              Floor(ib)%ix1=Floor(ib)%igx1/2**(-Floor(ib)%RefineX)
              Floor(ib)%iy0=Floor(ib)%igy0/2**(-Floor(ib)%RefineY)
              Floor(ib)%iy1=Floor(ib)%igy1/2**(-Floor(ib)%RefineY)
              Floor(ib)%iz0=Floor(ib)%igz0/2**(-Floor(ib)%RefineZ)
              Floor(ib)%iz1=Floor(ib)%igz1/2**(-Floor(ib)%RefineZ)
              Floor(ib)%xShift=mod(Floor(ib)%igx0,2**(-Floor(ib)%RefineX))
              Floor(ib)%yShift=mod(Floor(ib)%igy0,2**(-Floor(ib)%RefineY)) ! Hinneburg (vorher RefineX)
              Floor(ib)%zShift=mod(Floor(ib)%igz0,2**(-Floor(ib)%RefineZ)) ! Hinneburg (vorher RefineX)
            END IF
            Floor(ib)%ib=ib

            Floor(ib)%nx = Floor(ib)%ix1 - Floor(ib)%ix0
            Floor(ib)%ny = Floor(ib)%iy1 - Floor(ib)%iy0
            Floor(ib)%nz = Floor(ib)%iz1 - Floor(ib)%iz0
            Floor(ib)%nc = Floor(ib)%nx * Floor(ib)%ny * Floor(ib)%nz

            Floor(ib)%TypeW='iw'
            IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West/='Period') THEN
              Floor(ib)%TypeW='ow'
            END IF 
            Floor(ib)%TypeE='ie'
            IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East/='Period') THEN
              Floor(ib)%TypeE='oe'
            END IF 
            Floor(ib)%TypeS='is'
            IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South/='Period') THEN
              Floor(ib)%TypeS='os'
            END IF 
            Floor(ib)%TypeN='in'
            IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North/='Period') THEN
              Floor(ib)%TypeN='on'
            END IF 
            Floor(ib)%TypeB='ib'
            IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom/='Period') THEN
              Floor(ib)%TypeB='ob'
            END IF 
            Floor(ib)%TypeT='it'
            IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top/='Period') THEN
              Floor(ib)%TypeT='ot'
            END IF 
            ix0=ix1
            ib=ib+1
          END DO  
          iy0=iy1
        END DO  
        iz0=iz1
      END DO  
      EXIT
    ELSE IF (INDEX(Line,'#SelfMultiblock')>0.AND.SelfMultiBlock) THEN
      CALL SelfMultiBlockCompute
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)

 END SUBROUTINE read_grid2

SUBROUTINE inp_part(FileName)
!

!---------------------------------------------------------------
!---  Determine Multiblock Structure
!---------------------------------------------------------------

   CHARACTER(*) :: FileName

   INTEGER :: ibx, iby, ibz

!---------------------------------------------------------------
!  Read data from file

   CALL InputModelControl(FileName)
   CALL InputModelPos(FileName)
   CALL InputEnvPos(FileName)
   CALL InputModelPhysics(FileName)
!  CALL InputModelOutput(FileName)
!  NameList Input
   CALL InputModelBCVel(FileName)

!
!  Decision on read routine

   CALL read_grid2(FileName)

    
!---------------------------------------------------------------
!   Refinement of blocks

    IF (ref_glob /= 0) THEN
       CALL refine_grid
    END IF
    
!---------------------------------------------------------------
!   Determine neighbors

    CALL get_ngbrs
    CALL ngbr_bounds
    nbGlob=nb
! -- Anfangs-Zuordnung der Bloecke zu den Prozessoren --
    CALL assignment
! -- Vorbereitung des Rand-Austauschs --
    CALL getbufsize

!--------------------------------------------------------
  CALL SetLocalBlocks

END SUBROUTINE inp_part
SUBROUTINE get_ngbrs

!---------------------------------------------------------------
!---  Determine Neighbors w.r.t. Multiblock Structure
!---------------------------------------------------------------

  INTEGER :: ib1,in
  TYPE(NAchbar_T), ALLOCATABLE :: curr_nachb(:)
! TYPE(NAchbar_T), POINTER :: Nachbar

  MaxBlkNgbr = 0
  ALLOCATE(curr_nachb(nb+5))
  DO ib=1,nb

    ix0=Floor(ib)%ix0
    ix1=Floor(ib)%ix1
    iy0=Floor(ib)%iy0
    iy1=Floor(ib)%iy1
    iz0=Floor(ib)%iz0
    iz1=Floor(ib)%iz1

    AnzahlNachbar=0

!--------------------------------------------------------
!   West Boundary

    IF (Floor(ib)%igx0==domain%ix0) THEN
      IF (BCVel%West/='Period') THEN
        AnzahlNachbar=AnzahlNachbar+1
        curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
        curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
        curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
        curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
        curr_nachb(AnzahlNachbar)%ib=ib
        curr_nachb(AnzahlNachbar)%nTYPE='ow'

      ELSE
        DO ib1=1,nb
          IF (domain%ix1==Floor(ib1)%igx1) THEN
            IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
                Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                Floor(ib)%igz1>Floor(ib1)%igz0) THEN
              AnzahlNachbar=AnzahlNachbar+1
              curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
              curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
              curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
              curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
              curr_nachb(AnzahlNachbar)%ib=ib1
              curr_nachb(AnzahlNachbar)%nTYPE='iw'
            END IF
          END IF
        END DO
      END IF
    ELSE
      DO ib1=1,nb
        IF (Floor(ib)%igx0==Floor(ib1)%igx1) THEN
          IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
            Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
            Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
            Floor(ib)%igz1>Floor(ib1)%igz0) THEN
            AnzahlNachbar=AnzahlNachbar+1
            curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
            curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
            curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
            curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
            curr_nachb(AnzahlNachbar)%ib=ib1
            curr_nachb(AnzahlNachbar)%nTYPE='iw'
          END IF
        END IF
      END DO
    END IF         ! Floor(ib)%igx0==domain%ix0

!--------------------------------------------------------
!   East Boundary

    IF (Floor(ib)%igx1==domain%ix1) THEN
      IF (BCVel%East/='Period') THEN
        AnzahlNachbar=AnzahlNachbar+1
        curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
        curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
        curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
        curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
        curr_nachb(AnzahlNachbar)%ib=ib
        curr_nachb(AnzahlNachbar)%nTYPE='oe'

      ELSE
        DO ib1=1,nb
          IF (domain%ix0==Floor(ib1)%igx0) THEN
            IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
                Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                Floor(ib)%igz1>Floor(ib1)%igz0) THEN
               AnzahlNachbar=AnzahlNachbar+1
               curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
               curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
               curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
               curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
               curr_nachb(AnzahlNachbar)%ib=ib1
               curr_nachb(AnzahlNachbar)%nTYPE='ie'
             END IF
           END IF
         END DO
       END IF
     ELSE
       DO ib1=1,nb
         IF (Floor(ib)%igx1==Floor(ib1)%igx0) THEN
           IF (Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
               Floor(ib)%igy1>Floor(ib1)%igy0 .and. &
               Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
               Floor(ib)%igz1>Floor(ib1)%igz0) THEN
              AnzahlNachbar=AnzahlNachbar+1
              curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
              curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
              curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
              curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
              curr_nachb(AnzahlNachbar)%ib=ib1
              curr_nachb(AnzahlNachbar)%nTYPE='ie'
           END IF
         END IF
       END DO
     END IF         ! Floor(ib)%igx1==domain%ix1

!--------------------------------------------------------
!    South Boundary

     IF (Floor(ib)%igy0==domain%iy0) THEN
       IF (BCVel%South/='Period') THEN
         AnzahlNachbar=AnzahlNachbar+1
         curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
         curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
         curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
         curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
         curr_nachb(AnzahlNachbar)%ib=ib
         curr_nachb(AnzahlNachbar)%nTYPE='os'
       ELSE
         DO ib1=1,nb
           IF (domain%iy1==Floor(ib1)%igy1) THEN
             IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                 Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                 Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                 Floor(ib)%igz1>Floor(ib1)%igz0) THEN
               AnzahlNachbar=AnzahlNachbar+1
               curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
               curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
               curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
               curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
               curr_nachb(AnzahlNachbar)%ib=ib1
               curr_nachb(AnzahlNachbar)%nTYPE='is'
             END IF
           END IF
         END DO
       END IF
     ELSE
       DO ib1=1,nb
         IF (Floor(ib)%igy0==Floor(ib1)%igy1) THEN
           IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
               Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
               Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
               Floor(ib)%igz1>Floor(ib1)%igz0) THEN
             AnzahlNachbar=AnzahlNachbar+1
             curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
             curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
             curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
             curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
             curr_nachb(AnzahlNachbar)%ib=ib1
             curr_nachb(AnzahlNachbar)%nTYPE='is'
           END IF
         END IF
       END DO
     END IF         ! Floor(ib)%igy0==domain%iy0

!--------------------------------------------------------
!    North Boundary

     IF (Floor(ib)%igy1==domain%iy1) THEN
       IF (BCVel%North/='Period') THEN
         AnzahlNachbar=AnzahlNachbar+1
         curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
         curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
         curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
         curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
         curr_nachb(AnzahlNachbar)%ib=ib
         curr_nachb(AnzahlNachbar)%nTYPE='on'
       ELSE
         DO ib1=1,nb
           IF (domain%iy0==Floor(ib1)%igy0) THEN
             IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                 Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                 Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
                 Floor(ib)%igz1>Floor(ib1)%igz0) THEN
               AnzahlNachbar=AnzahlNachbar+1
               curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
               curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
               curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
               curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
               curr_nachb(AnzahlNachbar)%ib=ib1
               curr_nachb(AnzahlNachbar)%nTYPE='in'
             END IF
           END IF
         END DO
       END IF
     ELSE
       DO ib1=1,nb
         IF (Floor(ib)%igy1==Floor(ib1)%igy0) THEN
           IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
               Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
               Floor(ib)%igz0<Floor(ib1)%igz1 .and. &
               Floor(ib)%igz1>Floor(ib1)%igz0) THEN
             AnzahlNachbar=AnzahlNachbar+1
             curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
             curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
             curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
             curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
             curr_nachb(AnzahlNachbar)%ib=ib1
             curr_nachb(AnzahlNachbar)%nTYPE='in'
           END IF
         END IF
       END DO
     END IF         ! Floor(ib)%igy1==domain%iy1

!--------------------------------------------------------
!      Bottom Boundary

       IF (Floor(ib)%igz0==domain%iz0) THEN
         IF (BCVel%Bottom/='Period') THEN
           AnzahlNachbar=AnzahlNachbar+1
           curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
           curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
           curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
           curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
           curr_nachb(AnzahlNachbar)%ib=ib
           curr_nachb(AnzahlNachbar)%nTYPE='ob'
         ELSE
           DO ib1=1,nb
             IF (domain%iz1==Floor(ib1)%igz1) THEN
               IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                   Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                   Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                   Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                 AnzahlNachbar=AnzahlNachbar+1
                 curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                 curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                 curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                 curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                 curr_nachb(AnzahlNachbar)%ib=ib1
                 curr_nachb(AnzahlNachbar)%nTYPE='ib'
               END IF
             END IF
           END DO
         END IF
       ELSE
          DO ib1=1,nb
             IF (Floor(ib)%igz0==Floor(ib1)%igz1) THEN
                IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                    Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                    Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                    Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                   AnzahlNachbar=AnzahlNachbar+1
                   curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                   curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                   curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                   curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                   curr_nachb(AnzahlNachbar)%ib=ib1
                   curr_nachb(AnzahlNachbar)%nTYPE='ib'
                END IF
             END IF
           END DO
       END IF         ! Floor(ib)%igz0==domain%iz0

!--------------------------------------------------------
!      Top Boundary

       IF (Floor(ib)%igz1==domain%iz1) THEN
         IF (BCVel%Top/='Period') THEN
           AnzahlNachbar=AnzahlNachbar+1
           curr_nachb(AnzahlNachbar)%Refine=Floor(ib)%Refine
           curr_nachb(AnzahlNachbar)%RefineX=Floor(ib)%RefineX
           curr_nachb(AnzahlNachbar)%RefineY=Floor(ib)%RefineY
           curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib)%RefineZ
           curr_nachb(AnzahlNachbar)%ib=ib
           curr_nachb(AnzahlNachbar)%nTYPE='ot'
         ELSE
           DO ib1=1,nb
             IF (domain%iz0==Floor(ib1)%igz0) THEN
               IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                   Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                   Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                   Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                 AnzahlNachbar=AnzahlNachbar+1
                 curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                 curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                 curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                 curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                 curr_nachb(AnzahlNachbar)%ib=ib1
                 curr_nachb(AnzahlNachbar)%nTYPE='it'
               END IF
             END IF
           END DO
         END IF
       ELSE
          DO ib1=1,nb
             IF (Floor(ib)%igz1==Floor(ib1)%igz0) THEN
                IF (Floor(ib)%igx0<Floor(ib1)%igx1 .and. &
                    Floor(ib)%igx1>Floor(ib1)%igx0 .and. &
                    Floor(ib)%igy0<Floor(ib1)%igy1 .and. &
                    Floor(ib)%igy1>Floor(ib1)%igy0) THEN
                   AnzahlNachbar=AnzahlNachbar+1
                   curr_nachb(AnzahlNachbar)%Refine=Floor(ib1)%Refine
                   curr_nachb(AnzahlNachbar)%RefineX=Floor(ib1)%RefineX
                   curr_nachb(AnzahlNachbar)%RefineY=Floor(ib1)%RefineY
                   curr_nachb(AnzahlNachbar)%RefineZ=Floor(ib1)%RefineZ
                   curr_nachb(AnzahlNachbar)%ib=ib1
                   curr_nachb(AnzahlNachbar)%nTYPE='it'
                END IF
             END IF
           END DO
       END IF         ! Floor(ib)%igz1==domain%iz1

       MaxBlkNgbr = MAX(MaxBlkNgbr,AnzahlNachbar)

       Floor(ib)%AnzahlNachbar=AnzahlNachbar
       ALLOCATE(Floor(ib)%Nachbars(AnzahlNachbar))
       Floor(ib)%Nachbars(1:AnzahlNachbar)= curr_nachb(1:AnzahlNachbar)

  END DO         ! ib

  DEALLOCATE(curr_nachb)

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      Nachbar%IncrX=ABS(-RefineX+RefineNachbarX)+1
      Nachbar%IncrY=ABS(-RefineY+RefineNachbarY)+1
      Nachbar%IncrZ=ABS(-RefineZ+RefineNachbarZ)+1
      IF (nType=='iw'.OR.nType=='ie') THEN
        Nachbar%CopyCase=2*(RefineY-RefineNachbarY+1)+ &
                           (RefineZ-RefineNachbarZ+1)
      ELSE IF (nType=='is'.OR.nType=='in') THEN
        Nachbar%CopyCase=2*(RefineX-RefineNachbarX+1)+ &
                           (RefineZ-RefineNachbarZ+1)
      ELSE IF (nType=='ib'.OR.nType=='it') THEN
        Nachbar%CopyCase=2*(RefineX-RefineNachbarX+1)+ &
                           (RefineY-RefineNachbarY+1)
      END IF
    END DO
  END DO

END SUBROUTINE get_ngbrs

SUBROUTINE SetLocalBlocks

  INTEGER :: i,in,ip,jbLoc

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      Nachbar%ibLoc=0 
      DO jbLoc=1,nbLoc
        IF (ibn==LocGlob(jbLoc)) THEN
          Nachbar%ibLoc=jbLoc
          EXIT
        END IF
      END DO
    END DO
  END DO
  DO ip=0,NumProcs-1
    IF (chain(ip)%nlink > 0) THEN
      DO i=1,SIZE(chain(ip)%glied)
        CALL Set(chain(ip)%glied(i))
        chain(ip)%glied(i)%OwnBlockLoc=0
        DO jbLoc=1,nbLoc
          IF (ibC==LocGlob(jbLoc)) THEN
            chain(ip)%glied(i)%OwnBlockLoc=jbLoc
            EXIT
          END IF
        END DO
      END DO
    END IF  
  END DO

END SUBROUTINE SetLocalBlocks

SUBROUTINE SelfMultiblockCompute

  INTEGER :: NX,NY,NZ,BX,BY,BZ,xi,yi,zi
  INTEGER :: xRes,yRes,zRes
  INTEGER :: xIncr,yIncr,zIncr
  INTEGER :: ix,ix0,ix1
  INTEGER :: iy,iy0,iy1
  INTEGER :: iz,iz0,iz1
  CHARACTER :: lc

  TYPE Block_T
    INTEGER :: ix0,ix1,xC=1
    INTEGER :: iy0,iy1,yC=1
    INTEGER :: iz0,iz1,zC=1
    INTEGER :: Coarse=1
    LOGICAL :: Active=.TRUE.
  END TYPE Block_T

  TYPE(Block_T), ALLOCATABLE :: Blocks(:,:,:)
  INTEGER :: xCMax,yCMax,zCMax
  INTEGER :: xBMax,yBMax,zBMax

  INTEGER :: i,j,k
  INTEGER :: NumFineBlocks
  REAL(8) :: xC,yC,zC
  INTEGER :: ib,iF,iC
  INTEGER :: ii,jj,kk
  INTEGER :: iMaxC
  INTEGER :: NumberCells
  INTEGER :: iCoarse,jCoarse,kCoarse
  INTEGER :: iCoarseAct,jCoarseAct,kCoarseAct
  INTEGER, ALLOCATABLE :: xCLoc(:,:,:)
  INTEGER, ALLOCATABLE :: yCLoc(:,:,:)
  INTEGER, ALLOCATABLE :: zCLoc(:,:,:)

  READ(InputUnit,*) Domain%ix0,Domain%ix1,BX
  READ(InputUnit,*) Domain%iy0,Domain%iy1,BY
  READ(InputUnit,*) Domain%iz0,Domain%iz1,BZ
  ix0=Domain%ix0
  ix1=Domain%ix1
  iy0=Domain%iy0
  iy1=Domain%iy1
  iz0=Domain%iz0
  iz1=Domain%iz1
  Domain%igx0=Domain%ix0
  Domain%igx1=Domain%ix1
  Domain%igy0=Domain%iy0
  Domain%igy1=Domain%iy1
  Domain%igz0=Domain%iz0
  Domain%igz1=Domain%iz1
  NX=ix1-ix0
  NY=iy1-iy0
  NZ=iz1-iz0
  ALLOCATE(Blocks(BX,BY,BZ))
  READ(InputUnit,*) NumFineBlocks
  READ(InputUnit,*) xCMax,yCMax,zCMax
  READ(InputUnit,*) xBMax,yBMax,zBMax
  DO iF=1,NumFineBlocks
    READ(InputUnit,*) i,j,k,xC,yC,zC
    Blocks(i,j,k)%xC=xC
    Blocks(i,j,k)%yC=yC
    Blocks(i,j,k)%zC=zC
  END DO

  xRes=MOD(NX,BX)
  xIncr=NX/BX
  yRes=MOD(NY,BY)
  yIncr=NY/BY
  zRes=MOD(NZ,BZ)
  zIncr=NZ/BZ

  iz0=0
  DO iz=1,BZ
    IF (iz<=zRes) THEN
      iz1=iz0+zIncr+1
    ELSE  
      iz1=iz0+zIncr
    END IF  
    iy0=0
    DO iy=1,BY
      IF (iy<=yRes) THEN
        iy1=iy0+yIncr+1
      ELSE  
        iy1=iy0+yIncr
      END IF  
      ix0=0 
      DO ix=1,BX
        IF (ix<=xRes) THEN
          ix1=ix0+xIncr+1
        ELSE  
          ix1=ix0+xIncr
        END IF  
        Blocks(ix,iy,iz)%ix0=ix0
        Blocks(ix,iy,iz)%ix1=ix1
        Blocks(ix,iy,iz)%iy0=iy0
        Blocks(ix,iy,iz)%iy1=iy1
        Blocks(ix,iy,iz)%iz0=iz0
        Blocks(ix,iy,iz)%iz1=iz1
        ix0=ix1
      END DO  
      iy0=iy1
    END DO  
    iz0=iz1
  END DO  


  ALLOCATE(xCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                 LBOUND(Blocks,2):UBOUND(Blocks,2), &
                 LBOUND(Blocks,3):UBOUND(Blocks,3)))  
  ALLOCATE(yCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                 LBOUND(Blocks,2):UBOUND(Blocks,2), &
                 LBOUND(Blocks,3):UBOUND(Blocks,3)))  
  ALLOCATE(zCLoc(LBOUND(Blocks,1):UBOUND(Blocks,1), &
                 LBOUND(Blocks,2):UBOUND(Blocks,2), &
                 LBOUND(Blocks,3):UBOUND(Blocks,3)))  

  DO iC=1,xCMax
    xCLoc=1
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (Blocks(i,j,k)%xC==1) THEN
            IF (i>LBOUND(Blocks,1)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i-1,j,k)%xC)
            END IF  
            IF (j>LBOUND(Blocks,2)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j-1,k)%xC)
            END IF  
            IF (k>LBOUND(Blocks,3)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j,k-1)%xC)
            END IF  
            IF (i<UBOUND(Blocks,1)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i+1,j,k)%xC)
            END IF  
            IF (j<UBOUND(Blocks,2)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j+1,k)%xC)
            END IF  
            IF (k<UBOUND(Blocks,3)) THEN
              xCLoc(i,j,k)=MIN(xCLoc(i,j,k),Blocks(i,j,k+1)%xC)
            END IF  
          END IF  
        END DO
      END DO
    END DO
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (xCLoc(i,j,k)<1) THEN
            Blocks(i,j,k)%xC=xCLoc(i,j,k)-1
            Blocks(i,j,k)%Coarse=MAX(-xCLoc(i,j,k),Blocks(i,j,k)%Coarse)
          END IF  
        END DO
      END DO
    END DO
  END DO
  DO iC=1,yCMax
    yCLoc=1
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (Blocks(i,j,k)%yC==1) THEN
            IF (i>LBOUND(Blocks,1)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i-1,j,k)%yC)
            END IF  
            IF (j>LBOUND(Blocks,2)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j-1,k)%yC)
            END IF  
            IF (k>LBOUND(Blocks,3)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j,k-1)%yC)
            END IF  
            IF (i<UBOUND(Blocks,1)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i+1,j,k)%yC)
            END IF  
            IF (j<UBOUND(Blocks,2)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j+1,k)%yC)
            END IF  
            IF (k<UBOUND(Blocks,3)) THEN
              yCLoc(i,j,k)=MIN(yCLoc(i,j,k),Blocks(i,j,k+1)%yC)
            END IF  
          END IF  
        END DO
      END DO
    END DO
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (yCLoc(i,j,k)<1) THEN
            Blocks(i,j,k)%yC=yCLoc(i,j,k)-1
            Blocks(i,j,k)%Coarse=MAX(-yCLoc(i,j,k),Blocks(i,j,k)%Coarse)
          END IF  
        END DO
      END DO
    END DO
  END DO
  DO iC=1,zCMax
    zCLoc=1
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (Blocks(i,j,k)%zC==1) THEN
            IF (i>LBOUND(Blocks,1)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i-1,j,k)%zC)
            END IF  
            IF (j>LBOUND(Blocks,2)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j-1,k)%zC)
            END IF  
            IF (k>LBOUND(Blocks,3)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j,k-1)%zC)
            END IF  
            IF (i<UBOUND(Blocks,1)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i+1,j,k)%zC)
            END IF  
            IF (j<UBOUND(Blocks,2)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j+1,k)%zC)
            END IF  
            IF (k<UBOUND(Blocks,3)) THEN
              zCLoc(i,j,k)=MIN(zCLoc(i,j,k),Blocks(i,j,k+1)%zC)
            END IF  
          END IF  
        END DO
      END DO
    END DO
    DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
      DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
        DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
          IF (zCLoc(i,j,k)<1) THEN
            Blocks(i,j,k)%zC=zCLoc(i,j,k)-1
          END IF  
        END DO
      END DO
    END DO
  END DO
  NumberCells=0
  nb=0
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
        IF (Blocks(i,j,k)%xC==1) Blocks(i,j,k)%xC=-xCMax
        IF (Blocks(i,j,k)%yC==1) Blocks(i,j,k)%yC=-yCMax
        IF (Blocks(i,j,k)%zC==1) Blocks(i,j,k)%zC=-zCMax
        Blocks(i,j,k)%Coarse=MAX(-Blocks(i,j,k)%xC,-Blocks(i,j,k)%yC,-Blocks(i,j,k)%zC)
        nb=nb+2**MAX(xBMax+Blocks(i,j,k)%xC,0) &
             *2**MAX(yBMax+Blocks(i,j,k)%yC,0) &
             *2**MAX(zBMax+Blocks(i,j,k)%zC,0) 
!       NumberCells=NumberCells+(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/(2**(-Blocks(i,j,k)%xC))* &
!                               (Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/(2**(-Blocks(i,j,k)%yC))* &
!                               (Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/(2**(-Blocks(i,j,k)%zC))
      END DO
    END DO
  END DO
! Coarsening 
  
  iCoarse=2 
  jCoarse=1 
  kCoarse=2 
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1),iCoarse
    IF (i+iCoarse-1>UBOUND(Blocks,1)) THEN
      iCoarseAct=UBOUND(Blocks,1)-i+1
    ELSE
      iCoarseAct=iCoarse
    END IF  
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2),jCoarse
      IF (j+jCoarse-1>UBOUND(Blocks,2)) THEN
        jCoarseAct=UBOUND(Blocks,2)-j+1
      ELSE
        jCoarseAct=jCoarse
      END IF  
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3),kCoarse
        IF (k+kCoarse-1>UBOUND(Blocks,3)) THEN
          kCoarseAct=UBOUND(Blocks,3)-k+1
        ELSE
          kCoarseAct=kCoarse
        END IF  
        IF (MAXVAL(Blocks(i:i+iCoarseAct-1,j:j+jCoarseAct-1,k:k+kCoarseAct-1)%Coarse)>0.AND. &
          MAXVAL(Blocks(i:i+iCoarseAct-1,j:j+jCoarseAct-1,k:k+kCoarseAct-1)%Coarse)==          &
          MINVAL(Blocks(i:i+iCoarseAct-1,j:j+jCoarseAct-1,k:k+kCoarseAct-1)%Coarse)) THEN
          DO ii=1,iCoarseAct
            DO jj=1,jCoarseAct
              DO kk=1,kCoarseAct
                Blocks(i+ii-1,j+jj-1,k+kk-1)%Active=.FALSE.
              END DO  
            END DO  
          END DO  
          Blocks(i,j,k)%Active=.TRUE.
          nb=nb+1-iCoarseAct*jCoarseAct*kCoarseAct
          Blocks(i,j,k)%ix1=Blocks(i+iCoarseAct-1,j+jCoarseAct-1,k+kCoarseAct-1)%ix1
          Blocks(i,j,k)%iy1=Blocks(i+iCoarseAct-1,j+jCoarseAct-1,k+kCoarseAct-1)%iy1
          Blocks(i,j,k)%iz1=Blocks(i+iCoarseAct-1,j+jCoarseAct-1,k+kCoarseAct-1)%iz1
        END IF
      END DO
    END DO
  END DO

  ALLOCATE(Floor(nb))
  ib=0
  DO i=LBOUND(Blocks,1),UBOUND(Blocks,1)
    DO j=LBOUND(Blocks,2),UBOUND(Blocks,2)
      DO k=LBOUND(Blocks,3),UBOUND(Blocks,3)
        IF (Blocks(i,j,k)%Active) THEN
          DO ii=1,2**MAX(xBMax+Blocks(i,j,k)%xC,0)
            DO jj=1,2**MAX(yBMax+Blocks(i,j,k)%yC,0)
              DO kk=1,2**MAX(xBMax+Blocks(i,j,k)%zC,0)
                ib=ib+1
                Floor(ib)%igx0=Blocks(i,j,k)%ix0+(ii-1)*(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/2**MAX(xBMax+Blocks(i,j,k)%xC,0) 
                Floor(ib)%igx1=Blocks(i,j,k)%ix0+(ii  )*(Blocks(i,j,k)%ix1-Blocks(i,j,k)%ix0)/2**MAX(xBMax+Blocks(i,j,k)%xC,0)  
                Floor(ib)%igy0=Blocks(i,j,k)%iy0+(jj-1)*(Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/2**MAX(yBMax+Blocks(i,j,k)%yC,0)
                Floor(ib)%igy1=Blocks(i,j,k)%iy0+(jj  )*(Blocks(i,j,k)%iy1-Blocks(i,j,k)%iy0)/2**MAX(yBMax+Blocks(i,j,k)%yC,0)  
                Floor(ib)%igz0=Blocks(i,j,k)%iz0+(kk-1)*(Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/2**MAX(zBMax+Blocks(i,j,k)%zC,0)
                Floor(ib)%igz1=Blocks(i,j,k)%iz0+(kk  )*(Blocks(i,j,k)%iz1-Blocks(i,j,k)%iz0)/2**MAX(zBMax+Blocks(i,j,k)%zC,0)  
                iMaxC=MIN(Blocks(i,j,k)%xC,Blocks(i,j,k)%yC,Blocks(i,j,k)%zC)
                Floor(ib)%Refine=iMaxC
                Floor(ib)%RefineX=Blocks(i,j,k)%xC
                Floor(ib)%RefineY=Blocks(i,j,k)%yC
                Floor(ib)%RefineZ=Blocks(i,j,k)%zC
                Floor(ib)%ix0=Floor(ib)%igx0/2**(-Floor(ib)%RefineX)
                Floor(ib)%ix1=Floor(ib)%igx1/2**(-Floor(ib)%RefineX)
                Floor(ib)%iy0=Floor(ib)%igy0/2**(-Floor(ib)%RefineY)
                Floor(ib)%iy1=Floor(ib)%igy1/2**(-Floor(ib)%RefineY)
                Floor(ib)%iz0=Floor(ib)%igz0/2**(-Floor(ib)%RefineZ)
                Floor(ib)%iz1=Floor(ib)%igz1/2**(-Floor(ib)%RefineZ)
                Floor(ib)%xShift=mod(Floor(ib)%igx0,2**(-Floor(ib)%RefineX))
                Floor(ib)%yShift=mod(Floor(ib)%igy0,2**(-Floor(ib)%RefineY)) 
                Floor(ib)%zShift=mod(Floor(ib)%igz0,2**(-Floor(ib)%RefineZ)) 
                Floor(ib)%ib=ib
                Floor(ib)%nx = Floor(ib)%ix1 - Floor(ib)%ix0
                Floor(ib)%ny = Floor(ib)%iy1 - Floor(ib)%iy0
                Floor(ib)%nz = Floor(ib)%iz1 - Floor(ib)%iz0
                Floor(ib)%RefLevel=2

                IF (Floor(ib)%Refine<0) THEN
                  Floor(ib)%JacAdvection=.FALSE.
                  Floor(ib)%JacDiffusion=.FALSE.
                ELSE
                  Floor(ib)%JacAdvection=.TRUE.
                  Floor(ib)%JacDiffusion=.TRUE.
                END IF

                Floor(ib)%nc = Floor(ib)%nx * Floor(ib)%ny * Floor(ib)%nz

                Floor(ib)%TypeW='iw'
                IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West/='Period') THEN
                  Floor(ib)%TypeW='ow'
                ELSE IF (Floor(ib)%igx0==domain%ix0.AND.BCVel%West=='Period') THEN
                   Floor(ib)%TypeW='pw'
                END IF
                Floor(ib)%TypeE='ie'
                IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East/='Period') THEN
                  Floor(ib)%TypeE='oe'
                ELSE IF (Floor(ib)%igx1==domain%ix1.AND.BCVel%East=='Period') THEN
                 Floor(ib)%TypeE='pe'
                END IF
                Floor(ib)%TypeS='is'
                IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South/='Period') THEN
                  Floor(ib)%TypeS='os'
                ELSE IF (Floor(ib)%igy0==domain%iy0.AND.BCVel%South=='Period') THEN
                   Floor(ib)%TypeS='ps'
                END IF
                Floor(ib)%TypeN='in'
                IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North/='Period') THEN
                  Floor(ib)%TypeN='on'
                ELSE IF (Floor(ib)%igy1==domain%iy1.AND.BCVel%North=='Period') THEN
                 Floor(ib)%TypeN='pn'
                END IF
                Floor(ib)%TypeB='ib'
                IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom/='Period') THEN
                  Floor(ib)%TypeB='ob'
                ELSE IF (Floor(ib)%igz0==domain%iz0.AND.BCVel%Bottom=='Period') THEN
                  Floor(ib)%TypeB='pb'
                END IF
                Floor(ib)%TypeT='it'
                IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top/='Period') THEN
                  Floor(ib)%TypeT='ot'
                ELSE IF (Floor(ib)%igz1==domain%iz1.AND.BCVel%Top=='Period') THEN
                  Floor(ib)%TypeT='pt'
                END IF
              END DO
            END DO
          END DO
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE SelfMultiblockCompute
END MODULE Grid_Mod

