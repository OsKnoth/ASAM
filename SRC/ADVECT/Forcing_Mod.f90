MODULE Forcing_Mod

  USE Kind_Mod
  USE ReadProfile_Mod
  USE DataType_Mod
  USE Floor_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE Example_Mod
  IMPLICIT NONE

  INTEGER, PRIVATE :: InputUnitNudging=20
  INTEGER, PRIVATE :: InputUnitTendency=21
  INTEGER, PRIVATE :: InputUnitSurface=22
  INTEGER, PRIVATE :: nzData=1
  INTEGER, PRIVATE :: nzNudge=1
  INTEGER, PRIVATE :: NumDataNudging=1
  INTEGER, PRIVATE :: NumDataTendency=1
  INTEGER, PRIVATE :: NumDataSurface=1
  REAL(RealKind) :: dtNudging
  REAL(RealKind) :: dtTendency
  REAL(RealKind) :: dtSurface
  INTEGER :: uPosN,vPosN,ThetaPosN,RhoVPosN,RhoPosN
  TYPE(Vector1Cell_T), POINTER :: RhoProfile(:)
  TYPE(Vector1Cell_T), POINTER :: RhoProfile1(:)
  TYPE(Vector1Cell_T), POINTER :: RhoProfile2(:)
  TYPE(Vector1Cell_T), POINTER :: uG1(:)
  TYPE(Vector1Cell_T), POINTER :: uG2(:)
  TYPE(Vector1Cell_T), POINTER :: uG(:)
  TYPE(Vector1Cell_T), POINTER :: TendAdv1(:)
  TYPE(Vector1Cell_T), POINTER :: TendAdv2(:)
  TYPE(Vector1Cell_T), POINTER :: TendAdv(:)
  TYPE(Vector1Face_T), POINTER :: Subs1(:)
  TYPE(Vector1Face_T), POINTER :: Subs2(:)
  TYPE(Vector1Face_T), POINTER :: Subs(:)
  REAL(RealKind) :: TendTime1,TendTime2
  TYPE(Vector1Cell_T), POINTER :: NudgingProfile(:)
  TYPE(Vector1Cell_T), POINTER :: NudgingProfile1(:)
  TYPE(Vector1Cell_T), POINTER :: NudgingProfile2(:)
  TYPE(Vector1Cell_T), POINTER :: DampProfile(:)
  REAL(RealKind) :: NudgeTime1,NudgeTime2
  REAL(RealKind), ALLOCATABLE :: Surface(:)
  REAL(RealKind), ALLOCATABLE :: Surface1(:)
  REAL(RealKind), ALLOCATABLE :: Surface2(:)
  REAL(RealKind) :: SurfTime1,SurfTime2
  TYPE(Vector1Cell_T), POINTER :: ProfileMean(:)
  TYPE(Vector1Cell_T), POINTER :: RhoProfileMean(:)
  REAL(RealKind) :: RelaxTime1
  REAL(RealKind) :: RelaxTime2
  REAL(RealKind) :: RelaxTime
 

CONTAINS

SUBROUTINE PrepareForcing(VectorCell,Time)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  REAL(RealKind) :: Time

  INTEGER :: i,ic
  INTEGER :: ix,iy,iz

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    RhoProfile(ibLoc)%Vec(1)%c=((Time-NudgeTime1)*RhoProfile2(ibLoc)%Vec(1)%c &
                        +(NudgeTime2-Time)*RhoProfile1(ibLoc)%Vec(1)%c) &
                        /(NudgeTime2-NudgeTime1)
    IF (ForcingExternTendency) THEN                    
      uG(ibLoc)%Vec(1)%c=((Time-TendTime1)*uG2(ibLoc)%Vec(1)%c &
                        +(TendTime2-Time)*uG1(ibLoc)%Vec(1)%c) &
                      /(TendTime2-TendTime1)
      uG(ibLoc)%Vec(2)%c=((Time-TendTime1)*uG2(ibLoc)%Vec(2)%c &
                        +(TendTime2-Time)*uG1(ibLoc)%Vec(2)%c) &
                        /(TendTime2-TendTime1)
      TendAdv(ibLoc)%Vec(1)%c=((Time-TendTime1)*TendAdv2(ibLoc)%Vec(1)%c &
                             +(TendTime2-Time)*TendAdv1(ibLoc)%Vec(1)%c) &
                             /(TendTime2-TendTime1)
      TendAdv(ibLoc)%Vec(2)%c=((Time-TendTime1)*TendAdv2(ibLoc)%Vec(2)%c &
                             +(TendTime2-Time)*TendAdv1(ibLoc)%Vec(2)%c) &
                             /(TendTime2-TendTime1)
      Subs(ibLoc)%Vec(1)%c=((Time-TendTime1)*Subs2(ibLoc)%Vec(1)%c &
                          +(TendTime2-Time)*Subs1(ibLoc)%Vec(1)%c) &
                          /(TendTime2-TendTime1)
      TendAdv(ibLoc)%Vec(1)%c=TendAdv(ibLoc)%Vec(1)%c*RhoProfile(ibLoc)%Vec(1)%c
      TendAdv(ibLoc)%Vec(2)%c=TendAdv(ibLoc)%Vec(2)%c*RhoProfile(ibLoc)%Vec(1)%c
    END IF  
    
    DO ic=1,SIZE(NudgingProfile(ibLoc)%Vec)
      NudgingProfile(ibLoc)%Vec(ic)%c=((Time-NudgeTime1)*NudgingProfile2(ibLoc)%Vec(ic)%c &
                                     +(NudgeTime2-Time)*NudgingProfile1(ibLoc)%Vec(ic)%c) &
                                     /(NudgeTime2-NudgeTime1)
    ! IF (MyID==0) WRITE(*,*) 'Nudging: ',ic,NudgingProfile(ibLoc)%Vec(ic)%c(40)
    END DO  
    IF (ForcingExternSurface) THEN
      Surface=((Time-SurfTime1)*Surface2+(SurfTime2-Time)*Surface1) &
              /(SurfTime2-SurfTime1)
      DO i=1,NumBoundCell
        ix=Floor(ib)%BoundCell(i)%ix
        iy=Floor(ib)%BoundCell(i)%iy
        iz=Floor(ib)%BoundCell(i)%iz
        Floor(ib)%BoundCell(i)%TeS=Surface(1)
        Floor(ib)%BoundCell(i)%qv=Surface(2)*VectorCell(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1)
      END DO  
    END IF  
  END DO  
  RelaxTime=((Time-NudgeTime1)*RelaxTime2+(NudgeTime2-Time)*RelaxTime1) &
            /(NudgeTime2-NudgeTime1)
  RelaxTime=1.0d0/RelaxTime          

END SUBROUTINE PrepareForcing

SUBROUTINE DampProfileCompute

  INTEGER :: iz,ic
  REAL(RealKind) :: zM
  LOGICAL, SAVE :: Init=.TRUE.
  IF (Init) THEN
    CALL Allocate(DampProfile,NumDataNudging)
    Init=.FALSE.
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        zM=0.5d0*(zP(iz-1)+zP(iz))
        DO ic=1,SIZE(NudgingProfile2(ibLoc)%Vec)
          IF (PosE2Pos(ic)==uPosL.OR.PosE2Pos(ic)==vPosL) THEN
            DampProfile(ibLoc)%Vec(ic)%c(iz)=DampFun(zM,'VelProf')
          ELSE IF (PosE2Pos(ic)==RhoVPos) THEN
            DampProfile(ibLoc)%Vec(ic)%c(iz)=DampFun(zM,'RhoVProf')
          ELSE IF (PosE2Pos(ic)==thPos) THEN  
            DampProfile(ibLoc)%Vec(ic)%c(iz)=DampFun(zM,'ThProf')
          END IF  
        END DO  
      END DO  
    END DO  
  END IF
END SUBROUTINE DampProfileCompute

SUBROUTINE InputNudgingProfile(Time,FileName,Finish)
  REAL(RealKind) :: Time
  CHARACTER(*) :: FileName
  LOGICAL :: Finish

  INTEGER :: i,j
  INTEGER :: ic,iz
  REAL(RealKind) :: zM
  LOGICAL, SAVE :: Init=.TRUE.
  REAL(RealKind) :: cProfile(nzNudge,NumDataNudging)
  REAL(RealKind) :: RhoProf(nzNudge)
  REAL(RealKind) :: z(nzNudge)
  CHARACTER(300) :: Line
  REAL(RealKind) :: Dummy
  REAL(RealKind) :: RelHumLoc,TLoc,pLoc,RhoDLoc,RhoLLoc,RhoVLoc,RhoPotLoc


  IF (Init) THEN
    OPEN(UNIT=InputUnitNudging,FILE=FileName,STATUS='UNKNOWN')
    READ(InputUnitNudging,*) 
    READ(InputUnitNudging,*) 
    READ(InputUnitNudging,*) nzNudge,NumDataNudging
    WRITE(*,*) 'Nudge numbers',nzNudge,NumDataNudging
    Init=.FALSE.
!   dtNudging=3.0d0*3600.0d0
    CALL Allocate(ProfileMean,VectorComponentsME+2)
    CALL Allocate(RhoProfileMean,1)
    CALL Allocate(NudgingProfile,NumDataNudging)
    CALL Allocate(NudgingProfile1,NumDataNudging)
    CALL Allocate(NudgingProfile2,NumDataNudging)
    CALL Allocate(RhoProfile,1)
    CALL Allocate(RhoProfile1,1)
    CALL Allocate(RhoProfile2,1)
    DO 
      READ(InputUnitNudging,'(a300)',END=1) Line
      IF (Line(1:1)=='#') THEN 
        BACKSPACE InputUnitNudging
        BACKSPACE InputUnitNudging
        BACKSPACE InputUnitNudging
        EXIT
      END IF  
    END DO  
1   CONTINUE    
    Init=.FALSE.
  ELSE
    CALL Copy(NudgingProfile2,NudgingProfile1)
    CALL Copy(RhoProfile2,RhoProfile1)
    READ(InputUnitNudging,*)
    READ(InputUnitNudging,*)
    READ(InputUnitNudging,'(a300)') Line
    NudgeTime1=NudgeTime2
    RelaxTime1=RelaxTime2
    READ(Line(2:),*) NudgeTime2
    WRITE(*,*) 'NudgeTime1,NudgeTime2',NudgeTime1,NudgeTime2
    DO i=1,nzNudge
      READ(InputUnitNudging,*) z(i),RelaxTime2,cProfile(i,uPosEnv),cProfile(i,vPosEnv),Dummy,cProfile(i,ThPosEnv) &
                                  ,cProfile(i,RhoVPosEnv),RhoProf(i)
    END DO                              
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        zM=0.5d0*(zP(iz-1)+zP(iz))
        DO ic=1,SIZE(NudgingProfile2(ibLoc)%Vec)
          NudgingProfile2(ibLoc)%Vec(ic)%c(iz)=ProfileEqual(cProfile(:,ic),z,zM)
        END DO  
        RhoProfile2(ibLoc)%Vec(1)%c(iz)=ProfileEqual(RhoProf,z,zM)
      END DO  
    END DO  
    IF (Time>=NudgeTime1.AND.Time<NudgeTime2) THEN
      Finish=.TRUE.
      RETURN
    END IF
  END IF
END SUBROUTINE InputNudgingProfile

SUBROUTINE InputSurface(Time,FileName)
  REAL(RealKind) :: Time
  CHARACTER(*) :: FileName

  REAL(RealKind) :: Dummy
  LOGICAL, SAVE :: Init=.TRUE.

  IF (Init) THEN
    OPEN(UNIT=InputUnitSurface,FILE=FileName,STATUS='UNKNOWN')
    READ(InputUnitSurface,*) 
    READ(InputUnitSurface,*) 
    READ(InputUnitSurface,*) 
    READ(InputUnitSurface,*) NumDataSurface
    Init=.FALSE.
    dtSurface=3.0d0*3600.0d0
    ALLOCATE(Surface(NumDataSurface))
    ALLOCATE(Surface1(NumDataSurface))
    ALLOCATE(Surface2(NumDataSurface))
  ELSE
    Surface1=Surface2
    SurfTime1=SurfTime2
    READ(InputUnitSurface,*) SurfTime2,Dummy,Dummy,Surface2
  END IF
END SUBROUTINE InputSurface

SUBROUTINE InputTendencyProfile(Time,FileName,Finish)
  REAL(RealKind) :: Time
  CHARACTER(*) :: FileName
  LOGICAL :: Finish

  INTEGER :: i,j
  INTEGER :: iz
  REAL(RealKind) :: zM
  LOGICAL, SAVE :: Init=.TRUE.
  REAL(RealKind) :: z(nzData)
  REAL(RealKind) :: cProfile(nzData,NumDataTendency)
  CHARACTER(300) :: Line

  IF (Init) THEN
    OPEN(UNIT=InputUnitTendency,FILE=FileName,STATUS='UNKNOWN')
    READ(InputUnitTendency,*) 
    READ(InputUnitTendency,*) nzData,NumDataTendency
    CALL Allocate(uG1,2)
    CALL Allocate(uG2,2)
    CALL Allocate(uG,2)
    CALL Allocate(TendAdv1,2)
    CALL Allocate(TendAdv2,2)
    CALL Allocate(TendAdv,2)
    CALL Allocate(Subs1,1)
    CALL Allocate(Subs2,1)
    CALL Allocate(Subs,1)
    DO 
      READ(InputUnitTendency,'(a300)',END=1) Line
      IF (Line(1:1)=='#') THEN 
        BACKSPACE InputUnitTendency
        BACKSPACE InputUnitTendency
        BACKSPACE InputUnitTendency
        EXIT
      END IF  
    END DO  
1   CONTINUE    
    Init=.FALSE.
  ! dtTendency=3.0d0*3600.0d0
  ELSE  
    CALL Copy(uG2,uG1)
    CALL Copy(TendAdv2,TendAdv1)
    CALL Copy(Subs2,Subs1)
    READ(InputUnitTendency,*)
    READ(InputUnitTendency,*)
    READ(InputUnitTendency,'(a300)') Line
    TendTime1=TendTime2
    READ(Line(2:),*) TendTime2
    WRITE(*,*) 'TendTime1,TendTime2',TendTime1,TendTime2
    DO i=1,nzData
      READ(InputUnitTendency,*) z(i),cProfile(i,:)
    END DO  
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        zM=0.5d0*(zP(iz-1)+zP(iz))
        uG2(ibLoc)%Vec(1)%c(iz)=ProfileEqual(cProfile(:,1),z,zM)
        uG2(ibLoc)%Vec(2)%c(iz)=ProfileEqual(cProfile(:,2),z,zM)
        TendAdv2(ibLoc)%Vec(1)%c(iz)=ProfileEqual(cProfile(:,4),z,zM)
        TendAdv2(ibLoc)%Vec(2)%c(iz)=ProfileEqual(cProfile(:,5),z,zM)
   !    WRITE (*,*) 'TendAdv2 1',iz,TendAdv2(ibLoc)%Vec(1)%c(iz)
      END DO  
!      Subs2(ibLoc)%Vec(1)%c(iz0)=Zero
!      DO iz=iz0+1,iz1
      DO iz=iz0,iz1
        Subs2(ibLoc)%Vec(1)%c(iz)=ProfileEqual(cProfile(:,3),z,zP(iz))
      END DO  
    END DO  
    IF (Time>=TendTime1.AND.Time<TendTime2) THEN
      Finish=.TRUE.
      RETURN
    END IF
  END IF
END SUBROUTINE InputTendencyProfile

SUBROUTINE MeanVerticalProfile(VectorCell,UVec)
  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(Vector4Cell_T), POINTER, OPTIONAL :: UVec(:)

  INTEGER :: ix,iy,iz
  INTEGER :: ic0,ic
  INTEGER :: igz
  REAL(RealKind) :: ProfileLoc(Domain%nz,SIZE(ProfileMean(1)%Vec)+1)
  REAL(RealKind) :: ProfileSingle(Domain%nz,SIZE(ProfileMean(1)%Vec)+1)
  REAL(RealKind) :: VolLoc(Domain%nz)
  REAL(RealKind) :: Vol(Domain%nz)
  REAL(RealKind) :: dzLoc
  REAL(RealKind), POINTER :: c(:,:,:,:)

  ProfileLoc=Zero
  ProfileSingle=Zero
  VolLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ic=1
    CALL SingleVol
    DO ic0=1,VectorComponentsME
      IF (PosE2Pos(ic0)==uPosL) THEN
        IF (PRESENT(UVec)) THEN
          c=>UVec(ibLoc)%Vec(uPosL)%c
          CALL SingleProfile
          ic=ic+1
          c=>UVec(ibLoc)%Vec(uPosR)%c
          CALL SingleProfile
          ic=ic+1
        ELSE  
          c=>VectorCell(ibLoc)%Vec(uPosL)%c
          CALL SingleProfile
          ic=ic+1
          c=>VectorCell(ibLoc)%Vec(uPosR)%c
          CALL SingleProfile
          ic=ic+1
        END IF  
      ELSE IF (PosE2Pos(ic0)==vPosL) THEN
        IF (PRESENT(UVec)) THEN
          c=>UVec(ibLoc)%Vec(vPosL)%c
          CALL SingleProfile
          ic=ic+1
          c=>UVec(ibLoc)%Vec(vPosR)%c
          CALL SingleProfile
          ic=ic+1
        ELSE  
          c=>VectorCell(ibLoc)%Vec(vPosL)%c
          CALL SingleProfile
          ic=ic+1
          c=>VectorCell(ibLoc)%Vec(vPosR)%c
          CALL SingleProfile
          ic=ic+1
        END IF  
      ELSE IF (PosE2Pos(ic0)==RhoVPos) THEN
        c=>VectorCell(ibLoc)%Vec(RhoVPos)%c
        CALL SingleProfile
        IF (RhoCPos>0) THEN
          c=>VectorCell(ibLoc)%Vec(RhoCPos)%c
          CALL SingleProfile
        END IF  
        IF (RhoRPos>0) THEN
          c=>VectorCell(ibLoc)%Vec(RhoRPos)%c
          CALL SingleProfile
        END IF
        IF (RhoIPos>0) THEN
          c=>VectorCell(ibLoc)%Vec(RhoIPos)%c
          CALL SingleProfile
        END IF
        ic=ic+1
      ELSE  
        c=>VectorCell(ibLoc)%Vec(PosE2Pos(ic0))%c
        CALL SingleProfile
        ic=ic+1
      END IF  
    END DO
    c=>VectorCell(ibLoc)%Vec(RhoPos)%c
    CALL SingleProfile
  END DO
  CALL MPI_Allreduce(ProfileLoc,ProfileSingle,SIZE(ProfileSingle),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  CALL MPI_Allreduce(VolLoc,Vol,SIZE(Vol),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  DO ic=1,SIZE(ProfileSingle,2)
    ProfileSingle(:,ic)=ProfileSingle(:,ic)/(Vol+Eps)
  END DO  
  ProfileMean=Zero
  RhoProfileMean=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iz=iz0+1,iz1
      dzLoc=Zero
      DO igz=iz**(-RefineZ+1),(iz+1)**(-RefineZ+1)-1 
        dzLoc=dzLoc+Domain%dz(igz)
        DO ic=1,SIZE(ProfileMean(ibLoc)%Vec)
          ProfileMean(ibloc)%Vec(ic)%c(iz)=ProfileMean(iBloc)%Vec(ic)%c(iz)+ProfileSingle(igz,ic)*Domain%dz(igz)
        END DO  
        RhoProfileMean(ibloc)%Vec(1)%c(iz)=RhoProfileMean(iBloc)%Vec(1)%c(iz) &
                                          +ProfileSingle(igz,SIZE(ProfileMean(ibLoc)%Vec)+1)*Domain%dz(igz)
      END DO  
      DO ic=1,SIZE(ProfileMean(ibLoc)%Vec)
        ProfileMean(ibloc)%Vec(ic)%c(iz)=ProfileMean(iBloc)%Vec(ic)%c(iz)/dzLoc
      END DO  
      RhoProfileMean(ibloc)%Vec(1)%c(iz)=RhoProfileMean(iBloc)%Vec(1)%c(iz)/dzLoc
    END DO  
  END DO  
CONTAINS
SUBROUTINE SingleProfile
  DO iz=iz0+1,iz1
    DO igz=iz**(-RefineZ+1),(iz+1)**(-RefineZ+1)-1 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Fac=VolC(ix,iy,iz)*Domain%dz(igz)/dz(iz)
          ProfileLoc(igz,ic)=ProfileLoc(igz,ic)+Fac*c(ix,iy,iz,1)
        END DO  
      END DO  
    END DO  
  END DO
END SUBROUTINE SingleProfile

SUBROUTINE SingleVol
  DO iz=iz0+1,iz1
    DO igz=iz**(-RefineZ+1),(iz+1)**(-RefineZ+1)-1 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Fac=VolC(ix,iy,iz)*Domain%dz(igz)/dz(iz)
          VolLoc(igz)=VolLoc(igz)+Fac
        END DO  
      END DO  
    END DO  
  END DO
END SUBROUTINE SingleVol

END SUBROUTINE MeanVerticalProfile
END MODULE Forcing_Mod
