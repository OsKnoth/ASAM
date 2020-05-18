MODULE WindFarm_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Floor_Mod
  USE InputTool_Mod
  USE Physics_Mod

  IMPLICIT NONE
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vRhsR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wRhsL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wRhsR(:,:,:,:)

  TYPE WindTurbine_T 
    
    REAL(RealKind) :: xPos=0.0d0
    REAL(RealKind) :: yPos=0.0d0
    REAL(RealKind) :: zPos=0.0d0
    REAL(RealKind) :: Rad=0.0d0
    REAL(RealKind) :: Alpha=0.0d0
    
    REAL(RealKind) :: WindAcc(2)=0
    REAL(RealKind) :: LastTime=0
    REAL(RealKind) :: RotorPosition=0
    REAL(RealKind) :: LastTurn=0

    
  END TYPE WindTurbine_T 

  TYPE WindTurbineInfluence_T
    INTEGER :: Number=0
    INTEGER, POINTER :: Index(:,:)
   REAL(RealKind) , POINTER :: RadiusVec(:,:),Theta(:),Lambda(:),d(:)
   REAL(Realkind) , POINTER :: x(:),y(:),z(:)
  END TYPE WindTurbineInfluence_T

  TYPE(WindTurbine_T), ALLOCATABLE :: WindFarm(:)
  TYPE(WindTurbineInfluence_T), ALLOCATABLE :: WindTurbineInfluence(:,:)

  REAL(Realkind), ALLOCATABLE :: MeanWindGlob(:,:),MeanWindLoc(:,:),MeanWindNumberLoc(:),MeanWindNumberGlob(:),&
  &MeanRhoLoc(:),MeanRhoGlob(:)
  
  INTEGER :: NumberW

  REAL(Realkind), PARAMETER :: B1=3.0d0
  REAL(Realkind), PARAMETER :: B2=1.5d0 
  REAL(Realkind), PARAMETER :: CD0=2.5d0



CONTAINS

FUNCTION WindDragCoefficient(R,Rad,d,Theta)

  

  REAL(Realkind) :: WindDragCoefficient
  REAL(Realkind) :: R,Rad,d,Theta

  REAL(Realkind), PARAMETER :: B1=3.0d0
  REAL(Realkind), PARAMETER :: B2=1.5d0 
  REAL(Realkind), PARAMETER :: CD0=2.0d0 
  REAL(Realkind), PARAMETER :: WindRotCoeff=0.5d0
  

! 'Sphere1' -> all cells with r <= turbine radius have drag Cd(r) according to Heimann et al. 2011,
! 'Sphere2' -> all cells with r <= turbine radius have drag Cd(r) (modified)  according to Heimann et al. 2011,
! 'Disk' -> Cd(r) of 'Sphere2' + convolution to disk 
! 'Rot' -> 'Disk' + convolution rotor blades + rotation, rotation speed see  wind tip speed ratio in SUBROUTINE ROTATION


  IF (R>=5.0d0) THEN
    SELECT CASE(WindDragScheme)
      CASE('Sphere1')
        WindDragCoefficient=CD0*B1*3.0d0/(Pi*r)
      CASE('Sphere2')
        WindDragCoefficient=CD0*3.0d0*(B1+r*(B2-B1)/Rad)/(Pi*r)
      CASE('Disk')
        WindDragCoefficient=CD0*3.0d0*(B1+r*(B2-B1)/Rad)/(Pi*r)
        WindDragCoefficient=WindDragCoefficient*Kerneld(d)
      CASE('Rot')
        WindDragCoefficient=WindRotCoeff
        WindDragCoefficient=WindDragCoefficient*Kerneld(d)
        WindDragCoefficient=WindDragCoefficient*KernelTheta(Theta)
      END SELECT
  ElSE
     WindDragCoefficient=CD0*B1*3.0d0/(Pi*5.0d0)
  END IF
  WindDragCoefficient=WindDragCoefficient*MIN(TimeAct/10.0d0,1.0d0) !OSSI

END FUNCTION WindDragCoefficient

FUNCTION Kerneld(d)

  REAL(Realkind) :: Kerneld
  REAL(Realkind) :: d

  REAL(Realkind), PARAMETER :: Sigma=10, A=1/(SQRT(3.141)*5)

  Kerneld=A*EXP(-(d/Sigma)**2)

END FUNCTION Kerneld

FUNCTION KernelTheta(Theta)

  REAL(Realkind) :: KernelTheta
  REAL(Realkind) :: Theta

  REAL(Realkind), PARAMETER :: Sigma=0.2d0, A=6.2d0
  
  !blades position at 0,2/3*PI,4/3*PI
  
  KernelTheta=A*(EXP(-(Theta/Sigma)**2)+EXP(-((Theta-2*PI/3)/Sigma)**2)&
  &           +EXP(-((Theta-4*PI/3)/Sigma)**2)+EXP(-((Theta-2*PI)/Sigma)**2))

END FUNCTION KernelTheta

FUNCTION NormalWindMagnitude(u,v,w,Alpha) !computes WindMagnitude normal to wind turbine

  REAL(Realkind) :: u,v,w,Alpha,NormalWindMagnitude
  
 NormalWindMagnitude=(SIN(Alpha)*u-COS(Alpha)*v)**2/(SQRT(u*u+v*v+w*w))


END FUNCTION NormalWindMagnitude

SUBROUTINE InitWind(FileName)
  
  CHARACTER(*) :: FileName
  CHARACTER(20) :: S1,S2,End
  CHARACTER*300 :: Line
  LOGICAL :: Back
  REAL(RealKind) :: xW,yW
  REAL(RealKind) :: RadW
  REAL(RealKind) :: HeightW
  REAL(RealKind) :: Alpha,R,Lambda,RXY
  REAL(RealKind) :: xM,yM,zM
  REAL(RealKind) :: zPosLoc
  INTEGER :: ix,iy,iz,iw,k
  INTEGER, ALLOCATABLE :: ind(:,:,:,:)
  
  S1='BEGIN_WIND'
  End='END_WIND'

  CALL OpenFile(FileName)
  NumberW=0
  DO
    CALL LineFile(Back,S1,End=End,R1=xW,R2=yW)
    IF (Back) THEN
      EXIT
    END IF
    NumberW=NumberW+1
    CALL LineFile(Back,S1,End=End,R1=RadW)
    CALL LineFile(Back,S1,End=End,R1=HeightW)
    CALL LineFile(Back,S1,End=End,R1=Alpha)
  END DO
  CALL CloseFile

  ALLOCATE(WindFarm(NumberW))
  ALLOCATE(MeanWindGlob(NumberW,3))
  ALLOCATE(MeanWindLoc(NumberW,3))
  ALLOCATE(MeanWindNumberLoc(NumberW))
  ALLOCATE(MeanWindNumberGlob(NumberW))
  ALLOCATE(MeanRhoLoc(NumberW))
  ALLOCATE(MeanRhoGlob(NumberW))

  MeanWindGlob=0
  MeanRhoGlob=0

  CALL OpenFile(FileName)
  NumberW=0
  DO
    CALL LineFile(Back,S1,End=End,R1=xW,R2=yW)
    IF (Back) THEN
      EXIT
    END IF
    NumberW=NumberW+1
    CALL LineFile(Back,S1,End=End,R1=RadW)
    CALL LineFile(Back,S1,End=End,R1=HeightW)
    CALL LineFile(Back,S1,End=End,R1=Alpha)
    
    WindFarm(NumberW)%xPos=xW
    WindFarm(NumberW)%yPos=yW
    WindFarm(NumberW)%Rad=RadW
    WindFarm(NumberW)%Alpha=Alpha

    DO ib=1,nb
      IF (MyId==blMPI(ib)%Proc) THEN
        CALL Set(Floor(ib))
        IF (zP(iz0)==domain%zP(domain%igz0)) THEN
          W1:DO ix=ix0+1,ix1
            DO iy=iy0+1,iy1
              IF (xP(ix-1)<WindFarm(NumberW)%xPos.AND.WindFarm(NumberW)%xPos<=xP(ix).AND. &
&                 yP(iy-1)<WindFarm(NumberW)%yPos.AND.WindFarm(NumberW)%yPos<=yP(iy)) THEN
                WindFarm(NumberW)%zPos=HeightW+zH(ix,iy) !z= topo height + wind turbine height
                EXIT W1
              END IF
            END DO
          END DO W1
        END IF  
      END IF  
      CALL MPI_Bcast(WindFarm(NumberW)%zPos,1,MPI_RealKind,blMPI(ib)%Proc, &
 &                   MPI_COMM_WORLD,MPIErr)
    END DO
  END DO
  CALL CloseFile

  ALLOCATE(WindTurbineInfluence(NumberW,nbloc))
  DO ibloc=1,nbLoc
    ib=LocGlob(ibloc)
    CALL Set(Floor(ib))
    DO iW=1,NumberW
      WindTurbineInfluence(iW,ibloc)%Number=0
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          xM=0.5d0*(xP(ix-1)+xP(ix)) 
          yM=0.5d0*(yP(iy-1)+yP(iy)) 
          zM=0.5d0*(zP(iz-1)+zP(iz)) 
          DO iW=1,NumberW
            IF (SQRT((ym-WindFarm(iW)%yPos)**2+(xm-WindFarm(iW)%xPos)**2+(zm-WindFarm(iW)%zPos)**2)<WindFarm(iW)%Rad) THEN
              WindTurbineInfluence(iW,ibloc)%Number=WindTurbineInfluence(iW,ibloc)%Number+1
            END IF
          END DO
        END DO
      END DO
    END DO
    DO iW=1,NumberW
      ALLOCATE(WindTurbineInfluence(iW,ibloc)%Index(3,WindTurbineInfluence(iW,ibloc)%Number))
      ALLOCATE(WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%Number,3))
      ALLOCATE(WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%Number))
      ALLOCATE(WindTurbineInfluence(iW,ibloc)%Lambda(WindTurbineInfluence(iW,ibloc)%Number))
      ALLOCATE(WindTurbineInfluence(iW,ibloc)%d(WindTurbineInfluence(iW,ibloc)%Number))
    END DO
  END DO 
  DO ibloc=1,nbLoc
    ib=LocGlob(ibloc)
    CALL Set(Floor(ib))
    DO iW=1,NumberW
      WindTurbineInfluence(iW,ibloc)%Number=0
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          xM=0.5d0*(xP(ix-1)+xP(ix)) 
          yM=0.5d0*(yP(iy-1)+yP(iy)) 
          zM=0.5d0*(zP(iz-1)+zP(iz)) 
          DO iW=1,NumberW
            IF (SQRT((xm-WindFarm(iW)%xPos)**2 &
                    +(ym-WindFarm(iW)%yPos)**2 &
                    +(zm-WindFarm(iW)%zPos)**2)<WindFarm(iW)%Rad) THEN
              WindTurbineInfluence(iW,ibloc)%Number=WindTurbineInfluence(iW,ibloc)%Number+1
              WindTurbineInfluence(iW,ibloc)%Index(1,WindTurbineInfluence(iW,ibloc)%Number)=ix
              WindTurbineInfluence(iW,ibloc)%Index(2,WindTurbineInfluence(iW,ibloc)%Number)=iy
              WindTurbineInfluence(iW,ibloc)%Index(3,WindTurbineInfluence(iW,ibloc)%Number)=iz
             
              ! R = vector to wind turbine center
              WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,1)=xm-WindFarm(iW)%xPos 
              WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,2)=ym-WindFarm(iW)%yPos 
              WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,3)=zm-WindFarm(iW)%zPos 
              
              ! Lambda = angle R to x axes in xy-plane
              IF (xm-WindFarm(iW)%xPos/=0.0d0) THEN
                WindTurbineInfluence(iW,ibloc)%Lambda(WindTurbineInfluence(iW,ibloc)%number)=&
                &                             ATAN((ym-WindFarm(iW)%yPos)/(xm-WindFarm(iW)%xPos))
              ELSE
                WindTurbineInfluence(iW,ibloc)%Lambda(WindTurbineInfluence(iW,ibloc)%number)=PI/2
              END IF
              
              R=SQRT(WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,1)**2 &
                    +WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,2)**2 &
                    +WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,3)**2)   !explained in Ad_rho -> WindDragCompute 
              Lambda=WindTurbineInfluence(iW,ibloc)%Lambda(WindTurbineInfluence(iW,ibloc)%number) 
              RXY=SQRT(WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,1)**2 &
                      +WindTurbineInfluence(iW,ibloc)%RadiusVec(WindTurbineInfluence(iW,ibloc)%number,2)**2) 
              WindTurbineInfluence(iW,ibloc)%d(WindTurbineInfluence(iW,ibloc)%number)=SIN(WindFarm(iW)%Alpha-Lambda)*RXY  

              ! calculate Theta = angle R to z axes in yz-plane
              IF ((ym-WindFarm(iW)%yPos)>=0.AND.(zm-WindFarm(iW)%zPos)>0) THEN  
                WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%number)=&
                &ATAN((ym-WindFarm(iW)%yPos)/(zm-WindFarm(iW)%zPos)) 
              ELSE IF ((ym-WindFarm(iW)%yPos)>=0.AND.(zm-WindFarm(iW)%zPos)==0) THEN  
                WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%number)=PI/2 
              ELSE IF ((ym-WindFarm(iW)%yPos)>0.AND.(zm-WindFarm(iW)%zPos)<0) THEN  
                WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%number)=PI/2+&
                & ABS(ATAN((zm-WindFarm(iW)%zPos)/(ym-WindFarm(iW)%yPos))) 
              ELSE IF ((ym-WindFarm(iW)%yPos)==0.AND.(zm-WindFarm(iW)%zPos)<0) THEN  
                WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%number)=PI 
              ELSE IF ((ym-WindFarm(iW)%yPos)<0.AND.(zm-WindFarm(iW)%zPos)<0) THEN  
                WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%number)=PI+&
                & ABS(ATAN((ym-WindFarm(iW)%yPos)/(zm-WindFarm(iW)%zPos))) 
              ELSE IF ((ym-WindFarm(iW)%yPos)<0.AND.(zm-WindFarm(iW)%zPos)==0) THEN  
                WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%number)=3*PI/2 
              ELSE IF ((ym-WindFarm(iW)%yPos)<0.AND.(zm-WindFarm(iW)%zPos)>0) THEN  
                WindTurbineInfluence(iW,ibloc)%Theta(WindTurbineInfluence(iW,ibloc)%number)=3*PI/2+&
                &ABS(ATAN((zm-WindFarm(iW)%zPos)/(ym-WindFarm(iW)%yPos))) 
              END IF   
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE InitWind

SUBROUTINE ROTATION

  INTEGER :: Rotations,iW

  REAL(RealKind) :: MeanWindAcc(2)
  REAL(RealKind) :: AngVel

  REAL(RealKind), PARAMETER :: TSR=7.0d0
  
  DO iW=1,NumberW
    IF (TimeAct/=WindFarm(iW)%LastTime) THEN
     
      !tip speed ratio formula: AngVel=TSR*WindSpeed/Rad
      
      AngVel=TSR*(SIN(WindFarm(iW)%Alpha)*MeanWindGlob(iW,1)-COS(WindFarm(iW)%Alpha)*MeanWindGlob(iW,2))/WindFarm(iW)%Rad 
      
      !Compute CurrentRotorPosition
      
      Rotations=AINT((AngVel*(TimeAct-WindFarm(iW)%LastTime)+WindFarm(iW)%RotorPosition)/(2*Pi))
      WindFarm(iW)%RotorPosition=WindFarm(iW)%RotorPosition+AngVel*(TimeAct-Windfarm(iW)%LastTime)-2*Pi*Rotations
      WindFarm(iW)%LastTime=TimeAct
      
    END IF
  END DO
END SUBROUTINE Rotation    

SUBROUTINE UpdateWind(VelF,VectorCell,dtAct)
  
  TYPE(VelocityFace_T) :: VelF(:)
  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)

  REAL(RealKind) :: dtAct

  INTEGER :: ix,iy,iz
  INTEGER :: i,iW,n
  REAL(RealKind) :: TempU,TempV,TempW
  REAL(RealKind) :: MeanWindAcc(2),Radius(3)
  REAL(RealKind) :: TurnStep=60.0d0 
  REAL(RealKind), POINTER :: uF(:,:,:)
  REAL(RealKind), POINTER :: vF(:,:,:)
  REAL(RealKind), POINTER :: wF(:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)

    MeanWindLoc=0
    MeanWindGlob=0
    MeanWindNumberLoc=0
    MeanWindNumberGlob=0
    MeanRhoLoc=0
    MeanRhoGlob=0
    DO ibloc=1,nbLoc
      ib=LocGlob(ibloc)
      CALL Set(Floor(ib))
      uF=>VelF(ibLoc)%uF
      vF=>VelF(ibLoc)%vF
      wF=>VelF(ibLoc)%wF
      Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
      DO iW=1,NumberW
        n=WindTurbineInfluence(iW,ibloc)%Number
        DO i=1,n !MeanWindGlob is mean wind of half sphere in upstream direction (where wind comes from) in front of wind turbine
          Radius=WindTurbineInfluence(iW,ibloc)%RadiusVec(i,:)
          IF (Radius(1)*SIN(WindFarm(iW)%Alpha)-Radius(2)*COS(WindFarm(iW)%Alpha)<=0) THEN 
            ix=WindTurbineInfluence(iW,ibloc)%index(1,i) 
            iy=WindTurbineInfluence(iW,ibloc)%index(2,i) 
            iz=WindTurbineInfluence(iW,ibloc)%index(3,i) 
            
            TempU=(uF(ix,iy,iz)*FU(ix,iy,iz)+uF(ix-1,iy,iz)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
            TempV=(vF(ix,iy,iz)*FV(ix,iy,iz)+vF(ix,iy-1,iz)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
            TempW=(wF(ix,iy,iz)*FW(ix,iy,iz)+wF(ix,iy,iz-1)*FW(ix,iy,iz-1))/(FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)

            MeanWindLoc(iW,1)=MeanWindLoc(iW,1)+TempU
            MeanWindLoc(iW,2)=MeanWindLoc(iW,2)+TempV
            MeanWindLoc(iW,3)=MeanWindLoc(iW,3)+TempW
            MeanRhoLoc(iW)=MeanRhoLoc(iW)+Rho(ix,iy,iz,1)
            MeanWindNumberLoc(iW)=MeanWindNumberLoc(iW)+1
          END IF
        END DO
      END DO
    END DO
    CALL MPI_Allreduce(MeanWindLoc,MeanWindGlob,3*NumberW,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
    CALL MPI_Allreduce(MeanRhoLoc,MeanRhoGlob,NumberW,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
    CALL MPI_Allreduce(MeanWindNumberLoc,MeanWindNumberGlob,NumberW,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
    DO iW=1,NumberW
      MeanWindGlob(iW,:)=MeanWindGlob(iW,:)/MeanWindNumberGlob(iW)
      MeanRhoGlob(iW)=MeanRhoGlob(iW)/MeanWindNumberGlob(iW)
      WindFarm(iW)%WindAcc(1)=WindFarm(iW)%WindAcc(1)+MeanWindGlob(iW,1)
      WindFarm(iW)%WindAcc(2)=WindFarm(iW)%WindAcc(2)+MeanWindGlob(iW,2)
    END DO  
    !WindAcc accumulates MeanWind to calculate new wind turbine turn at each Turnstep*dtAct time
      
    
  IF (WindTurbineTurn) THEN !Turn wind turbine in direction of mean wind to maintain maixmum normal wind component
    DO iW=1,NumberW
      IF (TimeAct>=WindFarm(iW)%LastTurn+TurnStep) THEN
        MeanWindAcc(1)=WindFarm(iW)%WindAcc(1)*dtAct/TurnStep
        MeanWindAcc(2)=WindFarm(iW)%WindAcc(2)*dtAct/TurnStep
        IF (MeanWindAcc(1)>=0.AND.MeanWindAcc(2)<0) THEN
          WindFarm(iW)%Alpha=ATAN(-MeanWindAcc(1)/MeanWindAcc(2))
        ELSE IF (MeanWindAcc(1)>=0.AND.MeanWindAcc(2)==0) THEN
          WindFarm(iW)%Alpha=PI/2
        ELSE IF (MeanWindAcc(1)>0.AND.MeanWindAcc(2)>=0) THEN
          WindFarm(iW)%Alpha=ATAN(MeanWindAcc(2)/MeanWindAcc(1))+PI/2
        ELSE IF (MeanWindAcc(1)==0.AND.MeanWindAcc(2)>0) THEN
          WindFarm(iW)%Alpha=PI
        ELSE IF (MeanWindAcc(1)<=0.AND.MeanWindAcc(2)>0) THEN
          WindFarm(iW)%Alpha=ATAN(-MeanWindAcc(1)/MeanWindAcc(2))+PI
        ELSE IF (MeanWindAcc(1)<0.AND.MeanWindAcc(2)==0) THEN
          WindFarm(iW)%Alpha=3*PI/2
        ELSE IF (MeanWindAcc(1)<0.AND.MeanWindAcc(2)<=0) THEN
          WindFarm(iW)%Alpha=ATAN(MeanWindAcc(2)/MeanWindAcc(1))+3*PI/2
        END IF  
      END IF
    END DO
  END IF
END SUBROUTINE UpdateWind

SUBROUTINE WindDrag(Vector,Rhs,UVec)
  TYPE(Vector4Cell_T) :: Vector,Rhs
  TYPE (Vector4Cell_T), OPTIONAL :: UVec

  Rho=>RhoCell(ibLoc)%c
  IF (PRESENT(UVec)) THEN
    uCL=>UVec%Vec(uPosL)%c
    vCL=>UVec%Vec(vPosL)%c
    wCL=>UVec%Vec(wPosL)%c
    uCR=>UVec%Vec(uPosR)%c
    vCR=>UVec%Vec(vPosR)%c
    wCR=>UVec%Vec(wPosR)%c
  ELSE  
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
  END IF  
  uRhsL=>Rhs%Vec(uPosL)%c
  uRhsR=>Rhs%Vec(uPosR)%c
  vRhsL=>Rhs%Vec(vPosL)%c
  vRhsR=>Rhs%Vec(vPosR)%c
  wRhsL=>Rhs%Vec(wPosL)%c
  wRhsR=>Rhs%Vec(wPosR)%c

  CALL WindDragCompute
END SUBROUTINE WindDrag

SUBROUTINE WindDragCompute

  INTEGER :: i,ix,iy,iz,iW,n
  REAL(RealKind) :: TempU,TempV,TempW,TempVM ! velocity components and mean velocity
  REAL(RealKind) :: Cd ! drag coefficient of the wind turbine 
  
  REAL(RealKind) :: R,Rad,d,Theta,Lambda,RXY

  IF (WindDragScheme=='Rot') THEN
    CALL Rotation 
  END IF  

  DO iW=1,NumberW 
    n=WindTurbineInfluence(iW,ibloc)%Number
    DO i=1,n
      ix=WindTurbineInfluence(iW,ibloc)%index(1,i) 
      iy=WindTurbineInfluence(iW,ibloc)%index(2,i) 
      iz=WindTurbineInfluence(iW,ibloc)%index(3,i) 
     
      R=SQRT(WindTurbineInfluence(iW,ibloc)%RadiusVec(i,1)**2+WindTurbineInfluence(iW,ibloc)%RadiusVec(i,2)**2+ &
&       WindTurbineInfluence(iW,ibloc)%RadiusVec(i,3)**2)    !distance cell to wind turbine center, norm RVec
      Theta=WindTurbineInfluence(iW,ibloc)%Theta(i) !angle R to z-axes in yz plane
      Lambda=WindTurbineInfluence(iW,ibloc)%Lambda(i) !angle R to x axes in xy plane
      RXY=SQRT(WindTurbineInfluence(iW,ibloc)%RadiusVec(i,1)**2+WindTurbineInfluence(iW,ibloc)%RadiusVec(i,2)**2) !RVec proj xy plane
      d=SIN(WindFarm(iW)%Alpha-Lambda)*RXY !disntance to WindTurbine plane 
      WindTurbineInfluence(iW,ibloc)%d(i)=d

      Rad=Windfarm(iW)%Rad !Rotor radius
       
      Cd=WindDragCoefficient(R,Rad,d,ABS(Theta-WindFarm(iW)%RotorPosition))
      
      
      TempU=(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
      TempV=(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
      TempW=(wCR(ix,iy,iz,1)*FW(ix,iy,iz)+wCL(ix,iy,iz,1)*FW(ix,iy,iz-1))/(FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)
        
      TempVM=NormalWindMagnitude(TempU,TempV,TempW,WindFarm(iW)%Alpha)
        
      uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)-Cd*TempVM*TempU/(Rho(ix,iy,iz,1)+Eps)
      uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)-Cd*TempVM*TempU/(Rho(ix,iy,iz,1)+Eps)
      vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-Cd*TempVM*TempV/(Rho(ix,iy,iz,1)+Eps)
      vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)-Cd*TempVM*TempV/(Rho(ix,iy,iz,1)+Eps)
      wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)-Cd*TempVM*TempW/(Rho(ix,iy,iz,1)+Eps)
      wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)-Cd*TempVM*TempW/(Rho(ix,iy,iz,1)+Eps)
    END DO    
  END DO
END SUBROUTINE WindDragCompute
END MODULE WindFarm_Mod
