MODULE Rates_Mod

  USE Kind_Mod 
  USE DataType_Mod

  IMPLICIT NONE

  REAL(RealKind), PARAMETER :: InvRefTemp=1.0d0/298.15d0
  REAL(RealKind), PARAMETER :: RefTemp=298.15d0
  REAL(RealKind) :: Chi
  REAL(RealKind) :: Dust=0.5d0

CONTAINS

SUBROUTINE PhoABCCompute(Rate,Constants)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)

  INTEGER :: i
  REAL(RealKind) :: ChiZ,yChiZ,EyChiZ

  IF (Chi<PiHalf) THEN
    ChiZ=Chi*Constants(3)
    IF (ChiZ<PiHalf) THEN
      yChiZ=Constants(2)*(1.0d0-1.0d0/COS(ChiZ))
      IF (yChiZ>-30.0d0) THEN
        EyChiZ=EXP(yChiZ)
      ELSE
        EyChiZ=9.357d-14
      END IF
    ELSE
      EyChiZ=9.357d-14
    END IF
    Rate%c=Dust*Constants(1)*EyChiz
  ELSE
    Rate%c=0.0d0
  END IF
END SUBROUTINE PhoABCCompute

SUBROUTINE PhoABCompute(Rate,Constants)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)

  INTEGER :: i

  IF (Chi<PiHalf) THEN
    Rate%c=Dust*Constants(1)*EXP(-Constants(2)/COS(Chi))
  ELSE
    Rate%c=0.0d0
  END IF
END SUBROUTINE PhoABCompute

SUBROUTINE PhoMCMCompute(Rate,Constants)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)

  REAL(RealKind) :: chiz,ychiz

!---  MCM version
  IF (chi < PiHalf) then
    chiz=EXP(-Constants(3)*(one/COS(chi)))
    ychiz=(COS(chi))**(Constants(2))
    Rate%c=Dust*Constants(1)*ychiz*chiz
  ELSE
    Rate%c=0.0d0
  END IF

END SUBROUTINE PhoMCMCompute

SUBROUTINE PhoMCM1Compute(Rate,Constants,Temp)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  REAL(RealKind) :: chiz,ychiz

!---  MCM version
  IF (chi < PiHalf) then
    chiz=EXP(-Constants(3)*(one/COS(chi)))
    ychiz=(COS(chi))**(Constants(2))
    WHERE (Mask)
      Rate%Vec3%c(:,:,:,1)=Dust*Constants(1)*ychiz*chiz*Constants(4)*Constants(5) &
            /(Constants(6)*EXP(-Constants(7)/Temp%c(:,:,:,1)))
    END WHERE        
    Rate=Rate%Vec3
  ELSE
    Rate%c=0.0d0
  END IF

END SUBROUTINE PhoMCM1Compute

SUBROUTINE PhoMCMSCompute(Rate,Constants)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)

  INTEGER :: i,ix,iy,iz
  REAL(RealKind) :: chiz,ychiz

!---  MCM version
  IF (chi < PiHalf) then
    chiz=EXP(-Constants(3)*(one/COS(chi)))
    ychiz=(COS(chi))**(Constants(2))
    Rate%c=Zero
    DO i=1,NumBoundCell
      ix=BoundCell(i)%ix
      iy=BoundCell(i)%iy
      iz=BoundCell(i)%iz
      Rate%c(ix,iy,iz,1)=Dust*Constants(1)*ychiz*chiz
    END DO  
  ELSE
    Rate%c=0.0d0
  END IF

END SUBROUTINE PhoMCMSCompute
SUBROUTINE ConstCompute(Rate,Constants)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)

  Rate%c=Constants(1)

END SUBROUTINE ConstCompute

SUBROUTINE ConstSCompute(Rate,Constants)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)

  INTEGER :: i,ix,iy,iz
  
  Rate%c=Zero
  DO i=1,NumBoundCell
    ix=BoundCell(i)%ix
    iy=BoundCell(i)%iy
    iz=BoundCell(i)%iz
    Rate%c(ix,iy,iz,1)=Constants(1)
  END DO  

END SUBROUTINE ConstSCompute

SUBROUTINE Temp1Compute(Rate,Constants,Temp)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind) :: Constants(:)
  TYPE(Vec4_T) :: Temp

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=Constants(1)*EXP(-Constants(2)/Temp%c(:,:,:,1))
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE Temp1Compute

SUBROUTINE Temp2Compute(Rate,Constants,Temp)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=Constants(1)*Temp%c(:,:,:,1)*Temp%c(:,:,:,1) &
                        *EXP(-Constants(2)/Temp%c(:,:,:,1))
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE Temp2Compute

SUBROUTINE Temp3Compute(Rate,Constants,Temp)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=Constants(1)*EXP(Constants(2) &
          *(1.0d0/Temp%c(:,:,:,1)-InvRefTemp))
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE Temp3Compute

SUBROUTINE TroeCompute(Rate,Constants,Temp,mAir)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=(mAir%c(:,:,:,1)*Constants(1)/Constants(3))*(Temp%c(:,:,:,1)/3.0d2)**(Constants(4)-Constants(2))
    Rate%Vec3%c(:,:,:,1)= &
               mAir%c(:,:,:,1)*Constants(1)*(Temp%c(:,:,:,1)/3.0d2)**(-Constants(2)) &
               /(1.0d0+Rate%Vec3%c(:,:,:,1)) &
               *0.6d0**(1.0d0/(1.0d0+LOG10(Rate%Vec3%c(:,:,:,1))**2))
  END WHERE
  Rate=Rate%Vec3
END SUBROUTINE TroeCompute

SUBROUTINE TroeFCompute(Rate,Constants,Temp,mAir) ! Barthel

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
  Rate%Vec3%c(:,:,:,1)=(mAir%c(:,:,:,1)*Constants(1)/Constants(3))*(Temp%c(:,:,:,1)/3.0d2)**(Constants(4)-Constants(2))
  Rate%Vec3%c(:,:,:,1)= &
             mAir%c(:,:,:,1)*Constants(1)*(Temp%c(:,:,:,1)/3.0d2)**(-Constants(2)) &
             /(1.0d0+Rate%Vec3%c(:,:,:,1)) &
             *Constants(5)**(1.0d0/(1.0d0+LOG10(Rate%Vec3%c(:,:,:,1))**2))
  END WHERE
  Rate=Rate%Vec3
END SUBROUTINE TroeFCompute

SUBROUTINE TroeEqCompute(Rate,Constants,Temp,mAir)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=(mAir%c(:,:,:,1)*Constants(1)/Constants(3))*(Temp%c(:,:,:,1)/3.0d2)**(Constants(4)-Constants(2))
    Rate%Vec3%c(:,:,:,1)= &
               mAir%c(:,:,:,1)*Constants(1)*(Temp%c(:,:,:,1)/3.0d2)**(-Constants(2)) &
               /(1.0d0+Rate%Vec3%c(:,:,:,1)) &
               *0.6d0**(1.0d0/(1.0d0+LOG10(Rate%Vec3%c(:,:,:,1))**2)) &
               /(Constants(5)*EXP(Constants(6)/Temp%c(:,:,:,1)))  
  END WHERE
  Rate=Rate%Vec3
END SUBROUTINE TroeEqCompute

SUBROUTINE TroeEqfCompute(Rate,Constants,Temp,mAir) ! Barthel

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=(mAir%c(:,:,:,1)*Constants(1)/Constants(3))*(Temp%c(:,:,:,1)/3.0d2)**(Constants(4)-Constants(2))
    Rate%Vec3%c(:,:,:,1)= &
               mAir%c(:,:,:,1)*Constants(1)*(Temp%c(:,:,:,1)/3.0d2)**(-Constants(2)) &
               /(1.0d0+Rate%Vec3%c(:,:,:,1)) &
               *Constants(7)**(1.0d0/(1.0d0+LOG10(Rate%Vec3%c(:,:,:,1))**2)) &
               /(Constants(5)*EXP(Constants(6)/Temp%c(:,:,:,1)))  
  END WHERE
  Rate=Rate%Vec3
END SUBROUTINE TroeEqfCompute

SUBROUTINE TroeXPCompute(Rate,Constants,Temp,mAir) ! Barthel

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=mAir%c(:,:,:,1)*Constants(1)*exp(-Constants(2)/Temp%c(:,:,:,1))
    Rate%Vec3%c(:,:,:,1)=Rate%Vec3%c(:,:,:,1)/(1.0d0+Rate%Vec3%c(:,:,:,1)/(Constants(3)*exp(-Constants(4)/Temp%c(:,:,:,1)))) &
  *Constants(5)**(1.0d0/(1.0d0+LOG10(Rate%Vec3%c(:,:,:,1)/(Constants(3)*exp(-Constants(4)/Temp%c(:,:,:,1)))))**2.0d0)
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE TroeXPCompute

SUBROUTINE Spec1Compute(Rate,Constants,mAir)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=Constants(1)*(1.0d0+mAir%c(:,:,:,1)*Constants(2))
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE Spec1Compute

SUBROUTINE Spec2Compute(Rate,Constants,Temp,mAir)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir


  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=mAir%c(:,:,:,1)*Constants(1)*(Temp%c(:,:,:,1)/300.0d0)**Constants(2)
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE Spec2Compute

SUBROUTINE Spec3Compute(Rate,Constants,Temp,mAir)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)= &
                Constants(1)*EXP(Constants(2)/Temp%c(:,:,:,1)) &
               +1.0d0/(EXP(-Constants(4)/Temp%c(:,:,:,1))/Constants(3) &
                      +EXP(-Constants(6)/Temp%c(:,:,:,1))/(Constants(5)*mAir%c(:,:,:,1)) &
                      )
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE Spec3Compute

SUBROUTINE Spec4Compute(Rate,Constants,Temp,mAir)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=Constants(1)*EXP(Constants(2)/Temp%c(:,:,:,1)) &
          +mAir%c(:,:,:,1)*Constants(3)*EXP(Constants(4)/Temp%c(:,:,:,1))
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE Spec4Compute

SUBROUTINE T1H2OCompute(Rate,Constants,Temp,rHum)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: rHum

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=rHum%c(:,:,:,1)*Constants(1)*EXP(-Constants(2)/Temp%c(:,:,:,1)) 
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE T1H2OCompute

SUBROUTINE S4H2OCompute(Rate,Constants,Temp,mAir,rHum)

  TYPE(Vec4_T) :: Rate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp
  TYPE(Vec4_T) :: mAir
  TYPE(Vec4_T) :: rHum

  WHERE (Mask)
    Rate%Vec3%c(:,:,:,1)=Constants(1)*EXP(Constants(2)/Temp%c(:,:,:,1)) &
           +mAir%c(:,:,:,1)*Constants(3)*EXP(Constants(4)/Temp%c(:,:,:,1))
  END WHERE
  Rate=Rate%Vec3

END SUBROUTINE S4H2OCompute

SUBROUTINE DConstCompute(EquiRate,BackRate,Constants)

  TYPE(Vec4_T) :: EquiRate
  TYPE(Vec4_T) :: BackRate
  REAL(RealKind), POINTER :: Constants(:)

  EquiRate%c=Constants(1)
  BackRate%c=Constants(2)

END SUBROUTINE DConstCompute

SUBROUTINE DTempCompute(EquiRate,BackRate,Constants,Temp)

  TYPE(Vec4_T) :: EquiRate
  TYPE(Vec4_T) :: BackRate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  EquiRate%Vec3%c=Constants(1)*EXP(Constants(2)*(1.0d0/Temp%c-InvRefTemp))
  Equirate=EquiRate%Vec3
  BackRate%c=Constants(3)

END SUBROUTINE DTempCompute

SUBROUTINE DTemp2Compute(EquiRate,BackRate,Constants,Temp)

  TYPE(Vec4_T) :: EquiRate
  TYPE(Vec4_T) :: BackRate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  EquiRate%Vec3%c=Constants(1)*EXP(Constants(2)*(1.0d0/Temp%c-InvRefTemp))
  Equirate=EquiRate%Vec3
  BackRate%Vec3%c=Constants(3)*EXP(Constants(4)*(1.0d0/Temp%c-InvRefTemp))
  Backrate=BackRate%Vec3

END SUBROUTINE DTemp2Compute

SUBROUTINE DTemp3Compute(EquiRate,BackRate,Constants,Temp)

  TYPE(Vec4_T) :: EquiRate
  TYPE(Vec4_T) :: BackRate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  EquiRate%Vec3%c=Constants(1) &
                  *(Temp%c/RefTemp)**Constants(2) &
                  *EXP(Constants(3)*(One/Temp%c-One/RefTemp)) 
  Equirate=EquiRate%Vec3
  BackRate%Vec3%c=1.d0
  Backrate=BackRate%Vec3

END SUBROUTINE DTemp3Compute

SUBROUTINE MeskhidzeCompute(EquiRate,BackRate,Constants,Temp)

  TYPE(Vec4_T) :: EquiRate
  TYPE(Vec4_T) :: BackRate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  EquiRate%Vec3%c=Constants(1) &
                  *(Temp%c/RefTemp)**Constants(2) &
                  *EXP(Constants(3)*(One/Temp%c-One/RefTemp)) 
  Equirate=EquiRate%Vec3
  BackRate%Vec3%c=Constants(4)*EXP(Constants(5)*(One/Temp%c-One/RefTemp))*Constants(7)!*(ActivityHp)**m
  Backrate=BackRate%Vec3

END SUBROUTINE MeskhidzeCompute

SUBROUTINE DTemp4Compute(EquiRate,BackRate,Constants,Temp)

  TYPE(Vec4_T) :: EquiRate
  TYPE(Vec4_T) :: BackRate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  EquiRate%Vec3%c=Constants(1) &
                  *EXP(Constants(2)*(Temp%c*InvRefTemp-One) &
                      +Constants(3)*(One+LOG(Temp%c*InvRefTemp)-Temp%c*InvRefTemp))
  Equirate=EquiRate%Vec3
  BackRate%Vec3%c=1.d0  ! OSSI
  Backrate=BackRate%Vec3

END SUBROUTINE DTemp4Compute


SUBROUTINE DTemp5Compute(EquiRate,BackRate,Constants,Temp)

  TYPE(Vec4_T) :: EquiRate
  TYPE(Vec4_T) :: BackRate
  REAL(RealKind), POINTER :: Constants(:)
  TYPE(Vec4_T) :: Temp

  EquiRate%Vec3%c=Constants(1) &
                 *(Temp%c*InvRefTemp)**Constants(2) &
                 *EXP(Constants(3)*(1.0d0/Temp%c-InvRefTemp))
  Equirate=EquiRate%Vec3
  BackRate%Vec3%c=1.0d0  ! OSSI
  Backrate=BackRate%Vec3

END SUBROUTINE DTemp5Compute
END MODULE Rates_Mod
