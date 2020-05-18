MODULE Matrix_Mod
  USE Kind_Mod

  TYPE A22_T
    REAL(RealKind) :: a(1:2,1:2)=RESHAPE((/0.0d0,0.0d0,0.0d0,0.0d0/),(/2,2/))
  END TYPE A22_T

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE CopyScalar22 
  END INTERFACE
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE Add22 &
                    ,Add22_2
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE Sub22 & 
                    ,Sub22_2 & 
                    ,Neg22  
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE Mult22_22 &
                    ,Mult22_2 
  END INTERFACE
  INTERFACE OPERATOR(/)
    MODULE PROCEDURE Div22 &
                    ,Div2
  END INTERFACE

CONTAINS

SUBROUTINE CopyScalar22(a,Scalar)
  TYPE(A22_T), INTENT(OUT) :: a
  REAL(RealKind), INTENT(IN) :: Scalar

  a%a=Scalar

END SUBROUTINE CopyScalar22

FUNCTION Mult22_22(b,c)
  TYPE(A22_T) :: Mult22_22
  TYPE(A22_T), INTENT(IN) :: b,c

  TYPE(A22_T) :: Temp

  temp%a(1,1)=b%a(1,1)*c%a(1,1)+b%a(1,2)*c%a(2,1)
  temp%a(1,2)=b%a(1,1)*c%a(1,2)+b%a(1,2)*c%a(2,2)
  temp%a(2,1)=b%a(2,1)*c%a(1,1)+b%a(2,2)*c%a(2,1)
  temp%a(2,2)=b%a(2,1)*c%a(1,2)+b%a(2,2)*c%a(2,2)

  Mult22_22=Temp

END FUNCTION Mult22_22

FUNCTION Mult22_2(b,c)
  REAL(RealKind) :: Mult22_2(1:2)
  TYPE(A22_T), INTENT(IN) :: b
  REAL(RealKind), INTENT(IN) :: c(1:2)

  REAL(RealKind) :: Temp(1:2)

  temp(1)=b%a(1,1)*c(1)+b%a(1,2)*c(2)
  temp(2)=b%a(2,1)*c(1)+b%a(2,2)*c(2)

  Mult22_2=Temp
END FUNCTION Mult22_2

FUNCTION Mult2_2(b,c)
  REAL(RealKind) :: Mult2_2(1:2)
  REAL(RealKind), INTENT(IN) :: b(1:2)
  REAL(RealKind), INTENT(IN) :: c(1:2)

  Mult2_2=b*c(1)
END FUNCTION Mult2_2


FUNCTION Div22(b,c)
  TYPE(A22_T) :: Div22
  TYPE(A22_T), INTENT(IN) :: b,c

  TYPE(A22_T) :: Temp,Div
  REAL(RealKind) :: Det

  Det=c%a(1,1)*c%a(2,2)-c%a(1,2)*c%a(2,1)
  Div%a(1,1)=c%a(2,2)/Det
  Div%a(1,2)=-c%a(1,2)/Det
  Div%a(2,1)=-c%a(2,1)/Det
  Div%a(2,2)=c%a(1,1)/Det

  temp%a(1,1)=Div%a(1,1)*b%a(1,1)+Div%a(1,2)*b%a(2,1)
  temp%a(1,2)=Div%a(1,1)*b%a(1,2)+Div%a(1,2)*b%a(2,2)
  temp%a(2,1)=Div%a(2,1)*b%a(1,1)+Div%a(2,2)*b%a(2,1)
  temp%a(2,2)=Div%a(2,1)*b%a(1,2)+Div%a(2,2)*b%a(2,2)

  Div22=Temp

END FUNCTION Div22

FUNCTION Div2(b,c)
  REAL(RealKind) :: Div2(1:2)
  REAL(RealKind), INTENT(IN) :: b(1:2)
  TYPE(A22_T), INTENT(IN) ::c

  TYPE(A22_T) :: Div
  REAL(RealKind) :: Temp(1:2)
  REAL(RealKind) :: Det

  Det=c%a(1,1)*c%a(2,2)-c%a(1,2)*c%a(2,1)
  Div%a(1,1)=c%a(2,2)/Det
  Div%a(1,2)=-c%a(1,2)/Det
  Div%a(2,1)=-c%a(2,1)/Det
  Div%a(2,2)=c%a(1,1)/Det

  temp(1)=Div%a(1,1)*b(1)+Div%a(1,2)*b(2)
  temp(2)=Div%a(2,1)*b(1)+Div%a(2,2)*b(2)

  Div2=Temp

END FUNCTION Div2

FUNCTION Add22(b,c)
  TYPE(A22_T) :: Add22
  TYPE(A22_T), INTENT(IN) :: b
  TYPE(A22_T), INTENT(IN) :: c

  Add22%a=b%a+c%a

END FUNCTION Add22

FUNCTION Add22_2(b,c)
  TYPE(A22_T) :: Add22_2
  TYPE(A22_T), INTENT(IN) :: b
  REAL(RealKind), INTENT(IN) :: c(1:2)

  Add22_2%a(1,1)=b%a(1,1)+c(1)
  Add22_2%a(2,1)=b%a(2,1)+c(2)
  Add22_2%a(1,2)=b%a(1,2)
  Add22_2%a(2,2)=b%a(2,2)

END FUNCTION Add22_2

FUNCTION Sub22(b,c)
  TYPE(A22_T) :: Sub22
  TYPE(A22_T), INTENT(IN) :: b
  TYPE(A22_T), INTENT(IN) :: c

  Sub22%a=b%a-c%a

END FUNCTION Sub22

FUNCTION Sub22_2(b,c)
  TYPE(A22_T) :: Sub22_2
  TYPE(A22_T), INTENT(IN) :: b
  REAL(RealKind), INTENT(IN) :: c(1:2)

  Sub22_2%a(1,1)=b%a(1,1)-c(1)
  Sub22_2%a(2,1)=b%a(2,1)-c(2)
  Sub22_2%a(1,2)=b%a(1,2)
  Sub22_2%a(2,2)=b%a(2,2)

END FUNCTION Sub22_2


FUNCTION Neg22(b)
  TYPE(A22_T) :: Neg22
  TYPE(A22_T), INTENT(IN) :: b

  Neg22%a=-b%a

END FUNCTION Neg22

END MODULE Matrix_Mod
