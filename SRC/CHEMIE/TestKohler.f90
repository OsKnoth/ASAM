PROGRAM TestKohler

  REAL(8) :: x
  REAL(8) :: FacX
  REAL(8) :: A,B,C


  x=1.d-9
  FacX=1.014d0
  A=1.0d-5
  B=1.0d0
  C=1.0d-2
  DO i=1,1000
    WRITE(10,*) x,x/(A+B*x)*exp(C/x**(1.0d0/3.0d0)) &
                 ,(A/((A+B*x)*(A+B*x))+x/(A+B*x)*(-1.0d0/3.0d0*C*x**(-4.0d0/3.0d0)))&
                 ,(A/((A+B*x)*(A+B*x))+x/(A+B*x)*(-1.0d0/3.0d0*C*x**(-4.0d0/3.0d0)))*exp(C/x**(1.0d0/3.0d0)) &
                 ,1.0d0
    x=x*FacX
  END DO    


END PROGRAM TestKohler
