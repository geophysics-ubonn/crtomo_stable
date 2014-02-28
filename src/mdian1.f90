SUBROUTINE MDIAN1(X2,N,XMED)
use alloci, only: prec
  IMPLICIT none
  INTEGER :: N,N2
  REAL(prec),DIMENSION(N) :: X,X2
  REAL(prec) :: XMED

  X = X2 ! save a copy
  CALL SORT(N,X)
  N2=N/2
  IF(2*N2.EQ.N)THEN
     XMED=0.5*(X(N2)+X(N2+1))
  ELSE
     XMED=X(N2+1)
  ENDIF

END SUBROUTINE MDIAN1
