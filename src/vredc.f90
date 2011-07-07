subroutine vredc()

!!!$     Fuehrt das Vorwaerts- und Rueckwaertseinsetzen mit der Cholesky-Links-
!!!$     dreiecksmatrix aus;
!!!$     'bdc' bleibt unveraendert, 'pot' ist Loesungsvektor.

!!!$     ( Vgl. Subroutine 'VRBNDN' in Schwarz (1991) )

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   14-Nov-1997

!!!$.....................................................................

  USE alloci
  USE femmod
  USE elemmod
  USE errmod

  IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE :: potdc
  INTEGER (KIND=4)  :: idi,i0
  INTEGER (KIND=4)  :: m1,jlow
  REAL(KIND(0D0))   ::   s

!!!$     Indexvariablen
  INTEGER (KIND=4)  ::  i,j

!!!$.....................................................................

  ALLOCATE (potdc(sanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation potdc failed'
     errnr = 94
     RETURN
  END IF


  m1 = mb+1

  do i=1,sanz
     idi  = i*m1
     s    = bdc(i)
     i0   = idi-i
     jlow = max0(1,i-mb)

     do j=jlow,i-1
        s = s - adc(i0+j)*potdc(j)
     END do

     potdc(i) = s / adc(idi)
  END do

  do i=sanz,1,-1
     potdc(i) = - potdc(i) / adc(idi)

     jlow = max0(1,i-mb)
     i0   = idi-i

     do j=jlow,i-1
        potdc(j) = potdc(j) + adc(i0+j)*potdc(i)
     END do

     idi = idi-m1
  END do

  do i=1,sanz
     pot(i) = dcmplx(potdc(i))
  end do

  DEALLOCATE (potdc)
end subroutine vredc
