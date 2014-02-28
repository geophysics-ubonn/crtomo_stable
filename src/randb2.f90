subroutine randb2(my_a,my_b)

!!!$     Unterprogramm modifiziert die Matrix 'a' (Bandbreite 'mb') und den
!!!$     Konstantenvektor 'b' zur Beruecksichtigung der Dirichletschen Rand-
!!!$     bedingungen ('rwdanz' Randwerte 'rwd(rwMAX)' mit zugeh. Knotennummern
!!!$     'rwdnr(smax)').

!!!$     ( Vgl. Subroutine 'RBSTBNDN' in Schwarz (1991) )

!!!$     Andreas Kemna                                            12-Feb-1993
!!!$     Letzte Aenderung   15-Jul-2007

!!!$.....................................................................
use alloci, only: prec
  USE elemmod , ONLY : sanz, mb
  USE randbmod

  IMPLICIT none

!!!$.....................................................................
!!!$     EIN-/AUSGABEPARAMETER:

  COMPLEX (prec),DIMENSION ((mb+1)*sanz) :: my_a
  COMPLEX (prec),DIMENSION (sanz)        :: my_b

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  INTEGER(KIND = 4)   ::     m1
  COMPLEX(prec)  ::     rwert

!!!$     Indexvariablen
  INTEGER(KIND = 4)   ::     ir,i,j,k,idk,ia,ki,ji

!!!$.....................................................................

  if (rwdanz.eq.0) return

  m1 = mb + 1

  do ir=1,rwdanz

     k      = rwdnr(ir)

     rwert  = rwd(ir)
     my_b(k)   = -rwert
     idk    = k*m1
     my_a(idk) = CMPLX(1d0)

     IF (k /= 1) THEN

        ia = max0(1,mb+2-k)

        do i=ia,mb
           j  = k+i-m1
           ki = idk+i-m1

           my_b(j)  = my_b(j) + rwert * my_a(ki)
           my_a(ki) = CMPLX(0d0)

        END do
     END IF

     if (k.eq.sanz) CYCLE

     ia = max0(1,k-sanz+m1)

     do i=ia,mb
        j  = k-i+m1
        ji = (j-1)*m1+i

        my_b(j)  = my_b(j) + rwert * my_a(ji)
        my_a(ji) = CMPLX(0d0)

     END do

  END do

  return
end subroutine randb2
