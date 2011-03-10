subroutine blam0()

!!!$     Unterprogramm zum Bestimmen des Start-Regularisierungsparameters.

!!!$     Andreas Kemna                                            20-Feb-1997
!!!$     Letzte Aenderung   07-Mar-2003

!!!$.....................................................................

  USE alloci
  USE femmod
  USE datmod
  USE invmod
  USE modelmod
  USE konvmod

  IMPLICIT none


!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Hilfsvariablen
  COMPLEX (KIND(0D0)) ::  cdum
  REAL (KIND(0D0))    ::   dum

!!!$     Indexvariablen
  INTEGER (KIND = 4)  ::  i,j,k,i_count

!!!$.....................................................................

!!!$     Start-Regularisierungsparameter bestimmen
  IF (nz<0) THEN
     IF (nz<-1) lammax = -REAL(nz)
     IF (nz==-1) lammax = MAX(REAL(manz),REAL(nanz))
     WRITE (*,'(t5,a,F12.1)')'taking easy lam_0 ',lammax
     RETURN
  END IF

  lammax = 0d0
  i_count = 0
  !$OMP PARALLEL DEFAULT (SHARED) PRIVATE(i,k,dum)
  if (ldc) then
     !$OMP DO
     do j=1,manz
        i_count = i_count + 1
        IF (lverb) write(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
             'blam0/ ',REAL( i_count * (100./manz)),'%'
        dum = 0d0

        do i=1,nanz
           do k=1,manz
              dum = dum + sensdc(i,j) * sensdc(i,k) * &
                   wmatd(i)*dble(wdfak(i))
           end do
        end do

        lammax = lammax + dabs(dum)
     end do
     !$OMP END DO
  else if (lip) then
     !$OMP DO
     do j=1,manz
        i_count = i_count + 1
        IF (lverb) write(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
             'blam0/ ',REAL( i_count * (100./manz)),'%'
        dum = 0d0

        do i=1,nanz
           do k=1,manz
              dum = dum + dble(sens(i,j)) * dble(sens(i,k)) * &
                   wmatd(i)*dble(wdfak(i))
           end do
        end do

        lammax = lammax + dabs(dum)
     end do
     !$OMP END DO
  else
     !$OMP DO
     do j=1,manz
        i_count = i_count + 1
        IF (lverb) write(*,'(a,t70,F6.2,A)',advance='no')ACHAR(13)//&
             'blam0/ ',REAL( i_count * (100./manz)),'%'
        cdum = dcmplx(0d0)

        do i=1,nanz
           do k=1,manz
              cdum = cdum + dconjg(sens(i,j)) * sens(i,k) * &
                   dcmplx(wmatd(i)*dble(wdfak(i)))
           end do
        end do

        lammax = lammax + cdabs(cdum)
     end do

  end if
  !$OMP END PARALLEL

  lammax = lammax/dble(manz)

  lammax = lammax * 2d0/(alfx+alfz)
!!!$     ak Default
  lammax = lammax * 5d0

  WRITE (*,'(t5,a,F12.1)')'found lam_0 ',lammax
!!!$     ak Synthetic Example (JoH)
!!!$     ak        lammax = lammax * 1d1

!!!$     ak MinFrac
!!!$     ak        lammax = lammax * 5d1

!!!$     ak Test
!!!$     ak        lammax = lammax * 1d1

!!!$     ak AAC
!!!$     ak        lammax = lammax * 5d0
  RETURN
end subroutine blam0
