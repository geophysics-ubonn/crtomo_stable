      SUBROUTINE bmcm(kanal)
c     
c     Unterprogramm berechnet (einfache) Modell Kovarianz Matrix
c     MCM = (A^TC_d^-1A + C_m^-1)^-1
c     Fuer beliebige Triangulierung
c
c     Copyright Andreas Kemna
c     Andreas Kemna / Roland Martin                            02-Nov-2009
c     
c     Letzte Aenderung    RM                                   20-Feb-2010
c
c.........................................................................
      USE alloci
      
      IMPLICIT none

      INCLUDE 'parmax.fin'
      INCLUDE 'konv.fin'
      INCLUDE 'model.fin'
      INCLUDE 'elem.fin'
      INCLUDE 'err.fin'
!.....................................................................
!     PROGRAMMINTERNE PARAMETER:
!     Hilfsvariablen 
      INTEGER                                       :: i,kanal
      COMPLEX(KIND(0D0)),DIMENSION(:,:),ALLOCATABLE :: work
      COMPLEX(KIND(0D0)),DIMENSION(:),ALLOCATABLE   :: ipiv
      REAL(KIND(0D0)),DIMENSION(:),ALLOCATABLE      :: dig,dig2
      REAL(KIND(0D0))                               :: dig_min,dig_max
!.....................................................................

c$$$  invert A^TC_d^-1A + C_m^-1

      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

      ALLOCATE (work(manz,manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem WORK in bmcm'
         errnr = 97
         RETURN
      END IF
      ALLOCATE (ipiv(manz),STAT=errnr)
      IF (errnr/=0) THEN
         WRITE (*,'(/a/)')'Allocation problem IPIV in bmcmdc'
         errnr = 97
         RETURN
      END IF
      
      work = ata_reg ! work is replaced by PLU decomposition
      
cc$$$  building Right Hand Side (unit matrix)
      DO i=1,manz
         cov_m(i,i) = CMPLX(1.d0,1.d0)
      END DO

cc$$$  Solving Linear System Ax=B -> B=A^-1
      WRITE (*,'(a)')'Solving Ax=B (ZGESV)'
      CALL ZGESV(manz,manz,work,manz,ipiv,cov_m,manz,errnr)

      IF (errnr /= 0) THEN
         PRINT*,'Zeile::',cov_m(abs(errnr),:)
         PRINT*,'Spalte::',cov_m(:,abs(errnr))
         errnr = 108
         RETURN
      END IF

      work = MATMUL(cov_m,ata_reg)

c$$$      cov_m = ata_reg
c$$$      CALL gauss_cmplx(cov_m,manz,errnr)
c$$$      IF (errnr /= 0) THEN
c$$$         PRINT*,'Zeile::',cov_m(abs(errnr),:)
c$$$         PRINT*,'Spalte::',cov_m(:,abs(errnr))
c$$$         errnr = 108
c$$$         RETURN
c$$$      END IF

      ALLOCATE (dig(manz),dig2(manz)) !prepare to write out main diagonal
      DO i=1,manz
         dig(i) = DBLE(cov_m(i,i))
         dig2(i) = DBLE(work(i,i))
      END DO
      
      dig_min = MINVAL(dig)
      dig_max = MAXVAL(dig)
      
      WRITE (kanal,*)manz
      DO i=1,manz
         WRITE (kanal,*)LOG10(SQRT(ABS(dig(i)))),dig2(i)
      END DO

      WRITE (kanal,*)'Max/Min:',dig_max,'/',dig_min
      WRITE (*,*)'Max/Min:',dig_max,'/',dig_min

      CLOSE(kanal)

      DEALLOCATE (dig,dig2,work,ipiv)

      errnr = 0
 999  RETURN
      END
