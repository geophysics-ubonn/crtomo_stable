SUBROUTINE relem(kanal,datei)

!!!$     Unterprogramm zum Einlesen der FEM-Parameter aus 'datei'.

!!!$     Andreas Kemna                                            11-Oct-1993
!!!$     Letzte Aenderung   24-Oct-1996

!!!$.....................................................................

  USE elemmod
  USE errmod
  USE konvmod

  IMPLICIT NONE


!!!$.....................................................................

!!!$     EIN-/AUSGABEPARAMETER:

!!!$     Kanalnummer
  INTEGER (KIND = 4) ::    kanal

!!!$     Datei
  CHARACTER (80)     ::    datei

!!!$.....................................................................

!!!$     PROGRAMMINTERNE PARAMETER:

!!!$     Indexvariablen
  INTEGER (KIND =4)  ::    i,j,k

!!!$     Hilfsvariable
  INTEGER (KIND =4)  ::    idum,ifln,iflnr
  LOGICAL            ::    my_check, failed

!!!$ To check for border to element connection (rnr)
  INTEGER            :: ik1,ik2,jk1,jk2,ic,l
!!!$ NEW rnr
  !  INTEGER (KIND = 4),ALLOCATABLE,DIMENSION(:) :: my_rnr

!!!$.....................................................................

!!!$     'datei' oeffnen
  fetxt = datei

  errnr = 1
  OPEN(kanal,file=TRIM(fetxt),status='old',err=999)

  errnr = 3

!!!$     Anzahl der Knoten (bzw. Knotenvariablen), Anzahl der Elementtypen
!!!$     sowie Bandbreite der Gesamtsteifigkeitsmatrix einlesen
  READ(kanal,*,END=1001,err=1000) sanz,typanz,mb

!!$ now get some memory for the fields..
!!$ first the sanz fields
  ALLOCATE (sx(sanz),sy(sanz),snr(sanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation sx failed'
     errnr = 94
     GOTO 999
  END IF

  ALLOCATE (typ(typanz),nelanz(typanz),selanz(typanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation selanz failed'
     errnr = 94
     GOTO 999
  END IF

!!!$     Elementtypen, Anzahl der Elemente eines bestimmten Typs sowie
!!!$     Anzahl der Knoten in einem Elementtyp einlesen
  READ(kanal,*,END=1001,err=1000)(typ(i),nelanz(i),selanz(i),i=1,typanz)

!!$ set number of node points for regular elements
  smaxs = MAXVAL(selanz)

!!!$     Anzahl der Elemente (ohne Randelemente) und Anzahl der Randelemente
!!!$     bestimmen
  relanz = 0
  elanz  = 0

  my_check = .FALSE.

  DO i=1,typanz
     IF (typ(i).GT.10) THEN
        relanz = relanz + nelanz(i)
     ELSE
        elanz  = elanz  + nelanz(i)
     END IF
     my_check = my_check .OR. (typ(i) == 11)
  END DO

!!!$ if all no flow boundaries, we do not have to search for a 
!!!$ average sy top...
!!$  lsytop = .NOT. my_check

!!$ get memory for the element integer field      
  ALLOCATE (nrel(elanz+relanz,smaxs),rnr(relanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation nrel,rnr failed'
     errnr = 94
     GOTO 999
  END IF
!!$ get memory for the regular element midpoint coordinates
  ALLOCATE (espx(elanz),espy(elanz),stat=errnr)
  IF (errnr /= 0) THEN
     fetxt = 'Error memory allocation espx failed'
     errnr = 94
     GOTO 999
  END IF
  espx = 0.;espy = 0.
!!!$     Zeiger auf Koordinaten, x-Koordinaten sowie y-Koordinaten der Knoten
!!!$     einlesen
  READ(kanal,*,END=1001,err=1000) (snr(i),sx(i),sy(i),i=1,sanz)
!!!$     Knotennummern der Elemente einlesen
  idum = 0;ifln = 0;iflnr = 0
  DO i=1,typanz
     DO j=1,nelanz(i)
        READ(kanal,*,END=1001,err=1000)(nrel(idum+j,k),k=1,selanz(i))

        IF (typ(i) < 10) THEN ! set midpoints

           ifln = ifln + 1

           DO k = 1,selanz(i)
              espx(ifln) = espx(ifln) + sx(snr(nrel(idum+j,k)))
              espy(ifln) = espy(ifln) + sy(snr(nrel(idum+j,k)))
           END DO

           espx(ifln) = espx(ifln) / selanz(i)
           espy(ifln) = espy(ifln) / selanz(i)

        END IF
     END DO
     idum = idum + nelanz(i)

  END DO

  failed = .FALSE.
  ik1 = nrel(elanz + 1,1)
  ik2 = nrel(elanz + 1,2)
  ic = 1
  DO i=2,relanz
!!!$ consistency check of node ordering, assume ascending?

     jk1 = nrel(elanz + i,1)
     jk2 = nrel(elanz + i,2)
     IF (ik2 /= jk1) THEN
        failed = .TRUE.
        ic = ic + 1
     END IF
  END DO
  IF (ic == relanz) THEN
     WRITE(*,'(/a/)',ADVANCE='no')'+++ consistency check '//&
          'border-line-elements failed -> swapping elements'
  END IF

  DO i=1,relanz

     ik1 = nrel(elanz + i,2)
     nrel(elanz + i,2) = nrel(elanz + i,1)
     nrel(elanz + i,1) = ik1
  END DO

!!!$     Zeiger auf Werte der Randelemente einlesen not anymore
!!$!  READ(kanal,*,END=1001,err=1000) (rnr(i),i=1,relanz)

!!!$ >> RM
!!!$ IF THIS IS NOT a pointer to the
!!!$ REGULAR ELEMENT adjacent to the border line
!!!$ THIS MAY CAUSE THE MIXED BOUNDARY TO BLOW UP
!!$  ic = 0
!!$  DO i=1,relanz
!!$     IF (rnr(i) > elanz.OR.rnr(i)<1) ic = ic + 1
!!$  END DO

!!$  IF (ic > 0) THEN
!!!$ CHECK where the BORDER ELE begin in nrel

  ic = 0
  DO i=1,relanz
!!!$ define the node points of the border-line
     ik1 = nrel(elanz + i,1)
     ik2 = nrel(elanz + i,2)

!!!$ now search for the corresponding element
!!!$ how to do so?
!!!$ suppose we have a ordered input file, where the
!!!$ big elements come first, than we might have something like
     DO k = 1,elanz
        jk1 = 0;jk2 = 0
        DO  l = 1, selanz(1) ! TODO : fix this for multi FE
           IF (nrel(k,l) == ik1) jk1 = 1
           IF (nrel(k,l) == ik2) jk2 = 1
        END DO
        IF (jk1 == 1 .AND. jk2 == 1.AND.rnr(i) /= k) THEN
           ic = ic + 1
           WRITE (*,'(a,I7,1x,a,I7)',ADVANCE='no')ACHAR(13)//&
                'wrong element to border connection: ',ic,'/',relanz
           rnr(i) = k
           failed = .TRUE.
        END IF
     END DO

  END DO

  ic = 0
  DO i=1,relanz
     IF (rnr(i) > elanz.OR.rnr(i) < 1) ic = ic + 1
  END DO

  IF (ic > 0) THEN
     PRINT*,'Kann nich sein'
     fetxt = 'border to element context is wrong'
     errnr = 111
     RETURN
  ELSE
     PRINT*,'done!'
  END IF

  !  DEALLOCATE (my_rnr)

!!!$ << RM
!!!$     'datei' schliessen
  CLOSE(kanal)


!!!$ 
  IF (failed) THEN

101  FORMAT (I9)
102  FORMAT (2(I9,1X))
103  FORMAT (3(I9,1X))
104  FORMAT (4(I9,1X))
105  FORMAT (I9,2(1X,G12.4))

     fetxt = TRIM(datei)//'_new'
     PRINT*,'+++ WRITING OUT IMPROVED GRID --> Writing',TRIM(fetxt)
     errnr = 1
     OPEN(kanal,file=TRIM(fetxt),status='replace')
     
     errnr = 3

!!!$     Anzahl der Knoten (bzw. Knotenvariablen), Anzahl der Elementtypen
!!!$     sowie Bandbreite der Gesamtsteifigkeitsmatrix einlesen
     WRITE(kanal,103) sanz,typanz,mb
     WRITE(kanal,103) (typ(i),nelanz(i),selanz(i),i=1,typanz)
     DO i=1,sanz
        WRITE(kanal,105) snr(i),sx(i),sy(i)
     END DO
!     STOP
!!!$     Knotennummern der Elemente einlesen
     idum = 0;ifln = 0;iflnr = 0
     DO i=1,typanz
        DO j=1,nelanz(i)
           IF (typ(i) == 8)  THEN
              WRITE(kanal,104)(nrel(idum+j,k),k=1,selanz(i))
           ELSE IF (typ(i) == 4)  THEN
              WRITE(kanal,103)(nrel(idum+j,k),k=1,selanz(i))
           ELSE
              WRITE(kanal,102)(nrel(idum+j,k),k=1,selanz(i))
           END IF
        END DO

        idum = idum + nelanz(i)
        
     END DO
     
     WRITE(kanal,101) (rnr(i),i=1,relanz)

     CLOSE (kanal)

  END IF


  errnr = 0
  RETURN

!!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!!!$     Fehlermeldungen

999 RETURN

1000 CLOSE(kanal)
  RETURN

1001 CLOSE(kanal)
  errnr = 2
  RETURN

END SUBROUTINE relem
