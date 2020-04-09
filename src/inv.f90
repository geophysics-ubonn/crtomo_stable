! Copyright © 1990-2020 Andreas Kemna <kemna@geo.uni-bonn.de>
! Copyright © 2008-2020 CRTomo development team (see AUTHORS file)
!
!
! Permission is hereby granted, free of charge, to any person obtaining a copy of
! this software and associated documentation files (the “Software”), to deal in
! the Software without restriction, including without limitation the rights to
! use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
! of the Software, and to permit persons to whom the Software is furnished to do
! so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.

!> \file inv.f90
!! \brief Main inversion algorithm
!> \author Andreas Kemna, Roland Martin, Johannes Kenkel
!> \date 11/4/2014
!> \date {original version} 1990
!> @param aggr information
PROGRAM inv

    !>   Hauptprogramm zur Complex-Resistivity-2.5D-Inversion.

    !!!$   Belegte Kanaele:
    !!!$   9 - error.dat   -> fperr
    !!!$   10 - run.ctr    -> fprun
    !!!$   11 - in-/output -> kanal
    !!!$   12 - crtomo.cfg -> fpcfg
    !!!$   13 - inv.ctr    -> fpinv
    !!!$   14 - cjg.ctr    -> fpcjg
    !!!$   15 - eps.ctr    -> fpeps
    !!!$
    !!!$   Andreas Kemna                                        02-May-1995
    !!!$   Letzte Aenderung                                     Jul-2010

    !!!$.....................................................................

    USE alloci
    USE tic_toc
    USE femmod
    USE datmod
    USE invmod
    USE cjgmod
    USE sigmamod
    USE electrmod
    USE modelmod
    USE elemmod
    USE wavenmod
    USE randbmod
    USE konvmod
    USE errmod
    USE pathmod
    USE bsmatm_mod
    USE bmcm_mod
    USE brough_mod
    USE invhpmod
    USE omp_lib
    USE ompmod
    USE get_ver
    !!!$   USE portlib

    IMPLICIT NONE

    CHARACTER(256)         :: ftext
    INTEGER                :: c1,i,count,mythreads,maxthreads
    REAL(KIND(0D0))        :: lamalt,bdalt
    LOGICAL                :: converged,l_bsmat
    INTEGER,PARAMETER      :: clrln_len=50
    CHARACTER(clrln_len)   :: clrln ! clear line

    INTEGER :: getpid,pid,myerr
    !!!$:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    !!!$   SETUP UND INPUT
    !!!$   'crtomo.cfg' oeffnen
    errnr = 1
    !!!$   Kanal nummern belegen.. die variablen sind global!
    fperr = 9
    fprun = 10
    kanal = 11
    fpcfg = 12
    fpinv = 13
    fpcjg = 14
    fpeps = 15

    DO i=1,clrln_len-1
        clrln(i:i+1) = ' ' ! fill clear line CHARACTER array
    END DO

    pid = getpid()
    fetxt = 'crtomo.pid'
    PRINT*,'######### CRTomo ############'
    WRITE(*,"(a)") 'Licence:', &
    'Copyright © 1990-2020 Andreas Kemna <kemna@geo.uni-bonn.de>', &
    'Copyright © 2008-2020 CRTomo development team (see AUTHORS file)',&
    'Permission is hereby granted, free of charge, to any person obtaining a copy of',&
    'this software and associated documentation files (the “Software”), to deal in',&
    'the Software without restriction, including without limitation the rights to',&
    'use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies',&
    'of the Software, and to permit persons to whom the Software is furnished to do',&
    'so, subject to the following conditions:',&
    '',&
    'The above copyright notice and this permission notice shall be included in all',&
    'copies or substantial portions of the Software.',&
    '',&
    'THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR',&
    'IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,',&
    'FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE',&
    'AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER',&
    'LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,',&
    'OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE',&
    'SOFTWARE.'

    PRINT*, ''
    PRINT*, ''

    PRINT*,'Process_ID ::',pid
    OPEN (fprun,FILE=TRIM(fetxt),STATUS='replace',err=999)
    WRITE (fprun,*)pid
    CLOSE (fprun)
    maxthreads = OMP_GET_MAX_THREADS()
    WRITE(6,"(a, i5)") " OpenMP max threads: ", maxthreads

    PRINT*,'Version control of binary:'
    CALL get_git_ver

    fetxt = 'crtomo.cfg'
    OPEN(fpcfg,file=TRIM(fetxt),status='old',err=999)

    lagain=.TRUE. ! is set afterwards by user input file to false

    !!!$   Allgemeine Parameter setzen
    DO WHILE (lagain) ! this loop exits after all files are processed
        errnr2 = 0
        !!!$   Benoetigte Variablen einlesen
        CALL rall(kanal,delem,delectr,dstrom,drandb,&
            dsigma,dvolt,dsens,dstart,dd0,dm0,dfm0,lagain)
        !!!$   diff+<
        !!!$   diff-     1            dsigma,dvolt,dsens,dstart,lsens,lagain)
        !!!$   diff+>
        IF (errnr.NE.0) GOTO 999

        NTHREADS = maxthreads
        !!!$ now that we know nf and kwnanz, we can adjust the OMP environment..
        IF (maxthreads > 2) THEN ! single or double processor machines don't need scheduling..
            mythreads = MAX(kwnanz,2)
            PRINT*,'Rescheduling..'
            IF (mythreads <= maxthreads) THEN ! best case,
                !!!$ the number of processors is greater or equal the assumed
                !!!$ workload
                PRINT*,'perfect match'
            ELSE
                !!!$ is smaller than the minimum workload.. now we have to devide a bit..
                PRINT*,'less nodes than wavenumbers'
                DO i = 1, INT(kwnanz/2)
                    mythreads = INT(kwnanz / i) + 1
                    IF (mythreads < maxthreads) EXIT
                END DO
            END IF
            NTHREADS = mythreads
        END IF
        CALL OMP_SET_NUM_THREADS ( NTHREADS )
        ! recheck ..
        i = OMP_GET_MAX_THREADS()
        WRITE(6,'(2(a, i3),a)') " OpenMP threads: ",i,'(',maxthreads,')'
        !!!$

        !!!$   Element- und Randelementbeitraege sowie ggf. Konfigurationsfaktoren
        !!!$   zur Berechnung der gemischten Randbedingung bestimmen
        CALL precal()

        IF (errnr.NE.0) GOTO 999

        IF (.NOT.lbeta) THEN
            lsr = .FALSE.

            !!!$   Ggf. Fehlermeldungen
            IF (.NOT.lsink) THEN
                fetxt = 'no mixed boundary specify sink node'
                errnr = 102
            END IF
            !!!$ RM this is handeled in rall..
            !!$        if (swrtr.eq.0.and..not.lrandb2) then
            !!$           fetxt = ' '
            !!$           errnr = 103
            !!$        end if

            IF (errnr.NE.0) GOTO 999
        ELSE

            !!!$   Ggf. Fehlermeldungen
            IF (swrtr.EQ.0) THEN
                fetxt = ' '
                errnr = 106
            END IF
        END IF

        !!!$   getting dynamic memory
        errnr = 94
        !!!$ physical model
        fetxt = 'allocation problem sigma'
        ALLOCATE (sigma(elanz),STAT=myerr)
        IF (myerr /= 0) GOTO 999
        fetxt = 'allocation problem sigma2'
        ALLOCATE (sigma2(elanz),STAT=myerr)
        IF (myerr /= 0) GOTO 999
        !!!$  model parameters
        fetxt = 'allocation problem par'
        ALLOCATE (par(manz),STAT=myerr)
        IF (myerr /= 0) GOTO 999
        fetxt = 'allocation problem dpar'
        ALLOCATE (dpar(manz),STAT=myerr)
        IF (myerr /= 0) GOTO 999
        fetxt = 'allocation problem dpar2'
        ALLOCATE (dpar2(manz),STAT=myerr)
        IF (myerr /= 0) GOTO 999
        fetxt = 'allocation problem pot'
        ALLOCATE(pot(sanz),STAT=myerr)
        IF (myerr /= 0) GOTO 999
        fetxt = 'allocation problem pota'
        ALLOCATE (pota(sanz),STAT=myerr)
        IF (myerr /= 0) GOTO 999
        !!!$ now the big array are coming..
        fetxt = 'allocation problem fak'
        ALLOCATE (fak(sanz),STAT=myerr) ! fak for modeling
        IF (myerr /= 0) GOTO 999
        IF (ldc) THEN
            fetxt = 'allocation problem sensdc'
            ALLOCATE (sensdc(nanz,manz),STAT=myerr)
            IF (myerr /= 0) GOTO 999
            fetxt = 'allocation problem kpotdc'
            ALLOCATE (kpotdc(sanz,eanz,kwnanz),STAT=myerr)
        ELSE
            fetxt = 'allocation problem sens'
            ALLOCATE (sens(nanz,manz),STAT=myerr)
            IF (myerr /= 0) GOTO 999
            fetxt = 'allocation problem kpot'
            ALLOCATE (kpot(sanz,eanz,kwnanz),STAT=myerr)
        END IF
        IF (myerr /= 0) GOTO 999

        !!!$ get CG data storage of residuums and bvec, which is global
        CALL con_cjgmod (1,fetxt,myerr)
        IF (myerr /= 0) GOTO 999

        !!!$ >> RM ref model regu
        !!!$ assign memory to global variables
        fetxt = 'allocation problem reference model'
        ALLOCATE (w_ref_re(manz),w_ref_im(manz),m_ref(manz),STAT=myerr)
        w_ref_re = 0d0;m_ref = DCMPLX(0d0);w_ref_im = 0d0
        IF (myerr /= 0) GOTO 999
        !!!$ << RM ref model regu

        !!!!$ INITIALIZE
        !!!$   Startparameter setzen
        bdalt = 1d0
        it     = 0;itr    = 0
        rmsalt = 0d0; lamalt = 1d0; bdpar = 1d0
        !     IF (lamnull_cri > 0d0) llamalt = lamnull_cri
        IF (BTEST(llamf,0)) lamalt = lamfix
        betrms = 0d0; pharms = 0d0
        lsetup = .TRUE.; lsetip = .FALSE.; lfpi    = .FALSE.
        llam   = .FALSE.; ldlami = .TRUE.; lstep  = .FALSE.
        lfstep = .FALSE.; l_bsmat = .TRUE.
        step   = 1d0; stpalt = 1d0; alam   = 0d0

        !!!$   Kontrolldateien oeffnen
        errnr = 1
        !!!$ OPEN CONTRL FILES
        fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
        OPEN(fpinv,file=TRIM(fetxt),status='replace',err=999)
        CLOSE(fpinv)
        fetxt = ramd(1:lnramd)//slash(1:1)//'run.ctr'
        OPEN(fprun,file=TRIM(fetxt),status='replace',err=999)
        !!!$  close(fprun) muss geoeffnet bleiben da sie staendig beschrieben wird
        fetxt = ramd(1:lnramd)//slash(1:1)//'cjg.ctr'
        OPEN(fpcjg,file=TRIM(fetxt),status='replace',err=999)
        CLOSE(fpcjg)
        fetxt = ramd(1:lnramd)//slash(1:1)//'eps.ctr'
        OPEN(fpeps,file=TRIM(fetxt),status='replace',err=999)

        !!!$  SET ERRORS for all measurements and write control to fpeps
        !!!$ >> RM
        IF (ldc) THEN
            WRITE (fpeps,'(a)')'1/eps_r eps_r datum'
            WRITE (fpeps,'(G12.3,2x,G12.3,2X,G14.5)')( &
                SQRT(wmatdr(i)),&
                1/SQRT(wmatdr(i)), &
                REAL(dat(i)), i=1,nanz)
        ELSE
            WRITE (fpeps,'(t5,a,t14,a,t27,a,t38,a,t50,a,t62,a,t71,a,t87,a)')'1/eps_r','1/eps_p',&
                '1/eps','eps_r','eps_p','eps','-log(|R|)', '-Phase (rad)'
            WRITE (fpeps,'(3F10.1,2x,3G12.3,2G15.7)')&
                (SQRT(wmatdr(i)),SQRT(wmatdp(i)),SQRT(wmatd(i)),1/SQRT(wmatdr(i)),&
                1/SQRT(wmatdp(i)),1/SQRT(wmatd(i)),REAL(dat(i)),AIMAG(dat(i)),i=1,nanz)
        END IF
        CLOSE(fpeps)
        errnr = 4

        !!!$ set starting model
        CALL bsigm0(kanal,dstart)
        IF (errnr.NE.0) GOTO 999

        !!!$   Kontrolldateien initialisieren
        !!!$   diff-        call kont1(delem,delectr,dstrom,drandb)
        !!!$   diff+<
        CALL kont1(delem,delectr,dstrom,drandb,dd0,dm0,dfm0,lagain)
        !!!$   diff+>
        IF (errnr.NE.0) GOTO 999

        !!$     write(6,"(a, i3)") " OpenMP max threads: ", OMP_GET_MAX_THREADS()
!!$     !$OMP PARALLEL
!!$     write(6,"(2(a,i3))") " OpenMP: N_threads = ",&
!!$          OMP_GET_NUM_THREADS()," thread = ", OMP_GET_THREAD_NUM()
!!$     !$OMP END PARALLEL

!!!$-------------
!!!$   get current time
     CALL tic(c1)
!!!$.................................................
     converged = .FALSE.

     DO WHILE (.NOT.converged) ! optimization loop

!!!$   Control output
        WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//clrln
        WRITE(*,'(a,i3,a,i3,a)',ADVANCE='no')ACHAR(13)//&
             ' Iteration ',it,', ',itr,' : Calculating Potentials'
        WRITE(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,&
             ' : Calculating Potentials'

!!!$   MODELLING
        count = 0
        IF (ldc) THEN
           fetxt = 'allocation problem adc'
           ALLOCATE (adc((mb+1)*sanz),STAT=myerr)
           IF (myerr /= 0) GOTO 999
           fetxt = 'allocation problem hpotdc'
           ALLOCATE (hpotdc(sanz,eanz),STAT=myerr)
           IF (myerr /= 0) GOTO 999
           ALLOCATE (bdc(sanz),STAT=myerr)
           fetxt = 'allocation problem adc'
           IF (myerr /= 0) GOTO 999

           !$OMP PARALLEL DEFAULT (none) &
           !$OMP FIRSTPRIVATE (pota,fak,pot,adc,bdc,fetxt) &
           !$OMP PRIVATE (j,l) &
           !$OMP SHARED (kwnanz,lverb,eanz,lsr,lbeta,lrandb,&
           !$OMP  lrandb2,sanz,kpotdc,swrtr,hpotdc,elbg,count,mb)
           !$OMP DO
!!!$   DC CASE

           DO k=1,kwnanz
              !$OMP ATOMIC
              count = count + 1
              fetxt = 'DC-Calculation wavenumber'
!!$              IF (lverb) WRITE (*,'(a,t35,I4,t100,a)',ADVANCE='no')&
!!$                   ACHAR(13)//TRIM(fetxt),count,''
              DO l=1,eanz
                 IF (lsr.OR.lbeta.OR.l.EQ.1) THEN
!!!$   Evtl calculation of analytical potentials
                    IF (lsr) CALL potana(l,k,pota)

!!!$   COMPilation of the linear system
                    fetxt = 'kompadc'
                    CALL kompadc(l,k,adc,bdc)
                    !                    if (errnr.ne.0) goto 999

!!!$   Evtl take Dirichlet boundary values into account
                    IF (lrandb) CALL randdc(adc,bdc)
                    IF (lrandb2) CALL randbdc2(adc,bdc)

!!!$   Scale the linear system (preconditioning stores fak)
                    fetxt = 'scaldc'
                    CALL scaldc(adc,bdc,fak)
                    !                    if (errnr.ne.0) goto 999
!!!$   Cholesky-Factorization of the Matrix
                    fetxt = 'choldc'
                    CALL choldc(adc)
                    !                    if (errnr.ne.0) goto 999
                 ELSE
                    fetxt = 'kompbdc'
!!!$   Modification of the current vector (Right Hand Side)
                    CALL kompbdc(l,bdc,fak)
                 END IF

!!!$   Solve linear system
                 fetxt = 'vredc'

                 CALL vredc(adc,bdc,pot)

!!!$   Scale back the potentials, save them and
!!!$   eventually add the analytical response
                 DO j=1,sanz
                    kpotdc(j,l,k) = DBLE(pot(j)) * fak(j)
                    IF (lsr) kpotdc(j,l,k) = kpotdc(j,l,k) + &
                         DBLE(pota(j))
                    IF (swrtr.EQ.0) hpotdc(j,l) = kpotdc(j,l,k)
                 END DO
              END DO

           END DO
           !$OMP END DO
           !$OMP END PARALLEL

        ELSE

           fetxt = 'allocation problem a'
           ALLOCATE (a((mb+1)*sanz),STAT=myerr)
           IF (myerr /= 0) GOTO 999
           fetxt = 'allocation problem hpot'
           ALLOCATE (hpot(sanz,eanz),STAT=myerr)
           IF (myerr /= 0) GOTO 999
           fetxt = 'allocation problem b'
           ALLOCATE (b(sanz),STAT=myerr)
           IF (myerr /= 0) GOTO 999
           !$OMP PARALLEL DEFAULT (none) &
           !$OMP FIRSTPRIVATE (pota,fak,pot,a,b,fetxt) &
           !$OMP PRIVATE (j,l,k) &
           !$OMP SHARED (kwnanz,lverb,eanz,lsr,lbeta,lrandb,lrandb2,sanz,kpot,swrtr,hpot,count)
           !$OMP DO
!!!$   COMPLEX CASE
           DO k=1,kwnanz
              !$OMP ATOMIC
              count = count + 1
              fetxt = 'IP-Calculation wavenumber'
              IF (lverb) WRITE (*,'(a,t35,I4,t100,a)',ADVANCE='no')&
                   ACHAR(13)//TRIM(fetxt),count,''
              DO l=1,eanz
                 IF (lsr.OR.lbeta.OR.l.EQ.1) THEN

!!!$   Ggf. Potentialwerte fuer homogenen Fall analytisch berechnen
                    IF (lsr) CALL potana(l,k,pota)

!!!$   Kompilation des Gleichungssystems (fuer Einheitsstrom !)
                    fetxt = 'kompab'
                    CALL kompab(l,k,a,b)
                    !                    if (errnr.ne.0) goto 999

!!!$   Ggf. Randbedingung beruecksichtigen
                    IF (lrandb) CALL randb(a,b)
                    IF (lrandb2) CALL randb2(a,b)

!!!$   Gleichungssystem skalieren
                    fetxt = 'scalab'
                    CALL scalab(a,b,fak)
                    !                    if (errnr.ne.0) goto 999

!!!$   Cholesky-Zerlegung der Matrix
                    fetxt = 'chol'
                    CALL chol(a)
                    !                    if (errnr.ne.0) goto 999
                 ELSE

!!!$   Stromvektor modifizieren
                    fetxt = 'kompb'
                    CALL kompb(l,b,fak)

                 END IF

!!!$   Gleichungssystem loesen
                 fetxt = 'vre'
                 CALL vre(a,b,pot)

!!!$   Potentialwerte zurueckskalieren und umspeichern sowie ggf.
!!!$   analytische Loesung addieren
                 DO j=1,sanz
                    kpot(j,l,k) = pot(j) * dcmplx(fak(j))
                    IF (lsr) kpot(j,l,k) = kpot(j,l,k) + pota(j)
                    IF (swrtr.EQ.0) hpot(j,l) = kpot(j,l,k)
                 END DO
              END DO
           END DO
           !$OMP END DO
           !$OMP END PARALLEL

        END IF


!!!$   Ggf. Ruecktransformation der Potentialwerte
        IF (swrtr.EQ.1) THEN
           CALL rtrafo(errnr)
           IF (errnr.NE.0) GOTO 999
        END IF
!!!$   Spannungswerte berechnen
        CALL bvolti()
        IF (errnr.NE.0) THEN
!!!$   reset model and data
           sigma = sigma2

           sigmaa = sgmaa2

           EXIT
        END IF
!!$  free some memory..
        IF (ldc) THEN
           DEALLOCATE(adc,hpotdc,bdc)
        ELSE
           DEALLOCATE(a,hpot,b)
        END IF

        IF (lsetup.OR.lsetip) THEN
!!!$   Ggf. background auf ratio-Daten "multiplizieren"
           IF (lratio) THEN
              DO j=1,nanz
                 dat(j) = dat(j) + sigmaa(j)
              END DO
           END IF

!!!$   Polaritaeten checken
           CALL chkpol(lsetup.OR.lsetip)
        END IF

!!!$   Daten-CHI berechnen
        CALL dmisft(lsetup.OR.lsetip)
        !        print*,nrmsd,betrms,pharms,lrobust,l1rat
        IF (errnr.NE.0) GOTO 999

        IF ((llam.AND..NOT.lstep).OR.lsetup) THEN
           IF (itmax == 0) THEN
              WRITE (*,'(t20,a,G14.4/)')'++ Calculated Fit',nrmsd
           ELSE
              IF (it == 0)  THEN
                 WRITE (*,'(a,G14.4/)')'++ Starting Fit',nrmsd
              ELSE
                 WRITE (*,'(t10,a,G14.4/)',ADVANCE='no')'++ Actual Fit',nrmsd
              END IF
           END IF
        ELSE
           WRITE (*,'(t10,a,G14.4)',ADVANCE='no')'-- Update Fit',nrmsd
        END IF

!!!$   'nrmsd=0' ausschliessen
        IF (nrmsd.LT.1d-12) nrmsd=nrmsdm*(1d0-mqrms)

!!!$   tst
!!!$   tst        if (lfphai) then
!!!$   tst            llam = .true.
!!!$   tst            if (.not.lfpi) nrmsd = 1d0
!!!$   tst        end if

!!!$.............................
        IF (it == 0) THEN
           WRITE (*,'(/a,t100/)')ACHAR(13)//&
                'WRITING STARTING MODEL'
           CALL wout(kanal,dsigma,dvolt)
        END IF
!!!$   Kontrollvariablen ausgeben
        CALL kont2(lsetup.OR.lsetip)
        IF (errnr.NE.0) GOTO 999

!!!$   ABBRUCHBEDINGUNGEN
        IF (llam.AND..NOT.lstep) THEN

!!!$   Polaritaeten checken
           CALL chkpol(lsetup.OR.lsetip)

!!!$   Wiederholt minimale step-length ?
           IF (stpalt.EQ.0d0) THEN
              errnr2 = 92
              fetxt = 'repeated step length'
           END IF
           WRITE (*,'(/a,G12.4,a/)')'+++ Convergence check (CHI (old/new)) ',&
                100.0*(1d0-rmsalt/nrmsd),' %'
!!!$   Keine Verbesserung des Daten-CHI ?
           IF (dabs(1d0-rmsalt/nrmsd).LE.mqrms) THEN
              errnr2 = 81
              WRITE (fetxt,*)'No further CHI approvement ',&
                   REAL(ABS(1d0-rmsalt/nrmsd))
           END IF
!!!$   Minimaler Daten-CHI erreicht ?
!!!$   tst            if (dabs(1d0-nrmsd/nrmsdm).le.mqrms) errnr2=80
           IF (dabs(1d0-nrmsd/nrmsdm).LE.mqrms.AND.ldlamf) THEN
              errnr2 = 80
              WRITE (fetxt,*)'Optimal RMS ',REAL(nrmsd),' reached'
           END IF

!           IF (llam) THEN
!              WRITE (6,'(a)',ADVANCE='no')'convergence '
!              IF (nrmsd < 1d0 ) errnr2 = 94
!              IF (nrmsd > rmsalt) errnr2 = 95
!              IF (dabs(1d0-nrmsd/rmsalt) < mqrms) errnr2 = 93
!              PRINT*,errnr2
!           END IF

!!!$   Maximale Anzahl an Iterationen ?
           IF (it.GE.itmax) THEN
              errnr2 = 79
              WRITE (fetxt,*)'Reached max number of iterations ',itmax
           END IF
!!!$   Minimal stepsize erreicht ?
           IF (dabs(1d0-bdpar/bdalt) < bdmin) THEN
              errnr2 = 109
              WRITE (fetxt,*)' Min model update reached ',dabs(1d0-bdpar/bdalt)
           END IF

!!!$   Ggf. abbrechen oder "final phase improvement"
           IF (errnr2.NE.0) THEN

              IF (lfphai.AND.errnr2.NE.79) THEN
                 PRINT*,'CRI termination '//TRIM(fetxt),errnr2
!!!$   ak
!!!$   Widerstandsverteilung und modellierte Daten ausgeben
                 CALL wout(kanal,dsigma,dvolt)
                 IF (errnr.NE.0) GOTO 999

!!!$   Kontrollausgaben
                 WRITE(*,'(/a/)')&
                      '**** Final phase improvement ****'

                 WRITE(fprun,'(a24)',err=999)&
                      ' Final phase improvement'

                 fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
                 OPEN(fpinv,file=TRIM(fetxt),status='old',&
                      POSITION='append',err=999)
                 WRITE(fpinv,'(/a/)',err=999)&
                      '------------------------------------------------'//&
                      '------------------------------------------------'//&
                      '-----------------'
                 CLOSE(fpinv)

!!!$   Wichtungsfeld umspeichern
                 wmatd = wmatdp
                 lam_cri = lamalt

                 WRITE (*,'(/a,g12.4/)')'++ (FPI) setting phase error '//&
                      'and saving lam_cri: ',REAL(lam_cri)
                 WRITE (fprun,'(/a,g12.4/)')'++ (FPI) setting phase error '//&
                      'and saving lam_cri: ',REAL(lam_cri)

                 lfpi    = .TRUE.
                 lsetip = .TRUE. !
                 lfphai = .FALSE.
                 llam   = .FALSE.
                 ldlami = .TRUE.
                 lfstep = .TRUE.
                 step   = 1d0
                 errnr2 = 0 ! reset convergence case "error"

!!!$   ak
                 fetxt = 'cp -f inv.lastmod inv.lastmod_rho'
                 CALL SYSTEM (TRIM(fetxt))
                 IF (lffhom) THEN
                    WRITE(*,*)&
                         ' ******* Restarting phase model ********'
                    WRITE(fprun,*)&
                         ' ******* Restarting phase model ********'
                    fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
                    OPEN (fpinv,file=TRIM(fetxt),status='old',err=999,&
                         position='append')
                    WRITE (fpinv,*)&
                         ' ******* Resetting phase model ********'
                    CLOSE (fpinv)
                    DO j=1,elanz
                       sigma(j) = dcmplx(&
                            dcos(pha0/1d3)*cdabs(sigma(j)) ,&
                            -dsin(pha0/1d3)*cdabs(sigma(j)) )
                    END DO
                    lsetup = .TRUE. ! ensure proper misfit and kont2 output
                    CYCLE       ! neues calc
                 ELSE
                    lsetup = .FALSE.
                 END IF
!!!$   ak

!!!$   Daten-CHI berechnen
                 CALL dmisft(lsetip)
                 IF (errnr.NE.0) GOTO 999
                 WRITE (*,'(a,t45,a,t78,F14.4)') &
                      ACHAR(13),'++ Phase data fit',nrmsd

!!!$   Kontrollvariablen ausgeben
                 CALL kont2(lsetip)
                 IF (errnr.NE.0) GOTO 999
              ELSE
                 EXIT
              END IF
           ELSE
!!!$ >> RM
              lamalt = lam ! save the lambda of the previous iteration
!!!$ if, and only if the iterate was successful...
!!!$ << RM
!!!$   ak
!!!$   Widerstandsverteilung und modellierte Daten ausgeben
!!$              WRITE (*,'(a,t30,I4,t100,a)')ACHAR(13)//&
!!$                   'WRITING MODEL ITERATE',it,''
              CALL wout(kanal,dsigma,dvolt)
              IF (errnr.NE.0) GOTO 999
           END IF
        END IF

        IF ((llam.AND..NOT.lstep).OR.lsetup.OR.lsetip) THEN

!!!$   Iterationsindex hochzaehlen
           it = it+1

!!!$   Parameter zuruecksetzen
           itr    = 0
           rmsreg = 0d0
           dlam   = 1d0
           dlalt  = 1d0
           ldlamf = .FALSE.
           bdpar = 1d0
           bdalt = 1d0

!!!$   Daten-CHI speichern
           rmsalt = nrmsd

!!!$   Lambda speichen
           IF (it>1) lamalt = lam ! mindestens einmal konvergiert

!!!$   Felder speichern
           sigma2 = sigma

           sgmaa2 = sigmaa

           IF (lrobust) wmatd2 = wmatd

!!!$   Kontrollausgaben
           WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//clrln
           WRITE (*,'(a,i3,a,i3,a)',ADVANCE='no') &
                ACHAR(13)//' Iteration ',it,', ',itr,&
                ' : Calculating Sensitivities'

           WRITE(fprun,'(a,i3,a,i3,a)')' Iteration ',it,', ',itr,&
                ' : Calculating Sensitivities'

!!!$   SENSITIVITAETEN berechnen
           IF (ldc) THEN
              CALL bsendc((it==0))
           ELSE
              CALL bsensi((it==0))
           END IF

!!!$   evtl   Rauhigkeitsmatrix oder CmS-Matrix belegen
           IF (l_bsmat) CALL bsmatm(it,l_bsmat)

        ELSE

           IF (lverb) THEN
              CALL wout_up(kanal,it,itr,.FALSE.)
              IF (errnr.NE.0) GOTO 999
           END IF

!!!$   Felder zuruecksetzen
           sigma = sigma2

           sigmaa = sgmaa2

           IF (lrobust) wmatd = wmatd2

        END IF
        IF (itmax == 0) THEN ! only precalcs..
           errnr2 = 109
           PRINT*,'Only precalcs'
           EXIT
        END IF

        IF (BTEST(llamf,0)) THEN ! for fixed lambda we do not want any parabola fitting?
           lam = lamfix
           IF (BTEST(llamf,1).AND.it>1) lam = lamfix / (2d0 * DBLE(it - 1))
           llam = .FALSE. ! in order to calculate any update this needs to be false
           IF (lsetup.OR.lsetip) THEN
              lsetup = .FALSE.
              lsetip = .FALSE.
           END IF

!!!$   REGULARISIERUNG / STEP-LENGTH einstellen
        ELSE
           IF (.NOT.lstep) THEN

              IF (llam) THEN
!!!$   "Regularisierungsschleife" initialisieren und step-length zuruecksetzen
                 llam = .FALSE.
                 step = 1d0
              ELSE

!!!$   Regularisierungsindex hochzaehlen
                 itr = itr+1
                 IF ((((nrmsd.LT.rmsreg.AND.itr.LE.nlam).OR. &
                      (dlam.GT.1d0.AND.itr.LE.nlam)).AND.&
                      (.NOT.ldlamf.OR.dlalt.LE.1d0).AND.&
                   (dabs(1d0-bdpar/bdalt) > bdmin).AND.&
                      (dabs(1d0-rmsreg/nrmsdm).GT.mqrms)).OR.&
                      (rmsreg.EQ.0d0)) THEN

                    IF (rmsreg > 0d0) THEN
                       WRITE (fprun,'(/a,G12.4,a)')'Chi increase:',&
                            100.0*(1d0-rmsalt/nrmsd),' %'
                       WRITE (fprun,'(a,G12.4,a)')'Stepsize :',bdpar
                       WRITE (fprun,'(a,G12.4/)')'nrmsd/rmsreg :',nrmsd/rmsreg
                    END IF
!!!$   Regularisierungsparameter bestimmen
                    IF (lsetup.OR.lsetip) THEN ! general initialization, lam0

!!!$   Kontrollausgabe
                       WRITE(*,'(a,i3,a,i3,a,t100,a)',ADVANCE='no')&
                            ACHAR(13)//' Iteration ',it,', ',itr,&
                            ' : Calculating 1st regularization parameter',''
                       WRITE(fprun,'(a,i3,a,i3,a)',ADVANCE='no')&
                            ' Iteration ',it,', ',itr,&
                            ' : Calculating 1st regularization parameter'
                       CALL blam0()
                       WRITE (*,'(a,G10.2)',ADVANCE='no')'lam_0:: ',lammax
                       WRITE (fprun,'(a,G10.2)')'lam_0 ',lammax
                       lam = lammax
                       lamalt = lammax
!!!$   ak Model EGS2003, ERT2003                        call blam0()
!!!$   ak Model EGS2003, ERT2003                        lam = lammax
!!!$   ak                        lam = 1d4
!!!$
!!!$ GENERAL REMARK ON OUR REGULARISATION PARAMETER
!!!$ Standard normal equations are
!!!$ (A^h * Cd^-1 A + Cm^-1) dm = A_q^h * C_d^-1 * (d - f(m_q)) + C_m^-1 * ({m_q,(m_q-m_0)})
!!!$ Where we identify lam as inverse a-priori variance:
!!!$ Cm^-1 = \lam R^T * R.
!!!$ If we solve the normal equation
!!!$  Cm * (A^h * Cd^-1 A + I) dm = Cm * (A_q^h * C_d^-1 * (d - f(m_q)) - ({m_q,(m_q-m_0)}))
!!!$ instead, we use \lam as prior variance.
!!!$ This also results in the fundamental problem on how to adjust the search direction which depends on the
!!!$ CHI decrease and thus should also be reciprocal.
!!!$
!!!$ nrmsdm and mqrms, fstart and fstop are set in rall.f90.
!!!$ Defaults are
!!!$ nrmsdm := 1d0; mqrms := 2d-2; fstart := 0.5; fstop := 0.9
!!!$
                    ELSE ! for lambda search we go here
                       dlalt = dlam
                       IF (ldlami) THEN ! initialize search
                          ldlami = .FALSE.
                          alam   = dmax1(dabs(dlog(nrmsd/nrmsdm)),&
                               dlog(1d0+mqrms))
!!!$ alam = MAX(log(actual chi),log(1+0.02))
                          dlam   = fstart ! sets first dlam (0.5 default, which should be fine)
                       ELSE
                          alam = dmax1(alam,dabs(dlog(nrmsd/nrmsdm)))
!!!$ CHI dependend partial fraction of lam variation
!!!$ alam = MAX(alam,log(actual chi))
                          dlam = dlog(fstop)*&
                               SIGN(1d0,dlog(nrmsd/nrmsdm))+&
                               dlog(fstart/fstop)*&
                               dlog(nrmsd/nrmsdm)/alam
!!!$
!!!$ dlam = ln(0.9) * sign(1,ln(act chi)) + (ln(0.5/0.9)) * ln(act chi)/alam
!!!$ this makes mostly the same and dlam = exp(ln(0.9) + ln(0.5/0.9)) = 0.5
!!!$ This holds until the act chi value drops below 1. Than sign gives -1 and
!!!$ also alam is now not identical to log(act chi) but log(1.02), which results
!!!$ (for example act chi = 0,98)
!!!$ dlam = exp(-ln(0.9) + ln(0.5/0.9) * ln(0.98)/log(1.02)) = 2
!!!$
                          dlam = dexp(dlam)

                       END IF
                       lam = lam*dlam
                       IF (dlalt.GT.1d0.AND.dlam.LT.1d0) ldlamf=.TRUE.
!!!$   tst                        if (dlam.gt.1d0) lfstep=.true.
!!!$   ak Model EGS2003
                       IF (dlam.GT.1d0) lrobust=.FALSE.

                    END IF
                 ELSE

!!!$   Regularisierungsparameter zuruecksetzen und step-length verkleinern
!!!$ if no Chi decrease found
                    llam = .TRUE.
                    lam  = lam/dlam

                    IF (lfstep) THEN
                       lfstep = .FALSE.
                    ELSE
                       lstep = .TRUE.
                       step  = 5d-1
                    END IF
                 END IF

!!!$   Ggf. Daten-CHI speichern
                 IF (lsetup.OR.lsetip) THEN
                    lsetup = .FALSE.
                    lsetip = .FALSE.
                 ELSE
                    IF (.NOT.lstep) rmsreg=nrmsd
                 END IF
              END IF
           ELSE
              lstep = .FALSE.

!!!$   Parabolische Interpolation zur Bestimmung der optimalen step-length
              CALL parfit(rmsalt,nrmsd,rmsreg,nrmsdm,stpmin)

              IF (step.EQ.stpmin.AND.stpalt.EQ.stpmin)THEN

!!!$   Nach naechstem Modelling abbrechen
                 stpalt = 0d0
              ELSE

!!!$   Step-length speichern
                 stpalt = step
              END IF

           END IF

        END IF
!!!$   Kontrollausgaben
        WRITE (*,'(a)',ADVANCE='no')ACHAR(13)//clrln
        WRITE(*,'(a,i3,a,i3,a)',ADVANCE='no')&
             ACHAR(13)//' Iteration ',it,', ',itr,&
             ' : Updating'

        WRITE(fprun,*)' Iteration ',it,', ',itr,&
             ' : Updating'

!!!$ Modell parameter mit aktuellen Leitfaehigkeiten belegen

        CALL bpar
        IF (errnr.NE.0) GOTO 999

!!!$   UPDATE anbringen
        bdalt = bdpar
        CALL update
        IF (errnr.NE.0) GOTO 999

!!!$ Leitfaehigkeiten mit verbessertem Modell belegen
        CALL bsigma
        IF (errnr.NE.0) GOTO 999

        IF (lverb) THEN
           CALL wout_up(kanal,it,itr,.TRUE.)
           IF (errnr.NE.0) GOTO 999
        END IF
!!!$   Roughness bestimmen
        CALL brough

!!!$   Ggf. Referenzleitfaehigkeit bestimmen
!!!$ >>> RM
        !! $THIS has now a general meaning with the
!!!$ mixed boundary, since we like to set sigma0
!!!$ as reference sigma as "mean" boundary value
!!!$ <<< RM
        IF (lbeta) CALL refsig()

        IF (BTEST(llamf,0)) THEN
           llam = .TRUE.
        END IF

!!!$   Neues Modelling
     END DO ! DO WHILE (.not. converged)

!!!$ RESET FPI status variable to proceed with full COMPLEX calculus
     ! lfpi = .FALSE.
!!!$
!!!$.................................................

!!!$   OUTPUT
     WRITE (*,'(a,t25,I4,t35,a,t100,a)')ACHAR(13)//&
          'MODEL ESTIMATE AFTER',it,'ITERATIONS',''


!!!$   Kontrollausgaben
!!!$ errnr2 is more a status variable than a pure error number
!!!$ it shows which case was the reason for termination of the algorithm
     SELECT CASE (errnr2)
     CASE (95)

        WRITE(*,'(a22,a31)') ' Iteration terminated:',&
             ' no CHI decrease'
        WRITE(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' no CHI decrease'
!!!$ reset model to previous state
        IF (BTEST(llamf,0)) THEN
           WRITE (*,'(a)',ADVANCE='no')'Taking model state of previous iteration... '
           sigma = sigma2
           sigmaa = sgmaa2
           nrmsd = rmsalt
           PRINT*,'CHI = ',REAL(nrmsd)
        END IF

     CASE (94)

        WRITE(*,'(a22,a31)') ' Iteration terminated:',&
             ' CHI < 1'
        WRITE(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' CHI < 1'

     CASE (93)

        WRITE(*,'(a22,a31)') ' Iteration terminated:',&
             ' CHI decrease sufficiently'
        WRITE(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' CHI decreasse sufficiently'

     CASE (92)

        WRITE(*,'(a22,a31)') ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'

        WRITE(fprun,'(a22,a31)',err=999) ' Iteration terminated:',&
             ' Min. step-length for 2nd time.'
     CASE (80)
        WRITE(*,'(a22,a10)') ' Iteration terminated:', &
             ' Min. CHI.'

        WRITE(fprun,'(a22,a10)',err=999) ' Iteration terminated:',&
             ' Min. CHI.'
     CASE (81)
        WRITE(*,'(a22,a24)') ' Iteration terminated:',&
             ' Min. rel. CHI decrease.'

        WRITE(fprun,'(a22,a24)',err=999) ' Iteration terminated:',&
             ' Min. rel. CHI decrease.'
     CASE (79)
        WRITE(*,'(a22,a19)') ' Iteration terminated:',&
             ' Max. # iterations.'

        WRITE(fprun,'(a22,a19)',err=999) ' Iteration terminated:',&
             ' Max. # iterations.'

     CASE (109)
        WRITE(*,'(a)') ' Iteration terminated:'//&
             ' Min. model changes reached'

        WRITE(fprun,'(a)',err=999) ' Iteration terminated:'//&
             ' Min. model changes reached'

     END SELECT

     CALL wout(kanal,dsigma,dvolt)
     IF (errnr.NE.0 .AND..NOT. errnr == 82) GOTO 999

!!!$   Kontrollausgaben

!!!$   Run-time abfragen und ausgeben
     fetxt = ' CPU time: '
     CALL toc(c1,fetxt)
!!$     PRINT*,''
!!$     PRINT*,''
     WRITE(fprun,'(a)',err=999)TRIM(fetxt)

!!!$   Kontrolldateien schliessen
     CLOSE(fpinv)
     CLOSE(fpcjg)
     CLOSE(fpeps)


!!!$   Ggf. Summe der Sensitivitaeten aller Messungen ausgeben
     IF (lsens) THEN
        CALL BBSENS(kanal,dsens)
        IF (errnr.NE.0) GOTO 999
     END IF

     IF (lvario) THEN
        !        IF (lsens) CALL bvariogram_s ! calculate experimental variogram
        CALL bvariogram
     END IF

     IF (lcov1) CALL buncert (kanal,lamalt)

     CALL des_cjgmod(1,fetxt,myerr) ! call cjgmod destructor
     IF (myerr /= 0) GOTO 999

!!!$   'sens' und 'pot' freigeben
     IF (ldc) THEN
        fetxt = 'allocation sensdc'
        DEALLOCATE (sensdc)
        fetxt = 'allocation koptdc'
        DEALLOCATE (kpotdc)
     ELSE
        DEALLOCATE(sens,kpot)
     END IF

     fetxt = 'allocation smatm'
     IF (ALLOCATED (smatm)) DEALLOCATE (smatm)

     fetxt = 'allocation pot,pota,fak'
     IF (ALLOCATED (pot)) DEALLOCATE (pot,pota,fak)

     fetxt = 'allocation snr,sx,sy'
     IF (ALLOCATED (snr)) DEALLOCATE (snr,sx,sy)

     fetxt = 'allocation typ'
     IF (ALLOCATED (typ)) DEALLOCATE (typ,nelanz,selanz)

     fetxt = 'allocation nrel'
     IF (ALLOCATED (nrel)) DEALLOCATE (nrel,rnr)

     fetxt = 'allocation esp'
     IF (ALLOCATED (espx)) DEALLOCATE (espx,espy)

     fetxt = 'allocation kwn'
     IF (ALLOCATED (kwn)) DEALLOCATE (kwn)

     fetxt = 'allocation kwni'
     IF (ALLOCATED (kwnwi)) DEALLOCATE (kwnwi)

     fetxt = 'allocation elbg'
     IF (ALLOCATED (elbg)) DEALLOCATE (elbg,relbg,kg)

     fetxt = 'allocation enr'
     IF (ALLOCATED (enr)) DEALLOCATE (enr)

     fetxt = 'allocation mnr'
     IF (ALLOCATED (mnr)) DEALLOCATE (mnr)

     fetxt = 'allocation strnr,strom,volt,etc'
     IF (ALLOCATED (strnr)) DEALLOCATE (strnr,strom,volt,sigmaa,&
          kfak,wmatdr,wmatdp,vnr,dat,wmatd,wmatd2,sgmaa2,wdfak,&
          wmatd_cri) !!! these are allocated in rdati!!!

     fetxt = 'allocation par,dpar,dpar2'

     IF (ALLOCATED (par)) DEALLOCATE (par,dpar,dpar2)

     fetxt = 'allocation sigma'

     IF (ALLOCATED (sigma)) DEALLOCATE (sigma,sigma2)

     IF (ALLOCATED (d0)) DEALLOCATE (d0,fm0)
     IF (ALLOCATED (m0)) DEALLOCATE (m0)

     IF (ALLOCATED (rwddc)) DEALLOCATE (rwddc)
     IF (ALLOCATED (rwndc)) DEALLOCATE (rwndc)
     IF (ALLOCATED (rwd)) DEALLOCATE (rwd)
     IF (ALLOCATED (rwn)) DEALLOCATE (rwn)
     IF (ALLOCATED (rwdnr)) DEALLOCATE (rwdnr)

     IF (ALLOCATED (w_ref_re)) DEALLOCATE (w_ref_re,w_ref_im,m_ref)

    !!! write final statement in run.ctr
     WRITE(fprun,'(a)') '***finished***'
     CLOSE(fprun)

!!!$   Ggf. weiteren Datensatz invertieren
  END DO

!!!$   'crtomo.cfg' schliessen
  CLOSE (fpcfg)

!!!$ crtomo.pid löschen
  fetxt = 'crtomo.pid'
  OPEN (fprun,FILE=TRIM(fetxt),STATUS='old',err=999)
  CLOSE (fprun,STATUS='delete')

!!!$ finished-string an inv.ctr anhängen
  fetxt = ramd(1:lnramd)//slash(1:1)//'inv.ctr'
  OPEN(fpinv,file=TRIM(fetxt),status='old',POSITION='append',err=999)
  WRITE (fpinv,'(a)')'***finished***'
  CLOSE(fpinv)

  !!! also write the final statement to STDOU
  WRITE(*,*) '***finished***'
  STOP '0'

!!!$.....................................................................

!!!$   Fehlermeldung
999 OPEN(fperr,file='error.dat',status='replace')
  errflag = 2
  CALL get_error(ftext,errnr,errflag,fetxt)
  WRITE(fperr,*) 'CRTomo PID ',pid,' exited abnormally'
  WRITE(fperr,'(a80,i3,i1)') fetxt,errnr,errflag
  WRITE(fperr,*)ftext
  CLOSE(fperr)

  STOP '-1'

END PROGRAM inv
