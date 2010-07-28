      subroutine wpot(kanal,datei,np)

c     Unterprogramm zum Schreiben der Potentialwerte.

c     Andreas Kemna                                            17-Dec-1993
c     Letzte Aenderung   10-Mar-2007

c.....................................................................

      USE femmod
      USE datmod
      USE elemmod

      IMPLICIT none

      INCLUDE 'err.fin'

c.....................................................................

c     EIN-/AUSGABEPARAMETER:

c     Kanalnummer
      integer         * 4     kanal

c     Datei
      character       * 80    datei

c     Nummer der Projektion
      integer         * 4     np

c.....................................................................

c     PROGRAMMINTERNE PARAMETER:

c     Indexvariable
      integer         * 4     i

c     Hilfsvariablen
      integer         * 4     idum,idum2,lnanz
      character       * 80    htxt
      character       * 12    htxt2

c     Hilfsfunctions
      character       * 12    intcha
      character       * 80    filpat

c     (Back-) Slash
      character       * 1     slash

c     tst        real            * 8     dum_re,dum_im,dum_mag,dum_pha

c.....................................................................

c     Slash
      slash = '/'

c     'datei' modifizieren
      lnanz = int(log10(real(nanz)))+1
      htxt  = filpat(datei,idum2,1,slash(1:1))
      idum  = idum2+index(datei(idum2+1:80),'.')-1

      if ((idum-idum2-1).gt.(8-lnanz)) then
         fetxt = datei
         errnr = 15
         goto 999
      end if

      htxt2 = intcha(np,lnanz)
      htxt  = datei(1:idum)//htxt2(1:lnanz)//datei(idum+1:idum+4)

c     'datei' oeffnen
      fetxt = htxt
      errnr = 1
      open(kanal,file=fetxt,status='replace',err=999)
      errnr = 4

c     Koordinaten und Potentialwerte (Real- und Imaginaerteil) der
c     Knotenpunkte schreiben (in Reihenfolge der urspruenglichen
c     Numerierung)
      do i=1,sanz
         write(kanal,*,err=1000)
     1        real(sx(snr(i))),real(sy(snr(i))),
     1        real(dble(pot(i))),real(dimag(pot(i)))
c     ak     1                real(cdabs(pot(i))),
c     ak     1                real(1d3*datan2(dimag(pot(i)),dble(pot(i))))

c     tst                dum_re  = dble (pot(i))
c     tst                dum_im  = dimag(pot(i))
c     tst                dum_mag = cdabs(pot(i))
c     tst                dum_pha = datan2(dum_im,dum_re)
c     tst                if (dum_re.lt.0d0) then
c     tst                    dum_mag = -dum_mag
c     tst                    if (dum_im.lt.0d0) then
c     tst                        dum_pha = dum_pha + dacos(-1d0)
c     tst                    else
c     tst                        dum_pha = dum_pha - dacos(-1d0)
c     tst                    end if
c     tst                end if
c     tst                write(kanal,'(2(f6.2,2x),4(g12.5,2x))',err=1000)
c     tst     1                sx(snr(i)),sy(snr(i)),
c     tst     1                dum_re,dum_im,dum_mag,1d3*dum_pha
      end do

c     'datei' schliessen
      close(kanal)

      errnr = 0
      return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     Fehlermeldungen

 999  return

 1000 close(kanal)
      return

      end
