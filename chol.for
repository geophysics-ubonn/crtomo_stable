        subroutine chol()
 
c Cholesky-Zerlegung der positiv definiten Matrix 'a'; erfolgt auf dem
c Platz von 'a', d.h. bei Auftreten eines Fehlers ist gegebene Matrix
c 'a' zerstoert.

c ( Vgl. Subroutine 'CHOBNDN' in Schwarz (1991) )

c Andreas Kemna                                            11-Oct-1993
c                                       Letzte Aenderung   07-Mar-2003

c.....................................................................

        USE alloci

        INCLUDE 'parmax.fin'
        INCLUDE 'err.fin'
        INCLUDE 'elem.fin'

c.....................................................................

c PROGRAMMINTERNE PARAMETER:

c Hilfsvariablen
        integer         * 4     idi,i0,ij,j0
        integer         * 4     m1,fi
        complex         * 16    s

c Indexvariablen
        integer         * 4     i,j,k

c.....................................................................

        m1 = mb+1

        do 30 i=1,sanz

            idi = i*m1
            fi  = max0(1,i-mb)
            i0  = idi-i

            do 20 j=fi,i

                ij = i0+j
                j0 = j*mb
                s  = a(ij)

                do 10 k=fi,j-1
                    s = s - a(i0+k)*a(j0+k)
10              continue

                if (j.lt.i) then

                    a(ij) = s / a(j*m1)

                else

                    if (cdabs(s).le.0d0) then
                        fetxt = ' '
                        errnr = 28
                        goto 1000
                    end if

                    a(idi) = cdsqrt(s)

                end if

20          continue

30      continue

        errnr = 0
        return

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c Fehlermeldungen

1000    return

        end
