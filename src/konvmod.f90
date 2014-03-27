MODULE konvmod
!!!$'konv.fin'
!!!$ Andreas Kemna                                          01-Mar-1995
!!!$                                      
!!!$ Last changed       RM                                  Jul-2010
!!!$
!!!$.....................................................................
use alloci
!!!$Regularisierungsparameter
  REAL(prec), PUBLIC ::     lam

!!!$Maximaler Regularisierungsparameter
  REAL(prec), PUBLIC ::     lammax

!!$ >> RM
!!!$ Post CRI regularization parameter
  REAL(prec), PUBLIC ::     lam_cri
!!$ Starting regularization parameter (CRI)
  REAL(prec), PUBLIC ::     lamnull_cri
!!$ Starting regularization parameter (FPI)
  REAL(prec), PUBLIC ::     lamnull_fpi
!!$ << RM

!!$ >> RM
!!!$ Reference model regularization lambda (currently its a factor for lam)
  REAL(prec), PUBLIC ::     lam_ref
!!!$ switch if we want the absolute or just the gradient regularized
  INTEGER, PUBLIC         ::     lam_ref_sw
!!! << RM

!!!$Fester Regularisierungsparameter
  REAL(prec), PUBLIC ::     lamfix

!!!$Schrittweitenfaktor, maximale Anzahl der Schritte
  REAL(prec), PUBLIC ::     dlam,dlalt
  INTEGER,PUBLIC ::     nlam

!!!$Parameter "a", "initial/final step factor"
  REAL(prec), PUBLIC ::     alam
  REAL(prec), PUBLIC ::     fstart
  REAL(prec), PUBLIC ::     fstop

!!!$Step-length der Modellverbesserung
  REAL(prec), PUBLIC ::     step

!!!$Step-length der vorherigen Iteration
  REAL(prec), PUBLIC ::     stpalt

!!!$Mindest-step-length
  REAL(prec), PUBLIC ::     stpmin

!!!$Schalter zur Steuerung der Regularisierung und der step-length-Wahl
  LOGICAL,PUBLIC  ::     llam,lstep,ldlami,ldlamf
!!!$wenn llamf >0 ist wird ein fixes Lambda gesetzt, llamf = BTEST(llamf,1) fuer lambda 
!!!$ cooling 
  INTEGER,PUBLIC  ::     llamf
!!!$Schalter ob volle step-length angebracht werden soll
  LOGICAL,PUBLIC  ::     lfstep

!!!$Daten-RMS
  REAL(prec), PUBLIC ::     nrmsd

!!!$Daten-RMS der vorherigen Iteration
  REAL(prec), PUBLIC ::     rmsalt

!!!$Daten-RMS des vorherigen "Regularisierungsschritts"
  REAL(prec), PUBLIC ::     rmsreg

!!!$Daten-Misfit insgesamt (in der objective function),
!!!$Misfit der Betraege, Misfit der Phasen
  REAL(prec), PUBLIC ::     rmssum,betrms,pharms

!!!$Roughness
  REAL(prec), PUBLIC ::     rough

!!!$Minimaler Daten-RMS
  REAL(prec), PUBLIC ::     nrmsdm

!!!$Minimaler Quotient zweier aufeinanderfolgender Daten-RMS-Werte
!!!$bzw. zwischen tatsaechlichem und minimalem Daten-RMS
  REAL(prec), PUBLIC ::     mqrms

!!!$Minimale "L1-ratio" (Grenze der "robust inversion")
  REAL(prec), PUBLIC ::     l1min

!!!$"L1-ratio"
  REAL(prec), PUBLIC ::     l1rat

!!!$Norm des Verbesserungsvektors (stepsize)
  REAL(prec), PUBLIC ::     bdpar

!!!$Minimale step size 
  REAL(prec), PUBLIC ::     bMIN

!!!$Felddimensionen in x- und z-Richtung
  INTEGER,PUBLIC ::     nx,nz

!!!$Maximale Anzahl der Iterationen
  INTEGER,PUBLIC ::     itmax

!!!$Iterationsindex
  INTEGER,PUBLIC ::     it

!!!$Index des Regularisierungsschritts
  INTEGER,PUBLIC ::     itr

!!!$Parameter zur Glaettung in x- und z-Richtung
!!!$(=0 bedeutet keine Glaettung in der jeweiligen Richtung)
  REAL(prec), PUBLIC ::     alfx,alfz
!!!$MGS beta
  REAL(prec), PUBLIC ::     betamgs

!!!$Schalter ob "robust inversion" durchgefuehrt werden soll
  LOGICAL,PUBLIC  ::     lrobust

!!$cdiff+<
!!!$Schalter ob "difference inversion" durchgefuehrt werden soll
  LOGICAL,PUBLIC  ::     ldiff
!!$cdiff+>
!!!$Schalter ob positive Phasen auf Null korrigiert werden sollen
  LOGICAL,PUBLIC  ::     lphi0
!!!$Schalter ob "final phase improvement" durchgefuehrt werden soll
  LOGICAL,PUBLIC  ::     lfphai
!!!$Schalter ob "final phase improvement" mit homogenem phasenwert startet
  LOGICAL,PUBLIC  ::     lffhom
!!$c	triang>
!!!$Schalter für die Triangularisierungsinversion
  INTEGER,PUBLIC :: ltri
!!$c	triang<	
!!!$Schalter ob prior model vorhanden
  LOGICAL,PUBLIC  :: lprior
!!!$ >> RM ref model regu
!!!$switch is true for reference model regularization (DAMPING)
  LOGICAL,PUBLIC  :: lw_ref
!!!$ << RM ref model regu
!!!$Schalter ob prior model verrauschen
  LOGICAL,PUBLIC  :: lnsepri
!!!$Schalter ob Summe der Sensitivitaeten aller Messungen ausgegeben
!!!$werden soll
  LOGICAL,PUBLIC  :: lsens
!!!$Schalter fuer Resolution matrix
  LOGICAL,PUBLIC  :: lres
!!!$Schalter fuer posterior cov
  LOGICAL,PUBLIC  :: lcov1
!!!$Schalter fuer posterior cov
  LOGICAL,PUBLIC  :: lcov2
!!!$Mega schalter, wird binaer getestet
  INTEGER,PUBLIC :: mswitch
!!!$Schalter für ols Loeser
  LOGICAL,PUBLIC  :: lgauss
!!!$Schalter ob experimentelles variogram berechnet werden soll
  LOGICAL,PUBLIC  :: lvario
!!!$ Verbose output switch
  LOGICAL,PUBLIC  :: lverb
!!!$ Verbose output of full resoultion, cm0, cm0_inv, full covariance, etc..
  LOGICAL,PUBLIC  :: lverb_dat
!!!$ lsytop = .NOT.BTEST (mswitch,8)  disables sy top check of 
!!$ no flow boundary electrodes for enhanced beta calculation (bsytop). 
!!$  This is useful for including topographical effects and should be used
  LOGICAL,PUBLIC  ::     lsytop
!!!$ switch to regard error ellipses in complex plane 
!!!$ is srt false as default if we have FPI but may be overwritten with
!!!$ mswitch + 32 setting it to true
  LOGICAL, PUBLIC ::     lelerr

END MODULE konvmod
