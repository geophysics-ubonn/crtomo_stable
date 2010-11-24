# $Id: Makefile,v 1.0 2009/07/02 12:43:39 roland Exp $
#
#      Makefile for CRTomo (default)
#       

RM		= rm -f
CP		= cp -f
MV		= mv -f
WPATH 		= ~/bin

F90		= gfortran
FFLAG90         = -O4 -march=native -ftree-vectorize
FFLAG90         = -g -fbounds-check -Wuninitialized -O -ftrapv \
		-fimplicit-none -fno-automatic
#F90		= ifort
#FFLAG90		= -O3 -fast
#FFLAG90         = -C -g -debug all -check all -implicitnone \
		-warn unused -fp-stack-check -heap-arrays -ftrapuv \
		-check pointers -check bounds

FFLAGMPI        = -I/usr/include/lam
FFLAGMPI        = 
FLIBMPI         = -L/usr/lib/lam/lib -llammpio -llamf77mpi -lmpi -llam -lutil -ldl -lnsl
FLIBMPI         = 
FLIB            = -lm
FLIBF77         = -lm

# definition der default targets..
#  BLAS und LAPACK tools
LALIB		= -llapack -lblas
LALIB		= 
# das hier chek obs ein bin im home gibt
C1		= cbn
# macht CRTomo
PR1		= crt
# macht CRMod
PR2		= crm
# macht CutMckee
PR3		= ctm

MACHINE		= $(shell uname -n)
################################################################
# default
all:		$(C1) $(PR1) $(PR2) $(PR3) install
################################################################
# this is for evry one here
ferr		= get_error.o
# CRTomo objects
f90crt		= alloci.o femmod.o \
		  datmod.o invmod.o cjgmod.o sigmamod.o \
		  electrmod.o modelmod.o elemmod.o wavenmod.o \
		  randbmod.o errmod.o konvmod.o pathmod.o \
		  invhpmod.o

f90crtsub	= bbsedc.o bbsens.o besp_elem.o \
		  bessi0.o bessi1.o bessk0.o bessk1.o \
		  beta.o bkfak.o blam0.o bnachbar.o bpar.o bpot.o \
		  bsendc.o bsens.o bsensi.o \
		  bsigm0.o bsigma.o bsytop.o bvolt.o bvolti.o chareal.o \
		  chkpol.o choldc.o chol.o dmisft.o \
		  elem1.o elem3.o elem4.o \
		  elem5.o elem8.o filpat.o \
		  chold.o cholz.o linvd.o linvz.o \
		  mdian1.o parfit.o potana.o precal.o rall.o \
		  gammln.o gaulag.o gauleg.o intcha.o kompab.o \
		  kompadc.o kompbdc.o kompb.o kont1.o kont2.o \
		  randb2.o randb.o randdc.o rdati.o rdatm.o \
		  refsig.o relectr.o relem.o rrandb.o rsigma.o \
		  rtrafo.o rwaven.o scalab.o scaldc.o sort.o \
		  update.o vredc.o vre.o wdatm.o \
		  wkpot.o wout.o wpot.o wsens.o \
		  gauss_dble.o gauss_cmplx.o get_unit.o \
		  make_noise.o tic_toc.o variomodel.o bvariogram.o \
		  cg_mod.o bsmatm_mod.o bmcm_mod.o brough_mod.o \

fcrt		= inv.o

forcrt		= 
# CRMod objects
f90crm		= alloci.o femmod.o datmod.o \
		  invmod.o sigmamod.o electrmod.o modelmod.o \
		  elemmod.o wavenmod.o randbmod.o errmod.o konvmod.o \
		  pathmod.o

fcrm		= fem.o

forcrm		= 

f90crmsub	= bbsens.o bessi0.o bessi1.o bessk0.o bessk1.o \
		  beta.o bkfak.o bpot.o \
		  bsendc.o bsens.o bsensi.o \
		  bsytop.o bvolt.o bvolti.o chareal.o \
		  chkpol.o choldc.o chol.o \
		  elem1.o elem3.o elem4.o \
		  elem5.o elem8.o filpat.o \
		  kompadc.o kompbdc.o kompb.o \
		  potana.o precal.o \
		  randb2.o randb.o randdc.o rdatm.o \
		  relectr.o relem.o rrandb.o rsigma.o refsig.o \
		  rtrafo.o rwaven.o scalab.o scaldc.o sort.o \
		  vredc.o vre.o wdatm.o \
		  wkpot.o wout.o wpot.o wsens.o \
		  gammln.o gaulag.o gauleg.o intcha.o kompab.o \
		  tic_toc.o make_noise.o get_unit.o 
################################################################
# rules
%.o:		%.for
		$(F90) $(FFLAG90) -c $<

#$(forcrt):	%.o : %.for
#		$(F90) $(FFLAG90) -c $<

$(fcrt):	%.o : %.f90
		$(F90) $(FFLAG90) -c $<

$(f90crt):	%.o : %.f90		
		$(F90) $(FFLAG90) -c $<

$(f90crtsub):	%.o : %.f90		
		$(F90) $(FFLAG90) -c $<

$(fcrm):	%.o : %.f		
		$(F90) $(FFLAG90) -c $<

$(ferr):	error.txt get_error.f90
		./make_crerr.sh
		$(F90) $(FFLAG90) -c get_error.f90

################################################################
# Dependencies

besp_elem.o:	alloci.o elemmod.o

bvariogram.o:	invmod.o

bmcm_mod.o:	tic_toc.o alloci.o femmod.o elemmod.o invmod.o \
		errmod.o konvmod.o modelmod.o datmod.o sigmamod.o \
		pathmod.o

bnachbar.o:	alloci.o modelmod.o elemmod.o errmod.o konvmod.o

bsmatm_mod.o:	tic_toc.o alloci.o femmod.o elemmod.o invmod.o \
		errmod.o konvmod.o modelmod.o datmod.o sigmamod.o \
		pathmod.o variomodel.o

bvariogram.o:	invmod.o variomodel.o sigmamod.o modelmod.o elemmod.o \
		errmod.o konvmod.o

cjgmod.o:	modelmod.o datmod.o

cg_mod.o:	cjgmod.o alloci.o femmod.o elemmod.o invmod.o errmod.o \
		konvmod.o modelmod.o datmod.o

brough_mod.o:	alloci.o invmod.o konvmod.o modelmod.o elemmod.o \
		errmod.o datmod.o

kont1.o:	variomodel.o

rall.o:		make_noise.o variomodel.o

update.o:	cg_mod.o

inv.o:		$(f90crt) $(forcrt) $(f90crtsub)
###############################################################
.SILENT:	cbn
###################################
LALIB:		./libla/%.f	
		make -C libla

cbn:		
		echo "Pruefe ~/bin"
		if [ -d ~/bin ]; then \
			echo "ok"; \
		else \
			echo "Du hast kein bin in deinem home.--"; \
			mkdir ~/bin; \
		fi

crt:		$(C1) $(f90crt) $(f90crtsub) $(forcrt) $(fcrt) $(ferr)
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRTomo \
		$(f90crt) $(f90crtsub) $(forcrt) $(fcrt) $(ferr) $(LALIB)
		$(CP) CRTomo $(WPATH)/CRTomo_$(MACHINE) 

crm:		$(C1) $(f90crm) $(f90crmsub) $(forcrm) $(fcrm) $(ferr)
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRMod \
		$(f90crm) $(f90crmsub) $(forcrm) $(fcrm) $(ferr) $(LALIB)
		$(CP) CRMod $(WPATH)/CRMod_$(MACHINE)

ctm:		
		cd ./cutmckee ; make

install:	$(C1) $(crt) $(crm)				
		$(CP) CRTomo $(WPATH)/CRTomo_$(MACHINE)
		$(CP) CRMod $(WPATH)/CRMod_$(MACHINE)
		cd ./cutmckee ; make install

clean:		
		$(RM) CRTomo CRMod *~ *.mod *.o
		cd ./cutmckee ; make clean
