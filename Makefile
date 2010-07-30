# $Id: Makefile,v 1.0 2009/07/02 12:43:39 roland Exp $
#
#      Makefile for CRTomo (default)
#       

RM		= rm -f
CP		= cp -f
MV		= mv -f
WPATH 		= ~/bin

F90		= gfortran
F77		= gfortran
FFLAG90         = -O3 -march=native -ftree-vectorize -fexpensive-optimizations -ffast-math
FFLAG90         = -Wconversion -Wunderflow -fbacktrace
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

# kopiert die matlab tools
PRM		= mtools
MACHINE		= $(shell uname -n)
################################################################
# default
all:		$(C1) $(PR1) $(PR2) $(PR3) $(PRM) install
################################################################
# this is for evry one here
ferr		= get_error.o
# CRTomo objects
f90crt		= alloci.o gauss_dble.o gauss_cmplx.o get_unit.o \
		  make_noise.o tic_toc.o variomodel.o \
		  chold.o cholz.o linvd.o linvz.o femmod.o \
		  datmod.o invmod.o cjgmod.o sigmamod.o \
		  electrmod.o modelmod.o elemmod.o wavenmod.o \
		  randbmod.o errmod.o konvmod.o cg_mod.o

fcrt		= inv.o
forcrt		= bbsedc.o bbsens.o besp_elem.o bessi0.o bessi1.o \
		  bessk0.o bessk1.o beta.o bkfak.o blam0.o bnachbar.o \
		  bpot.o brough.o broughsto.o broughtri.o bsendc.o \
		  bsens.o bsensi.o bsigm0.o bsmatm.o bsmatmsto.o \
		  bsmatmtri.o bvolt.o bvolti.o chareal.o chkpol.o \
		  choldc.o chol.o dmisft.o elem1.o \
		  elem3.o elem4.o elem5.o elem8.o filpat.o findinv.o \
		  gammln.o gaulag.o gauleg.o intcha.o kompab.o \
		  kompadc.o kompbdc.o kompb.o kont1.o kont2.o \
		  mdian1.o parfit.o potana.o precal.o rall.o \
		  randb2.o randb.o randdc.o rdati.o rdatm.o \
		  refsig.o relectr.o relem.o rrandb.o rsigma.o \
		  rtrafo.o rwaven.o scalab.o scaldc.o sort.o \
		  update.o vredc.o vre.o wdatm.o wkpot.o wout.o \
		  wpot.o wsens.o bsmatmmgs.o bsytop.o bsmatmtv.o \
		  bata_dc.o bata_reg_dc.o bmcm_dc.o bmcm2_dc.o \
		  bres_dc.o bata.o bata_reg.o bmcm.o bmcm2.o bres.o \
		  bsmatmlma.o broughlma.o \
		  bvariogram.o bpar.o bsigma.o
# CRMod objects
f90crm		= alloci.o tic_toc.o make_noise.o femmod.o datmod.o \
		  invmod.o sigmamod.o electrmod.o modelmod.o \
		  elemmod.o wavenmod.o randbmod.o errmod.o konvmod.o
fcrm		= fem.o
forcrm		= bbsens.o bessi0.o bessi1.o \
		  bessk0.o bessk1.o bkfak.o beta.o bpot.o \
		  bsendc.o bsens.o bsensi.o \
		  bvolt.o bvolti.o chareal.o chkpol.o \
		  choldc.o chol.o elem1.o bsytop.o \
		  elem3.o elem4.o elem5.o elem8.o filpat.o \
		  gammln.o gaulag.o gauleg.o intcha.o kompab.o \
		  kompadc.o kompbdc.o kompb.o \
		  potana.o precal.o \
		  randb2.o randb.o randdc.o rdatm.o \
		  relectr.o relem.o rrandb.o rsigma.o refsig.o \
		  rtrafo.o rwaven.o scalab.o scaldc.o sort.o \
		  vredc.o vre.o wdatm.o wkpot.o wout.o \
		  wpot.o wsens.o get_unit.o
################################################################
# rules
%.o:		%.for
		$(F90) $(FFLAG90) -c $<

#$(forcrt):	%.o : %.for
#		$(F90) $(FFLAG90) -c $<

$(fcrt):	%.o : %.f
		$(F90) $(FFLAG90) -c $<
$(f90crt):	%.o : %.f90		
		$(F90) $(FFLAG90) -c $<
$(fcrm):	%.o : %.f		
		$(F90) $(FFLAG90) -c $<
$(ferr):	error.txt get_error.f90
		./make_crerr.sh
		$(F90) $(FFLAG90) -c get_error.f90
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

crt:		$(C1) $(f90crt) $(forcrt) $(fcrt) $(ferr)
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRTomo \
		$(f90crt) $(forcrt) $(fcrt) $(ferr) $(LALIB)
		$(CP) CRTomo $(WPATH)/CRTomo_$(MACHINE) 

crm:		$(C1) $(f90crm) $(forcrm) $(fcrm) $(ferr)
		$(F90) $(FFLAG90) $(FFLAGMPI) -o CRMod \
		$(f90crm) $(forcrm) $(fcrm) $(ferr) $(LALIB)
		$(CP) CRMod $(WPATH)/CRMod_$(MACHINE)

mtools:
		cd ./m_tools ; make

ctm:		
		cd ./cutmckee ; make

install:	$(C1) $(crt) $(crm)				
		$(CP) CRTomo $(WPATH)/CRTomo_$(MACHINE)
		$(CP) CRMod $(WPATH)/CRMod_$(MACHINE)
		cd ./m_tools ; make install
		cd ./cutmckee ; make install

clean:		
		$(RM) CRTomo CRMod *~ *.mod *.o
		cd ./m_tools ; make clean
		cd ./cutmckee ; make clean
