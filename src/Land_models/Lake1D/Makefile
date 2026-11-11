# Makefile
#
# 1-D lake model (module_sf_lake) and its WRF-Hydro driver (module_lake1d),
# built with the offline Noah-MP model (compile_offline_NoahMP.sh).
.SUFFIXES:
.SUFFIXES: .o .F

include ../NoahMP/user_build_options

OBJS = \
	module_sf_lake.o \
	module_lake1d.o

CPPHRLDAS = -D_HRLDAS_OFFLINE_

all:	$(OBJS)

.F.o:
	@echo ""
	$(COMPILERF90) $(CPPINVOKE) $(CPPFLAGS) $(CPPHRLDAS) -o $(@) -c $(F90FLAGS) $(FREESOURCE) -I. \
	-I../NoahMP/Utility_routines -I../../mod -I../../MPP $(NETCDFMOD) $(*).F
	@echo ""

#
# Dependencies:
#
module_lake1d.o: module_sf_lake.o

#
# This command cleans up object (etc) files:
#
clean:
	rm -f *.o *.mod *.stb *~
