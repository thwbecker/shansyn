# output directories
ARCH=x86_64
#F77 = ifort 
#FFLAGS = -O2 -extend_source -cpp

F77 = /usr/bin/gfortran
FFLAGS = -O -cpp -ffixed-line-length-132
#FFLAGS = -g -cpp -ffixed-line-length-132
#FFLAGS =  -cpp -ffixed-line-length-132 -static
LDFLAGS = -llapack -lblas
#LDFLAGS = 

BDIR = ../bin/
ODIR = ./objects/

HDR_FLS = prec.h


# binaries
BINS = $(BDIR)/xyz2aniso $(BDIR)/xyz2gsh $(BDIR)/aniso2xyz $(BDIR)/gsh2xyz

# library
LIBS = $(ODIR)/libgsh.a
GSHL_LIB_OBJS = $(ODIR)/ylmv4.o $(ODIR)/rotmx2.o $(ODIR)/choles.o



all: obj_dirs $(LIBS) $(BINS)


clean:
	rm $(ODIR)/*.o
dist_clean:
	rm $(BINS)

# archive with auxiliary routines
$(ODIR)/libgsh.a: $(GSHL_LIB_OBJS)
	$(AR) -rv $(ODIR)/libgsh.a $(GSHL_LIB_OBJS)


# directories
obj_dirs:
	if [ ! -s ./$(ARCH)/ ]; then\
		mkdir $(ARCH);\
	fi;

# object rule
$(ODIR)/%.o: %.f   $(HDR_FLS)
	$(F77) $(FFLAGS) -c $< -o $(ODIR)/$*.o

# binary rule
$(BDIR)/%: %.f   $(HDR_FLS) $(LIBS)
	$(F77) $(FFLAGS) $< -o $(BDIR)/$* -L$(ODIR)/ -lgsh $(LDFLAGS)

