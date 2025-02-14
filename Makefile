#
# this file holds all machine dependent compiler flags and 
# directories for GMT source code and libraries, most modifications
# should be done in this file or in environment variables such as CC
#

ARCH=$(shell uname -m | awk '{print(tolower($$1))}')
-include machine_dependent.$(ARCH).m

#
# these are settings for GMT version < 5
#
#GMT_VERSION_OPTION = -DUSE_GMT4
#GMT_INC = -I$(GMT4HOME)/include/ -I$(NETCDFDIR)/include 
#GMT_LIB = -L$(GMT4HOME)/lib -lgmt -lpsl
#
# these should be defined for newer version >= 5
#
GMT_VERSION_OPTION =
GMT_INC = -I$(GMTHOME)/include/gmt/  -I/usr/include/gdal/  #-I$(NETCDFDIR)/include # gmt-config --cflags
GMT_LIB = -L$(GMTHOME)/lib/ -lgmt # gmt-config --libes
#
# netcdf should work with both
#
#NETCDF_LIB = -L$(NETCDFDIR)/lib/ -lnetcdf
NETCDF_LIB = -L/usr/lib64/mpich/lib/ -lnetcdf

#

#
#
#
#CFLAGS = -g -DDEBUG
#FFLAGS = -g -DDEBUG
#
#CFLAGS = -fast
#FFLAGS = -fast 
#  makefile for the shansyn package: spherical harmonic analysis
# 
#  $Id: tmp.dat,v 1.1 2009/11/25 21:34:48 twb Exp twb $
#
#  for license, see COPYRIGHT
#
FORMATS =
# options for the programs
#
OPTIONS = -DBE_VERBOSE

#
SHELL = /bin/sh

#
# object files and libraries
ODIR = objects/$(ARCH)/
# binaries
BDIR = bin/$(ARCH)/
#
INCLUDES = $(GMT_INC) -I$(NETCDFDIR)/include $(GMT_VERSION_OPTION)
HFILES = abconvert.h   shana.h spear.h determine_coeff.h  precision.h  shsyn.h  shansyn.h auto_proto.h
#	vectmath.h


#
# spherical hamornic synthesis and analysis
#

SHANA_OBJS =   $(ODIR)/shana.o   $(ODIR)/gauleg.o   $(ODIR)/plgndr.o      \
	$(ODIR)/spline.o $(ODIR)/misc.o 

SHSYN_OBJS =   $(ODIR)/shsyn.o   $(ODIR)/gauleg.o   $(ODIR)/plgndr.o  \
	$(ODIR)/powercorr.o   $(ODIR)/dfour1.o   $(ODIR)/mygrdio.o $(ODIR)/misc.o \
	$(ODIR)/spear.o $(ODIR)/nr_utils.o

# broken, need to fix again
SHSYN_XS_OBJS = $(ODIR)/shsyn-xsection.o $(ODIR)/gauleg.o $(ODIR)/plgndr.o  \
	$(ODIR)/dfour1.o 

#
# ab model format stuff
#
ABMODEL_OBJ =   $(ODIR)/powercorr.o         $(ODIR)/write_coefficients.o      $(ODIR)/chebyshev.o  \
		$(ODIR)/numrec_svd.o        $(ODIR)/interpolate_she_model.o   $(ODIR)/read_she_model.o    \
                $(ODIR)/determine_coeff.o   $(ODIR)/splinesc.o                $(ODIR)/splinesf.o \
		$(ODIR)/spear.o 	    $(ODIR)/nr_utils.o  	      $(ODIR)/select_lms.o \
		$(ODIR)/gsh_handling.o

#
# double precision binaries
PROGRAMS = $(BDIR)/shana $(BDIR)/shsyn $(BDIR)/plotlegendre \
	$(BDIR)/abconvert $(BDIR)/cmodellinreg $(BDIR)/calc_pt_corr \
	$(BDIR)/cmodelcorr $(BDIR)/cmodelcorr_gsh  \
	$(BDIR)/abadd $(BDIR)/model2scatter $(BDIR)/cmodelmeancorr \
	$(BDIR)/calc_radial_corr \
	$(BDIR)/extract_layer $(BDIR)/extract_layer_gsh \
	$(BDIR)/cmodelpower $(BDIR)/cmodelpower_gsh \
	$(BDIR)/ab2centroid $(BDIR)/modmodellmax $(BDIR)/ab2scatter \
	$(BDIR)/extract_model_depths 	$(BDIR)/extract_model_depths_gsh \
	$(BDIR)/cradialcorr $(BDIR)/cradialcorr_gsh \
	$(BDIR)/scale_model $(BDIR)/mod_modelbase \
	$(BDIR)/convert_she_model $(BDIR)/cmradialcorr 

# broken, need to fix
#	$(BDIR)/shsyn-xsection


# 
# targets
#
#
all: obj_directories  $(PROGRAMS)


install: all 

auto_proto.h:
	rm auto_proto.h;\
	cproto $(INCLUDES) -f2 *.c | grep -v main > auto_proto.h

obj_directories:
	if [ ! -s ./objects/ ]; then\
		mkdir objects;\
	fi;
	if [ ! -s $(ODIR) ];then \
		mkdir $(ODIR);\
	fi;\
	if [ ! -s ./bin/ ];then\
		mkdir bin;\
	fi;\
	if [ ! -s bin/$(ARCH)/ ];then \
		mkdir bin/$(ARCH);\
	fi;

clean: 
	rm $(ODIR)/*.o $(ODIR)/*.a



#
# individual programs
#

$(BDIR)/shana: $(SHANA_OBJS) 
	$(CC) $(SHANA_OBJS)  $(INCLUDES) \
		$(GMT_LIB) $(NETCDF_LIB) -o $@  $(MATHLIBS)  $(LDFLAGS) 


$(BDIR)/shana_s: $(SHANA_S_OBJS) 
	$(CC)  $(SHANA_OBJS)  $(INCLUDES) \
		$(GMT_LIB) $(NETCDF_LIB) -o $@  $(LDFLAGS)   $(MATHLIBS)

$(BDIR)/shsyn: $(SHSYN_OBJS) 
	$(CC)  $(SHSYN_OBJS)  $(INCLUDES) \
		$(GMT_LIB) $(NETCDF_LIB)  -o $@ $(LDFLAGS) 

$(BDIR)/shsyn-xsection: $(SHSYN_XS_OBJS) $(ABMODEL_OBJ)
	$(F77) $(SHSYN_XS_OBJS) $(ABMODEL_OBJ) $(INCLUDES) \
		$(GMT_LIB) $(NETCDF_LIB)  -o $@ $(FTRN_LIB) $(LDFLAGS) 


$(BDIR)/test: $(ODIR)/test.o $(ODIR)/mygrdio.o
	$(CC)  $(ODIR)/test.o $(ODIR)/mygrdio.o \
		$(INCLUDES) \
		$(GMT_LIB) $(NETCDF_LIB)  -o $@ $(LDFLAGS) 

$(BDIR)/shsyn_s: $(SHSYN_S_OBJS) 
	$(CC)  $(SHSYN_OBJS)  $(INCLUDES) \
		$(GMT_LIB) $(NETCDF_LIB)  -o $@ $(LDFLAGS) 

$(BDIR)/plotlegendre: $(ODIR)/plotlegendre.o $(ODIR)/plgndr.o $(ODIR)/gauleg.o
	$(CC) $(ODIR)/plotlegendre.o  $(ODIR)/plgndr.o $(ODIR)/gauleg.o \
	$(INCLUDES) $(LDFLAGS) -o $@ 


$(BDIR)/abconvert: $(ODIR)/abconvert.o $(ODIR)/select_lms.o \
	$(ODIR)/powercorr.o  	$(ODIR)/gsh_handling.o  \
	$(ODIR)/nr_utils.o \
	$(ODIR)/write_coefficients.o $(ODIR)/rand.o $(ODIR)/spear.o
	$(CC) $(INCLUDES)  $(ODIR)/rand.o $(ODIR)/write_coefficients.o $(ODIR)/powercorr.o  \
	$(ODIR)/abconvert.o $(ODIR)/nr_utils.o  \
	$(ODIR)/select_lms.o $(ODIR)/gsh_handling.o \
	$(ODIR)/spear.o -o $@ $(LDFLAGS)

$(BDIR)/calc_pt_corr: $(ODIR)/calc_pt_corr.o $(ODIR)/select_lms.o \
	$(ODIR)/powercorr.o  	$(ODIR)/gsh_handling.o  $(ODIR)/nr_utils.o \
	$(ODIR)/write_coefficients.o $(ODIR)/rand.o $(ODIR)/spear.o
	$(CC) $(INCLUDES)  $(ODIR)/powercorr.o  $(ODIR)/nr_utils.o \
	$(ODIR)/calc_pt_corr.o  $(ODIR)/select_lms.o $(ODIR)/gsh_handling.o \
	$(ODIR)/spear.o -o $@ $(LDFLAGS) -lz

$(BDIR)/calc_radial_corr: $(ODIR)/calc_radial_corr.o $(ABMODEL_OBJ) 
	$(CC) $(INCLUDES)  $(ABMODEL_OBJ)   \
	$(ODIR)/calc_radial_corr.o -o $@ $(FTRN_LIB) $(LDFLAGS)


$(BDIR)/abadd: abadd.c
	$(CC)  abadd.c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		 -o $@

$(BDIR)/cmodelcorr: $(ABMODEL_OBJ) $(ODIR)/cmodelcorr.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/cmodelcorr.o -o $@ $(ODIR)/myopen.o $(FTRN_LIB) $(LDFLAGS) 

$(BDIR)/cmodelcorr_gsh: $(ABMODEL_OBJ) $(ODIR)/cmodelcorr_gsh.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/cmodelcorr_gsh.o -o $@ $(ODIR)/myopen.o $(FTRN_LIB) $(LDFLAGS) 

$(BDIR)/cmodelmeancorr: $(ABMODEL_OBJ) $(ODIR)/cmodelmeancorr.o
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/cmodelmeancorr.o -o $@ $(FTRN_LIB) $(LDFLAGS)  $(ODIR)/myopen.o 

$(BDIR)/model2scatter: $(ABMODEL_OBJ) $(ODIR)/model2scatter.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/model2scatter.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o 


$(BDIR)/mod_modelbase: $(ABMODEL_OBJ) $(ODIR)/mod_modelbase.o  $(ODIR)/myopen.o
	$(CC) $(INCLUDES)  \
	$(ABMODEL_OBJ)  $(ODIR)/mod_modelbase.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/cmodellinreg: $(ABMODEL_OBJ) $(ODIR)/cmodellinreg.o  $(ODIR)/fitxy.o  \
	$(ODIR)/sphex_lin_reg.o $(ODIR)/myopen.o
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/fitxy.o $(ODIR)/sphex_lin_reg.o $(ODIR)/cmodellinreg.o \
	-o $@ $(ODIR)/myopen.o $(FTRN_LIB) $(LDFLAGS) 

$(BDIR)/scale_model: $(ABMODEL_OBJ) $(ODIR)/scale_model.o
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/scale_model.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/cradialcorr: $(ABMODEL_OBJ) $(ODIR)/cradialcorr.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/cradialcorr.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/cradialcorr_gsh: $(ABMODEL_OBJ) $(ODIR)/cradialcorr_gsh.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/cradialcorr_gsh.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/cmradialcorr: $(ABMODEL_OBJ) $(ODIR)/cmradialcorr.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/cmradialcorr.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/extract_model_depths: $(ABMODEL_OBJ) $(ODIR)/extract_model_depths.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/extract_model_depths.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/extract_model_depths_gsh: $(ABMODEL_OBJ) $(ODIR)/extract_model_depths_gsh.o  
	$(CC) $(INCLUDES)   \
	$(ABMODEL_OBJ) $(ODIR)/extract_model_depths_gsh.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/ab2centroid: $(ODIR)/ab2centroid.o $(ODIR)/plgndr.o $(ODIR)/gauleg.o 
	$(CC) $(INCLUDES) $(ODIR)/plgndr.o  $(ODIR)/gauleg.o $(ODIR)/ab2centroid.o  \
	-o $@ $(LDFLAGS) 

$(BDIR)/ab2scatter: ab2scatter.c
	$(CC) $(CFLAGS)  $(INCLUDES)  ab2scatter.c  -o $@ $(LDFLAGS)

$(BDIR)/modmodellmax: $(ODIR)/modmodellmax.o $(ABMODEL_OBJ) $(ODIR)/modmodellmax.o
	$(CC) $(INCLUDES)  $(ODIR)/modmodellmax.o \
	$(ABMODEL_OBJ) -o $@ $(ODIR)/myopen.o $(FTRN_LIB) $(LDFLAGS)

$(BDIR)/cmodelpower: $(ODIR)/cmodelpower.o  $(ABMODEL_OBJ) 
	$(CC)   $(INCLUDES)  $(ABMODEL_OBJ) \
	$(ODIR)/cmodelpower.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/cmodelpower_gsh: $(ODIR)/cmodelpower_gsh.o  $(ABMODEL_OBJ) 
	$(CC)   $(INCLUDES)  $(ABMODEL_OBJ) \
	$(ODIR)/cmodelpower_gsh.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

$(BDIR)/extract_layer: $(ODIR)/extract_layer.o $(ABMODEL_OBJ) 
	$(CC) $(INCLUDES)   \
		$(ABMODEL_OBJ) $(ODIR)/extract_layer.o -o $@ $(FTRN_LIB) $(LDFLAGS)  $(ODIR)/myopen.o

$(BDIR)/extract_layer_gsh: $(ODIR)/extract_layer_gsh.o $(ABMODEL_OBJ) 
	$(CC) $(INCLUDES)   \
		$(ABMODEL_OBJ) $(ODIR)/extract_layer_gsh.o -o $@ $(FTRN_LIB) $(LDFLAGS)  $(ODIR)/myopen.o

$(BDIR)/convert_she_model: $(ODIR)/convert_she_model.o $(ABMODEL_OBJ) 
	$(CC) $(INCLUDES)   \
		$(ABMODEL_OBJ) $(ODIR)/convert_she_model.o -o $@ $(FTRN_LIB) $(LDFLAGS) $(ODIR)/myopen.o

#
#
# object files
#
#
$(ODIR)/shana.o: shana.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(LAPACK_DEFINE) $(GMT_INC) $(OPTIONS) $(FORMATS) shana.c -o $@


$(ODIR)/shsyn.o: shsyn.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(GMT_INC) $(OPTIONS) $(FORMATS) shsyn.c -o $@

$(ODIR)/shsyn-xsection.o: shsyn-xsection.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(GMT_INC) $(OPTIONS) $(FORMATS) shsyn-xsection.c -o $@


$(ODIR)/test.o: test.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(GMT_INC) $(OPTIONS) $(FORMATS) test.c -o $@


$(ODIR)/plotlegendre.o: plotlegendre.c  $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) plotlegendre.c -o $@

$(ODIR)/plotlegendre_s.o: plotlegendre.c  $(HFILES)
	$(CC) -c plotlegendre.c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/plotlegendre_s.o

$(ODIR)/abconvert.o: abconvert.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) abconvert.c -o $@

$(ODIR)/calc_pt_corr.o: calc_pt_corr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) calc_pt_corr.c -o $@

$(ODIR)/calc_radial_corr.o: calc_radial_corr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) calc_radial_corr.c -o $@ 


$(ODIR)/abconvert_s.o: abconvert.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) abconvert.c -o $(ODIR)/abconvert_s.o

$(ODIR)/cradialcorr.o: cradialcorr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  \
		$(FORMATS) cradialcorr.c -o $@

$(ODIR)/cradialcorr_gsh.o: cradialcorr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  \
		$(FORMATS) cradialcorr.c -DSHANA_EXPECT_GSH -o $@

$(ODIR)/cmradialcorr.o: cmradialcorr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  \
		$(FORMATS) cmradialcorr.c -o $@

$(ODIR)/cmodelcorr.o: cmodelcorr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  \
		$(FORMATS) cmodelcorr.c -o $@ 

$(ODIR)/cmodelcorr_gsh.o: cmodelcorr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  \
		$(FORMATS) cmodelcorr.c -DSHANA_EXPECT_GSH -o $@ 

$(ODIR)/model2scatter.o: model2scatter.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  \
		$(FORMATS) model2scatter.c -o $@

$(ODIR)/cmodelcorr_s.o: cmodelcorr.c $(HFILES)
	$(CC) -c cmodelcorr.c $(CFLAGS) $(INCLUDES)  \
	-DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/cmodelcorr_s.o

$(ODIR)/cmodellinreg.o: cmodellinreg.c $(HFILES)
	$(CC) -c cmodellinreg.c $(CFLAGS) $(INCLUDES)   $(OPTIONS) $(FORMATS) -o $@

$(ODIR)/cmodelpower.o: cmodelpower.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) cmodelpower.c -o $@
$(ODIR)/cmodelpower_gsh.o: cmodelpower.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) -DSHANA_EXPECT_GSH cmodelpower.c -o $@

$(ODIR)/mod_modelbase.o: mod_modelbase.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) mod_modelbase.c -o $@

$(ODIR)/cmodelpower_s.o: cmodelpower.c $(HFILES)
	$(CC) -c cmodelpower.c $(CFLAGS) $(INCLUDES)  -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/cmodelpower_s.o

$(ODIR)/plgndr.o: plgndr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) plgndr.c -o $@

$(ODIR)/plgndr_s.o: plgndr.c $(HFILES)
	$(CC) -c plgndr.c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/plgndr_s.o

$(ODIR)/sphex_lin_reg.o: sphex_lin_reg.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) sphex_lin_reg.c -o $@

$(ODIR)/sphex_lin_reg_s.o: sphex_lin_reg.c $(HFILES)
	$(CC) -c sphex_lin_reg.c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/sphex_lin_reg_s.o

$(ODIR)/fitxy.o: fitxy.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) fitxy.c -o $@

$(ODIR)/fitxy_s.o: fitxy.c $(HFILES)
	$(CC) -c fitxy.c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/fitxy_s.o

$(ODIR)/powercorr.o: powercorr.c	$(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) \
	$(FORMATS) powercorr.c -o $@

$(ODIR)/powercorr_s.o: powercorr.c	$(HFILES)
	$(CC) -c powercorr.c $(CFLAGS) $(INCLUDES) \
	-DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/powercorr_s.o

$(ODIR)/dfour1.o: dfour1.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) dfour1.c -o $@

$(ODIR)/dfour1_s.o: dfour1.c $(HFILES)
	$(CC) -c dfour1.c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/dfour1_s.o

$(ODIR)/mygrdio.o: mygrdio.c $(HFILES) 
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) $(GMT_INC) mygrdio.c -o $@

$(ODIR)/mygrdio_s.o: mygrdio.c $(HFILES)
	$(CC) -c mygrdio.c $(CFLAGS) $(INCLUDES) $(GMT_INC) \
	-DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/mygrdio_s.o

$(ODIR)/nr_utils.o: nr_utils.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) nr_utils.c -o $@

$(ODIR)/nr_utils_s.o: nr_utils.c $(HFILES)
	$(CC) -c nr_utils.c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/nr_utils_s.o

$(ODIR)/gauleg.o: gauleg.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) gauleg.c -o $@

$(ODIR)/gauleg_s.o: gauleg.c $(HFILES)
	$(CC) -c gauleg.c $(CFLAGS) $(INCLUDES) -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/gauleg_s.o

$(ODIR)/spline.o: spline.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(GMT_INC) $(FORMATS) spline.c -o $@

$(ODIR)/spline_s.o: spline.c $(HFILES)
	$(CC) -c spline.c $(CFLAGS) $(INCLUDES) $(GMT_INC)  -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/spline_s.o

$(ODIR)/spear.o: spear.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(GMT_INC) $(FORMATS) spear.c -o $@

$(ODIR)/spear_s.o: spear.c $(HFILES)
	$(CC) -c spear.c $(CFLAGS) $(INCLUDES) $(GMT_INC)  -DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/spear_s.o

$(ODIR)/interpolate_she_model.o: interpolate_she_model.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) interpolate_she_model.c -o $@

$(ODIR)/interpolate_she_model_s.o: interpolate_she_model.c $(HFILES)
	$(CC) -c interpolate_she_model.c $(CFLAGS) $(INCLUDES) \
	-DSINGLE_PREC $(OPTIONS) $(FORMATS) -o $(ODIR)/interpolate_she_model_s.o

$(ODIR)/read_she_model.o: read_she_model.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) read_she_model.c -o $@

$(ODIR)/scale_model.o: scale_model.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) scale_model.c -o $@

$(ODIR)/modmodellmax.o: modmodellmax.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) modmodellmax.c -o $@

$(ODIR)/cmodelmeancorr.o: cmodelmeancorr.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) cmodelmeancorr.c -o $@

$(ODIR)/read_she_model_s.o: read_she_model.c $(HFILES)
	$(CC) -c read_she_model.c $(CFLAGS) $(INCLUDES) \
	-DSINGLE_PREC $(OPTIONS) $(FORMATS)  -o $(ODIR)/read_she_model_s.o

$(ODIR)/extract_layer.o: extract_layer.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) extract_layer.c -o $@

$(ODIR)/extract_layer_gsh.o: extract_layer.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) -DSHANA_EXPECT_GSH extract_layer.c -o $@

$(ODIR)/convert_she_model.o: convert_she_model.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS)  $(FORMATS) convert_she_model.c -o $@

$(ODIR)/extract_layer_s.o: extract_layer.c $(HFILES)
	$(CC) -c extract_layer.c $(CFLAGS) $(INCLUDES) \
	-DSINGLE_PREC $(OPTIONS) $(FORMATS)  -o $(ODIR)/extract_layer_s.o

$(ODIR)/ab2centroid.o: ab2centroid.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		ab2centroid.c	-o $(ODIR)/ab2centroid.o

$(ODIR)/ab2centroid_s.o: ab2centroid.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		-DSINGLE_PREC ab2centroid.c -o $(ODIR)/ab2centroid_s.o

$(ODIR)/write_coefficients.o: write_coefficients.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		write_coefficients.c -o $(ODIR)/write_coefficients.o

$(ODIR)/write_coefficients_s.o: write_coefficients.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		-DSINGLE_PREC write_coefficients.c \
	-o $(ODIR)/write_coefficients_s.o

$(ODIR)/extract_model_depths.o: extract_model_depths.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		extract_model_depths.c -o $@

$(ODIR)/extract_model_depths_gsh.o: extract_model_depths.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		extract_model_depths.c -DSHANA_EXPECT_GSH -o $@

$(ODIR)/chebyshev.o: chebyshev.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		chebyshev.c -o $@

$(ODIR)/chebyshev_s.o: chebyshev.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		-DSINGLE_PREC chebyshev.c -o $(ODIR)/chebyshev_s.o

$(ODIR)/splinesc.o: splinesc.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		splinesc.c -o $@

$(ODIR)/splinesc_s.o: splinesc.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		-DSINGLE_PREC splinesc.c -o $(ODIR)/splinesc_s.o

$(ODIR)/misc.o: misc.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		misc.c -o $@

$(ODIR)/misc_s.o: misc.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		-DSINGLE_PREC misc.c -o $(ODIR)/misc_s.o

$(ODIR)/determine_coeff.o: determine_coeff.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		determine_coeff.c -o $@

$(ODIR)/determine_coeff_s.o: determine_coeff.c $(HFILES)
	$(CC) -c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		-DSINGLE_PREC determine_coeff.c -o $(ODIR)/determine_coeff_s.o

$(ODIR)/rand_s.o: rand.c $(HFILES)
	$(CC) -c rand.c -DSINGLE_PREC $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		 -o $(ODIR)/rand_s.o

$(ODIR)/myopen.o: myopen.c $(HFILES)
	$(CC) -c myopen.c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		 -o $(ODIR)/myopen.o

$(ODIR)/rand.o: rand.c $(HFILES)
	$(CC) -c rand.c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		 -o $(ODIR)/rand.o

$(ODIR)/select_lms.o: select_lms.c $(HFILES)
	$(CC) -c select_lms.c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		 -o $(ODIR)/select_lms.o

$(ODIR)/gsh_handling.o: gsh_handling.c $(HFILES)
	$(CC) -c gsh_handling.c $(CFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		 -o $(ODIR)/gsh_handling.o


$(ODIR)/numrec_svd.o: numrec_svd.F 
	$(F77) -c $(FFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		numrec_svd.F -o $@

$(ODIR)/numrec_svd_s.o: numrec_svd.F
	$(F77) -c -DSINGLE_PREC $(FFLAGS) $(INCLUDES) \
		$(OPTIONS) $(FORMATS) \
		numrec_svd.F -o $(ODIR)/numrec_svd_s.o 

$(ODIR)/splinesf.o: splinesf.F
	$(F77) -c $(FFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		splinesf.F -o $@

$(ODIR)/splinesf_s.o: splinesf.F
	$(F77) -c -DSINGLE_PREC $(FFLAGS) $(INCLUDES) $(OPTIONS) $(FORMATS) \
		splinesf.F -o $(ODIR)/splinesf_s.o 



