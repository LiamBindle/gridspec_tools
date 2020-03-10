# The following three directory may need to set.
NETCDFPATH = /usr/local/netcdf-3.6.2
SHAREDIR = $(PWD)/../../shared/mosaic
COREDIR  = $(PWD)
TOOLSHAREDIR = $(PWD)/../shared
CFLAGS = -O2 -fast -I$(TOOLSHAREDIR) -I$(COREDIR) -I$(SHAREDIR) -I${NETCDFPATH}/include -I/usr/include
LDFLAGS =  -L${NETCDFPATH}/lib -lm -lnetcdf 
LNFLAGS = -v
CC      = icc

OBJS = mosaic_util.o mpp.o mpp_domain.o mpp_io.o tool_util.o make_boundary_contact.o  \
       make_xgrid_contact.o  interp.o create_xgrid.o make_mosaic.o

HEADERS = Makefile $(TOOLSHAREDIR)/mpp.h  $(TOOLSHAREDIR)/mpp_domain.h   \
          $(TOOLSHAREDIR)/mpp_io.h $(TOOLSHAREDIR)/tool_util.h               \
          $(COREDIR)/make_boundary_contact.h $(SHAREDIR)/interp.h  $(SHAREDIR)/create_xgrid.h \
          $(COREDIR)/make_xgrid_contact.h

make_mosaic: $(OBJS)
	$(CC) $(LNFLAGS) -o $@ $(OBJS) $(LDFLAGS)

make_mosaic.o: $(COREDIR)/make_mosaic.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(COREDIR)/make_mosaic.c 

mpp.o: $(TOOLSHAREDIR)/mpp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp.c 

mpp_domain.o: $(TOOLSHAREDIR)/mpp_domain.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_domain.c 

mpp_io.o: $(TOOLSHAREDIR)/mpp_io.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/mpp_io.c

mosaic_util.o: $(SHAREDIR)/mosaic_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/mosaic_util.c

tool_util.o: $(TOOLSHAREDIR)/tool_util.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(TOOLSHAREDIR)/tool_util.c

make_boundary_contact.o: $(COREDIR)/make_boundary_contact.c 
	$(CC) $(CFLAGS) -c $(COREDIR)/make_boundary_contact.c 

make_xgrid_contact.o: $(COREDIR)/make_xgrid_contact.c 
	$(CC) $(CFLAGS) -c $(COREDIR)/make_xgrid_contact.c 

interp.o: $(SHAREDIR)/interp.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/interp.c 

create_xgrid.o: $(SHAREDIR)/create_xgrid.c $(HEADERS)
	$(CC) $(CFLAGS) -c $(SHAREDIR)/create_xgrid.c
