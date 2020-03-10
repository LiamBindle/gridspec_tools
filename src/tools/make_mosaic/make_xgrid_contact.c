#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "mpp_io.h"
#include "mpp.h"
#include "constant.h"
#include "make_xgrid_contact.h"
#include "create_xgrid.h"
#include "tool_util.h"
#define  MAX_XGRID_FILE (100)
 

typedef struct {
  int nx, ny;
  double *x, *y;
  double *mask;
  char tile_name[STRING];
} grid_type;


/* This routine will get all the tile grid information */
void get_mosaic_grid(const char *mosaic, int ntiles, grid_type *Grid)
{
  int nx, ny, n, i, j;
  int m_fid, g_fid, vid, vid_childfile;
  char dir[STRING], grid_file[STRING], filename[STRING];
  double *tmpx=NULL, *tmpy=NULL;
  size_t start[4], nread[4];
  
  m_fid = mpp_open(mosaic, MPP_READ);
  vid_childfile = mpp_get_varid(m_fid, TILE_FILES_NAME);
  get_file_path(mosaic , dir);
  for(n=0; n<ntiles ; n++) {
    start[0] = n; start[1] = 0; nread[0] = 1; nread[1] = STRING;
    mpp_get_var_value_block(m_fid, vid_childfile, start, nread, filename);
    sprintf(grid_file, "%s/%s", dir, filename);
    g_fid = mpp_open(grid_file, MPP_READ);
    vid = mpp_get_varid(g_fid, "tile");
    mpp_get_var_value(g_fid, vid, Grid[n].tile_name);
  
    nx = mpp_get_dimlen(g_fid, "nx");
    ny = mpp_get_dimlen(g_fid, "ny");
    if(nx%2) mpp_error("make_xgrid_mosaic: the size of dimension nx of mosaic should be even (on supergrid)");
    if(ny%2) mpp_error("make_xgrid_mosaic: the size of dimension ny of mosaic should be even (on supergrid)");
    tmpx = (double *)malloc((nx+1)*(ny+1)*sizeof(double));
    tmpy = (double *)malloc((nx+1)*(ny+1)*sizeof(double));
    nx /= 2;
    ny /= 2;
    Grid[n].nx =nx;
    Grid[n].ny =ny;
    Grid[n].x = (double *)malloc((nx+1)*(ny+1)*sizeof(double));
    Grid[n].y = (double *)malloc((nx+1)*(ny+1)*sizeof(double));
    Grid[n].mask = (double *)malloc(nx*ny*sizeof(double));
    vid = mpp_get_varid(g_fid, "x");
    mpp_get_var_value(g_fid, vid, tmpx);
    vid = mpp_get_varid(g_fid, "y");
    mpp_get_var_value(g_fid, vid, tmpy);
    for(j=0; j<=ny; j++) for(i=0; i<=nx; i++) {
      Grid[n].x[j*(nx+1)+i] = tmpx[2*j*(2*nx+1)+2*i];
      Grid[n].y[j*(nx+1)+i] = tmpy[2*j*(2*nx+1)+2*i];
    }
    
    for(i=0; i<nx*ny; i++) Grid[n].mask[i] = 1.0;
    free(tmpx);
    free(tmpy);
    mpp_close(g_fid);
  }  

  mpp_close(m_fid);
  
}
  


/****************************************************************************
  This routine will compute the exchange grid information between two mosaic
  and write out those information.

****************************************************************************/
int make_xgrid_contact(const char* mosaic1, const char* mosaic2, const char* contact_file,
		       const char* history)
{
  int fid, vid, n;
  int ntiles1, ntiles2, n1, n2, nxgrid, nxfile=0;
  grid_type *Grid1=NULL, *Grid2=NULL;
  int *i1=NULL, *j1=NULL, *i2=NULL, *j2=NULL;
  double *di1=NULL, *dj1=NULL, *di2=NULL, *dj2=NULL, *xarea=NULL;
  double *clon=NULL, *clat=NULL;
  
  char xfile[MAX_XGRID_FILE][STRING];
  char mosaic1_name[STRING], mosaic2_name[STRING];  
  char tagname[] = "$Name: quebec_gridV3_z1l $";
  
  /* First read the grid information of mosaic1*/
  fid = mpp_open(mosaic1, MPP_READ);
  ntiles1 = mpp_get_dimlen(fid, NTILES_NAME);
  vid = mpp_get_varid(fid, MOSAIC_NAME);    
  mpp_get_var_value(fid, vid, mosaic1_name);
  mpp_close(fid);
  fid = mpp_open(mosaic2, MPP_READ);
  vid = mpp_get_varid(fid, MOSAIC_NAME);    
  mpp_get_var_value(fid, vid, mosaic2_name);
  ntiles2 = mpp_get_dimlen(fid, NTILES_NAME);
  mpp_close(fid);
  Grid1 = (grid_type *)malloc(ntiles1*sizeof(grid_type));
  Grid2 = (grid_type *)malloc(ntiles2*sizeof(grid_type));
  get_mosaic_grid(mosaic1, ntiles1, Grid1);
  get_mosaic_grid(mosaic2, ntiles2, Grid2);
  i1    = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  j1    = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  i2    = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  j2    = (int    *)malloc(MAXXGRID   * sizeof(int   ));
  clon  = (double *)malloc(MAXXGRID   * sizeof(double));
  clat  = (double *)malloc(MAXXGRID   * sizeof(double));
  di1   = (double *)malloc(MAXXGRID   * sizeof(double));
  dj1   = (double *)malloc(MAXXGRID   * sizeof(double));
  di2   = (double *)malloc(MAXXGRID   * sizeof(double));
  dj2   = (double *)malloc(MAXXGRID   * sizeof(double));
  xarea = (double *)malloc(MAXXGRID   * sizeof(double));

  nxfile = 0;
  for(n1=0; n1<ntiles1; n1++) {
    for(n2=0;n2<ntiles2;n2++) {
      nxgrid = create_xgrid_2dx2d_order2(&(Grid1[n1].nx), &(Grid1[n1].ny), &(Grid2[n2].nx), &(Grid2[n2].ny),
					 Grid1[n1].x, Grid1[n1].y, Grid1[n2].x, Grid1[n2].y, Grid1[n1].mask,
					 i1, j1, i2, j2, xarea, clon, clat);
      if(nxgrid>0) {
	int dim_ncontact, dim_two, dims[2], i;
	int id_tile1_cell, id_tile2_cell, id_xgrid_area, id_tile1_dist, id_tile2_dist;
	size_t start[4], nwrite[4];
	if(nxfile==MAX_XGRID_FILE) mpp_error("make_xgrid_contact: number of exchange grid file is greater "
					     "than MAX_XGRID_FILE, increase MAX_XGRID_FILE");
	sprintf(xfile[nxfile], "%s_%sXlnd_%s_%s.nc", mosaic1_name, Grid1[n1].tile_name, mosaic2_name, Grid2[n2].tile_name);
	fid = mpp_open(xfile[nxfile], MPP_WRITE);
	dim_ncontact = mpp_def_dim(fid, "ncontact", nxgrid);
	dim_two    = mpp_def_dim(fid, "two", 2);
	dims[0] = dim_ncontact;
	dims[1] = dim_two;
	id_tile1_cell = mpp_def_var(fid, "tile1_cell", NC_INT, 2, dims, 1,
				    "standard_name", "parent_cell_indices_in_mosaic1");
	id_tile2_cell = mpp_def_var(fid, "tile2_cell", NC_INT, 2, dims, 1,
				    "standard_name", "parent_cell_indices_in_mosaic2");
	id_xgrid_area = mpp_def_var(fid, "xgrid_area", NC_DOUBLE, 1, &dim_ncontact, 2,
				    "standard_name", "exchange_grid_area", "units", "m2");
	id_tile1_dist = mpp_def_var(fid, "tile1_distance", NC_DOUBLE, 2, dims, 1,
				    "standard_name", "distance_from_parent1_cell_centroid");
	id_tile2_dist = mpp_def_var(fid, "tile2_distance", NC_DOUBLE, 2, dims, 1,
				    "standard_name", "distance_from_parent2_cell_centroid");

	mpp_end_def(fid);
	for(i=0; i<4; i++) {
	  start[i] = 0; nwrite[i] = 1;
	}
	   
	mpp_put_var_value(fid, id_tile1_cell, i1);
	mpp_put_var_value(fid, id_tile2_cell, i2);
	mpp_put_var_value(fid, id_xgrid_area, xarea);
	/*	mpp_put_var_value(fid, id_tile1_dist, di1);
	  mpp_put_var_value(fid, id_tile2_dist, di2); */
	start[1] = 1;
	mpp_put_var_value(fid, id_tile1_cell, j1);
	mpp_put_var_value(fid, id_tile2_cell, j2);
	/*	mpp_put_var_value(fid, id_tile1_dist, dj1);
	  mpp_put_var_value(fid, id_tile2_dist, dj2); */
	mpp_close(fid);
	nxfile ++;
      }
    }
  }

  if(nxfile >0) {
    int dim_string, dim_nchild, dim[2];
    int id_childfile, i;
    size_t start[4], nwrite[4];
    
    fid = mpp_open(contact_file, MPP_WRITE);
    dim_string = mpp_def_dim(fid, STRING_NAME, STRING);  
    dim_nchild = mpp_def_dim(fid, NCONTACT_NAME, nxfile);
    dim[0] = dim_nchild; dim[1] = dim_string;
    id_childfile = mpp_def_var(fid, TILE_FILES_NAME, MPP_CHAR, 2, dim, 0);
    mpp_def_global_att(fid, GRID_VERSION_NAME, grid_version);
    mpp_def_global_att(fid, CODE_VERSION_NAME, tagname);
    mpp_def_global_att(fid, HISTORY_NAME, history);
    mpp_end_def(fid);    
    /* write out data */
    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1;
    }
    for(n=0; n<nxfile; n++) {
      start[0] = n;
       nwrite[1] = strlen(xfile[n]);
       mpp_put_var_value_block(fid, id_childfile, start, nwrite, xfile[n]);
    }
    mpp_close(fid);
  }

  free(i1);
  free(j1);
  free(i2);
  free(j2);  
  free(di1);
  free(dj1);
  free(di2);
  free(dj2);
  free(xarea);

  for(n=0; n<ntiles1; n++) {
    free(Grid1[n].x);
    free(Grid1[n].y);
  }
  for(n=0; n<ntiles2; n++) {
    free(Grid2[n].x);
    free(Grid2[n].y);
  }
  
  free(Grid1);
  free(Grid2);
  
  return nxfile;
    
}; /* make_xgrid_contact */





  
