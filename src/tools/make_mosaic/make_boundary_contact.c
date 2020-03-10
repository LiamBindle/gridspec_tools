#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpp_io.h"
#include "mpp.h"
#include "make_boundary_contact.h"
#include "mosaic_util.h"
#include "tool_util.h"

const int MAXCONTACT = 100;
double* east_bound(const double *data, int nx, int ny);
double* west_bound(const double *data, int nx, int ny);
double* south_bound(const double *data, int nx, int ny);
double* north_bound(const double *data, int nx, int ny);
#define EPSLN (1.0e-10)

/************************************************************************************************
  This routine will create the boundary condition of the mosaic, whose tile grid files
  are tilefile. Then save the boundary condition (line contact) into file file
  contact_file. Return ncontact.
  This routine assume the starting and ending points of the  contact line are coincidence with
  the grid points of both tiles. lrg_rectangle tiles are assumed. 
************************************************************************************************/
int get_contact_index( int size1, int size2, double *x1, double *y1, double *x2, double *y2, double periodx,
		       double periody, int *start1, int *end1, int *start2, int *end2);
int make_boundary_contact(const char *mosaic_name, int ntile, char tilefile[][STRING],
			  const char *contact_file, double periodx, double periody,
			  const char *history, char *grid_type, char *congruence)
{
  char tagname[] = "$Name: quebec_gridV3_z1l $";

  double *xb1, *yb1, *xb2, *yb2;
  int ncontact, start1, end1, start2, end2;
  int nx1, ny1, nx2, ny2, n, m, i;
  int *nxp=NULL, *nyp=NULL;
  double **x=NULL, **y=NULL;
  double *x1=NULL, *x2=NULL, *y1=NULL, *y2=NULL;
  int tile1[MAXCONTACT], tile2[MAXCONTACT];
  int istart1[MAXCONTACT], iend1[MAXCONTACT], jstart1[MAXCONTACT], jend1[MAXCONTACT];
  int istart2[MAXCONTACT], iend2[MAXCONTACT], jstart2[MAXCONTACT], jend2[MAXCONTACT];
  char tile_history[512];
  char **tile_name=NULL;
  
  /* First read the grid */
  nxp = (int *)malloc(ntile*sizeof(int));
  nyp = (int *)malloc(ntile*sizeof(int));
  x = (double **)malloc(ntile*sizeof(double *));
  y = (double **)malloc(ntile*sizeof(double *));
  tile_name = (char **)malloc(ntile*sizeof(char *));
  strcpy(grid_type, "not-specified");
  for(n=0; n<ntile; n++) {
    int fid, vid;
    fid = mpp_open(tilefile[n], MPP_READ);
    nxp[n] = mpp_get_dimlen(fid, "nxp");
    nyp[n] = mpp_get_dimlen(fid, "nyp");
    x[n] = (double *)malloc(nxp[n]*nyp[n]*sizeof(double));
    y[n] = (double *)malloc(nxp[n]*nyp[n]*sizeof(double));
    tile_name[n] = (char *)malloc(STRING*sizeof(char));
    vid = mpp_get_varid(fid, "tile");
    mpp_get_var_value(fid, vid, tile_name[n]);
    vid = mpp_get_varid(fid, "x");
    mpp_get_var_value(fid, vid, x[n]);
    vid = mpp_get_varid(fid, "y");
    mpp_get_var_value(fid, vid, y[n]);
    if(n==0) {
      if(mpp_global_att_exist(fid, "history")) {
	char *pch;
	mpp_get_global_att(fid, "history", tile_history);
	pch = strtok (tile_history," ");
	while (pch != NULL) {
	  pch = strtok (NULL, " ");
	  if( !strcmp(pch, "--grid_type") ) {
	    pch = strtok (NULL, " ");
	    if(!strcmp(pch, "gnomonic_ed"))
	      strcpy(grid_type, "gnomonic_cubic_sphere");
	    else
	      strcpy(grid_type, pch);
	    break;
	  }
	}
      }
    }
    
    mpp_close(fid);
  }  

  strcpy(congruence, "true");
  for(n=1; n<ntile; n++) {
    if(nxp[n] != nxp[n-1] || nyp[n] != nyp[n-1] ) {
      strcpy(congruence, "false");
      break;
    }
  }
  
  ncontact = 0;
  for(n=0; n<ntile; n++) {
    for(m=n; m<ntile; m++) {
      x1 = x[n];
      x2 = x[m];
      y1 = y[n];
      y2 = y[m];
      nx1 = nxp[n]; ny1 = nyp[n];
      nx2 = nxp[m]; ny2 = nyp[m];
      /* East bound of tile1 and west bound of tile2 */
      xb1 = east_bound(x1, nx1, ny1);
      yb1 = east_bound(y1, nx1, ny1);
      xb2 = west_bound(x2, nx2, ny2);
      yb2 = west_bound(y2, nx2, ny2);
      if( get_contact_index( ny1, ny2, xb1, yb1, xb2, yb2, periodx, 0.0, &start1, &end1, &start2, &end2) ) {
	if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 1");
	istart1[ncontact] = nx1-1;
	iend1[ncontact]   = nx1-1;
	istart2[ncontact] = 1;
	iend2[ncontact]   = 1;    
	jstart1[ncontact] = start1;
	jend1[ncontact]   = end1;
	jstart2[ncontact] = start2;
	jend2[ncontact]   = end2;
	tile1[ncontact]   = n;
	tile2[ncontact]   = m;
	ncontact++;   
      }

      /* East bound of tile1 and SOUTH bound of tile2, tile1 and tile must be different tile */
      if(n != m) {
	free(xb2);
	free(yb2);
	xb2 = south_bound(x2, nx2, ny2);
	yb2 = south_bound(y2, nx2, ny2);
	if( get_contact_index( ny1, nx2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
	  if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 2");
	  istart1[ncontact] = nx1-1;
	  iend1[ncontact]   = nx1-1;
	  istart2[ncontact] = start2;
	  iend2[ncontact]   = end2;    
	  jstart1[ncontact] = start1;
	  jend1[ncontact]   = end1;
	  jstart2[ncontact] = 1;
	  jend2[ncontact]   = 1;
	  tile1[ncontact]   = n;
	  tile2[ncontact]   = m;  
	  ncontact++;
	}
      }
      free(xb1);
      free(yb1);   
      free(xb2);
      free(yb2);  
  
      /* South bound of tile1 and NORTH bound of tile2 */
      xb1 = south_bound(x1, nx1, ny1);
      yb1 = south_bound(y1, nx1, ny1);
      xb2 = north_bound(x2, nx2, ny2);
      yb2 = north_bound(y2, nx2, ny2);
      if( get_contact_index( nx1, nx2, xb1, yb1, xb2, yb2, 0.0, periody, &start1, &end1, &start2, &end2) ) {
	if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 3");
	istart1[ncontact] = start1;
	iend1[ncontact]   = end1;
	istart2[ncontact] = start2;
	iend2[ncontact]   = end2;    
	jstart1[ncontact] = 1;
	jend1[ncontact]   = 1;
	jstart2[ncontact] = ny2-1;
	jend2[ncontact]   = ny2-1;
	tile1[ncontact]   = n;
	tile2[ncontact]   = m;
	ncontact++;
      }

      /* South bound of tile1 and East bound of tile2, tile1 and tile must be different tile*/  
      if(n != m ) {
	free(xb2);
	free(yb2);
	xb2 = east_bound(x2, nx2, ny2);
	yb2 = east_bound(y2, nx2, ny2);
	if( get_contact_index( nx1, ny2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
	  if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 4");
	  istart1[ncontact] = start1;
	  iend1[ncontact]   = end1;
	  istart2[ncontact] = nx2-1;
	  iend2[ncontact]   = nx2-1;    
	  jstart1[ncontact] = 1;
	  jend1[ncontact]   = 1;
	  jstart2[ncontact] = start2;
	  jend2[ncontact]   = end2;
	  tile1[ncontact]   = n;
	  tile2[ncontact]   = m;
	  ncontact++;
	}
      }
      free(xb1);
      free(yb1);   
      free(xb2);
      free(yb2);  

      /* to avoid duplicate, the following will be done only when tile1 not equal to tile2 */
      if(n != m) {
	/* West bound of tile1 and east bound of tile2*/
	xb1 = west_bound(x1, nx1, ny1);
	yb1 = west_bound(y1, nx1, ny1);
	xb2 = east_bound(x2, nx2, ny2);
	yb2 = east_bound(y2, nx2, ny2);
	if( get_contact_index( ny1, ny2, xb1, yb1, xb2, yb2, periodx, 0.0, &start1, &end1, &start2, &end2) ) {
	  if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 5");
	  istart1[ncontact] = 1;
	  iend1[ncontact]   = 1;
	  istart2[ncontact] = nx2-1;
	  iend2[ncontact]   = nx2-1;    
	  jstart1[ncontact] = start1;
	  jend1[ncontact]   = end1;
	  jstart2[ncontact] = start2;
	  jend2[ncontact]   = end2;
	  tile1[ncontact]   = n;
	  tile2[ncontact]   = m;
	  ncontact++;
	}
	free(xb2);
	free(yb2);

	/* West bound of tile1 and North bound of tile2 */  
	xb2 = north_bound(x2, nx2, ny2);
	yb2 = north_bound(y2, nx2, ny2);
	if( get_contact_index( ny1, nx2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
	  if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 6");
	  istart1[ncontact] = 1;
	  iend1[ncontact]   = 1;
	  istart2[ncontact] = start2;
	  iend2[ncontact]   = end2;    
	  jstart1[ncontact] = start1;
	  jend1[ncontact]   = end1;
	  jstart2[ncontact] = ny2-1;
	  jend2[ncontact]   = ny2-1;
	  tile1[ncontact]   = n;
	  tile2[ncontact]   = m;
	  ncontact++;
	}
	free(xb1);
	free(yb1);   
	free(xb2);
	free(yb2);  

  
	/* North bound of tile1 and South bound of tile2 */
	xb1 = north_bound(x1, nx1, ny1);
	yb1 = north_bound(y1, nx1, ny1);
	xb2 = south_bound(x2, nx2, ny2);
	yb2 = south_bound(y2, nx2, ny2);
	if( get_contact_index( nx1, nx2, xb1, yb1, xb2, yb2, 0.0, periody, &start1, &end1, &start2, &end2) ) {
	  if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 7");
	  istart1[ncontact] = start1;
	  iend1[ncontact]   = end1;
	  istart2[ncontact] = start2;
	  iend2[ncontact]   = end2;    
	  jstart1[ncontact] = ny1-1;
	  jend1[ncontact]   = ny1-1;
	  jstart2[ncontact] = 1;
	  jend2[ncontact]   = 1;
	  tile1[ncontact]   = n;
	  tile2[ncontact]   = m;
	  ncontact++;
	}
	free(xb2);
	free(yb2);

	/* North bound of tile1 and West bound of tile2 */  
	xb2 = west_bound(x2, nx2, ny2);
	yb2 = west_bound(y2, nx2, ny2);
	if( get_contact_index( nx1, ny2, xb1, yb1, xb2, yb2, 0.0, 0.0, &start1, &end1, &start2, &end2) ) {
	  if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 8");
	  istart1[ncontact] = start1;
	  iend1[ncontact]   = end1;
	  istart2[ncontact] = 1;
	  iend2[ncontact]   = 1;    
	  jstart1[ncontact] = ny1-1;
	  jend1[ncontact]   = ny1-1;
	  jstart2[ncontact] = start2;
	  jend2[ncontact]   = end2;
	  tile1[ncontact]   = n;
	  tile2[ncontact]   = m;
	  ncontact++;
	}
	free(xb1);
	free(yb1);   
	free(xb2);
	free(yb2);  
      }

      /* when tile1 = tile2, we need to consider about folded. Only foled north is considered here */
      if(n == m) {
	int i, folded = 1;
	double dx;
	xb1 = north_bound(x1, nx1, ny1);
	yb1 = north_bound(y1, nx1, ny1);
	for(i=0; i<nx1/2; i++) {
	  if( yb1[i] != yb1[nx1-i-1] ) {
	    folded = 0;
	    break;
	  }
	  dx = fabs(xb1[i] - xb1[nx1-i-1]);
	  if( dx !=0 && dx != 360 ) {
	    folded = 0;
	    break;
	  }
	}
	if(folded) {
	  if(ncontact == MAXCONTACT) mpp_error("make_boundary_contact: number of contacts is more than MAXCONTACT 9");
	  istart1[ncontact] = 1;
	  iend1[ncontact]   = nx1/2;
	  istart2[ncontact] = nx1-1;
	  iend2[ncontact]   = nx1/2+1;    
	  jstart1[ncontact] = ny1-1;
	  jend1[ncontact]   = ny1-1;
	  jstart2[ncontact] = ny1-1;
	  jend2[ncontact]   = ny1-1;
	  tile1[ncontact]   = n;
	  tile2[ncontact]   = m;
	  ncontact++;
	}
	free(xb1);
	free(yb1);
      }
    }
  }
  
  /* write out data */
  if(ncontact>0) {
    char str[STRING], outfile[STRING];
    int fid, dim_ncontact, dim_string, id_mosaic, id_contacts;
    int id_contact_index, dim[2];
    size_t start[4], nwrite[4];

    fid = mpp_open(contact_file, MPP_WRITE);
    /* define dimenison */
    dim_ncontact = mpp_def_dim(fid, "ncontact", ncontact);
    dim_string = mpp_def_dim(fid, "string", STRING);    
    /* define variable */
    dim[0] = dim_ncontact; dim[1] = dim_string;
    id_contacts = mpp_def_var(fid, "contacts", MPP_CHAR, 2, dim, 5, "standard_name", "grid_contact_spec",
			      "contact_type", "boundary", "alignment", "true",
			      "contact_index", "contact_index", "orientation", "orient");
    id_contact_index = mpp_def_var(fid, "contact_index", MPP_CHAR, 2, dim, 1, "standard_name",
				   "starting_ending_point_index_of_contact");

    mpp_def_global_att(fid, "grid_version", grid_version);
    mpp_def_global_att(fid, "code_version", tagname);
    mpp_def_global_att(fid, "history", history);
    mpp_end_def(fid);

    /* write out data */
    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1;
    }

    for(n=0; n<ncontact; n++) {
      sprintf(str,"%s:%s::%s:%s", mosaic_name, tile_name[tile1[n]], mosaic_name,
	      tile_name[tile2[n]]);
      start[0] = n; nwrite[1] = strlen(str);
      mpp_put_var_value_block(fid, id_contacts, start, nwrite, str);
      sprintf(str,"%d:%d,%d:%d::%d:%d,%d:%d", istart1[n], iend1[n], jstart1[n], jend1[n],
	      istart2[n], iend2[n], jstart2[n], jend2[n] );
      nwrite[1] = strlen(str);
      mpp_put_var_value_block(fid, id_contact_index, start, nwrite, str);
    }
    mpp_close(fid);    
  }  


  /* release memory */
  for(n=0; n<ntile; n++) {
    free(tile_name[n]);
    free(x[n]);
    free(y[n]);
  }
  free(tile_name);
  free(nxp);
  free(nyp);
  free(x);
  free(y);
  
  return ncontact;
}; /*make_boundary_contact*/


/*************************************************************************************************/
int get_contact_index( int size1, int size2, double *x1, double *y1, double *x2, double *y2, double periodx,
		       double periody, int *start1, int *end1, int *start2, int *end2)
{
  int i1, i2;
  double dx, dy;

  /* Find the first point in tile 1 cocindent with a point in tile2  */
  *start1 = -1;
  *start2 = -1;
  for(i1=0; i1<size1; i1++) {
    for(i2=0; i2<size2; i2++) {
      dx = fabs(x1[i1]- x2[i2]);
      dx = min(dx, fabs(dx-periodx));
      dy = fabs(y1[i1]- y2[i2]);
      dy = min(dy, fabs(dy-periody));
      if( dx < EPSLN && dy <EPSLN ) {
	*start1 = i1+1;
	*start2 = i2+1;
	goto foundstart;
      }
    }
  }

  return 0;
    
  foundstart:

  /* Find the last point in tile 1 cocindent with a point in tile2 */
  *end1 = -1;
  *end2 = -1;
  for(i1=size1-1; i1>=0; i1--) {
    for(i2=size2-1; i2>=0; i2--) {
      dx = fabs(x1[i1]- x2[i2]);
      dx = min(dx, fabs(dx-periodx));
      dy = fabs(y1[i1]- y2[i2]);
      dy = min(dy, fabs(dy-periody));
      if( dx < EPSLN && dy <EPSLN ) {
	*end1 = i1+1;
	*end2 = i2+1;
	goto foundend;
      }
    }
  }

  return 0;
    
 foundend: if( *start1 == *end1 || *start2 == *end2 ) return 0;

  if(*start1 > *end1 )
    (*start1)--;
  else
    (*end1)--;

  if(*start2 > *end2 )
    (*start2)--;
  else
    (*end2)--;  

  return 1;
    
};


double* west_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(ny*sizeof(double));
  for(i=0; i<ny; i++) bnd[i] = data[i*nx];
  return bnd;
}

double* east_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(ny*sizeof(double));
  for(i=0; i<ny; i++) bnd[i] = data[i*nx+nx-1];
  return bnd;
}

double* south_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(nx*sizeof(double));
  for(i=0; i<nx; i++) bnd[i] = data[i];
  return bnd;
}

double* north_bound(const double *data, int nx, int ny)
{
  int i;
  double *bnd;
  bnd = (double *)malloc(nx*sizeof(double));
  for(i=0; i<nx; i++) bnd[i] = data[(ny-1)*nx+i];
  return bnd;
}
