#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "make_boundary_contact.h"
#include "make_xgrid_contact.h"
#include "constant.h"
#include "mpp_io.h"
#include "mpp.h"
#include "tool_util.h"

char *usage[] = {
  "",
  " make_mosaic --num_tiles num_tiles --tile_file tile_file [--dir directory]       ",
  "             [--num_mosaic num_mosaic] [--mosaic_file mosaic_file]               ",
  "             [--mosaic_name mosaic_name] [--periodx #] [--periody #]             ",
  "             [--generate_contact]                                                ",
  "                                                                                 ",
  "make_mosaic generates Mosaic of tile grid files. If generate_overlap is          ",
  "specified, it also contains a link to the mosaic file that contains contact      ",
  "information of child files, specified by index, contact type.                    ",
  "                                                                                 ",
  "make_solo_mosaic takes the following flags:                                      ",
  "                                                                                 ",
  "REQUIRED:                                                                        ",
  "                                                                                 ",
  "--num_tiles num_tiles    Number of grid tiles.                                   ",
  "                                                                                 ",
  "--tile_file tile_file    File name of all grid tiles. The number of entry must   ",
  "be ntiles (The file name should be full file name) or 1 (The file name will be   ",
  "base file name, which does not include tile name).                               ",
  "                                                                                 ",    
  "OPTIONAL FLAGS                                                                   ",
  "                                                                                 ",  
  "--num_mosaic #     Number of child mosaic file.                                  ",
  "                                                                                 ",
  "--mosaic_file mosaic_file File name of child mosaic file. The number of entry    ",
  "                          must be num_mosaic.                                    ",
  "                                                                                 ",    
  "--dir directory          The directory that contains all the child files. If dir ",
  "                         is not specified, the file path should be contained in  ",
  "                         child_file.                                             ",
  "                                                                                 ",  
  "--mosaic_name name       mosaic name. The output file will be mosaic_name.nc.    ",
  "                         default is 'mosaic'.                                    ",
  "                                                                                 ",
  "--periodx #              Specify the period in x-direction of mosaic. Default    ",
  "                         value is 0 (not periodic). This is needed only when     ",
  "                         child_file are tile grid files.                         ",
  "                                                                                 ",
  "--periody #              Specify the period in y-direction of mosaic. Default    ",
  "                         value is 0 (not periodic).  This is needed only when    ",
  "                         child_file are tile grid files.                         ",
  "                                                                                 ",
  "--generate_contact       Specify whether to generate the contact between mosaic  ",
  "                         when num_mosaic is greater than 1                       ",
  "                                                                                 ",
  NULL};
  
const int MAXCHILD = 100;
const int SHORTSTRING = 32;
char tagname[] = "$Name: quebec_gridV3_z1l $";

main (int argc, char *argv[])
{

  extern char *optarg;
  char *pch=NULL, *tile_dir=NULL, history[512], entry[1280];
  char tilefile[MAXCHILD][STRING];
  char gridfile[MAXCHILD][STRING];
  char mosaicfile[MAXCHILD][STRING];
  char basefile[STRING];
  char contact_file[STRING], overlap_file[STRING];
  char congruence[STRING];
  char grid_type[STRING];
  int  ntile=0, ntile2=0, ncontact=0, noverlap=0;
  int  nmosaic=0, nmosaic2=0;
  int  generate_contact = 0;

  double periodx=0, periody=0;

  char mosaic_name[128] = "mosaic";
  char grid_descriptor[128] = "";
  int c, i, n, m, l, errflg;
  
  int option_index = 0;
  static struct option long_options[] = {
    {"mosaic_name",       required_argument, NULL, 'm'},
    {"num_tiles",         required_argument, NULL, 'n'},
    {"grid_descriptor",   required_argument, NULL, 'g'},
    {"tile_file",         required_argument, NULL, 'f'},
    {"num_mosaic",  required_argument, NULL, 'l'},
    {"mosaic_file", required_argument, NULL, 'k'},
    {"periodx",           required_argument, NULL, 'x'},
    {"periody",           required_argument, NULL, 'y'},
    {"tile_directory",    required_argument, NULL, 'd'},
    {"generate_contact",  no_argument,       NULL, 'c'},
    {NULL, 0, NULL, 0}
  };

  mpp_init(&argc, &argv);
  /* this tool must be run one processor */
  if(mpp_npes()>1) mpp_error("make_mosaic: this tool must be run on one processor");
  
  errflg = (argc == 1);
  /* First read command line arguments. */

  while ((c = getopt_long(argc, argv, "h", long_options, &option_index)) != -1) {
    switch (c) {
    case 'n':
      ntile = atoi(optarg);
      break;
    case 'm':
      strcpy(mosaic_name, optarg);
      break;
    case 'g':
      strcpy(grid_descriptor, optarg);
      break;
    case 'f':
      strcpy(entry, optarg); 
      pch = strtok(entry, ", ");
      while( pch != NULL) {
        strcpy(tilefile[ntile2++], pch);
	pch = strtok(NULL, ", ");
      }
      break;
    case 'l':
      nmosaic = atoi(optarg);
      break;
    case 'k':
      strcpy(entry, optarg); 
      pch = strtok(entry, ", ");
      while( pch != NULL) {
        strcpy(mosaicfile[nmosaic2++], pch);
	pch = strtok(NULL, ", ");
      }
      break;
    case 'x':
      periodx = atof(optarg);
      break;
    case 'y':
      periody = atof(optarg);
      break; 
    case 'd':  // path of the simple grid file.
      tile_dir = optarg;
      break;
    case 'c':
      generate_contact = 1;
      break;
    case '?':
      errflg++;
    }
  }

  if (errflg || ntile < 1) {
    char **u = usage;
    while (*u) { fprintf(stderr, "%s\n", *u); u++; }
    exit(2);      
  }

  strcpy(history,argv[0]);

  for(i=1;i<argc;i++) {
    strcat(history, " ");
    strcat(history, argv[i]);
  }

  if(ntile > MAXCHILD) {
    mpp_error("make_mosaic: number of tiles is greater than MAXCHILD.");
  }

  /* ntile2 should equal ntile or ntile2 should be 1
     if ntile2 is 1 and ntile is greater than 1, the tile file name is only the base file name,
     need to add .tile#.nc at the end of the file
  */
  if(ntile2 == ntile) {
    /* do nothing */
  }
  else if(ntile2 == 1) {
    strcpy(basefile, tilefile[0]);
    for(n=0; n<ntile; n++) sprintf(tilefile[n], "%s.tile%d.nc", basefile, n+1);
  }
  else {
    mpp_error("make_mosaic: number of files (>1) specified through --child_file does not equal to num_tiles");
  }
  
  if(tile_dir) {
    n = strlen(tile_dir);
    if( tile_dir[n-1] != '/') strcat(tile_dir, "/");
    for(n=0; n<ntile; n++) sprintf(gridfile[n], "%s%s", tile_dir, tilefile[n]);
  }
  else {
    for(n=0; n<ntile; n++) strcpy(gridfile[n], tilefile[n]);
  }

  /*
    if num_mosaic is greater than 0, then number of entry specified through --mosaic_file
    must equal to num_mosaic
  */
  if(nmosaic>0) {
    if(nmosaic > MAXCHILD) mpp_error("make_mosaic: number of child mosaic is greater than MAXCHILD.");
    if(nmosaic != nmosaic2) mpp_error("make_mosaic: then number of entry specified through --mosaic_file "
				      "must equal to num_mosaic" );
  } 
  
  ncontact = 0;

  sprintf(contact_file, "%s.boundary_contact.nc", mosaic_name);
  ncontact = make_boundary_contact(mosaic_name, ntile, gridfile, contact_file, periodx,
				   periody, history, grid_type, congruence);
  if(nmosaic > 0 && generate_contact ) {
    int fid, vid;
    char mosaic1_name[STRING], mosaic2_name[STRING];
    /* when generating the overlap between mosaic, the number of mosaic must be two. */
    if(nmosaic != 2) mpp_error("make_mosaic: num_mosaic must be 2 when generate_contact is specified");
    fid = mpp_open(mosaicfile[0], MPP_READ);
    vid = mpp_get_varid(fid, "mosaic");    
    mpp_get_var_value(fid, vid, mosaic1_name);
    mpp_close(fid);
    fid = mpp_open(mosaicfile[1], MPP_READ);
    vid = mpp_get_varid(fid, "mosaic");    
    mpp_get_var_value(fid, vid, mosaic2_name);
    mpp_close(fid);
      
    sprintf(contact_file, "%sX%s.nc", mosaic1_name, mosaic2_name);
    noverlap = make_xgrid_contact(mosaicfile[0], mosaicfile[1], overlap_file, history);
  }
  
  /* write out mosaic */
  {
    char outfile[STRING];
    size_t start[4], nwrite[4];
    int dim_ntile, dim_nmosaic, dim_ncontact, dim_string, id_mosaic, id_contact;
    int fid, dim[2], id_tilefile, id_mosaicfile, id_overlap;
    
    sprintf(outfile, "%s.nc", mosaic_name);
    fid = mpp_open(outfile, MPP_WRITE);
    /* define dimenison */
    dim_ntile = mpp_def_dim(fid, NTILES_NAME, ntile);
    if(nmosaic>0)  dim_nmosaic = mpp_def_dim(fid, NMOSAIC_NAME, nmosaic);
    dim_string = mpp_def_dim(fid, STRING_NAME, STRING);    
    /* define variable */
    id_mosaic = mpp_def_var(fid, MOSAIC_NAME, MPP_CHAR, 1, &dim_string, 6, "standard_name",
			    "grid_mosaic_spec", "grid_type", grid_type, "children", "gridtiles",
			    "contact_regions", "contacts", "grid_descriptor", grid_descriptor,
			    CONGRUCENCE_NAME, congruence);
    dim[0] = dim_ntile; dim[1] = dim_string;
    id_tilefile = mpp_def_var(fid, TILE_FILES_NAME, MPP_CHAR, 2, dim, 0);

    if(nmosaic>0) {
      dim[0] = dim_nmosaic;
      id_mosaicfile = mpp_def_var(fid, MOSAIC_FILES_NAME, MPP_CHAR, 2, dim, 0);
    }
    
    if(ncontact>0) {
      id_contact = mpp_def_var(fid, CONTACT_FILES_NAME, MPP_CHAR, 1, &dim_string, 2, "standard_name", "mosaic boundary contact",
			       "contact_type", "boundary_contact");
    }
    if(noverlap>0) {
      id_overlap = mpp_def_var(fid, OVERLAP_FILES_NAME, MPP_CHAR, 1, &dim_string, 2, "standard_name", "overlap between mosaic",
			       "contact_type", "overlap");
    }
    mpp_def_global_att(fid, GRID_VERSION_NAME, grid_version);
    mpp_def_global_att(fid, CODE_VERSION_NAME, tagname);
    mpp_def_global_att(fid, HISTORY_NAME, history);
    mpp_end_def(fid);

    /* write out data */
    for(i=0; i<4; i++) {
      start[i] = 0; nwrite[i] = 1;
    }
    nwrite[0] = strlen(mosaic_name);
    mpp_put_var_value_block(fid, id_mosaic, start, nwrite, mosaic_name);
    nwrite[0] = 1;
    if(ncontact>0) {
      nwrite[0]=strlen(contact_file); 
      mpp_put_var_value_block(fid, id_contact, start, nwrite, contact_file);
    }
    if(noverlap>0) {
      nwrite[0]=strlen(overlap_file); 
      mpp_put_var_value_block(fid, id_overlap, start, nwrite, overlap_file);
    }

    for(n=0; n<ntile; n++) {
      start[0]=n; nwrite[0] =1;
      nwrite[1] = strlen(tilefile[n]);
      mpp_put_var_value_block(fid, id_tilefile, start, nwrite, tilefile[n]);
    }

    for(n=0; n<nmosaic; n++) {
      start[0]=n;
      nwrite[1] = strlen(mosaicfile[n]);
      mpp_put_var_value_block(fid, id_mosaicfile, start, nwrite, mosaicfile[n]);
    }
    
    mpp_close(fid);
  }
    
  
  printf("congradulation: You have successfully run make_mosaic\n");

  return 0;
  
}; // end of main



