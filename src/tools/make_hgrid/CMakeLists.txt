add_executable(make_hgrid
   create_conformal_cubic_grid.c
   create_gnomonic_cubic_grid.c
   create_grid_from_file.c
   create_lonlat_grid.c
   make_hgrid.c
)

target_link_libraries(make_hgrid shared mpp m)