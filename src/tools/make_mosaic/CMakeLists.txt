add_executable(make_mosaic
   make_boundary_contact.c
   make_mosaic.c
   make_xgrid_contact.c
)
target_link_libraries(make_mosaic mpp shared m)