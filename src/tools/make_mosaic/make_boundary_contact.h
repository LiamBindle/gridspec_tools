#ifndef MAKE_BOUNDATRY_CONTACT_
#define MAKE_BOUNDATRY_CONTACT_
#include "constant.h"
int make_boundary_contact(const char *mosaic_name, int ntile, char tilefile[][STRING],
			  const char *contact_file, double periodx, double periody, 
                          const char *history, char *grid_type, char *congruence);
#endif
