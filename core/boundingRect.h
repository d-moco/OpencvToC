#ifndef BOUNDINGRECT_H
#define BOUNDINGRECT_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "types_c.h"




#ifdef __cplusplus
extern "C" {
#endif 

	box_rect* findbox(unsigned char* mask, int width, int height, int areaMax, int boxMax, int* ret);

#ifdef __cplusplus
}
#endif 

#endif