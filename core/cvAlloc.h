#ifndef CVALLOC_H
#define CVALLOC_H

#ifdef __cplusplus
extern "C" {
#endif
#include <stdlib.h>

#define ALLOC_ARRAY_SIZE (10240 * 1024)


#define CvMalloc 	malloc
#define CvFree 		free
#define CvCalloc 	calloc
#define CvRealloc 	realloc



#ifdef __cplusplus
}
#endif /* end of __cplusplus */

#endif /* CVALLOC_H */
