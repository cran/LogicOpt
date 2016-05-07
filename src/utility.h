#ifndef UTILITY_H
#define UTILITY_H

/*
 *  assumes the memory manager is libmm.a
 *	- allows malloc(0) or realloc(obj, 0)
 *	- catches out of memory (and calls MMout_of_memory())
 *	- catch free(0) and realloc(0, size) in the macros
 */
#define NIL(type)		((type *) 0)

#ifndef R_PACKAGE
#define ALLOC(type, num)	\
    ((type *) malloc(sizeof(type) * (num)))
#define REALLOC(type, obj, num)	\
    (obj) ? ((type *) realloc((char *) obj, sizeof(type) * (num))) : \
	    ((type *) malloc(sizeof(type) * (num)))
#define FREE(obj)		\
    if ((obj)) { (void) free((char *) (obj)); (obj) = 0; }

#define Strdup(X) strdup(X)
#else
#include <R.h>
#include <Rinternals.h>
#define ALLOC(type, num) Calloc(num,type)
#define REALLOC(type, obj, num) Realloc(obj, num, type)
#define FREE(obj) Free(obj)
#define Strdup(X) strcpy(Calloc(strlen(X)+1, char), X) 

#define printf(...) Rprintf (__VA_ARGS__)

#endif

extern long  util_cpu_time(void);
extern char *util_print_time (long t);

#ifndef NIL_FN
#define NIL_FN(type) ((type (*)()) 0)
#endif /* NIL_FN */

#ifndef MAX
#define MAX(a,b)	((a) > (b) ? (a) : (b))
#endif /* MAX */
#ifndef MIN
#define MIN(a,b)	((a) < (b) ? (a) : (b))
#endif /* MIN */
#ifndef ABS
#define ABS(a)		((a) > 0 ? (a) : -(a))
#endif /* ABS */

#endif

