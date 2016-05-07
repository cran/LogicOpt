#include "espresso.h"

/**
static void dump_irredundant(pset_family E, pset_family Rt, pset_family Rp, sm_matrix *table);
**/

/* assume pre-sorted covers */
static bool equal_cover(pcover A, pcover B)
{
    if (A->count != B->count) return FALSE;
    for(int i = 0; i < A->count; i++) {
       if (!setp_equal(GETSET(A,i),GETSET(B,i))) 
          return FALSE;
    }
    return TRUE;
}


static pcover do_minimize(pset_family F, pset_family D, pset_family R, int exact_cover, int weighted, int max_solutions, int full_pi_cover) ;


/*
 *  minimize_exact -- main entry point for exact minimization
 *
 *  Global flags which affect this routine are:
 *
 *      debug
 *      skip_make_sparse
 */

pcover
minimize_exact(pset_family F, pset_family D, pset_family R, int exact_cover, int max_solutions, int full_pi_cover)
{
    return do_minimize(F, D, R, exact_cover, /*weighted*/ 0, max_solutions, full_pi_cover);
}


pcover
minimize_exact_literals(pset_family F, pset_family D, pset_family R, int exact_cover, int full_pi_cover)
{
    return do_minimize(F, D, R, exact_cover, /*weighted*/ 1, 0, 0);
}



static pcover
do_minimize(pset_family F, pset_family D, pset_family R, int exact_cover, int weighted, int max_solutions, int full_pi_cover)
{
    pcover newFprev, newFhead, newF, E, Rt, Rp;
    pset p, last;
    int heur, level, *weights;
    sm_matrix *table;
    sm_row *cover, *cover_last;
    sm_element *pe;
    int debug_save = debug;
    int n_solutions=0; 

    newFhead = (pcover) NULL;
    newFprev = (pcover) NULL;
    newF = (pcover) NULL;

    if (debug & EXACT) {
	debug |= (IRRED | MINCOV);
    }
#if defined(sun) || defined(bsd4_2)			/* hack ... */
    if (debug & MINCOV) {
	setlinebuf(stdout);
    }
#endif
    level = (debug & MINCOV) ? 4 : 0;
    heur = ! exact_cover;

    /* Generate all prime implicants */
    EXEC(F = primes_consensus(cube2list(F, D)), "PRIMES     ", F);

    /* Setup the prime implicant table */
    EXEC(irred_split_cover(F, D, &E, &Rt, &Rp), "ESSENTIALS ", E);
    EXEC(table = irred_derive_table(D, E, Rp),  "PI-TABLE   ", Rp);

    /* Solve either a weighted or nonweighted covering problem */
    if (weighted) {
	/* correct only for all 2-valued variables */
	weights = ALLOC(int, F->count);
	foreach_set(Rp, last, p) {
	    weights[SIZE(p)] = cube.size - set_ord(p);
	}
    } else {
	weights = NIL(int);
    }

    EXEC(cover=sm_minimum_cover(table,weights,heur,level), "MINCOV     ", F);

/**
    if (debug & EXACT) {
	dump_irredundant(E, Rt, Rp, table);
    }

dump_irredundant(E, Rt, Rp, table);
**/

    int orig_rows = table->nrows;
    cover_last = NIL(sm_row);

    while ((table->ncols > 0) && (table->nrows > 0)) {
/**
       (void) printf("SOLUTION %d\n", n_solutions);
       sm_print(stdout, table);
       (void) printf("BEST COVER :");
       sm_row_print(stdout,cover);
       (void) printf("\n");
**/
       sm_element *re;
       for (re = cover->first_col; re != 0; re = re->next_col) {
           sm_delcol(table, re->col_num); 

           /* if removing element shrinks table, add F term to E */
           if (table->nrows < orig_rows) {
              if (full_pi_cover) {
   	         E = sf_addset(E, GETSET(F, re->col_num));
                 orig_rows = table->nrows;
/**
                 printf("Updating essentials with %s\n", pc1(GETSET(F,re->col_num)));
**/
              }
           }
       }
       if (cover_last != NIL(sm_row))
          sm_row_free(cover_last);
       cover_last= sm_row_dup(cover);
       sm_row_free(cover);

       EXEC(cover=sm_minimum_cover(table,weights,heur,level), "MINCOV    ", F);
       newF = new_cover(100);

       /* load up the essentials */
       foreach_set(E, last, p) { 
   	   newF = sf_addset(newF, p);
       }

       /* add the cubes from partials */
       sm_foreach_row_element(cover_last, pe) { 
   	   newF = sf_addset(newF, GETSET(F, pe->col_num));
       }

       /* Attempt to make the results more sparse */ 
       debug &= ~ (IRRED | SHARP | MINCOV);
       if (! skip_make_sparse && R != 0) {
   	   newF = make_sparse(newF, D, R);
       }

       if (newFhead == (pcover) NULL) {
           newFhead = newF;
           newFprev = newF;
       }
       else {
           /** check to see if newF is equal to newFprev **/
           /** make_sparse can cause this to happen **/
           if (equal_cover(newF, newFprev)) {
              free_cover(newF);
           }
           else {
              newFprev->sfnext = newF;
              newFprev = newF;
              n_solutions++;
           }
       }
       if (n_solutions >=  max_solutions) 
           break;
       if (!full_pi_cover && (table->nrows < orig_rows))
           break;
    }

    if (cover_last == NIL(sm_row))
       cover_last = cover;

    /* if no solutions from PI-table, load up any essentials */
    if (newFhead == (pcover) NULL) {
       newFhead = new_cover(100);
       foreach_set(E, last, p) { 
          newFhead = sf_addset(newFhead, p);
       }
    }

    free_cover(E);
    free_cover(Rt);
    free_cover(Rp);
    sm_free(table);
    sm_row_free(cover_last);
    free_cover(F);
    if (weights != 0) {
	FREE(weights);
    }
    debug = debug_save;

/*
pcover X = newFhead;
int k = 0;
printf("FSOLS:");
while (X != (pcover) NULL) {
    printf("%d ", k++);
    X = X->sfnext;
}
printf("\n");
*/
    return newFhead;
}

#ifndef R_PACKAGE
/* static void */
void
dump_irredundant(pset_family E, pset_family Rt, pset_family Rp, sm_matrix *table)
{
    FILE *fp_pi_table, *fp_primes;
    pPLA PLA;
    pset last, p;
    char *file;

    if (filename == 0 || strcmp(filename, "(stdin)") == 0) {
	fp_pi_table = fp_primes = stdout;
    } else {
	file = ALLOC(char, strlen(filename)+20);
	(void) sprintf(file, "%s.primes", filename);
	if ((fp_primes = fopen(file, "w")) == NULL) {
	    printf("espresso: Unable to open %s\n", file);
	    fp_primes = stdout;
	}
	(void) sprintf(file, "%s.pi", filename);
	if ((fp_pi_table = fopen(file, "w")) == NULL) {
	    printf("espresso: Unable to open %s\n", file);
	    fp_pi_table = stdout;
	}
	FREE(file);
    }

    PLA = new_PLA();
    PLA_labels(PLA);

    fpr_header(fp_primes, PLA, F_type);
    free_PLA(PLA);

    (void) fprintf(fp_primes, "# Essential primes are\n");
    foreach_set(E, last, p) {
	(void) fprintf(fp_primes, "%s\n", pc1(p));
    }
    fprintf(fp_primes, "# Totally redundant primes are\n");
    foreach_set(Rt, last, p) {
	(void) fprintf(fp_primes, "%s\n", pc1(p));
    }
    fprintf(fp_primes, "# Partially redundant primes are\n");
    foreach_set(Rp, last, p) {
	(void) fprintf(fp_primes, "%s\n", pc1(p));
    }
    if (fp_primes != stdout) {
	(void) fclose(fp_primes);
    }
	
    if ((table->nrows == 0) || (table->ncols ==  0))
       (void) fprintf(fp_pi_table, "EMPTY PI TABLE\n");
    else
       sm_print(fp_pi_table, table);
    if (fp_pi_table != stdout) {
	(void) fclose(fp_pi_table);
    }
}
#endif
