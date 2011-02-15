#ifndef __SUPERLU_INTERFACE /* allow multiple inclusions */
#define __SUPERLU_INTERFACE

#include "slu_zdefs.h"

/*! \brief Driver routines */
extern void
dgssv(superlu_options_t *, SuperMatrix *, int *, int *, SuperMatrix *,
      SuperMatrix *, SuperMatrix *, SuperLUStat_t *, int *);
extern void
dgssvx(superlu_options_t *, SuperMatrix *, int *, int *, int *,
       char *, double *, double *, SuperMatrix *, SuperMatrix *,
       void *, int, SuperMatrix *, SuperMatrix *,
       double *, double *, double *, double *,
       mem_usage_t *, SuperLUStat_t *, int *);

/*! \brief Supernodal LU factor related */
extern void
dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_CompRow_Matrix(SuperMatrix *, int, int, int, double *,
		       int *, int *, Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_CompCol_Matrix(SuperMatrix *, SuperMatrix *);
extern void
dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
		     Stype_t, Dtype_t, Mtype_t);
extern void
dCreate_SuperNode_Matrix(SuperMatrix *, int, int, int, double *, 
		         int *, int *, int *, int *, int *,
			 Stype_t, Dtype_t, Mtype_t);
extern void
dCopy_Dense_Matrix(int, int, double *, int, double *, int);

extern void    dgstrf (superlu_options_t*, SuperMatrix*, double, 
                       int, int, int*, void *, int, int *, int *, 
                       SuperMatrix *, SuperMatrix *, SuperLUStat_t*, int *);
extern int     dsnode_dfs (const int, const int, const int *, const int *,
			     const int *, int *, int *, GlobalLU_t *);
extern int     dsnode_bmod (const int, const int, const int, double *,
                              double *, GlobalLU_t *, SuperLUStat_t*);
extern void    dpanel_dfs (const int, const int, const int, SuperMatrix *,
			   int *, int *, double *, int *, int *, int *,
			   int *, int *, int *, int *, GlobalLU_t *);
extern void    dpanel_bmod (const int, const int, const int, const int,
                           double *, double *, int *, int *,
			   GlobalLU_t *, SuperLUStat_t*);
extern int     dcolumn_dfs (const int, const int, int *, int *, int *, int *,
			   int *, int *, int *, int *, int *, GlobalLU_t *);
extern int     dcolumn_bmod (const int, const int, double *,
			   double *, int *, int *, int,
                           GlobalLU_t *, SuperLUStat_t*);
extern int     dcopy_to_ucol (int, int, int *, int *, int *,
                              double *, GlobalLU_t *);         
extern int     dpivotL (const int, const double, int *, int *, 
                         int *, int *, int *, GlobalLU_t *, SuperLUStat_t*);
extern void    dpruneL (const int, const int *, const int, const int,
			  const int *, const int *, int *, GlobalLU_t *);
extern void    dreadmt (int *, int *, int *, double **, int **, int **);
extern void    dGenXtrue (int, int, double *, int);
extern void    dFillRHS (trans_t, int, double *, int, SuperMatrix *,
			  SuperMatrix *);
extern void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                        SuperMatrix *, SuperLUStat_t*, int *);


/*! \brief Driver related */

extern void    dgsequ (SuperMatrix *, double *, double *, double *,
			double *, double *, int *);
extern void    dlaqgs (SuperMatrix *, double *, double *, double,
                        double, double, char *);
extern void    dgscon (char *, SuperMatrix *, SuperMatrix *, 
		         double, double *, SuperLUStat_t*, int *);
extern double   dPivotGrowth(int, SuperMatrix *, int *, 
                            SuperMatrix *, SuperMatrix *);
extern void    dgsrfs (trans_t, SuperMatrix *, SuperMatrix *,
                       SuperMatrix *, int *, int *, char *, double *, 
                       double *, SuperMatrix *, SuperMatrix *,
                       double *, double *, SuperLUStat_t*, int *);

extern int     sp_dtrsv (char *, char *, char *, SuperMatrix *,
			SuperMatrix *, double *, SuperLUStat_t*, int *);
extern int     sp_dgemv (char *, double, SuperMatrix *, double *,
			int, double, double *, int);

extern int     sp_dgemm (char *, char *, int, int, int, double,
			SuperMatrix *, double *, int, double, 
			double *, int);

/*! \brief Memory-related */
extern int     dLUMemInit (fact_t, void *, int, int, int, int, int,
			     SuperMatrix *, SuperMatrix *,
			     GlobalLU_t *, int **, double **);
extern void    dSetRWork (int, int, double *, double **, double **);
extern void    dLUWorkFree (int *, double *, GlobalLU_t *);
extern int     dLUMemXpand (int, int, MemType, int *, GlobalLU_t *);

extern double  *doubleMalloc(int);
extern double  *doubleCalloc(int);
extern int     dmemory_usage(const int, const int, const int, const int);
extern int     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);

/*! \brief Auxiliary routines */
extern void    dreadhb(int *, int *, int *, double **, int **, int **);
extern void    dCompRow_to_CompCol(int, int, int, double*, int*, int*,
		                   double **, int **, int **);
extern void    dfill (double *, int, double);
extern void    dinf_norm_error (int, SuperMatrix *, double *);

/*! \brief Routines for debugging */
extern void    dPrint_CompCol_Matrix(char *, SuperMatrix *);
extern void    dPrint_SuperNode_Matrix(char *, SuperMatrix *);
extern void    dPrint_Dense_Matrix(char *, SuperMatrix *);

#endif /* __SUPERLU_INTERFACE */
