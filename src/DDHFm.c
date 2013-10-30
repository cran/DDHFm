#include <R.h>
#include <math.h>

/*
 * C source code for replacing R code in ddhft.np.2 and isotone
 *
 * Author: G P Nason, August 2005
 *
 */


/*
 * which_minvC
 *
 * Find index where vector v is minimum, store answer in ans
 *
 * Variables:
 *
 * v the vector to search for the minimum value (not changed by function)
 * n the length of the vector v (not changed by function)
 * ans: integer storage for the minimum index (changed)
 *
 * Input: v,n
 * Output: ans
 *
 */

void which_minvC(double *v, int *n, int *ans)
{
int i;

int whichix=0;
double minv = *v;

if (*n > 1)
	for(i=1; i < *n; ++i)
		if (*(v+i) < minv)	{
			minv = *(v+i);
			whichix = i;
		}

*ans = whichix;	/* Need +1 in conversion from C to R */

}

/*
 * Compute the absolute difference between a and lvect
 * Return the index of where the smallest difference is
 *
 * Variables:
 *
 * a: a single double number(not changed)
 * vect: the other vector (not changed)
 * lvect: the length of vect (not changed)
 * ans: single integer,
 * 	where the index of the min |a-vect| is (changed by this function)
 *
 * Input: a,vect,lvect
 * Output: ans 
 */

void which_min_diffC(double *a, double *vect, int *lvect, int *ans)
{
/*
int i;
double *fabsvect;
void which_minvC();

fabsvect = (double *)Calloc(*lvect, double);

for(i=0; i < *lvect; ++i)
	*(fabsvect+i) = fabs(*a - *(vect+i));

which_minvC(fabsvect, lvect, ans);

Free(fabsvect);
*/

register int i;
double minv;
double tmp;

minv = fabs(*a - *vect);	
*ans = 0;

for(i=1; i < *lvect; ++i)	{
	tmp = fabs(*a - *(vect+i));
	if (tmp < minv)	{
		minv = tmp;
		*ans = i;
		}
	}
}

void function_from_vectorC(double *x, double *y, int *lx, double *argvec,
		double *ansvec, int *lav)
{
int i;
int ans;

void which_min_diffC();

for(i=0; i < *lav; ++i)	{
	which_min_diffC(argvec+i, x, lx, &ans);
	*(ansvec+i) = *(y+ans);	/* Convert from R to C */
	}
}


void CentralDDHFT(double *sm,
		double *det,
		double *mu,
		double *sigma, 
		int *nhalf,
		double *hft,
		double *factors,
		int *n,
		int *J)
{
int i,j;
int musiglength;
void function_from_vectorC();
double *v;


musiglength = *nhalf;


v = (double *)Calloc(*n, double);

for(i=1; i <= *J; ++i)	{

	for(j=0; j < *nhalf; ++j)	{

		*(sm+j) = (*(hft + 2*j) + *(hft + 2*j+1))/2.0;


		*(det+j) = (*(hft + 2*j) - *(hft + 2*j+1))/2.0;
		}

	function_from_vectorC(mu, sigma, &musiglength, sm, v, nhalf); 

	for(j=0; j < *nhalf; ++j)	{
		if (*(v+j) > 0)	{
			*(det+j) = *(det+j)/ *(v+j);
			}
		}

	for(j=0; j < *nhalf; ++j)	{
		*(hft+j) = *(sm+j);
		*(hft+*nhalf+j) = *(det+j);
		*(factors+*nhalf+j) = *(v+j);
		}
	*n = *n/2;

	*nhalf = *n/2;
	}


*nhalf = 1;
*n = 2;

for(i=1; i <= *J; ++i)	{

	for(j=0; j < *nhalf; ++j)	{
		*(sm+j) = *(hft+j);
		*(det+j) = *(hft+*nhalf+j);
		}

	for(j=0; j < *nhalf; ++j)	{

		*(hft+2*j) = *(sm+j) + *(det+j);

	
		*(hft+2*j+1) = *(sm+j) - *(det+j);

		}
	*nhalf = *n;
	*n = *n * 2;
	}


Free(v);
}

void isotoneC(double *x, double *wt, int *nn, int *increasing) 
{
int i,j,jb;
int nx, ndx;
int indl,indr;
int ipcnt, nuxcnt;
int *jmax, *jmin, jmaxcnt=0, jmincnt=0;
int *ip;
double tmp;
double *Lx;
double *dx, wtn;
int *killx;

double minvC();

if (*nn == 1)
	return;

killx = (int *)Calloc( *nn, int);
ip = (int *)Calloc( *nn, int);

for(i=0; i < *nn; ++i)
	*(ip+i) = i;

if (!*increasing)	{
	for(i=0; i < *nn; ++i)
		*(x+i) = - *(x+i);
	}

ndx = *nn-1;
dx = (double *)Calloc(ndx, double);

/* Do diff */

for(i=0; i < *nn-1; ++i)
	*(dx+i) = *(x+i+1) - *(x+i);

/* Make nx length of x */

nx = *nn;

/* Get space for two indicator arrays used later */
jmax = (int *)Calloc(nx, int);
jmin = (int *)Calloc(nx, int);

while( (nx > 1) && (minvC(dx, &ndx) < 0))	{

	/*
	Rprintf("nx is %d\n", nx);

	Rprintf("ndx is %d\n", ndx);
	for(j=0; j < ndx; ++j)
		Rprintf("diff[%d] = %lf\n", j, *(dx+j));
	*/


	jmaxcnt = jmincnt = 0;

	for(j=0; j < nx; ++j)	{	

		/* jmax comp */

		if (j==0)	{
			if ((*(dx+j) <= 0.0))	{
				*(jmax + jmaxcnt++) = j;
				/*
				Rprintf("jmax[%d]=%d\n", jmaxcnt-1,*(jmax+jmaxcnt-1));
				*/
				}
			}
		else if (j != (nx-1))	{
			if ( (*(dx+j) <= 0) && (*(dx+j-1) > 0))	{
				*(jmax + jmaxcnt++) = j;
				/*
				Rprintf("jmax[%d]=%d\n", jmaxcnt-1,*(jmax+jmaxcnt-1));
				*/
				}
			}

		/* jmin comp */

		if (j == (nx-1))	{
			/*
			Rprintf("Got to first jmin bit and j is %d\n", j);
			*/
			if (*(dx+j-1) <= 0)	{
			       *(jmin + jmincnt++) = j;
			       /*
				Rprintf("jmin[%d]=%d\n", jmincnt-1,*(jmin+jmincnt-1));
				*/
			}
			}

	        else if (j != 0)	{
			/*
			Rprintf("Got to jmin bit and j is %d\n", j);
			*/
			if ( (*(dx+j) > 0) && (*(dx+j-1) <= 0))	{
				*(jmin + jmincnt++) = j;
				/*
				Rprintf("jmin[%d]=%d\n", jmincnt-1,*(jmin+jmincnt-1));
				*/
			}
			}
		}
	for(jb=0; jb < jmaxcnt; ++jb)	{

		indl = *(jmax+jb);
		indr = *(jmin+jb);

		/*
		Rprintf("indl is %d, indr is %d\n", indl, indr);
		*/

		wtn = 0.0;

		for(j=indl; j <= indr; ++j)
			wtn += *(wt+j);

		tmp = 0.0;

		for(j=indl; j <= indr; ++j)
			tmp += *(wt+j) * *(x+j);

		*(x+indl) = tmp/wtn;
		*(wt+indl) = wtn;

		for(j=indl+1; j <= indr; ++j)
			*(killx+j) = 1;
		}

	/* Now copy non-killed x,wt */

	nuxcnt = 0;

	for(j=0; j < nx; ++j)	{
		/*
		Rprintf("killx[%d] = %d\n", j, *(killx+j));
		*/
		if (*(killx+j) == 0) { /* Copy */
			*(x+nuxcnt) = *(x+j);
			*(ip+nuxcnt) = *(ip+j);
			*(wt+nuxcnt++) = *(wt+j);
			}
		}

	Free(killx);
	
	/* Now make memory for new Lx and copy across */
	nx = nuxcnt;
	killx = (int *)Calloc( nx, double);


	/* Now make memory for new dx and make them */

	Free(dx);	
	ndx = nx-1;
	dx = Calloc(ndx, double);

	/* Do diff */

	for(i=0; i < ndx; ++i)
		*(dx+i) = *(x+i+1) - *(x+i);

			
	/* Free jmax and jmin */

	}

Lx = (double *)Calloc( *nn, double);

nuxcnt = 0;

if (nx==1)	{
	for(j=0; j < *nn; ++j)
		*(Lx+j) = *x;
	}
else	{

	ipcnt = 1;

	for(j = 0; j < *nn; ++j)	{

		if (ipcnt < nx && *(ip+ipcnt) == j)	{
			++ipcnt;
			++nuxcnt;
			}
		*(Lx+j) = *(x+nuxcnt);
		/*
		*(x+j) = *(Lx+nuxcnt);
		*/
		}
	}

if (!*increasing)
	for(j=0; j < *nn; ++j)
		*(x+j) = - *(Lx+j);
else
	for(j=0; j < *nn; ++j)
		*(x+j) = *(Lx+j);

/* Free stuff */ 

Free(ip);
Free(Lx);
Free(dx);
Free(killx);
Free(jmax);
Free(jmin);
}

/* Find minimum value of vector, R min function */

double minvC(double *v, int *n)
{
int i;

double themin=*v;

for(i=1; i < *n; ++i)
	if ( *(v+i) < themin)
		themin = *(v+i);

return(themin);
}
