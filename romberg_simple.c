// r o m b e r g _ s i m p l e . c
//
//    Robert Thacker
//    Romberg Integration
//    integrate a function - in this case natural log of 10.
//    John Tiller wrote this
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern double y_nat_log (double x);

int main ()
{
    int keep_irow, keep_j;
	double check;
    double T [50+1][51+1];

    FILE *file10 = fopen ("natural_log.out", "w");

//==========================================================
//    Preamble    equation 4
//    set limits of integration and convergence value
//
    double a         =   1.0;
    double b         =  10.0;

    double h         = b - a;

    double tol       = 1.0E-15;
    int niter_max = 26;

    fprintf (file10, "a   = %f\n", a);
    fprintf (file10, "b   = %f\n", b);
    fprintf (file10, "h   = %f\n", h);
    fprintf (file10, "tol = %e\n", tol);

    fprintf (file10, "   Best solution values for each n" \
        "  (i - 1)/i are the rows and j is the column of T" \
        "  n','  i','  j','   T(i,j)/T(i,j) ','      difference  '\n");
//==========================================================

    for (int i = 0; i <= 50; i++)
        for (int j = 1; j <= 51; j++)
            T [i][j] = 0;

//=====================================================

//    equation 10 
    int n = 0;

    double begin  = y_nat_log (a);

    double end    = y_nat_log (b);

    T [n][1] = h * ((begin + end) / 2);
    printf ("n = %d\n", n); // write to screen

//-------------------------------------------------
//    begin romberg iter n = 1,niter_max first column of T
//
//    Trapezoidal Rule with interval halving

    for (n = 1; n <= niter_max; n++) {
        printf ("n = %d\n", n); // write to screen

        double first   = h / pow (2, n - 1);  // see equation 11
        int iend       = (int) pow (2, n) - 1;
        double denom   = pow (2, n);

        double value   = 0.0;

        for (int i = 1; i <= iend; i+= 2) {
            double x = a + h * i / denom;
            value = value + y_nat_log (x);
        } // 20

        T [n][1] = 0.5 * (T [n-1][1] + first * value);  // equation 11

//-----------------------------------------------------
//    perform richardson extrapolation
//    equation 12

        int nrow     = n;            // start with this row of T
        int jcol_end = nrow + 1;     // this is the last column of T

        for (int j = 2; j <= jcol_end; j++) { // j is which column -
            nrow = nrow - 1;
            T [nrow][j] = (pow (4.0, j-1) * T [nrow+1][j-1] - T [nrow][j-1]) / (pow (4.0, j-1) - 1);
        } // 30

//------------------------------------------------------
//
//    check for convergence    
//
        double check_converge = 1.0E10;
        int irow     = n;

        if (n > 3) {
            for (int j = 1; j <= n - 1; j++) {
                if (fabs (T [irow][j] - T [irow-1][j]) <= tol) { // we have converged
                    printf ("solution has converged: see natural_log.out for details\n");
					fprintf (file10, "the run has converged: " \
						 "n        = %5d, " \
						 "irow-1,j = %5d%5d T(%2d,%2d) = %20.15f, " \
						 "irow,j = %5d%5d T(%2d,%2d) = %20.15f" \
						 "           T(%2d,%2d) - T(%2d,%2d) = %20.15f" \
						 " = %20.15e\n",
						 n,
						 irow - 1, j, irow - 1, j, T [irow-1][j],
						 irow, j, irow, j, T [irow][j],
						 irow, j, irow - 1, j,
						 T [irow][j]- T [irow-1][j], T [irow][j] - T [irow-1][j]);

	  				fprintf (file10, "1 234567890123456789012345678901234567890" \
						"                 solution      = %20.15f = T(%2d,%2d)\n",
						T [irow][j], irow, j);
                    exit (0); // goto 799
                }

//==============================================================
//
//    look for best solution in T fr this value of n
//
                check = T [irow][j] - T [irow-1][j];

                if (fabs (check) < fabs (check_converge)) {
                    check_converge = check;
                    keep_irow      = irow;
                    keep_j         =    j;
                }

//==============================================================

                irow = irow - 1;
            } // 35

//===============================================================
//
//    print out best solution for this n

			fprintf (file10, "%3d%3d%3d %20.15f %20.15f %3d%3d %20.15f\n",
				n, keep_irow - 1, keep_j, T [keep_irow-1][keep_j],
				check,
				keep_irow, keep_j, T [keep_irow][keep_j]);

			fprintf (file10, "n first iend denom i j = %5d %20.15f %15d %10.0f %3d %3d %20.12e\n",
				n, first, iend, denom, keep_irow, keep_j, check_converge);
        }
    } // 10
}

//-----------------------------------------------------------------------

//================================
//
//    natural log
//

double y_nat_log (double x)
{
    double y = 1.0 / x;
    return y;
}
