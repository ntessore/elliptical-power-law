// ELLIPTICAL POWER LAW PROFILE LENS
// SAMPLE IMPLEMENTATION
// 
// author: Nicolas Tessore
// date:   22 Jan 2016
// 
// You are free to use this code in any way you like, but if you do,
// please cite Tessore & Metcalf (2015) in your resulting publications.

// necessary for cos, sin
#include <math.h>

// Compute the angular dependency `omega` of the deflection angle of an
// elliptical power law profile lens with axis ratio `q` and slope `t`
// up to order `N`. Takes the elliptical angle `phi` in radians and writes
// the two components to `omega[0]` and `omega[1]`.
void epl_omega(double phi, double q, double t, int N, double omega[])
{
    // utility variables
    const double T = 2 - t;
    const double f = (1 - q)/(1 + q);
    
    // cosines and sines of phi
    const double c = cos(phi);
    const double s = sin(phi);
    const double c2 = c*c - s*s;
    const double s2 = 2*c*s;
    
    // vectors for iteration
    double A[2];
    double B[2];
    
    // start of iteration
    omega[0] = A[0] = c;
    omega[1] = A[1] = s;
    
    // compute deflection up to N'th order
    for(int k = 1; k <= N; ++k)
    {
        // rotate
        B[0] = c2*A[0] - s2*A[1];
        B[1] = s2*A[0] + c2*A[1];
        
        // iterate
        omega[0] += A[0] = -f*(2*k - T)/(2*k + T)*B[0];
        omega[1] += A[1] = -f*(2*k - T)/(2*k + T)*B[1];
    }
}

//----8<----------------------------------------------------------------------
// simple main routine to make a table of (phi, omega)
//
#include <stdlib.h>
#include <stdio.h>
int main(int argc, char* argv[])
{
    double omega[2];
    
    double q, t;
    
    if(argc < 3)
    {
        fprintf(stderr, "usage: epl <q> <t>\n");
        return EXIT_FAILURE;
    }
    
    q = atof(argv[1]); 
    t = atof(argv[2]);
    
    printf("%20s %20s %20s\n", "phi [deg]", "omega_1", "omega_2");
    for(int d = 0; d <= 360; d += 5)
    {
        epl_omega(d*M_PI/180, q, t, 10, omega);
        printf("%20d %20f %20f\n", d, omega[0], omega[1]);
    }
    
    return EXIT_SUCCESS;
}
//
// end of main
//----------------------------------------------------------------------------
