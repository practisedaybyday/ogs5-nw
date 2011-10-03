//-------------------------------------------------------------------
// $Id: num_methods.h 705 2006-04-28 19:39:01Z gems $
//
// C/C++ Numerical Methods used in GEMS-PSI and GEMS3K
// (c) 2006,2011 S.Dmytriyeva, D.Kulik
//
// This file is part of a GEM-Selektor library for thermodynamic
// modelling by Gibbs energy minimization and of the GEMIPM2K code
//
// This file may be distributed under the terms of the GEMS-PSI
// QA Licence (GEMSPSI.QAL)
//
// See http://gems.web.psi.ch/ for more information
// E-mail: gems2.support@psi.ch
//-------------------------------------------------------------------

#ifndef _num_methods_h_
#define _num_methods_h_

#include <math.h>

// Calculate number of points from iterators
long int  getNpoints( double Tai[4] );
double    getStep( double *Tai, int nPoints );


// Lagrangian interpolation functions
double LagranInterp(float *y, float *x, double *d, float yoi,
		float xoi, int M, int N, int pp );
double LagranInterp(float *y, float *x, float *d, float yoi,
		float xoi, int M, int N, int pp );
double LagranInterp(double *y, double *x, double *d, double yoi,
		double xoi, long int M, long int N, long int pp );


// generic functions for calculating partial derivatives
double quot( double u, double v, double du, double dv );

double quot( double u, double v, double du, double dv, double d2u, double d2v );


double prod2( double u, double v, double du, double dv );

double prod2 ( double u, double v, double du, double dv, double d2u, double d2v );


double prod3 ( double u, double v, double w, double du, double dv, double dw );

double prod3 ( double u, double v, double w, double du, double dv, double dw,
		double d2u, double d2v, double d2w );


typedef double (*minFunction)(double x, double y );

struct GoldenSelectionData
{
    double Fa;
    double Fb;
    double Fx1;
    double Fx2;
    double a;
    double b;
    double x1;
    double x2;
    double Xtol;


    GoldenSelectionData( double xstart, double xend, double Xtol_)
    {
      Xtol = Xtol_;
      x1 = xstart;
      x2 = xend;
      a = min( x1, x2 );
      b = max( x1, x2 );
   }

 };


// Class for minimization of convex one parameter function ( f(x)=>0 )
// Method of Gold Selection
class GoldenSelection
{
protected:

  GoldenSelectionData dat1;
  double Ftol;
  minFunction minF;

public:

  // Golden Selection in interval x1 to x2, to minimize function f_proc
  // xtol, ftol tolerance for the parameter and function
  GoldenSelection( double x1, double x2, double xtol, double ftol,
                   double (f_proc)(double val, double val2 )):
   dat1(x1,x2,xtol),Ftol(ftol)
  {
    minF = f_proc;
  }

  virtual double calcFunction( double x, double y )
  {
      return minF( x, y );
  }

  virtual double getMinimum( double val2=0 )
  {
     return getMinimumDat( dat1, val2 );
  }

  virtual double getMinimumDat( GoldenSelectionData dat, double val2 );

};

// Class for minimization of convex two parameter function ( f(x,y)=>0 )
// Method of Gold Selection
class GoldenSelectionTwo : public GoldenSelection
{
  GoldenSelectionData dat2;
  int nOperand;   // minimization for first or second parameter
  double minX;
  double minY;

public:

  // Golden Selection in intervals x1 to x2, y1 to y2 to minimize function f_proc
  // xtol, ytol, ftol tolerance for the parameters and function
  GoldenSelectionTwo( double x1, double x2, double xtol,
                      double y1, double y2, double ytol,
                      double ftol, double (f_proc)(double val, double val2 )):
  GoldenSelection( x1, x2, xtol, ftol, f_proc),  dat2(y1,y2,ytol)
  {}

  double calcFunction( double x, double y );
  double getMinimum( double val2=0 );

  double getMinX()
  { return minX; }
  double getMinY()
  { return minY; }

};


// Golden Selection in interval x1 to x2, to minimize function f_proc
//double GoldenSelection( double x1, double x2, double xtol, double ftol,
//                        double val2, double (f_proc)(double val, double val2 ));
// Method of Gold Selection for two parameters
//   ystart, yend, ytol,  parameter y    - start, end, tolerance
//   xstart, xend, xtol,  parameter x    - start, end, tolerance
//   ftol  function tolerance
//   ymin and xmin return values
//   f_proc function to minimize ( f( y, x)=>0 )
//void GoldenSelection2( double ystart, double yend, double Ytol,
//                       double xstart, double xend, double Xtol,
//                       double Ftol, double& ymin, double& xmin,
//                       double (f_proc)(double valy, double valx ));



#endif   // _num_methods_h_

//-----------------------End of num_methods.h--------------------------

