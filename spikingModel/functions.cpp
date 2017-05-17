#include "functions.h"

double r8_uniform_01 ( int *seed )
{
  int i4_huge = 2147483647;
  int k;
  double r;
  k = *seed / 127773;
  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;
  if ( *seed < 0 )
  {
    *seed = *seed + i4_huge;
  }
  r = ( double ) ( *seed ) * 4.656612875E-10;
  return r;
}


double r8_normal_01 ( int *seed )
{
  double pi = 3.141592653589793;
  double r1;
  double r2;
  static int used = -1;
  double x;
  static double y = 0.0;

  if ( used == -1 )
  {
    used = 0;
  }
  if ( ( used % 2 )== 0 )
  {
    for ( ; ; )
    {
      r1 = r8_uniform_01 ( seed );
      if ( r1 != 0.0 )
      {
        break;
      }
    }
    r2 = r8_uniform_01 ( seed );
    x = sqrt ( -2.0 * log ( r1 ) ) * cos ( 2.0 * pi * r2 );
    y = sqrt ( -2.0 * log ( r1 ) ) * sin ( 2.0 * pi * r2 );
  }
  else
  {
    x = y;
  }
  used = used + 1;
  return x;
}

double r8_exp_01 ( int *seed , double lambda)
{
  return (-lambda) * log(r8_uniform_01(seed));
}