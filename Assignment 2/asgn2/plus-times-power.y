/*  yacc  file for PLUS-TIMES-POWER expressions, for copying and then
    modification (of the copy) to make a quaternionic calculator
    It already contains C functions for use by the quaternionic calculator.
    Graham Farr, Monash University
    Last updated:  11 September 2021
*/

        /*   declarations   */


%{
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "quat.h"
int yylex(void);
void yyerror(char *);
int yydebug=0;   /*  set to 1 if using with yacc's debug/verbose flags   */
int  printQuaternion(Quaternion);
Quaternion  newQuaternion(double, double, double, double);
Quaternion  realPart(Quaternion);
Quaternion  imaginaryPart(Quaternion);
double  normSquared(Quaternion);
double  norm(Quaternion);
Quaternion  conjugate(Quaternion);
Quaternion  scalarTimesQuaternion(double, Quaternion);
Quaternion  unitQuaternion(Quaternion);
Quaternion  inverse(Quaternion);
Quaternion  sum(Quaternion, Quaternion);
Quaternion  diff(Quaternion, Quaternion);
Quaternion  product(Quaternion, Quaternion);
Quaternion  quotient(Quaternion, Quaternion);
Quaternion  power(Quaternion, int);
Quaternion  rotation(double, Quaternion);
/*
int  total = 0;
*/
%}

%union {
    int  iValue;
    char *str;
    Quaternion qtn;
    double num;
};

%token  <num>  NUMBER
%token  <iValue>  WHOLENUMBER
%token  <str>  POWER

%left '+'
%left '*'

%type  <iValue>  start
%type  <iValue>  line
%type  <iValue>  expr
%type  <iValue>  int

%start  start



%%       /*   rules section   */



start    :    line  '\n'       {  }
         |    start  line  '\n'          {  }
         ;

line     :    expr   {  printf("%d\n", $1);   }
         |        /*  allow "empty" expression  */           {     }
         

expr     :    int    {  $$ = $1;  }
         |    POWER '(' expr ',' expr ')'    {  $$ = pow($3,$5);  }
         |    expr '*' expr    {  $$ = $1 * $3;  }
         |    expr '+' expr    {  $$ = $1 + $3;  }
         ;

int      :    NUMBER             {  $$ = $1;   }
         ;


%%      /*   programs   */




int  printQuaternion(Quaternion q)
/*  print the quaternion  q  in the standard form  w + x i + y j + z k  */
{
  return  printf("%f + %f i + %f j + %f k\n", q.w, q.x, q.y, q.z);
}


Quaternion  newQuaternion(double w, double x, double y, double z)
/*  given four real numbers  w, x, y, z,  returns the quaternion
    w + x i + y j + z k
*/
{
  Quaternion  q;

  q.w = w;
  q.x = x;
  q.y = y;
  q.z = z;

  return  q;
}


Quaternion  realPart(Quaternion q)
/*  returns the real part of the quaternion  q.
    If  q = w + xi + yj + zk, then its real part is just  w.
*/
{
  return  newQuaternion(q.w, 0, 0, 0);
}


Quaternion  imaginaryPart(Quaternion q)
/*  returns the imaginary part of the quaternion  q.
    If  q = w + xi + yj + zk, then its imaginary part is
    the pure quaternion  xi + yj + zk.
*/
{
  return  newQuaternion(0, q.x, q.y, q.z);
}


double  normSquared(Quaternion q)
/*  returns the sum of the squares of the components of the quaternion  q.
    This is the square of the norm (i.e., of the length) of  q.
*/
{
  return  q.w*q.w + q.x*q.x + q.y*q.y + q.z*q.z;
}


double  norm(Quaternion q)
/*  returns the norm (i.e., length) of the quaternion  q.  */
{
  return  sqrt(normSquared(q));
}


Quaternion  conjugate(Quaternion q)
/*  returns the conjugate of the quaternion  q.
    This is obtained by negating the imaginary part and leaving the
    real part unchanged.
*/
{
  return  newQuaternion(q.w, -q.x, -q.y, -q.z);
}


Quaternion  scalarTimesQuaternion(double s, Quaternion q)
/*  returns the quaternion obtained by multiplying  q  by the scalar  s.
    Think of  s  as a real number, and we just multiply each component
    of  q  by  s.
*/
{
  return  newQuaternion(s*q.w, s*q.x, s*q.y, s*q.z);
}


Quaternion  unitQuaternion(Quaternion q)
/*  returns the unit quaternion in the direction of  q, obtained by
    dividing  q  by its length.
*/
{
  return  scalarTimesQuaternion(1/norm(q),q);
}


Quaternion  inverse(Quaternion q)
/*  returns the multiplicative inverse  1/q  of the quaternion  q.  */
{
  return  scalarTimesQuaternion(1.0/normSquared(q),conjugate(q));
}


Quaternion  sum(Quaternion q1, Quaternion q2)
/*  returns the sum  q1 + q2  of the quaternions  q1  and  q2.
    This may be viewed as just a vector sum.
*/
{
  Quaternion  qsum;

  qsum.w = q1.w + q2.w;
  qsum.x = q1.x + q2.x;
  qsum.y = q1.y + q2.y;
  qsum.z = q1.z + q2.z;

  return  qsum;
}


Quaternion  diff(Quaternion q1, Quaternion q2)
/*  returns the difference  q1 - q2  of the quaternions  q1  and  q2.
    This may be viewed as just vector subtraction.
*/
{
  Quaternion  qdiff;

  qdiff.w = q1.w - q2.w;
  qdiff.x = q1.x - q2.x;
  qdiff.y = q1.y - q2.y;
  qdiff.z = q1.z - q2.z;

  return  qdiff;
}


Quaternion  product(Quaternion q1, Quaternion q2)
/*  returns the product  q1 * q2  of the quaternions  q1  and  q2.
    This is *not* the same as the dot or cross product of vectors,
    although the three products are closely related.
*/
{
  Quaternion  qprod;

  qprod.w = q1.w*q2.w - q1.x*q2.x - q1.y*q2.y - q1.z*q2.z;
  qprod.x = q1.w*q2.x + q1.x*q2.w + q1.y*q2.z - q1.z*q2.y;
  qprod.y = q1.w*q2.y + q1.y*q2.w + q1.z*q2.x - q1.x*q2.z;
  qprod.z = q1.w*q2.z + q1.z*q2.w + q1.x*q2.y - q1.y*q2.x;

  return  qprod;
}


Quaternion  quotient(Quaternion q1, Quaternion q2)
/*  returns the quotient  q1 / q2  of the quaternions  q1  and  q2,
    calculated as  q1 * (1/q2).
*/
{
  Quaternion  qquotient;

  return  product(q1,inverse(q2));
}


Quaternion  power(Quaternion q, int n)
/*  returns the result of raising quaternion  q  to power  n.
*/
{
  Quaternion  qpower;
  int  i;

  qpower = newQuaternion(1,0,0,0);
  for  (i=0; i < n; i++)
  {
    qpower = product(qpower, q);
  }
  return  qpower;
}



Quaternion  rotation(double angleInDegrees, Quaternion direction)
/*  Arguments:
       -  an angle in degrees,
       -  a quaternion.
    The imaginary part of the supplied quaternion is viewed as a vector
    indicating an axis in 3D space (the "axis vector").  Its real part
    is ignored.
    Value returned:
    a quaternion that represents a 3D rotation by the specified
    angle, clockwise (as viewed from the origin looking outwards along the
    vector), around the axis specified by the axis vector.
*/
{
  double  angleInRadians;
  Quaternion qHat;

  qHat = unitQuaternion(imaginaryPart(direction));

  angleInRadians = angleInDegrees*M_PI/180;
  return  newQuaternion(cos(angleInRadians/2),
                        sin(angleInRadians/2)*qHat.x,
                        sin(angleInRadians/2)*qHat.y,
                        sin(angleInRadians/2)*qHat.z
                       );
}



/*
void yyerror(char *s) {
      fprintf(stderr, "%s\n", s);
      fprintf(stderr, "line %d: %s\n", yylineno, s);
}
*/


int main(void) {
    yyparse();
    return 0;
}
