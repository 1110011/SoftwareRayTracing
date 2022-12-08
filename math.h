#ifndef MATH_H_INCLUDED
#define MATH_H_INCLUDED

#include <math.h>

bool solveQuadratic(double a, double b, double c, double &x1, double &x2)
{
    // divide out quadratic term
    b /= a;
    c /= a;

    double pHalf = b/2;

    // compute square of half distance between roots
    double radicand = pHalf*pHalf - c;

    // quadratic has no real solutions
    if(radicand < 0.0d)
        return false;

    // quadratic has just on solution, or the solutions are negligably close together
    if(std::fabs(radicand) < 0.000001d)
    {
        x1 = -pHalf;
        x2 = -pHalf;
        return true;
    }

    // compute solutions
    double root = sqrt(radicand);
    x1 = -pHalf + root;
    x2 = -pHalf - root;

    return true;
}

#endif // MATH_H_INCLUDED
