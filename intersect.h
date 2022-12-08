#ifndef INTERSECT_H_INCLUDED
#define INTERSECT_H_INCLUDED

#include <math.h>

bool raySphereIntersect(ray r, sphere &s, vec3 &intersectionPoint, double &fragDist)
{
    vec3 rOrigAdjust = r.origin - s.origin;

    double a = dot(r.dir,r.dir);
    double b = 2.0d * dot(r.dir, rOrigAdjust);
    double c = dot(rOrigAdjust,rOrigAdjust) - s.qRadius;

    double t1,t2;
    if(!solveQuadratic(a,b,c, t1,t2))
        return false;

    if(t1 < 0.0d)
    {
        if(t2 < 0.0d)
            return false;
        intersectionPoint = r.origin + r.dir*t2;
        fragDist = t2;
        return true;
    }
    if(t2 < 0.0d)
    {
        intersectionPoint = r.origin + r.dir*t1;
        fragDist = t1;
        return true;
    }
    if(t1 < t2)
    {
        intersectionPoint = r.origin + r.dir*t1;
        fragDist = t1;
        return true;
    }
    intersectionPoint = r.origin + r.dir*t2;
    fragDist = t2;
    return true;
}

bool rayLightSphereIntersect(ray r, light &s, vec3 &intersectionPoint, double &fragDist)
{
    vec3 rOrigAdjust = r.origin - s.pos;

    double a = dot(r.dir,r.dir);
    double b = 2.0d * dot(r.dir, rOrigAdjust);
    double c = dot(rOrigAdjust,rOrigAdjust) - s.qRadius;

    double t1,t2;
    if(!solveQuadratic(a,b,c, t1,t2))
        return false;

    if(t1 < 0.0d)
    {
        if(t2 < 0.0d)
            return false;
        intersectionPoint = r.origin + r.dir*t2;
        fragDist = t2;
        return true;
    }
    if(t2 < 0.0d)
    {
        intersectionPoint = r.origin + r.dir*t1;
        fragDist = t1;
        return true;
    }
    if(t1 < t2)
    {
        intersectionPoint = r.origin + r.dir*t1;
        fragDist = t1;
        return true;
    }
    intersectionPoint = r.origin + r.dir*t2;
    fragDist = t2;
    return true;
}

bool rayPlaneIntersect(ray r, plane &p, vec3 &intersectionPoint, double &fragDist)
{
    double denominator = dot(r.dir, p.normal);

    if(std::fabs(denominator) < 0.000001d)
        return false;

    fragDist = dot(p.origin-r.origin, p.normal)/denominator;

    if(fragDist < 0.0f)
        return false;

    intersectionPoint = r.origin + r.dir*fragDist;
    return true;
}

bool rayTriangleIntersect(ray r, triangle t, vec3 &intersectionPoint, double &fragDist)
{
    double denominator = dot(r.dir, t.normal);

    if(std::fabs(denominator) < 0.000001d)
        return false;

    fragDist = dot(t.v0-r.origin, t.normal)/denominator;

    if(fragDist < 0.0f)
        return false;

    intersectionPoint = r.origin + r.dir*fragDist;

    if(dot(intersectionPoint, t.e0) >= 0)
    {
        return true;
    }

    return false;
}

bool rayDiskIntersect(ray r, disk d, vec3 &intersectionPoint, double &fragDist)
{
    double denominator = dot(r.dir, d.normal);

    if(std::fabs(denominator) < 0.000001d)
        return false;

    fragDist = dot(d.origin-r.origin, d.normal)/denominator;

    if(fragDist < 0.0f)
        return false;

    intersectionPoint = r.origin + r.dir*fragDist;

    if((intersectionPoint-d.origin).qmag() < d.qRadius)
        return true;

    return false;
}

#endif // INTERSECT_H_INCLUDED
