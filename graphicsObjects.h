#ifndef GRAPHICSOBJECTS_H_INCLUDED
#define GRAPHICSOBJECTS_H_INCLUDED

struct ray
{
    vec3 origin;
    vec3 dir;
    ray(vec3 orig, vec3 direction) : origin(orig), dir(direction) {}
};

struct light
{
    vec3 pos, col;
    double qRadius;
};

struct material
{
    vec3 color=vec3();
    double albedo=1;
    double sqRoughness=0;
};

struct sphere
{
    vec3 origin;
    double qRadius;
    material mat;
};

struct plane
{
    vec3 origin;
    vec3 normal;
    material mat;
};

struct disk
{
    vec3 origin;
    vec3 normal;
    double qRadius;
    material mat;
};

struct triangle
{
    vec3 v0,v1,v2;
    material mat;
    vec3 normal=vec3();
    vec3 e0=vec3(),e1=vec3(),e2=vec3();

    void preCompute()
    {
        this->normal = cross((v1-v0).norm(), (v2-v0).norm());
        this->e0 = cross(this->normal, this->v1 - this->v0);
        this->e1 = cross(this->normal, this->v0 - this->v2);
        this->e2 = cross(this->normal, this->v2 - this->v1);
    }
};

#endif // GRAPHICSOBJECTS_H_INCLUDED
