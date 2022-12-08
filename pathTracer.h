#ifndef PATHTRACER_H_INCLUDED
#define PATHTRACER_H_INCLUDED

#include <vector>
#include "lightScatteringDistributionSolver.h"

bool occlusionTrace(vec3 pos, vec3 pathDest, std::vector<plane> &planes, disk &d, std::vector<sphere> &s, unsigned long int planeExclusion, unsigned long diskExclusion, unsigned long int sphereExclusion)
{
    ray occlusionpath = {pos, pathDest-pos};
    vec3 ip;
    double fragDist;
    double qmDist = (pathDest-pos).qmag();

    for(unsigned int i = 0; i < planes.size(); ++i)
    {
        if(planeExclusion != i)
        {
            if(rayPlaneIntersect(occlusionpath, planes[i], ip, fragDist))
            {
                if((ip-pos).qmag() < qmDist)
                    return true;
            }
        }
    }
    for(size_t i = 0; i < s.size(); ++i)
    {
        if(sphereExclusion != i)
        {
            if(raySphereIntersect(occlusionpath, s[i], ip, fragDist))
            {
                if((ip-pos).qmag() < qmDist)
                    return true;
            }
        }
    }

    return false;
}

vec3 reflectionTrace(vec3 rayOrigin, vec3 intersect, vec3 &rNormal, std::vector<plane> &planes, disk &d, std::vector<sphere> &spheres, light &l)
{
    vec3 vIncidence = intersect-rayOrigin;
    vec3 vReflect = vIncidence - rNormal*(2*dot(vIncidence,rNormal));

    ray rReflect = ray(intersect, vReflect);

    vec3 ip;
    double fragDepth = std::numeric_limits<double>::max(); // background fragment is maximally far away
    double depth;

    vec3 normal = vec3(0,1,0);
    vec3 lightDir = vec3(0,0,-1);
    vec3 color = vec3(0,0,0);
    double ldpMag = 0.0;
    double diff = 0.0d;

    // render all planes in the scene
    for(unsigned int i = 0; i < planes.size(); ++i)
    {
        if(rayPlaneIntersect(rReflect,planes[i],ip,depth))
        {
            if(depth < fragDepth) // depth sorting
            {
                fragDepth = depth;
                normal = planes[i].normal;
                lightDir = (l.pos-ip);
                if(occlusionTrace(ip, l.pos, planes, d, spheres, 0, -1, -1))
                    color = vec3();
                else
                {
                    ldpMag = lightDir.qmag();
                    diff = std::max(dot(normal, lightDir*(1.0d/sqrt(ldpMag))), 0.0d) * (1.0d/ldpMag);
                    color = planes[i].mat.color*diff*l.col;
                }
            }
        }
    }

    // render all disks in the scene
    /*if(rayDiskIntersect(rReflect,d,ip,depth))
    {
        if(depth < fragDepth) // depth sorting
        {
            fragDepth = depth;
            normal = d.normal;
            lightDir = (l.pos-ip);
            color = d.color * reflectionTrace(rReflect.origin, ip, normal, p, d, spheres, l);
        }
    }*/

    // render all spheres in the scene
    for(size_t i = 0; i < spheres.size(); ++i)
    {
        if(raySphereIntersect(rReflect,spheres[i],ip,depth))
        {
            if(depth < fragDepth) // depth sorting
            {
                fragDepth = depth;
                normal = (ip-spheres[i].origin).norm();
                lightDir = (l.pos-ip);
                if(occlusionTrace(ip, l.pos, planes, d, spheres, -1, -1, i))
                    color = vec3();
                else
                {
                    ldpMag = lightDir.qmag();
                    diff = std::max(dot(normal, lightDir*(1.0d/sqrt(ldpMag))), 0.0d) * (1.0d/ldpMag);
                    color = spheres[i].mat.color*diff*l.col;
                }
            }
        }
    }

    if(rayLightSphereIntersect(rReflect,l,ip,depth))
    {
        if(depth < fragDepth) // depth sorting
        {
            fragDepth = depth;
            color = l.col;
        }
    }

    return color;
}

vec3 sceneIntersect(ray testRay, std::vector<plane> &planes, disk &d, std::vector<sphere> &spheres, vec3 &ip, vec3 &col)
{
    double fragDepth = std::numeric_limits<double>::max(); // background fragment is maximally far away
    double depth;
    vec3 normal = vec3();

    // render all planes in the scene
    for(unsigned int i = 0; i < planes.size(); ++i)
    {
        if(rayPlaneIntersect(testRay,planes[i],ip,depth))
        {
            if(depth < fragDepth) // depth sorting
            {
                fragDepth = depth;
                normal = planes[i].normal;
                col = planes[i].mat.color;
                //lightDir = (l.pos-ip);
            }
        }
    }

    // render all spheres in the scene
    for(size_t i = 0; i < spheres.size(); ++i)
    {
        if(raySphereIntersect(testRay,spheres[i],ip,depth))
        {
            if(depth < fragDepth) // depth sorting
            {
                fragDepth = depth;
                normal = (ip-spheres[i].origin).norm();
                col = spheres[i].mat.color;
                //lightDir = (l.pos-ip);
            }
        }
    }
    return normal;
}

vec3 irradianceTrace(ray testRay, vec3 color, std::vector<plane> &planes, disk &d, std::vector<sphere> &spheres, sphere refractSphere, std::vector<light> lights, double &pathLength, int recursionLimit, unsigned int traceDepth=0)
{
    vec3 ip;

    if(recursionLimit == 0)
    {
        unsigned int rLight = rand()%lights.size();
        vec3 rDir = sampleSphereUniform();
        ray rRay = ray(lights[rLight].pos, rDir);

        vec3 col;
        vec3 sNormal = sceneIntersect(rRay, planes, d, spheres, ip, col);
        if(sNormal.qmag() > 0.5)
        {
            if(!occlusionTrace(testRay.origin, ip, planes, d, spheres, -1, -1, -1))
            {
                vec3 lightDir = (ip-rRay.origin);
                double ldpMag = lightDir.qmag();
                return color * col * std::max(dot(sNormal, lightDir*(1.0d/sqrt(ldpMag))), 0.0d) * (1/(ip-testRay.origin).qmag());
            }
        }

        return vec3();
    }
    if(recursionLimit > 0)
        --recursionLimit;

    double fragDepth = std::numeric_limits<double>::max(); // background fragment is maximally far away
    double depth;

    vec3 normal = vec3(0,1,0);
    //vec3 lightDir = vec3(0,0,-1);
    material mat;
    int lightIndex = -1;

    // render all planes in the scene
    for(unsigned int i = 0; i < planes.size(); ++i)
    {
        if(rayPlaneIntersect(testRay,planes[i],ip,depth))
        {
            if(depth < fragDepth) // depth sorting
            {
                fragDepth = depth;
                normal = planes[i].normal;
                //lightDir = (l.pos-ip);
                mat = planes[i].mat;
            }
        }
    }

    // render all disks in the scene
    /*if(rayDiskIntersect(testRay,d,ip,depth))
    {
        if(depth < fragDepth) // depth sorting
        {
            fragDepth = depth;
            normal = d.normal;
            //lightDir = (l.pos-ip);
            mat = d.mat;
        }
    }*/

    // render all spheres in the scene
    for(size_t i = 0; i < spheres.size(); ++i)
    {
        if(raySphereIntersect(testRay,spheres[i],ip,depth))
        {
            if(depth < fragDepth) // depth sorting
            {
                fragDepth = depth;
                normal = (ip-spheres[i].origin).norm();
                //lightDir = (l.pos-ip);
                mat = spheres[i].mat;
            }
        }
    }

    // refract light from refraction sphere
    if(raySphereIntersect(testRay,refractSphere,ip,depth))
    {
        if(depth < fragDepth)
        {
            fragDepth = depth;
            normal = (ip-refractSphere.origin).norm();
            double angle = dot(normal, testRay.dir.norm());
            if(angle < 0.0f)
            {
                normal *= -1.0d;
                angle = dot(normal, testRay.dir.norm());
            }
            vec3 refractVector = testRay.dir.norm() + normal*2*angle;
            //pathLength += (ip-testRay.origin).qmag();
            return irradianceTrace({ip+refractVector*0.1, refractVector}, (color*refractSphere.mat.albedo+refractSphere.mat.color), planes, d, spheres, refractSphere, lights, pathLength, recursionLimit, traceDepth+1);
        }
    }

    for(unsigned int i = 0; i < lights.size(); ++i)
    {
        if(rayLightSphereIntersect(testRay,lights[i],ip,depth))
        {
            if(depth < fragDepth) // depth sorting
            {
                fragDepth = depth;
                lightIndex = i;
            }
        }
    }

    if(fragDepth == std::numeric_limits<double>::max())
        return vec3();

    if(lightIndex != -1)
    {
        if(traceDepth != 1)
            return color * lights[lightIndex].col;
        return vec3();
    }

    // randomly sample reflection vector on hemisphere around reflection normal
    vec3 reflectVector = sampleHemisphereUniform(normal);
    /*vec3 vIncidence = (testRay.origin-ip).norm();

    vec3 sideVector = cross(normal,vIncidence).norm();

    vec3 incidenceTangentVector = cross(normal,sideVector);
    vec3 reflectionTangentVector = cross(normal,cross(reflectVector,normal).norm());
    double azimuth = acos(dot(incidenceTangentVector, sideVector)) + acos(dot(reflectionTangentVector, sideVector));

    double theta_r = acos(dot(reflectVector,normal));
    double irrFactor = Li(vIncidence, normal, mat, theta_r, azimuth);
    if(irrFactor < 0.0d)
        return vec3();*/
    double diff = std::max(dot(normal, reflectVector), 0.0d);
    pathLength += (ip-testRay.origin).qmag();
    return irradianceTrace({ip, reflectVector}, (color*mat.color*diff), planes, d, spheres, refractSphere, lights, pathLength, recursionLimit,traceDepth+1);
}

#endif // PATHTRACER_H_INCLUDED
