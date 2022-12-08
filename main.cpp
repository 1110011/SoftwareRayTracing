#include "fileExport.h"
#include "math.h"
#include "vectorMath.h"
#include "graphicsObjects.h"
#include "intersect.h"
#include "pathTracer.h"

#include <limits>
#include <vector>

#include <chrono>
#include <omp.h>

#define imgSize 2048
#define imgSizeHalf 1024.0d

#define numSamples 64

vec3 toneMap(vec3 x)
{
    return ((x*(x*0.15+0.05)+0.004)/(x*(x*0.15+0.50)+0.06))-0.0666666666d;
}

int main()
{

    material matPlane = {vec3(0.5,0.4,0.1), 0.9, 0.54};
    material graySphere = {vec3(0.9,0.8,0.8), 0.9, 0.54};
    material greenSphere = {vec3(0.1,1,0.1), 0.5, 0.54};
    material blueSphere = {vec3(0.1,0.1,1), 0.9, 0.54};
    material matMirror = {vec3(1,1,1), 0.9, 0.8};

    material matGlass = {vec3(0.05,0.05,0.05), 0.98, 0.0};

    std::vector<plane> planes = {{vec3(0,-1,0), vec3(0,1,0), matPlane},
                                 {vec3(0,0,-4), vec3(0,0,1), matPlane}};

    std::vector<sphere> spheres;
    spheres.push_back({vec3(0,0,-3), 1.0d, graySphere});
    spheres.push_back({vec3(-3.0,-0.2,-1.5), 0.25d, greenSphere});
    spheres.push_back({vec3(1.0,-0.2,2.0), 0.5d, blueSphere});

    sphere refractSphere = {vec3(-2.0,-0.5,0.0), 0.4d, matGlass};

    disk d = {vec3(-0.25,-0.6,2), vec3(0,1,0), 2.0, graySphere};

    /*triangle t = {vec3(-1,1,0), vec3(-1,-1,0), vec3(1,1,0), vec3(1,0,1)};
    t.preCompute();
    t.e0.print();*/

    std::vector<light> lights = {{vec3(-2.0, 3.0, -3.0), vec3(1,1,1)*0.8f, 0.75},
                                 {vec3(-2.0, -0.5, 2.5), vec3(0.3,1.0,0.3)*0.3f, 0.25}};

    auto tp1 = std::chrono::steady_clock::now();

    unsigned char* buf = new unsigned char[imgSize*imgSize*3];

    // openMP parallelization directive
    int maxThreads = omp_get_max_threads();
    #pragma omp parallel for num_threads(maxThreads)
    for(unsigned long int y = 0; y < imgSize; ++y)
    {
        for(unsigned long int x = 0; x < imgSize; ++x)
        {
            ray testRay = ray(vec3(0,0,4),
                              vec3((double)x/imgSizeHalf-1.0d,
                                   (double)y/imgSizeHalf-1.0d,
                                   -1.0d));

            vec3 hdrColor = vec3(0,0,0);
            vec3 ip;
            double fragDepth = std::numeric_limits<double>::max(); // background fragment is maximally far away
            double depth;

            vec3 normal = vec3(0,1,0);
            vec3 lightDir = vec3(0,0,-1);

            bool tmp = false;
            double ldpMag = 0.0d;
            double diff = 0.0d;
            double spec = 0.0f;

            // render all planes in the scene
            for(unsigned int i = 0; i < planes.size(); ++i)
            {
                if(rayPlaneIntersect(testRay,planes[i],ip,depth))
                {
                    if(depth < fragDepth) // depth sorting
                    {
                        fragDepth = depth;
                        normal = planes[i].normal;
                        hdrColor = vec3();
                        for(unsigned int li = 0; li < lights.size(); ++li)
                        {
                            lightDir = (lights[li].pos-ip);
                            if(!occlusionTrace(ip, lights[li].pos, planes, d, spheres, i, -1, -1))
                            {
                                ldpMag = lightDir.qmag();
                                diff = std::max(dot(normal, lightDir*(1.0d/sqrt(ldpMag))), 0.0d) * (1.0d/((lights[li].pos-ip).qmag() + (ip-testRay.origin).qmag()));
                                spec = 0.0f;
                                if(diff > 0.0d)
                                {
                                    vec3 halfwayDir = (testRay.dir.norm()*-1+lightDir).norm();
                                    spec = pow(std::max(dot(normal, halfwayDir), 0.0), 8.0f)*0.02d;
                                }

                                hdrColor += planes[i].mat.color*planes[i].mat.albedo*(diff+spec)*lights[li].col;
                            }
                        }
                    }
                }
            }

            /*if(rayTriangleIntersect(testRay, t, ip, depth))
            {
                if(depth < fragDepth)
                {
                    fragDepth = depth;
                    normal = t.normal;
                    //lightDir = (l.pos-ip);
                    //ldpMag = lightDir.qmag();
                    //diff = std::max(dot(normal, lightDir*(1.0d/sqrt(ldpMag))), 0.0d) * (1.0d/ldpMag);
                    //hdrColor = t.color*diff*l.col;
                    hdrColor = t.color;
                    //if(occlusionTrace(ip, l.pos, p, d, spheres, -1, -1, -1))
                    //    hdrColor = vec3();
                    //else
                    //{
                    //    ldpMag = lightDir.qmag();
                    //    diff = std::max(dot(normal, lightDir*(1.0d/sqrt(ldpMag))), 0.0d) * (1.0d/ldpMag);
                    //    hdrColor = t.color*diff*l.col;
                    //}
                }
            }*/

            /*// render all disks in the scene
            if(rayDiskIntersect(testRay,d,ip,depth))
            {
                if(depth < fragDepth) // depth sorting
                {
                    fragDepth = depth;
                    normal = d.normal;
                    lightDir = (l.pos-ip);
                    hdrColor = d.mat->color * reflectionTrace(testRay.origin, ip, normal, p, d, spheres, l);
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
                        hdrColor = vec3();
                        for(unsigned int li = 0; li < lights.size(); ++li)
                        {
                            lightDir = (lights[li].pos-ip);
                            if(!occlusionTrace(ip, lights[li].pos, planes, d, spheres, -1, -1, i))
                            {
                                ldpMag = lightDir.qmag();
                                diff = std::max(dot(normal, lightDir*(1.0d/sqrt(ldpMag))), 0.0d) * (1.0d/((lights[li].pos-ip).qmag() + (ip-testRay.origin).qmag()));
                                spec = 0.0f;
                                if(diff > 0.0d)
                                {
                                    vec3 halfwayDir = (testRay.dir.norm()*-1+lightDir).norm();
                                    spec = pow(std::max(dot(normal, halfwayDir), 0.0), 8.0f)*0.02d;
                                }
                                hdrColor += spheres[i].mat.color*spheres[i].mat.albedo*(diff+spec)*lights[li].col;
                            }
                        }
                    }
                }
            }


            vec3_precision hdrc_p = vec3_precision();

            for(unsigned int n = 0; n < numSamples; ++n)
            {
                double pathLength = 0.0d;
                vec3 irr = irradianceTrace(testRay, vec3(1,1,1), planes, d, spheres, refractSphere, lights, pathLength, 8) * 2;
                hdrc_p = hdrc_p + irr;
            }

            hdrc_p = hdrc_p / numSamples;

            if(hdrc_p.getDoubleVector() > vec3())
                hdrColor += hdrc_p.getDoubleVector();


            // hdr filmic tone-mapping
            vec3 mapped = toneMap(hdrColor*2.0);
            vec3 whitescale = vec3(1.371034966875795,1.371034966875795,1.371034966875795); // <- 1/toneMap(11.5)
            mapped = mapped*whitescale;

            // gamma correction
            mapped.x = std::pow(mapped.x,1.0/2.2);
            mapped.y = std::pow(mapped.y,1.0/2.2);
            mapped.z = std::pow(mapped.z,1.0/2.2);

            buf[x*3+y*3*imgSize]   = (unsigned char)(mapped.z*255.0d);
            buf[x*3+y*3*imgSize+1] = (unsigned char)(mapped.y*255.0d);
            buf[x*3+y*3*imgSize+2] = (unsigned char)(mapped.x*255.0d);
        }
    }

    auto tp2 = std::chrono::steady_clock::now();

    exportArrayAsBMP("image.bmp", &buf[0], imgSize, imgSize);

    std::chrono::duration<double> timeDelta = tp2-tp1;
    std::cout << "Rendering done!\n" << timeDelta.count() << "sec/frame\n";

    delete buf;
}
