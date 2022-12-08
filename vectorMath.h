#ifndef VECTORMATH_H_INCLUDED
#define VECTORMATH_H_INCLUDED

#define M_PI 3.14159265358979

#include <math.h>
#include <random>
#include <limits>
#include <chrono>

std::mt19937 mt(time(0));

struct vec2
{
    double x,y;

    vec2() : x(0), y(0) {}
    vec2(double X, double Y) : x(X), y(Y) {}

    double mag()                           { return sqrt(x*x + y*y);                               }
    double qmag()                          { return this->x*this->x + this->y*this->y;             }
    vec2 norm()                           { double r = 1 / mag(); return vec2(x*r, y*r);           }
    vec2 perp()                           { return vec2(-y, x);                                   }

    vec2 operator + (const vec2& rhs)     { return vec2(this -> x + rhs.x, this -> y + rhs.y);    }
    vec2 operator - (const vec2& rhs)     { return vec2(this -> x - rhs.x, this -> y - rhs.y);    }
    vec2 operator * (const double& rhs)    { return vec2(this -> x * rhs,   this -> y * rhs);      }
    vec2 operator / (const double& rhs)    { return vec2(this -> x / rhs,   this -> y / rhs);      }

    vec2& operator += (const vec2& rhs)   { this -> x += rhs.x; this -> y += rhs.y; return *this; }
    vec2& operator -= (const vec2& rhs)   { this -> x -= rhs.x; this -> y -= rhs.y; return *this; }
    vec2& operator *= (const double& rhs)  { this -> x *= rhs;   this -> y *= rhs;   return *this; }
    vec2& operator /= (const double& rhs)  { this -> x /= rhs;   this -> y /= rhs;   return *this; }

};

struct vec3
{
    double x,y,z;

    vec3() : x(0), y(0), z(0) {}
    vec3(double X, double Y,double Z) : x(X), y(Y), z(Z) {}

    double mag()                           { return sqrt(x*x + y*y + z*z );                                            }
    double qmag()                          { return x*x + y*y + z*z;                                                   }
    vec3 norm()                           { double r = 1 / mag(); return vec3(x*r, y*r, z*r);                          }
    vec3 perp()                           { return vec3( -z, y, x);                                                   }

    vec3 operator + (const vec3& rhs)     { return vec3( this -> x + rhs.x, this -> y + rhs.y, this -> z + rhs.z );   }
    vec3 operator + (const double& rhs)   { return vec3( this -> x + rhs, this -> y + rhs, this -> z + rhs );         }
    vec3 operator - (const vec3& rhs)     { return vec3( this -> x - rhs.x, this -> y - rhs.y, this -> z - rhs.z );   }
    vec3 operator - (const double& rhs)   { return vec3( this -> x - rhs, this -> y - rhs, this -> z - rhs );         }
    vec3 operator * (const double& rhs)    { return vec3( this -> x * rhs, this -> y * rhs, this -> z * rhs );         }
    vec3 operator * (const vec3& rhs)     { return vec3( this -> x * rhs.x, this -> y * rhs.y, this -> z * rhs.z );   }
    vec3 operator / (const vec3& rhs)     { return vec3( this -> x / rhs.x, this -> y / rhs.y, this -> z / rhs.z );   }
    vec3 operator / (const double& rhs)    { return vec3( this -> x / rhs, this -> y / rhs, this -> z / rhs );         }

    vec3& operator += (const vec3& rhs)   { this -> x += rhs.x; this -> y += rhs.y; this -> z += rhs.z; return *this; }
    vec3& operator -= (const vec3& rhs)   { this -> x -= rhs.x; this -> y -= rhs.y; this -> z -= rhs.z; return *this; }
    vec3& operator *= (const double& rhs)  { this -> x *= rhs;   this -> y *= rhs;   this -> z *= rhs;   return *this; }
    vec3& operator /= (const double& rhs)  { this -> x /= rhs;   this -> y /= rhs;   this -> z /= rhs;   return *this; }

    bool operator < (vec3 rhs)
    {
        return (this->x < rhs.x) || (this->y < rhs.y) || (this->z < rhs.z);
    }
    bool operator > (vec3 rhs)
    {
        return (this->x > rhs.x) && (this->y > rhs.y) && (this->z > rhs.z);
    }

    void print()
    {
        std::cout << this->x << ',' << this->y << ',' << this->z << '\n';
    }
};

double sign(double a)
{
    if(a > 0.0d)
        return 1.0d;
    else return -1.0d;
}

struct vec3_precision
{
    double xf=0.0d,yf=0.0d,zf=0.0d;
    int xi,yi,zi;

    vec3_precision() : xi(0), yi(0), zi(0) {}
    vec3_precision(double x, double y, double z) : xi(floor(x)), yi(floor(y)), zi(floor(z))
    {
        this->xf = x-this->xi;
        this->yf = y-this->yi;
        this->zf = z-this->zi;
    }

    void print()
    {
        std::cout << this->xi << ',' << this->yi << ',' << this->zi << ':' << this->xf << ',' << this->yf << ',' << this->zf << '\n';
    }

    vec3_precision operator + (const vec3 rhs)
    {
        vec3_precision tmp;

        tmp.xi = this->xi + floor(rhs.x);
        tmp.yi = this->yi + floor(rhs.y);
        tmp.zi = this->zi + floor(rhs.z);

        tmp.xf = this->xf + (rhs.x-floor(rhs.x));
        tmp.yf = this->yf + (rhs.y-floor(rhs.y));
        tmp.zf = this->zf + (rhs.z-floor(rhs.z));

        if(std::fabs(tmp.xf) > 1.0f)
        {
            while(tmp.xf > 1.0d)
            {
                tmp.xi++;
                tmp.xf -= 1.0d;
            }
        }

        if(std::fabs(tmp.xf) > 1.0f)
        {
            while(tmp.xf < 1.0d)
            {
                tmp.xi--;
                tmp.xf += 1.0d;
            }
        }

        if(std::fabs(tmp.yf) > 1.0f)
        {
            while(tmp.yf > 1.0d)
            {
                tmp.yi++;
                tmp.yf -= 1.0d;
            }
        }

        if(std::fabs(tmp.yf) > 1.0f)
        {
            while(tmp.yf < 1.0d)
            {
                tmp.yi--;
                tmp.yf += 1.0d;
            }
        }

        if(std::fabs(tmp.zf) > 1.0f)
        {
            while(tmp.zf > 1.0d)
            {
                tmp.zi++;
                tmp.zf -= 1.0d;
            }
        }

        if(std::fabs(tmp.zf) > 1.0f)
        {
            while(tmp.zf < 1.0d)
            {
                tmp.zi--;
                tmp.zf += 1.0d;
            }
        }

        return tmp;
    }

    vec3_precision operator / (unsigned int rhs)
    {
        vec3_precision tmp;

        tmp.xi = this->xi / (int)rhs;
        unsigned int remainder = this->xi % rhs;
        tmp.xf = (this->xf / rhs);
        tmp.xf += sign(this->xi) * ((double)remainder/(double)rhs);

        tmp.yi = this->yi / (int)rhs;
        remainder = this->yi % rhs;
        tmp.yf = (this->yf / rhs);
        tmp.yf += sign(this->yi) * ((double)remainder/(double)rhs);

        tmp.zi = this->zi / (int)rhs;
        remainder = this->zi % rhs;
        tmp.zf = (this->zf / rhs);
        tmp.zf += sign(this->zi) * ((double)remainder/(double)rhs);

        return tmp;
    }

    vec3 getDoubleVector()
    {
        vec3 tmp;
        tmp.x = (double)this->xi + this->xf;
        tmp.y = (double)this->yi + this->yf;
        tmp.z = (double)this->zi + this->zf;
        return tmp;
    }
};

vec3 cross(vec3 x, vec3 y)
{
    return vec3(x.y*y.z - y.y*x.z,
                x.z*y.x - y.z*x.x,
                x.x*y.y - y.x*x.y);
}
double dot(vec3 a, vec3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

double randDouble(double min, double max)
{
    double a = (double)mt()/(double)std::numeric_limits<uint32_t>::max();
    return a*(max-min)+min;
}

vec3 sampleSphereUniform()
{
    vec3 point;
    do point = vec3(randDouble(-1,1),randDouble(-1,1),randDouble(-1,1));
    while(point.qmag() > 1);

    return point.norm();
}

vec3 sampleHemisphereUniform(vec3 normal)
{
    vec3 point;
    do point = sampleSphereUniform();
    while(dot(point,normal) < 0);
    return point;
}

int hexToBin(char h)
{
    if(h < 65)
        return h-48;
    else if(h < 97)
        return h-55;
    else return h-87;
}

vec3 colorFromHex(std::string hexColor)
{
    vec3 tmp;
    tmp.x = float((hexToBin(hexColor[0])<<4) | hexToBin(hexColor[1])) / 255.0f;
    tmp.y = float((hexToBin(hexColor[2])<<4) | hexToBin(hexColor[3])) / 255.0f;
    tmp.z = float((hexToBin(hexColor[4])<<4) | hexToBin(hexColor[5])) / 255.0f;
    return tmp;
}

uint_fast8_t i3(uint_fast8_t index, uint_fast8_t r)
{
    return index>=r?index+1:index;
}
struct mat4
{
    float mat[4][4] = {};

    mat4() {}

    void loadIdentity()
    {
        for(uint16_t i = 0; i < 4; ++i)
        {
            for(uint16_t j = 0; j < 4; j++)
            {
                if(j == i)
                    this->mat[i][j] = 1.0f;
                else this->mat[i][j] = 0.0f;
            }
        }
    }

    void print()
    {
        std::cout << this->mat[0][0] << ',' << this->mat[1][0] << ',' << this->mat[2][0] << ',' << this->mat[3][0] << ",\n";
        std::cout << this->mat[0][1] << ',' << this->mat[1][1] << ',' << this->mat[2][1] << ',' << this->mat[3][1] << ",\n";
        std::cout << this->mat[0][2] << ',' << this->mat[1][2] << ',' << this->mat[2][2] << ',' << this->mat[3][2] << ",\n";
        std::cout << this->mat[0][3] << ',' << this->mat[1][3] << ',' << this->mat[2][3] << ',' << this->mat[3][3] << '\n';
    }

    const float* getGLMat()
    {
        return (const float*)&this->mat;
    }

    mat4 operator + (const mat4& rhs)
    {
        mat4 tmp = mat4();
        for(uint_fast8_t i = 0; i < 4; ++i)
        {
            for(uint_fast8_t j = 0; j < 4; j++)
                tmp.mat[i][j] = this->mat[i][j] + rhs.mat[i][j];
        }
        return tmp;
    }

    mat4 operator - (const mat4& rhs)
    {
        mat4 tmp = mat4();
        for(uint_fast8_t i = 0; i < 4; ++i)
        {
            for(uint_fast8_t j = 0; j < 4; j++)
                tmp.mat[i][j] = this->mat[i][j] - rhs.mat[i][j];
        }
        return tmp;
    }

    mat4 operator * (const mat4& rhs)
    {
        mat4 tmp = mat4();
        for(uint_fast8_t i = 0; i < 4; ++i)
        {
            for(uint_fast8_t j = 0; j < 4; j++)
            {
                for(uint_fast8_t k = 0; k < 4; k++)
                    tmp.mat[i][j] += this->mat[i][k] * rhs.mat[k][j];
            }
        }
        return tmp;
    }
    vec3 operator * (const vec3& rhs)
    {
        vec3 tmp;
        tmp.x = rhs.x*this->mat[0][0] + rhs.y*this->mat[1][0] + rhs.z*this->mat[2][0] + this->mat[3][0];
        tmp.y = rhs.x*this->mat[0][1] + rhs.y*this->mat[1][1] + rhs.z*this->mat[2][1] + this->mat[3][1];
        tmp.z = rhs.x*this->mat[0][2] + rhs.y*this->mat[1][2] + rhs.z*this->mat[2][2] + this->mat[3][2];
        return tmp;
    }

    float det()
    {
        mat4 A = *(this);
        float c, r = 1.0f;
        for(uint_fast8_t i = 0; i < 4; ++i)
        {
            for(uint_fast8_t k = i+1; k < 4; ++k)
            {
                c = A.mat[i][k] / A.mat[i][i];
                for(uint_fast8_t j = i; j < 4; ++j)
                    A.mat[j][k] = A.mat[j][k] - c*A.mat[j][i];
            }
        }
        for(uint_fast8_t i = 0; i < 4; ++i)
            r *= A.mat[i][i];
        return r;
    }

    float det3x3(uint_fast8_t ni, uint_fast8_t nj)
    {
        return this->mat[i3(0,ni)][i3(0,nj)]*this->mat[i3(1,ni)][i3(1,nj)]*this->mat[i3(2,ni)][i3(2,nj)]+   // a*e*i +
               this->mat[i3(0,ni)][i3(1,nj)]*this->mat[i3(1,ni)][i3(2,nj)]*this->mat[i3(2,ni)][i3(0,nj)]+   // b*f*g +
               this->mat[i3(0,ni)][i3(2,nj)]*this->mat[i3(1,ni)][i3(0,nj)]*this->mat[i3(2,ni)][i3(1,nj)]-   // c*d*h -
               this->mat[i3(0,ni)][i3(2,nj)]*this->mat[i3(1,ni)][i3(1,nj)]*this->mat[i3(2,ni)][i3(0,nj)]-   // c*e*g -
               this->mat[i3(0,ni)][i3(1,nj)]*this->mat[i3(1,ni)][i3(0,nj)]*this->mat[i3(2,ni)][i3(2,nj)]-   // b*d*i -
               this->mat[i3(0,ni)][i3(0,nj)]*this->mat[i3(1,ni)][i3(2,nj)]*this->mat[i3(2,ni)][i3(1,nj)];   // a*f*h
    }

    mat4 coFactor()
    {
        mat4 tmp = mat4();
        for(uint_fast8_t i = 0; i < 4; ++i)
        {
            for(uint_fast8_t j = 0; j < 4; j++)
                tmp.mat[i][j] = ((i+j)%2==0?1:-1) * this->det3x3(i,j);
        }
        return tmp;
    }

    mat4 transpose()
    {
        mat4 tmp;
        for(uint_fast8_t i = 0; i < 4; ++i)
        {
            for(uint_fast8_t j = 0; j < 4; j++)
                tmp.mat[i][j] = this->mat[j][i];
        }
        return tmp;
    }

    mat4 inverse()
    {
        mat4 tmp = *(this);
        float d = 1.0f/this->det();

        tmp = tmp.coFactor();
        tmp = tmp.transpose();

        for(uint_fast8_t i = 0; i < 4; ++i)
        {
            for(uint_fast8_t j = 0; j < 4; j++)
                tmp.mat[i][j] *= d;
        }
        return tmp;
    }

};

mat4 getIdentityMatrix()
{
    mat4 tmp = mat4();
    tmp.loadIdentity();
    return tmp;
}

mat4 getScalingMatrix(float vec[4])
{
    mat4 tmp = mat4();
    for(uint16_t i = 0; i < 4; ++i)
        tmp.mat[i][i] = vec[i];
    return tmp;
}
mat4 getScalingMatrix(vec3 vec)
{
    mat4 tmp = mat4();
    tmp.mat[0][0] = vec.x;
    tmp.mat[1][1] = vec.y;
    tmp.mat[2][2] = vec.z;
    tmp.mat[3][3] = 1.0f;
    return tmp;
}
mat4 getScalingMatrix(float s)
{
    mat4 tmp = mat4();
    tmp.mat[0][0] = s;
    tmp.mat[1][1] = s;
    tmp.mat[2][2] = s;
    tmp.mat[3][3] = 1.0f;
    return tmp;
}

mat4 getTranslationMatrix(float vec[4])
{
    mat4 tmp = mat4();
    tmp.loadIdentity();
    tmp.mat[3][0] = vec[0];
    tmp.mat[3][1] = vec[1];
    tmp.mat[3][2] = vec[2];
    return tmp;
}
mat4 getTranslationMatrix(vec3 vec)
{
    mat4 tmp = mat4();
    tmp.loadIdentity();
    tmp.mat[3][0] = vec.x;
    tmp.mat[3][1] = vec.y;
    tmp.mat[3][2] = vec.z;
    return tmp;
}

mat4 getRotationMatrix(float theta, vec3 r)
{
    float st = sin(theta);
    float ct = cos(theta);
    mat4 tmp = mat4();

    tmp.mat[0][0] = ct+(r.x*r.x)*(1-ct);
    tmp.mat[0][1] = r.x*r.y*(1-ct)-r.z*st;
    tmp.mat[0][2] = r.x*r.z*(1-ct)+r.y*st;
    tmp.mat[0][3] = 0.0f;

    tmp.mat[1][0] = r.y*r.x*(1-ct)+r.z*st;
    tmp.mat[1][1] = ct+(r.y*r.y)*(1-ct);
    tmp.mat[1][2] = r.y*r.z*(1-ct)-r.x*st;
    tmp.mat[1][3] = 0.0f;

    tmp.mat[2][0] = r.z*r.x*(1-ct)-r.y*st;
    tmp.mat[2][1] = r.z*r.y*(1-ct)+r.x*st;
    tmp.mat[2][2] = ct+(r.z*r.z)*(1-ct);
    tmp.mat[2][3] = 0.0f;

    tmp.mat[3][0] = 0.0f;
    tmp.mat[3][1] = 0.0f;
    tmp.mat[3][2] = 0.0f;
    tmp.mat[3][3] = 1.0f;

    return tmp;

}

mat4 getPerspectiveProjectionMatrix(float FOV, float aspect, float zNear, float zFar)
{
    mat4 tmp = mat4();
    float tanHalfFovy = tan(FOV/2);

    tmp.mat[0][0] = 1 / (aspect * tanHalfFovy);
    tmp.mat[1][1] = 1 / tanHalfFovy;
    tmp.mat[2][2] = -(zFar + zNear) / (zFar - zNear);
    tmp.mat[2][3] = -1.0f;
    tmp.mat[3][2] = -(2 * zFar * zNear) / (zFar - zNear);

    return tmp;
}

mat4 getOrthographicProjectionMatrix(float left, float right, float bottom, float top, float zNear, float zFar)
{
    mat4 tmp = mat4();
    tmp.loadIdentity();

    tmp.mat[0][0] =  2.0f/(right-left);
    tmp.mat[1][1] =  2.0f/(top-bottom);
    tmp.mat[2][2] = -2.0f/(zFar-zNear);
    tmp.mat[3][0] = -(right+left)/(right-left);
    tmp.mat[3][1] = -(top+bottom)/(top-bottom);
    tmp.mat[3][2] = -(zFar+zNear)/(zFar-zNear);

    return tmp;
}

mat4 lookAt(vec3 eye, vec3 center, vec3 up)
{
    mat4 c = mat4();

    vec3 d = eye-center; d = d.norm();
    vec3 r = cross(up,d); r = r.norm();
    vec3 u = cross(d,r);

    c.mat[0][0] = r.x;
    c.mat[1][0] = r.y;
    c.mat[2][0] = r.z;
    c.mat[3][0] = 0.0f;

    c.mat[0][1] = u.x;
    c.mat[1][1] = u.y;
    c.mat[2][1] = u.z;
    c.mat[3][1] = 0.0f;

    c.mat[0][2] = d.x;
    c.mat[1][2] = d.y;
    c.mat[2][2] = d.z;
    c.mat[3][2] = 0.0f;

    c.mat[0][3] = 0.0f;
    c.mat[1][3] = 0.0f;
    c.mat[2][3] = 0.0f;
    c.mat[3][3] = 1.0f;

    mat4 e = mat4();
    e.loadIdentity();

    e.mat[3][0] = -eye.x;
    e.mat[3][1] = -eye.y;
    e.mat[3][2] = -eye.z;

	return e*c;

}

#endif // VECTORMATH_H_INCLUDED
