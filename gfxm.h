#ifndef GRAPHICS_MATH_H
#define GRAPHICS_MATH_H

#include <math.h>

namespace gfxm{

const float pi = 3.14159265359f;
const double d_pi = 3.14159265359;
    
// ====== Types ==========

template<typename T>
struct tvec2
{
    T x;
    T y;
    
    tvec2() : x(0), y(0) {}
    tvec2(T x, T y) : x(x), y(y) {}
};

template<typename T>
struct tvec3
{
    T x;
    T y;
    T z;
    
    tvec3() : x(0), y(0), z(0) {}
    tvec3(T x, T y, T z) : x(x), y(y), z(z) {}
};

template<typename T>
struct tvec4
{
    T x;
    T y;
    T z;
    T w;
    
    tvec4() : x(0), y(0), z(0), w(0) {}
    tvec4(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) {}
};

template<typename T>
struct tquat
{
    T x;
    T y;
    T z;
    T w;
    
    tquat() : x(0), y(0), z(0), w(1) {}
    tquat(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) {}
};

template<typename T>
struct tmat3
{
    tmat3(){}
    explicit tmat3(T f)
    {
        col[0].x = f;
        col[1].y = f;
        col[2].z = f;
    }
    
    tvec3<T> operator[](const int &i) const {
        return col[i];
    }
    tvec3<T>& operator[](const int &i){
        return col[i];
    }
private:
    tvec3<T> col[3];
};

template<typename T>
struct tmat4
{
    tmat4(){}
    explicit tmat4(T f)
    {
        col[0].x = f;
        col[1].y = f;
        col[2].z = f;
        col[3].w = f;
    }
    
    tvec4<T> operator[](const int &i) const {
        return col[i];
    }
    tvec4<T>& operator[](const int &i){
        return col[i];
    }
private:
    tvec4<T> col[4];
};

typedef tvec2<float> vec2;
typedef tvec2<int> ivec2;
typedef tvec2<double> dvec2;

typedef tvec3<float> vec3;
typedef tvec3<int> ivec3;
typedef tvec3<double> dvec3;

typedef tvec4<float> vec4;
typedef tvec4<int> ivec4;
typedef tvec4<double> dvec4;

typedef tquat<float> quat;
typedef tquat<double> dquat;

typedef tmat3<float> mat3;
typedef tmat3<double> dmat3;

typedef tmat4<float> mat4;
typedef tmat4<double> dmat4;

// ====== Functions ======

template<typename T>
inline tvec2<T> operator+(const tvec2<T>& a, const tvec2<T>& b){
    return tvec2<T>(a.x + b.x, a.y + b.y);
}
template<typename T>
inline tvec3<T> operator+(const tvec3<T>& a, const tvec3<T>& b){
    return tvec3<T>(a.x + b.x, a.y + b.y, a.z + b.z);
}
template<typename T>
inline tvec4<T> operator+(const tvec4<T>& a, const tvec4<T>& b){
    return tvec4<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}

template<typename T>
inline tvec2<T> operator+=(tvec2<T>& a, const tvec2<T>& b){
    return a = a + b;
}
template<typename T>
inline tvec3<T> operator+=(tvec3<T>& a, const tvec3<T>& b){
    return a = a + b;
}
template<typename T>
inline tvec4<T> operator+=(tvec4<T>& a, const tvec4<T>& b){
    return a = a + b;
}

template<typename T>
inline tvec2<T> operator-(const tvec2<T>& a, const tvec2<T>& b){
    return tvec2<T>(a.x - b.x, a.y - b.y);
}
template<typename T>
inline tvec3<T> operator-(const tvec3<T>& a, const tvec3<T>& b){
    return tvec3<T>(a.x - b.x, a.y - b.y, a.z - b.z);
}
template<typename T>
inline tvec4<T> operator-(const tvec4<T>& a, const tvec4<T>& b){
    return tvec4<T>(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}

template<typename T>
inline tvec2<T> operator-=(tvec2<T>& a, const tvec2<T>& b){
    return a = a - b;
}
template<typename T>
inline tvec3<T> operator-=(tvec3<T>& a, const tvec3<T>& b){
    return a = a - b;
}
template<typename T>
inline tvec4<T> operator-=(tvec4<T>& a, const tvec4<T>& b){
    return a = a - b;
}

template<typename T>
inline tvec2<T> operator-(const tvec2<T>& v){
    return tvec2<T>(-v.x, -v.y);
}
template<typename T>
inline tvec3<T> operator-(const tvec3<T>& v){
    return tvec3<T>(-v.x, -v.y, -v.z);
}
template<typename T>
inline tvec4<T> operator-(const tvec4<T>& v){
    return tvec4<T>(-v.x, -v.y, -v.z, -v.w);
}

template<typename T, typename M>
inline tvec2<T> operator*(const tvec2<T>& a, const M& f){
    return tvec2<T>(a.x * f, a.y * f);
}
template<typename T, typename M>
inline tvec2<T> operator*(const M& f, const tvec2<T>& a){
    return tvec2<T>(a.x * f, a.y * f);
}
template<typename T, typename M>
inline tvec3<T> operator*(const tvec3<T>& a, const M& f){
    return tvec3<T>(a.x * f, a.y * f, a.z * f);
}
template<typename T, typename M>
inline tvec3<T> operator*(const M& f, const tvec3<T>& a){
    return tvec3<T>(a.x * f, a.y * f, a.z * f);
}
template<typename T, typename M>
inline tvec4<T> operator*(const tvec4<T>& a, const M& f){
    return tvec4<T>(a.x * f, a.y * f, a.z * f, a.w * f);
}
template<typename T, typename M>
inline tvec4<T> operator*(const M& f, const tvec4<T>& a){
    return tvec4<T>(a.x * f, a.y * f, a.z * f, a.w * f);
}

template<typename T, typename M>
inline tvec2<T> operator*=(tvec2<T>& a, const M& f){
    return a = a * f;
}
template<typename T, typename M>
inline tvec3<T> operator*=(tvec3<T>& a, const M& f){
    return a = a * f;
}
template<typename T, typename M>
inline tvec4<T> operator*=(tvec4<T>& a, const M& f){
    return a = a * f;
}

template<typename T, typename M>
inline tvec2<T> operator/(const tvec2<T>& a, const M& f){
    return tvec2<T>(a.x / f, a.y / f);
}
template<typename T, typename M>
inline tvec3<T> operator/(const tvec3<T>& a, const M& f){
    return tvec3<T>(a.x / f, a.y / f, a.z / f);
}
template<typename T, typename M>
inline tvec4<T> operator/(const tvec4<T>& a, const M& f){
    return tvec4<T>(a.x / f, a.y / f, a.z / f, a.w / f);
}

template<typename T, typename M>
inline tvec2<T> operator/=(tvec2<T>& a, const M& f){
    return a / f;
}
template<typename T, typename M>
inline tvec3<T> operator/=(tvec3<T>& a, const M& f){
    return a / f;
}
template<typename T, typename M>
inline tvec4<T> operator/=(tvec4<T>& a, const M& f){
    return a / f;
}

inline float qrsqrt(const float &n)
{
    long i;
    float x2, y;
    const float threehalves = 1.5f;
    x2 = n * 0.5f;
    y = n;
    i = *(long*)&y;
    i = 0x5f3759df - (i >> 1);
    y = *(float*)&i;
    y = y * (threehalves - (x2 * y * y));
    return y;
}

inline float sqrt(const float &n)
{
    return n * qrsqrt(n);
}

template<typename T>
inline T length(const tvec2<T>& v) { return sqrt(v.x*v.x + v.y*v.y); }
template<typename T>
inline T length(const tvec3<T>& v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z); }
template<typename T>
inline T length(const tvec4<T>& v) { return sqrt(v.x*v.x + v.y*v.y + v.z*v.z + v.w*v.w); }

template<typename T>
inline tvec2<T> normalize(const tvec2<T>& v) 
{
    T l = length(v);
    if(l == T())
        return v;
    return v / l;
}

template<typename T>
inline tvec3<T> normalize(const tvec3<T>& v) 
{
    T l = length(v);
    if(l == T())
        return v;
    return v / l;
}

template<typename T>
inline tvec4<T> normalize(const tvec4<T>& v) 
{
    T l = length(v);
    if(l == T())
        return v;
    return v / l;
}

template<typename T>
inline T dot(const tvec2<T>& a, const tvec2<T>& b)
{
    return a.x * b.x + a.y * b.y;
}
template<typename T>
inline T dot(const tvec3<T>& a, const tvec3<T>& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template<typename T>
inline T dot(const tvec4<T>& a, const tvec4<T>& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.z;
}

template<typename T>
inline tvec3<T> cross(const tvec3<T>& a, const tvec3<T>& b)
{
    return tvec3<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

template<typename T>
inline tmat4<T> inverse(const tmat4<T> &mat)
{
    const T* m;
    m = &mat;
    T det;
    int i;
    tmat4<T> inverse(1.0f);
    T* inv = &inverse;
    
    inv[0] = m[5]  * m[10] * m[15] - 
             m[5]  * m[11] * m[14] - 
             m[9]  * m[6]  * m[15] + 
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] - 
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] + 
              m[4]  * m[11] * m[14] + 
              m[8]  * m[6]  * m[15] - 
              m[8]  * m[7]  * m[14] - 
              m[12] * m[6]  * m[11] + 
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] - 
             m[4]  * m[11] * m[13] - 
             m[8]  * m[5] * m[15] + 
             m[8]  * m[7] * m[13] + 
             m[12] * m[5] * m[11] - 
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] + 
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] - 
               m[8]  * m[6] * m[13] - 
               m[12] * m[5] * m[10] + 
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] + 
              m[1]  * m[11] * m[14] + 
              m[9]  * m[2] * m[15] - 
              m[9]  * m[3] * m[14] - 
              m[13] * m[2] * m[11] + 
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] - 
             m[0]  * m[11] * m[14] - 
             m[8]  * m[2] * m[15] + 
             m[8]  * m[3] * m[14] + 
             m[12] * m[2] * m[11] - 
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] + 
              m[0]  * m[11] * m[13] + 
              m[8]  * m[1] * m[15] - 
              m[8]  * m[3] * m[13] - 
              m[12] * m[1] * m[11] + 
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] - 
              m[0]  * m[10] * m[13] - 
              m[8]  * m[1] * m[14] + 
              m[8]  * m[2] * m[13] + 
              m[12] * m[1] * m[10] - 
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] - 
             m[1]  * m[7] * m[14] - 
             m[5]  * m[2] * m[15] + 
             m[5]  * m[3] * m[14] + 
             m[13] * m[2] * m[7] - 
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] + 
              m[0]  * m[7] * m[14] + 
              m[4]  * m[2] * m[15] - 
              m[4]  * m[3] * m[14] - 
              m[12] * m[2] * m[7] + 
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] - 
              m[0]  * m[7] * m[13] - 
              m[4]  * m[1] * m[15] + 
              m[4]  * m[3] * m[13] + 
              m[12] * m[1] * m[7] - 
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] + 
               m[0]  * m[6] * m[13] + 
               m[4]  * m[1] * m[14] - 
               m[4]  * m[2] * m[13] - 
               m[12] * m[1] * m[6] + 
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] + 
              m[1] * m[7] * m[10] + 
              m[5] * m[2] * m[11] - 
              m[5] * m[3] * m[10] - 
              m[9] * m[2] * m[7] + 
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] - 
             m[0] * m[7] * m[10] - 
             m[4] * m[2] * m[11] + 
             m[4] * m[3] * m[10] + 
             m[8] * m[2] * m[7] - 
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] + 
               m[0] * m[7] * m[9] + 
               m[4] * m[1] * m[11] - 
               m[4] * m[3] * m[9] - 
               m[8] * m[1] * m[7] + 
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] - 
              m[0] * m[6] * m[9] - 
              m[4] * m[1] * m[10] + 
              m[4] * m[2] * m[9] + 
              m[8] * m[1] * m[6] - 
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    assert(det != 0);

    det = 1.0f / det;

    for (i = 0; i < 16; i++)
        inv[i] = inv[i] * det;

    return inverse;
}

}

#endif
