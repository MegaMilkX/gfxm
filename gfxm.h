#ifndef GRAPHICS_MATH_H
#define GRAPHICS_MATH_H

#include <math.h>

namespace gfxm{
    
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

}

#endif
