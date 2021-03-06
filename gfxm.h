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

    T operator[](const int &i) const {
        return *((&x) + i);
    }
    T& operator[](const int &i) {
        return *((&x) + i);
    }
};

template<typename T>
struct tvec3
{
    T x;
    T y;
    T z;
    
    tvec3() : x(0), y(0), z(0) {}
    tvec3(T x, T y, T z) : x(x), y(y), z(z) {}

    T operator[](const int &i) const {
        return *((&x) + i);
    }
    T& operator[](const int &i) {
        return *((&x) + i);
    }
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

    operator tvec3<T>() const { return tvec3<T>(x, y, z); }

    T operator[](const int &i) const {
        return *((&x) + i);
    }
    T& operator[](const int &i) {
        return *((&x) + i);
    }
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

    void operator=(const tmat3<T>& m)
    {
        (*this) = mat4(1.0f);
        col[0][0] = m[0][0]; col[0][1] = m[0][1]; col[0][2] = m[0][2];
        col[1][0] = m[1][0]; col[1][1] = m[1][1]; col[1][2] = m[1][2];
        col[2][0] = m[2][0]; col[2][1] = m[2][1]; col[2][2] = m[2][2];
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
T _min(T a, T b)
{
    if (a < b)
        return a;
    else
        return b;
}

template<typename T>
T _max(T a, T b)
{
    if (a > b)
        return a;
    else
        return b;
}

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
inline tquat<T> operator+(const tquat<T>& a, const tquat<T>& b) {
    return tquat<T>(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
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
inline tquat<T> operator+=(tquat<T>& a, const tquat<T>& b) {
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
inline tquat<T> operator-(const tquat<T>& a, const tquat<T>& b){
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
inline tquat<T> operator-=(tquat<T>& a, const tquat<T>& b) {
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
template<typename T>
inline tquat<T> operator-(const tquat<T>& v) {
    return tquat<T>(-v.x, -v.y, -v.z, -v.w);
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
inline tquat<T> operator/(const tquat<T>& a, const M& f) {
    return tquat<T>(a.x / f, a.y / f, a.z / f, a.w / f);
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
template<typename T, typename M>
inline tquat<T> operator/=(tquat<T>& a, const M& f) {
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
inline T length(const tquat<T>& q) { return sqrt(q.x*q.x + q.y*q.y + q.z*q.z + q.w*q.w); }

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
inline tquat<T> normalize(const tquat<T>& a) {
    if (length(a) == 0.0f)
        return a;
    return a / length(a);
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
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template<typename T>
inline T dot(const tquat<T>& a, const tquat<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

template<typename T>
inline tvec3<T> cross(const tvec3<T>& a, const tvec3<T>& b)
{
    return tvec3<T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

template<typename T>
inline tmat4<T> operator+(const tmat4<T>& m0, const tmat4<T>& m1) {
    tmat4<T> m;
    for (int i = 0; i < 4; ++i)
        m[i] = m0[i] + m1[i];
    return m;
}
template<typename T>
inline tmat3<T> operator*(const tmat3<T>& m0, const tmat3<T>& m1) {
    tmat3<T> m;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                m[i][j] += m0[k][j] * m1[i][k];
    return m;
}
template<typename T>
inline tmat4<T> operator*(const tmat4<T>& m0, const tmat4<T>& m1) {
    tmat4<T> m;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                m[i][j] += m0[k][j] * m1[i][k];
    return m;
}

template<typename T>
inline tvec4<T> operator*(const tmat4<T>& m, const tvec4<T>& v)
{
    tvec4<T> r;
    for (int i = 0; i < 4; ++i)
        for (int k = 0; k < 4; ++k)
            r[i] += m[k][i] * v[k];
    return r;
}

//Taking vec3 as vec4 and assuming v.w is zero
//so it transforms as a direction vector
template<typename T>
inline tvec3<T> operator*(const tmat4<T> &m, const tvec3<T> &v)
{
    tvec3<T> r;
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            r[i] += m[k][i] * v[k];
    return r;
}
template<typename T>
inline tvec3<T> operator*(const tmat3<T> &m, const tvec3<T> &v)
{
    tvec3<T> r;
    for (int i = 0; i < 3; ++i)
        for (int k = 0; k < 3; ++k)
            r[i] += m[k][i] * v[k];
    return r;
}
template<typename T>
inline tmat3<T> transpose(const tmat3<T>& m)
{
    tmat3<T> r(1.0f);
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            r[i][j] = m[j][i];
    return r;
}
template<typename T>
inline tmat4<T> transpose(const tmat4<T>& m)
{
    tmat4<T> r(1.0f);
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
            r[i][j] = m[j][i];
    return r;
}
template<typename T>
inline tmat4<T> scale(const tmat4<T>& m, const tvec3<T>& v)
{
    tmat4<T> r = m;
    r[0] *= v[0];
    r[1] *= v[1];
    r[2] *= v[2];
    return r;
}
template<typename T>
inline tmat4<T> translate(const tmat4<T> &m, const tvec3<T> &v)
{
    tmat4<T> r = m;
    r[3] = m[0] * v[0] + m[1] * v[1] + m[2] * v[2] + m[3];
    return r;
}

template<typename T>
inline tmat4<T> inverse(const tmat4<T> &mat)
{
    const T* m;
    m = (T*)&mat;
    T det;
    int i;
    tmat4<T> inverse(1.0f);
    T* inv = (T*)&inverse;
    
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

    //assert(det != 0);

    det = 1.0f / det;

    for (i = 0; i < 16; i++)
        inv[i] = inv[i] * det;

    return inverse;
}

template<typename T>
inline tquat<T> operator*(const tquat<T>& q0, const tquat<T>& q1)
{
    return normalize(tquat<T>((q0.w * q1.x + q1.w * q0.x) + (q0.y * q1.z - q1.y * q0.z),
        (q0.w * q1.y + q1.w * q0.y) + (q1.x * q0.z - q0.x * q1.z), //Inverted, since y axis rotation is inverted
        (q0.w * q1.z + q1.w * q0.z) + (q0.x * q1.y - q1.x * q0.y),
        (q1.w * q0.w) - (q1.x * q0.x) - (q1.y * q0.y) - (q1.z * q0.z)));
}
template<typename T>
inline tquat<T> operator*=(tquat<T>& q0, const tquat<T>& q1)
{
    return q0 = q0 * q1;
}

template<typename T, typename M>
inline tquat<T> operator*(const tquat<T>& a, const M& f) {
    return tquat<T>(a.x * f, a.y * f, a.z * f, a.w * f);
}
template<typename T, typename M>
inline tquat<T> operator*=(tquat<T>& a, const M& f) {
    return a = a * f;
}
template<typename T>
inline tquat<T> angle_axis(float a, const tvec3<T>& axis)
{
    float s = sinf(a * 0.5f);
    return normalize(tquat<T>(axis.x * s, axis.y * s, axis.z * s, cosf(a*0.5f)));
}
template<typename T>
inline tquat<T> inverse(const tquat<T>& q)
{
    tquat<T> i;
    float d = dot(q, q);
    i.x = -q.x / d;
    i.y = -q.y / d;
    i.z = -q.z / d;
    i.w = q.w / d;
    return i;
}
template<typename T>
inline tmat3<T> to_mat3(const tquat<T>& q)
{
    tmat3<T> m;
    m[0].x = 1 - 2 * q.y * q.y - 2 * q.z * q.z;
    m[0].y = q.z * 2 * q.w + 2 * q.x * q.y;
    m[0].z = -q.y * 2 * q.w + 2 * q.x * q.z;
    m[1].x = -q.z * 2 * q.w + 2 * q.x * q.y;
    m[1].y = 1 - 2 * q.x * q.x - 2 * q.z * q.z;
    m[1].z = q.x * 2 * q.w + 2 * q.y * q.z;
    m[2].x = q.y * 2 * q.w + 2 * q.x * q.z;
    m[2].y = -q.x * 2 * q.w + 2 * q.y * q.z;
    m[2].z = 1 - 2 * q.x * q.x - 2 * q.y * q.y;
    return m;
}
template<typename T>
inline tmat3<T> to_mat3(const tmat4<T>& m4)
{
    tmat3<T> m3;
    for (unsigned i = 0; i < 3; ++i)
        m3[i] = m4[i];
    return m3;
}
template<typename T>
inline tmat3<T> to_orient_mat3(const tmat4<T>& m)
{
    tmat3<T> mt = to_mat3(m);

    for (unsigned i = 0; i < 3; ++i)
    {
        tvec3<T> v3 = mt[i];
        v3 = normalize(v3);
        tvec4<T> v4 = tvec4<T>(v3.x, v3.y, v3.z, 0.0f);
        mt[i] = v4;
    }

    return mt;
}
template<typename T>
inline tmat4<T> to_mat4(const tquat<T>& q)
{
    tmat4<T> m(1.0f);
    m = to_mat3(q);
    return m;
}
template<typename T>
inline tquat<T> euler_to_quat(tvec3<T>& euler)
{
    tquat<T> qx = angle_axis(euler.x, tvec3<T>(1.0f, 0.0f, 0.0f));
    tquat<T> qy = angle_axis(euler.y, tvec3<T>(0.0f, 1.0f, 0.0f));
    tquat<T> qz = angle_axis(euler.z, tvec3<T>(0.0f, 0.0f, 1.0f));
    return normalize(qz * qy * qx);
}
template<typename T>
inline tquat<T> to_quat(tmat3<T>& m)
{
    tquat<T> q(0.0f, 0.0f, 0.0f, 1.0f);

    float t;
    if (m[2][2] < 0)
    {
        if (m[0][0] > m[1][1])
        {
            t = 1.0f + m[0][0] - m[1][1] - m[2][2];
            q = tquat<T>(t, m[0][1] + m[1][0], m[2][0] + m[0][2], m[1][2] - m[2][1]);
        }
        else
        {
            t = 1 - m[0][0] + m[1][1] - m[2][2];
            q = tquat<T>(m[0][1] + m[1][0], t, m[1][2] + m[2][1], m[2][0] - m[0][2]);
        }
    }
    else
    {
        if (m[0][0] < -m[1][1])
        {
            t = 1 - m[0][0] - m[1][1] + m[2][2];
            q = tquat<T>(m[2][0] + m[0][2], m[1][2] + m[2][1], t, m[0][1] - m[1][0]);
        }
        else
        {
            t = 1 + m[0][0] + m[1][1] + m[2][2];
            q = tquat<T>(m[1][2] - m[2][1], m[2][0] - m[0][2], m[0][1] - m[1][0], t);
        }
    }

    q *= 0.5f / sqrt(t);

    return normalize(q);
}

//

inline float clamp(float f, float a, float b)
{
    f = f < a ? a : (f > b ? b : f);
    return f;
}

inline float smoothstep(float a, float b, float x)
{
    x = clamp((x - a) / (b - a), 0.0f, 1.0f);
    return x * x * (3 - 2 * x);
}

inline float lerp(float a, float b, float x)
{
    return (a * (1.0f - x)) + (b * x);
}

template<typename T>
inline tvec3<T> lerp(tvec3<T>& a, tvec3<T>& b, float x)
{
    return tvec3<T>(lerp(a.x, b.x, x), lerp(a.y, b.y, x), lerp(a.z, b.z, x));
}

template<typename T>
inline tquat<T> lerp(tquat<T>& a, tquat<T>& b, float x)
{
    return normalize(a * (1.0f - x) + b * x);
}

template<typename T>
inline tquat<T> slerp(tquat<T>& a, tquat<T>& b, float x)
{
    tquat<T> r = tquat<T>(0.0f, 0.0f, 0.0f, 1.0f);
    float d = dot(a, b);

    if (d < 0.0f)
    {
        d = -d;
        r = -b;
    }
    else
    {
        r = b;
    }

    if (d < 0.95f)
    {
        float angle = acosf(d);
        return (a * sinf(angle * (1.0f - x)) + r * sinf(angle * x)) / sinf(angle);
    }
    else
    {
        return lerp(a, b, x);
    }
}

// Transform class (opengl coordinate system)

class transform
{
public:
    transform()
        : t(0.0f, 0.0f, 0.0f),
        r(0.0f, 0.0f, 0.0f, 1.0f),
        s(1.0f, 1.0f, 1.0f)
    {

    }
    void translate(float x, float y, float z)
    {
        translate(vec3(x, y, z));
    }
    void translate(const vec3& vec)
    {
        t = t + vec;
    }

    void rotate(float angle, float axisX, float axisY, float axisZ)
    {
	    rotate(angle, vec3(axisX, axisY, axisZ));
    }
    void rotate(float angle, const vec3& axis)
    {
	    rotate(angle_axis(angle, axis));
    }
    void rotate(const quat& q)
    {
	    r = normalize(q * r);
    }

    void position(float x, float y, float z)
    {
        position(vec3(x, y, z));
    }
    void position(const vec3& position)
    {
        t = position;
    }

    void rotation(float x, float y, float z)
    {
        r = euler_to_quat(vec3(x, y, z));
    }
    void rotation(float x, float y, float z, float w)
    {
        rotation(quat(x, y, z, w));
    }
    void rotation(const quat& rotation)
    {
        r = rotation;
    }

    void scale(float x)
    {
        scale(vec3(x, x, x));
    }
    void scale(float x, float y, float z)
    {
        scale(vec3(x, y, z));
    }
    void scale(const vec3& scale)
    {
        s = scale;
    }

    vec3 position()
    {
        return t;
    }
    quat rotation()
    {
        return r;
    }
    vec3 scale()
    {
        return s;
    }

    vec3 right(){ return matrix()[0]; }
    vec3 up(){ return matrix()[1]; }
    vec3 back() { return matrix()[2]; }
    vec3 left() { return -right(); }
    vec3 down() { return -up(); }
    vec3 forward() { return -back(); }

    void set_transform(mat4& mat)
    {
        t = vec3(mat[3].x, mat[3].y, mat[3].z);
        mat3 rotMat = to_orient_mat3(mat);
        r = to_quat(rotMat);
        vec3 right = mat[0];
        vec3 up = mat[1];
        vec3 back = mat[2];
        s = vec3(length(right), length(up), length(back));
    }

    mat4 matrix()
    {
        return 
            ::gfxm::translate(mat4(1.0f), t) * 
            to_mat4(r) * 
            ::gfxm::scale(mat4(1.0f), s);
    }

    void look_at(const vec3& target, const vec3& forward, const vec3& up = vec3(0.0f, 1.0f, 0.0f), float f = 1.0f)
    {
        f = _max(-1.0f, _min(f, 1.0f));

        mat4 mat = matrix();
        vec3 pos = mat[3];

        vec3 newFwdUnit = normalize(target - pos);
        vec3 rotAxis = normalize(cross(forward, newFwdUnit));

        quat q;
        float d = dot(forward, newFwdUnit);

        const float eps = 0.01f;
        if (fabs(d + 1.0f) <= eps)
        {
            q = angle_axis(pi * f, up);
        }/*
         else if(fabs(d - 1.0f) <= eps)
         {
         q = Au::Math::Quat(0.0f, 0.0f, 0.0f, 1.0f);
         }*/
        else
        {
            float rotAngle = acosf(_max(-1.0f, _min(d, 1.0f))) * f;
            q = angle_axis(rotAngle, rotAxis);
        }

        rotate(q);
    }
private:
    vec3 t;
    quat r;
    vec3 s;
	mat4 m;
};

}

#endif
