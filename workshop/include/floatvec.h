#ifndef FLOATVEC_H
#define FLOATVEC_H

#ifdef __AVX__
  #include <immintrin.h>
#else
#ifdef __SSE2__
  #include <emmintrin.h>
#endif
#endif

#include "workshop.h"

/** A simple class that provides a portable, optimised
    vector of floats
*/
class FloatVec
{
private:
    #ifdef __AVX__
        // if we have AVX, then we have a 8xfloat vector
        __m256 v;
        #define FLOATVEC_SIZE 8
    #else
    #ifdef __SSE2__
        // if we have SSE2, then we have a 4xfloat vector
        __m128 v;
        #define FLOATVEC_SIZE 4
    #else
        // in the general case, we will create a float array of 
        // size 8 (we can choose whatever size we want)
        float v[8];
        #define FLOATVEC_SIZE 8
    #endif
    #endif

public:
    typedef workshop::AlignedArrayN<FloatVec,4*FLOATVEC_SIZE> Array;

    /** Return the size of the vector (number of floats) */
    static int size()
    {
        return FLOATVEC_SIZE;
    }

    FloatVec();
    FloatVec(float val);

    #ifdef __AVX__
        FloatVec( __m256 _v ) : v(_v)
        {}
    #else
    #ifdef __SSE2__
        FloatVec( __m128 _v ) : v(_v)
        {}
    #endif
    #endif

    static FloatVec::Array fromArray(const workshop::Array<float> &values);

    static workshop::Array<float> toArray(const FloatVec::Array &values);

    FloatVec sqrt() const;

    FloatVec operator+(const FloatVec &other) const;
};

/** For speed, we don't initialise any values with the null constructor */
inline FloatVec::FloatVec()
{}

/** Create a vector where all elements equal 'val' */
inline FloatVec::FloatVec(float val)
{
    #ifdef __AVX__
        v = _mm256_set1_ps(val);
    #else
    #ifdef __SSE2__
        v = _mm_set1_ps(val);
    #else
        for (int i=0; i<size(); ++i)
        {
            v[i] = val;
        }
    #endif
    #endif
}

/** Return an array of vectors created from the passed
    array of floats */
inline FloatVec::Array FloatVec::fromArray(const workshop::Array<float> &v)
{
    if (v.size() == 0)
    {
        return FloatVec::Array();
    }

    //how many vectors do we need?
    int nvecs = v.size() / size();

    FloatVec::Array array(nvecs);

    for (int i=0; i<nvecs; ++i)
    {
        int j = i * size();

        #ifdef __AVX__
            array[i].v = _mm256_set_ps( v[j+7], v[j+6], v[j+5], v[j+4],
                                      v[j+3], v[j+2], v[j+1], v[j+0] );
        #else
        #ifdef __SSE2__
            array[i].v = _mm_set_ps( v[j+3], v[j+2], v[j+1], v[j+0] );
        #else
            for (int k=0; k<size(); ++k)
            {
                 array[i].v[k] = v[j+k];
            }
        #endif
        #endif
    }

    return array;
}

/** Return this array of vectors as an array of floats */
inline workshop::Array<float> FloatVec::toArray(const FloatVec::Array &v)
{
    if (v.size() == 0)
    {
        return workshop::Array<float>();
    }

    workshop::Array<float> array( size() * v.size() );

    for (int i=0; i<v.size(); ++i)
    {
        int j = i*size();

        float a[size()];

        #ifdef __AVX__
            _mm256_storeu_ps(a, v[i].v);
        #else
        #ifdef __SSE2__
            _mm_storeu_ps(a, v[i].v);
        #else
            for (int k=0; k<size(); ++k)
            {
                a[k] = v[i].v[k];
            }
        #endif
        #endif

        for (int k=0; k<size(); ++k)
        {
            array[j+k] = a[k];
        }
    }

    return array;
}

/** Add two vectors together */
inline FloatVec FloatVec::operator+(const FloatVec &other) const
{
    #ifdef __AVX__
        return FloatVec( _mm256_add_ps(v, other.v) );
    #else
    #ifdef __SSE2__
        return FloatVec( _mm_add_ps(v, other.v) );
    #else
        FloatVec result;

        #pragma omp simd
        for (int i=0; i<size(); ++i)
        {
            result.v[i] = v[i] + other.v[i];
        }

        return result;
    #endif
    #endif
}

/** Return the square root of the vector */
inline FloatVec FloatVec::sqrt() const
{
    #ifdef __AVX__
        return FloatVec( _mm256_sqrt_ps(v) );
    #else
    #ifdef __SSE2__
        return FloatVec( _mm_sqrt_ps(v) );
    #else
        FloatVec result;

        #pragma omp simd
        for (int i=0; i<size(); ++i)
        {
            result.v[i] = std::sqrt(v[i]);
        }

        return result;
    #endif
    #endif
}

/** Return the square root of the vector */
inline FloatVec sqrt(const FloatVec &v)
{
    return v.sqrt();
}

#endif

