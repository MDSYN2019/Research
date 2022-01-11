#ifndef workshop_h
#define workshop_h

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include <chrono>
#include <string>
#include <fstream>
#include <random>
#include <cstdint>

#if defined(__INTEL_COMPILER)
#include <malloc.h>
#else
#include <mm_malloc.h>
#endif // defined(__GNUC__)

namespace workshop
{

/**
 * Allocator for aligned data.
 *
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 *
 * Copied from https://gist.github.com/donny-dont/1471329
 */
template <typename T, std::size_t Alignment>
class aligned_allocator
{
	public:
 
		// The following will be the same for virtually all allocators.
		typedef T * pointer;
		typedef const T * const_pointer;
		typedef T& reference;
		typedef const T& const_reference;
		typedef T value_type;
		typedef std::size_t size_type;
		typedef std::ptrdiff_t difference_type;
 
		T * address(T& r) const
		{
			return &r;
		}
 
		const T * address(const T& s) const
		{
			return &s;
		}
 
		std::size_t max_size() const
		{
			// The following has been carefully written to be independent of
			// the definition of size_t and to avoid signed/unsigned warnings.
			return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
		}
 
 
		// The following must be the same for all allocators.
		template <typename U>
		struct rebind
		{
			typedef aligned_allocator<U, Alignment> other;
		} ;
 
		bool operator!=(const aligned_allocator& other) const
		{
			return !(*this == other);
		}
 
		void construct(T * const p, const T& t) const
		{
			void * const pv = static_cast<void *>(p);
 
			new (pv) T(t);
		}
 
		void destroy(T * const p) const
		{
			p->~T();
		}
 
		// Returns true if and only if storage allocated from *this
		// can be deallocated from other, and vice versa.
		// Always returns true for stateless allocators.
		bool operator==(const aligned_allocator& other) const
		{
			return true;
		}
 
 
		// Default constructor, copy constructor, rebinding constructor, and destructor.
		// Empty for stateless allocators.
		aligned_allocator() { }
 
		aligned_allocator(const aligned_allocator&) { }
 
		template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }
 
		~aligned_allocator() { }
 
 
		// The following will be different for each allocator.
		T * allocate(const std::size_t n) const
		{
			// The return value of allocate(0) is unspecified.
			// Mallocator returns NULL in order to avoid depending
			// on malloc(0)'s implementation-defined behavior
			// (the implementation can define malloc(0) to return NULL,
			// in which case the bad_alloc check below would fire).
			// All allocators can return NULL in this case.
			if (n == 0) {
				return NULL;
			}
 
			// All allocators should contain an integer overflow check.
			// The Standardization Committee recommends that std::length_error
			// be thrown in the case of integer overflow.
			if (n > max_size())
			{
				throw std::length_error("aligned_allocator<T>::allocate() - Integer overflow.");
			}
 
			// Mallocator wraps malloc().
			void * const pv = _mm_malloc(n * sizeof(T), Alignment);
 
			// Allocators should throw std::bad_alloc in the case of memory allocation failure.
			if (pv == NULL)
			{
				throw std::bad_alloc();
			}
 
			return static_cast<T *>(pv);
		}
 
		void deallocate(T * const p, const std::size_t n) const
		{
			_mm_free(p);
		}
 
 
		// The following will be the same for all allocators that ignore hints.
		template <typename U>
		T * allocate(const std::size_t n, const U * /* const hint */) const
		{
			return allocate(n);
		}
 
 
		// Allocators are not required to be assignable, so
		// all allocators should have a private unimplemented
		// assignment operator. Note that this will trigger the
		// off-by-default (enabled under /Wall) warning C4626
		// "assignment operator could not be generated because a
		// base class assignment operator is inaccessible" within
		// the STL headers, but that warning is useless.
	private:
		aligned_allocator& operator=(const aligned_allocator&);
};

/** This prints the elements of a vector to the screen */
template<class T>
void print_vector(const std::vector<T> &values)
{
    std::cout << "[";

    for (const T &value : values)
    {
        std::cout << " " << value;
    }

    std::cout << " ]" << std::endl;
}

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::time_point<Clock> TimePoint;

/** This returns a starter point for a timer */
TimePoint start_timer()
{
    return Clock::now();
}

/** Return the number of microseconds since the passed timer was started */
double get_duration(TimePoint timer)
{
    std::chrono::duration<double> secs = (Clock::now() - timer);
    return 1e6 * secs.count();
}

template<class T>
using Array = typename std::vector<T>;

template<class T>
using AlignedArray = typename std::vector< T, aligned_allocator<T,sizeof(T)> >;

template<class T, int N>
using AlignedArrayN = typename std::vector< T, aligned_allocator<T,N> >;

/** Simple function used to get a random height */
float getRandomHeight()
{
    return 1.5 + 2.0 * (0.5 - (std::rand()/(1.0*RAND_MAX)));
}

/** Simple function used to get a random weight */
float getRandomWeight()
{
    return 65.0 + 30 * (0.5 - (std::rand()/(1.0*RAND_MAX)));
}

/** A simple class that holds measurement data
    about a person */
class Person
{
public:
    Person(float _height=0, float _weight=0, float _girth=0)
      : height(_height), weight(_weight), girth(_girth)
    {}

    float height;
    float weight;
    float girth;
};

/** A simple class that holds the coordinates
    for a three-dimensional point */
class Point
{
public:
    Point(float _x=0, float _y=0, float _z=0)
      : x(_x), y(_y), z(_z)
    {}

    float x;
    float y;
    float z;
};

/** A simple class that holds the coordinates
    for an array of three-dimensional points, 
    automatically converting the "array of structs"
    into a "struct of arrays"
*/
class Points
{
public:
    Points(int size=0)
    {
        if (size > 0)
        {
            x = Array<float>(size, 0.0);
            y = Array<float>(size, 0.0);
            z = Array<float>(size, 0.0);
        }
    }

    Points(const Array<Point> &points)
    {
        x = Array<float>(points.size());
        y = Array<float>(points.size());
        z = Array<float>(points.size());

        for (int i=0; i<points.size(); ++i)
        {
            x[i] = points[i].x;
            y[i] = points[i].y;
            z[i] = points[i].z;
        }
    }

    Point operator[](int i) const
    {
        return Point( x[i], y[i], z[i] );
    }

    void set(int i, Point p)
    {
        x[i] = p.x;
        y[i] = p.y;
        z[i] = p.z;
    }

    Array<Point> toArray() const
    {
        Array<Point> p( x.size() );

        for (int i=0; i<x.size(); ++i)
        {
            p[i] = Point(x[i], y[i], z[i]);
        }

        return p;
    }

    Array<float> x;
    Array<float> y;
    Array<float> z;
};

} // end of namespace workshop

#endif

