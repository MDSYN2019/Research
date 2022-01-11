
#include "workshop.h"

int main(int argc, char **argv)
{
    const int size = 256;

    auto a = workshop::Array<float>(size);
    auto b = workshop::Array<float>(size);
    auto c = workshop::Array<float>(size);

    for (int i=0; i<size; ++i)
    {
        a[i] = 1.0*(i+1);
        b[i] = 2.5*(i+1);
        c[i] = 0.0;
    }

    auto timer = workshop::start_timer();

    for (int j=0; j<100000; ++j)
    {
        for (int i=0; i<size; ++i)
        {
            c[i] = a[i] + b[i];
        }
    }

    auto duration = workshop::get_duration(timer);

    timer = workshop::start_timer();

    for (int j=0; j<100000; ++j)
    {    
        #pragma omp simd
        for (int i=0; i<size; ++i)
        {
            c[i] = a[i] + b[i];
        }
    }

    auto duration_vec = workshop::get_duration(timer);

    if (duration_vec > 0.5*duration)
    {
        std::cout << "pragma omp simd has not sped up the loop!" << std::endl;
        std::cout << duration_vec << " microsecs versus " << duration << " microsecs"
                  << std::endl;

        std::cout << "This could be because either pragma omp simd doesn't work, "
                  << "or because your compiler is auto-vectorising the basic loop."
                  << std::endl;
    }
    else
    {
        std::cout << "Everything is OK :-)" << std::endl;
    }

    return 0;
}

