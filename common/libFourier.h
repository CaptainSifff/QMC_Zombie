#include <cmath>
#include <algorithm>
//See "Numerical Recipes 3rd Ed."
#ifndef LIBFOURIER_H
#define LIBFOURIER_H

template <class T>
inline void fourier1(T data[], const unsigned int N, const int isign)
{
    if ((N < 2) || (N & (N-1)) ) throw("N must be a power of two");
    unsigned int NN = N << 1;
    unsigned int j = 1;
    for (unsigned int i = 1; i < NN; i += 2)
    {
        if (j > i)
        {
            std::swap(data[j-1], data[i-1]);
            std::swap(data[j], data[i]);
        }
        unsigned int m = N;
        while ((m >= 2) && (j > m))
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    unsigned int mmax = 2;
    while (NN > mmax)
    {
        unsigned int istep = mmax << 1;
        T theta = isign * (M_PI * 2.0 / static_cast<T>(mmax));
        T wtemp = std::sin(0.5 * theta);
        T wpr = -2.0 * wtemp * wtemp;
        T wpi = std::sin(theta);
        T wr = 1.0;
        T wi = 0.0;
        for (unsigned int m = 1; m < mmax; m += 2)
        {
            for (unsigned int i = m; i <= NN; i += istep)
            {
                j = i + mmax;
                T tempr = wr * data[j-1] - wi * data[j];
                T tempi = wr * data[j] + wi * data[j-1];
                data[j-1] = data[i-1] - tempr;
                data[j] = data[i] - tempi;
                data[i-1] += tempr;
                data[i] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;//A Trigonometric Recurrence
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
    return;
}

/**
Calculates the fourier transform of a set of N real valued data points
This routine also calculates the inverse transform of a complex data array if it is the transform of real data(The result in this case must be scaled by 2/N). 
*/
template <class T>
void realFT(T* data, const unsigned int N, const unsigned int isign)
{
    T c2;
    T c1 = 0.5;
    T theta = M_PI / static_cast<T>(N >> 1);
    if (isign == 1)
    {
        c2 = -0.5;
        fourier1(data, N/2, 1);
    }
    else
    {
        c2 = 0.5;
        theta = -theta;
    }
    T wtemp = std::sin(0.5 * theta);
    T wpr = -2.0 * wtemp * wtemp;
    T wpi = std::sin(theta);
    T wr = 1.0 + wpr;
    T wi = wpi;
    for (unsigned int i = 1; i < (N >> 2); ++i)
    {
        unsigned int i1 = i + i;
        unsigned int i2 = 1 + i1;
        unsigned int i3 = N - i1;
        unsigned int i4 = 1 + i3;
        T h1r = c1 * (data[i1] + data[i3]);
        T h1i = c1 * (data[i2] - data[i4]);
        T h2r = -c2 * (data[i2] + data[i4]);
        T h2i = c2 * (data[i1] - data[i3]);
        data[i1] = h1r + wr * h2r - wi * h2i;
        data[i2] = h1i + wr * h2i + wi * h2r;
        data[i3] = h1r - wr * h2r + wi * h2i;
        data[i4] = -h1i+ wr * h2i + wi * h2r;
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    if (isign == 1)
    {
        T h1r = data[0];
        data[0] += data[1];
        data[1] = h1r - data[1];
    }
    else
    {
        T h1r = data[0];
        data[0] += data[1];
        data[0] *= c1;
        data[1] = c1 * (h1r - data[1]);
        fourier1(data, N/2, -1);
    }
    return;
}

template <class T>
void twoFFT(const T* const data1, const T* const data2, unsigned int len1, T* const fft1, T* const fft2, const unsigned int len2)
{
    //only forward transform!
    //UNTESTED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    unsigned int N(len1);
    unsigned int nn2 = N+N;
    unsigned int nn3 = 1 + nn2;
    for (unsigned int j = 0, jj = 0; j < N; ++j, jj += 2 )
    {
        fft1[jj] = data1[j];
        fft1[jj+1] = data2[j];
    }
    fourier1(fft1, N, 1);
    fft2[0] = fft1[1];
    fft1[1] = fft2[1] = 0.0;
    for (unsigned int j = 2; j < (N+1); j += 2)
    {
        T rep = 0.5 * (fft1[j] + fft1[nn2 - j]);
        T rem = 0.5 * (fft1[j] - fft1[nn2-j]);
        T aip = 0.5 * (fft1[j+1] + fft1[nn3 - j]);
        T aim = 0.5 * (fft1[j+1] - fft1[nn3-j]);
        fft1[j] = rep;
        fft1[j+1] = aim;
        fft1[nn2-j] = rep;
        fft1[nn3-j] = -aim;
        fft2[j] = aip;
        fft2[j+1] = -rem;
        fft2[nn2-j] = aip;
        fft2[nn3-j] = rem;
    }
    return;
}

template <class T>
void sinFFT(T* const data, unsigned int N)
{
    T theta = M_PI / static_cast<T>(N);
    T wtemp = std::sin(0.5 * theta);
    T wpr = -2.0 * wtemp * wtemp;
    T wpi = sin(theta);
    data[0] = 0.0;
    T wi = 0.0;
    T wr = 1.0;
    for (unsigned int j = 1; j < ((N>>1) + 1); ++j)
    {
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
        T y1 = wi * (data[j] + data[N-j]);
        T y2 = 0.5 * (data[j] - data[N-j]);
        data[j] = y1 + y2;
        data[N - j] = y1 - y2;
    }
    realFT(data, N, 1);
    data[0] *= 0.5;
    T sum = 0.0;
    data[1] = 0.0;
    for (unsigned int j = 0; j < (N - 1); j += 2)
    {
        sum += data[j];
        data[j] = data[j+1];
        data[j+1] = sum;
    }
    return;
}

template <class T>
inline void cosFFT2(T* const y, const unsigned int N, const unsigned int isign)
{
    T theta = 0.5 * M_PI / static_cast<T>(N);
    T wi = 0.0, wr = 1.0;
    T wtemp;
    T wr1 = cos(theta);
    T wi1 = sin(theta);
    T wpr = -2.0 * wi1 * wi1;
    T wpi = sin(2.0 * theta);
    if (isign == 1)
    {
        for (unsigned int i = 0; i < (N/2); ++i)
        {
            T y1 = 0.5 * (y[i] + y[N - 1 - i]);
            T y2 = wi1 * (y[i] - y[N - 1 - i]);
            y[i] = y1 + y2;
            y[N - 1 - i] = y1 - y2;
            wr1 = (wtemp = wr1) * wpr - wi1 * wpi + wr1;
            wi1 = wi1 * wpr + wtemp * wpi + wi1;
        }
        realFT(y, N, 1);
        for (unsigned int i = 2; i < N; i += 2)
        {
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
            T y1 = y[i] * wr - y[i+1] * wi;
            T y2 = y[i+1] * wr + y[i] * wi;
            y[i] = y1;
            y[i+1] = y2;
        }
        T sum = 0.5 * y[1];
        for (int i = N - 1; i > 0; i -= 2)
        {
            T sum1 = sum;
            sum += y[i];
            y[i] = sum1;
        }
    }
    else
    {
        T ytemp = y[N - 1];
        for (unsigned int i = N-1; i > 2; i -= 2)
            y[i] = y[i-2] - y[i];
        y[1] = 2.0 * ytemp;
        for (unsigned int i = 2; i < N; i += 2)
        {
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
            T y1 = y[i] * wr + y[i + 1] * wi;
            T y2 = y[i+1] * wr - y[i] * wi;
            y[i] = y1;
            y[i+1] = y2;
        }
        realFT(y, N , -1);
        for (unsigned int i = 0; i < N/2; ++i )
        {
            T y1 = y[i] + y[N - 1 - i];
            T y2 = (0.5/wi1) * (y[i] - y[N - 1 - i]);
            y[i] = 0.5 * (y1 + y2);
            y[N - 1 - i] = 0.5 * (y1 - y2);
            wr1 = (wtemp = wr1) * wpr - wi1 * wpi + wr1;
            wi1 = wi1 * wpr + wtemp * wpi + wi1;
        }
    }
    return;
}

#endif