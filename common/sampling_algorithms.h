#ifndef SAMPLING_ALGORITHMS_H
#define SAMPLING_ALGORITHMS_H
/**
@file sampling_algorithms.h This file contains the possible Classes for doing the Importance sampling
*/

/**
The classic Metropolis Algorithm
*/
class Metropolis
{
public:
template <typename T>
    static inline T F(T) throw();
private:
};

template <typename T>
T Metropolis::F(T z) throw()
{
    return z < 1.0 ? z : 1.0;
}

/**
The classic Heatbath Algorithm
*/
class HeatBath
{
public:
template <typename T>
    static inline T F(T) throw();
private:
};

template <typename T>
T HeatBath::F(T z) throw()
{
    return z/(1.0 + z);
}

#endif
