#ifndef ALPHA_H
#define ALPHA_H

/**
The Alpha function, which encapsulates the optimizing parameter delta for optimizing the QMC Algorithm.
The reason for stuffing this into its own object is, that alpha was needed for calculating the weights and for
calculating the full-particle Greensfunction
*/
template <class FPType>
class Alpha
{
public:
    /**
    Constructor of the alpha function
    @param d the parameter delta for optimizing the QMC simulation
    */
    inline Alpha(FPType d) : delta(d) {}
    inline Alpha() : delta(0xFFFFFFFF) {}
    /**
    determine the value of delta
    @param sigma which spin-sector
    @param s the spin of the vertex
    @return the value of alpha
    */
    inline FPType operator() (FPType sigma, FPType s) const
    {
        return 0.5 + delta * sigma * s;
    }
private:
    const FPType delta;///< here we store the parameter delta
};
#endif
