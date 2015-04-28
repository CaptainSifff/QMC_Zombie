# ifndef MONTE_CARLO_TOOLS_ARRAY_OBSERVABLE_H
# define MONTE_CARLO_TOOLS_ARRAY_OBSERVABLE_H
# include "errordata.h"
# include "analysis_functions.h"
# include "analysis_common.h"
namespace mc_analysis
{
	template <typename Scalar, class Cont = std::vector<std::valarray<Scalar> > >  class array_observable : public Cont
	{
		private:
			std::valarray< std::valarray<Scalar> > data;///< the binned data
			unsigned int num_bins;///< the number of bins
			unsigned int binsize;///< how many measurements make up a bin
                        using Cont::size;//make the size function of the base class inaccesible to the outside world
                        using Cont::operator[];

		public:
			inline array_observable(void);
                        inline virtual ~array_observable();
			inline std::size_t get_num_data_points(void){return Cont::size();};
			inline std::valarray<Scalar> mean(void);
			inline std::valarray<Scalar> sd(void); // This should be implemented by a general function, as mean
	//		std::valarray<Scalar> standard_error(void);
			inline errordata<std::valarray<Scalar> > jackknife(void);
                        /**
                        You call this function to rebin the measurements to the supplied binsize
                        @param binsize__ how many measurements make up a bin
                        */
			inline void rebin_data(unsigned int binsize__);
	//		void fractional_rebinning(unsigned int num_bins__);
			inline void print_data(void);

			template <typename Func> inline errordata<typename Func::res_t> jack_eval(Func & f);
                        using Cont::value_type;
                        /**
                        These access operators are handy for the jackknife analysis and provide access to the binned data
                        @param k the index of the Element
                        @return a reference to the element
                        */
                        inline std::valarray<Scalar>& operator[] (std::size_t k) {return data[k];}
                        inline const std::valarray<Scalar>& operator[] (std::size_t k) const {return data[k];}
                        inline std::size_t size() const {return data.size();}
	};

	template <typename Scalar, class Cont> 
        array_observable<Scalar, Cont>::~array_observable(void)
	{
            mc_analysis_common::CallClear<Cont>::destroy(static_cast<Cont&>(*this));
	}

	template <typename Scalar, class Cont> array_observable<Scalar, Cont>::array_observable(void)
	{
		num_bins=0;
		binsize=0;
		data.resize(0);	
	}


	template <typename Scalar, class Cont> std::valarray<Scalar> array_observable<Scalar, Cont>::mean(void)
	{
		return mc_analysis::mean(data);		
	}

	template <typename Scalar, class Cont> std::valarray<Scalar> array_observable<Scalar,Cont>::sd(void)
	{
		std::valarray<Scalar> sum(data[0].size());
		sum=0.;
		std::valarray<Scalar> mean(data[0].size());
		mean=this->mean();
		for(unsigned int i=0; i<num_bins; i++) sum += (data[i]-mean)*(data[i]-mean);
		return std::sqrt( sum/static_cast<Scalar> (num_bins -1.) );
	}



	template <typename Scalar, class Cont> errordata<std::valarray<Scalar> > array_observable<Scalar, Cont>::jackknife(void)
		/**
		 *	\brief this function performs a simple jackknife analysis for the data
		 *	it calculates the jackknife mean and the standard error.
		 */
	{
		// For this special case, the Jackknife mean is the same as the arithmetic average of the data
		errordata<std::valarray<Scalar> > errd;
		errd.set_mean(this->mean());
		// The Jackknife error for this case also reduces to the standard error:
		errd.set_error(this->sd()/std::sqrt( static_cast<Scalar>(num_bins) ));
		return errd;		
	}




	template <typename Scalar, class Cont> void array_observable<Scalar, Cont>::rebin_data(unsigned int binsize__)
		/**
		 *	\brief This function has to be called before using the data analysis functions. Rebinning is mandatory
		 *	as it fills the data vector. The original data remains intact in the class
		 *	\param binsize__ Specifies, how many data points should be averaged to give a new bin
		 */
	{
		this->binsize=binsize__;
		this->num_bins= Cont::size()/binsize__;//how many functions do we have
		data.resize(num_bins, Cont::operator[](0));//make place for num_bins functions
		for(unsigned int i=0; i<this->num_bins; i++)
		{
			std::valarray<Scalar> sum(Cont::operator[]( i * this->binsize  ));
			for(unsigned int j = 1; j < (this->binsize); j++)
			{
				sum += Cont::operator[]( i*this->binsize + j );
			}
			data[i] = sum;
			data[i] *= 1./static_cast<double>(this->binsize); 
		}
	}


	template <typename Scalar, class Cont> void  array_observable<Scalar, Cont>::print_data(void)
	{
		for(unsigned int i=0; i< data.size() ; i++)
		{
			for(unsigned int j=0; j< data[i].size(); j++)
			{
				std::cout << data[i][j] << "\t";			
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}



	template <typename Scalar, class Cont> template <typename Func> errordata<typename Func::res_t> array_observable<Scalar, Cont>::jack_eval(Func & f)
		/**
		 * \brief This function calculates the mean and error of the data after applying a function
		 * specified in a Function object.
		 * 
		 * \param F is a functor defined by the class Functor. You need to overload the operator().
		 * and to specify the typedef res_t for the result/return type.
		 *
		 * \warning Rebinning before calling this function is absolutely required!
		 */
	{
		// Generate Jackknife bins
		Scalar length = static_cast<Scalar>(this->num_bins); 
		std::valarray<Scalar> mean( this->mean() );
		typename Func::res_t fun_mean=f(mean);

		std::valarray<typename Func::res_t> jackdata(fun_mean,this->num_bins); // Init the jackknife bins with f(mean), to make them large enough


		for(unsigned int i=0; i < this->num_bins; i++)
		{
			std::valarray<Scalar> x_J(mean.size());
			x_J = (length * mean - data[i]);
			x_J *= 1./(length - 1.0);
			jackdata[i] = f(x_J);
		}

		typename Func::res_t jackmean(jackdata[0]);
		for(unsigned int i=1; i< this->num_bins;i++)
		{
			jackmean += jackdata[i];
		}
		jackmean *= 1.0/(length);
		

		typename Func::res_t fun_error= (jackmean - jackdata[0])*(jackmean - jackdata[0]);

		for(unsigned int i=1; i< this->num_bins; i++)
		{
			fun_error+= (jackmean - jackdata[i])*(jackmean - jackdata[i]);
		}
		fun_error *= (length-1.)/length;
		fun_error = std::sqrt(fun_error);
		
		errordata<typename Func::res_t> errd;
		errd.set_mean(fun_mean);
		errd.set_error(fun_error);
		return errd;
	}
}
# endif
