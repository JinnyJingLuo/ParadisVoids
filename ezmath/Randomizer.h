// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef RANDOMIZER_H_
#define RANDOMIZER_H_

#include "vector"

using namespace std;

namespace EZ
{
	class Randomizer
	{
	public:
		Randomizer();
		Randomizer(const Randomizer& oRNG);
		~Randomizer();
		Randomizer& operator=(Randomizer& oRNG);
		static double Random();
		static double Random(const double& dMin,const double& dMax);
		static int RandomInteger(const int& iMin,const int& iMax);
		static int RandomSign();
		static double RandomNormal();
		static double RandomNormal(const double& dMean,const double& dStandardDeviation);
		static double RandomExponential(const double& dMean);
		static vector<unsigned int> Shuffle(const unsigned int& iSize,const unsigned int& iPasses = 1);
	private:

	protected:
		static void Seed(unsigned int seed);
		static double MTRNG();
	};
}

#endif


