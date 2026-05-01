// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "iostream"
#include "time.h"
#include "math.h"
#include "Randomizer.h"

using namespace std;

namespace EZ
{
	#define N 624
	#define M 397
	#define MATRIX_A 0x9908b0df
	#define UPPER_MASK 0x80000000
	#define LOWER_MASK 0x7fffffff

	#define TEMPERING_MASK_B 0x9d2c5680
	#define TEMPERING_MASK_C 0xefc60000
	#define TEMPERING_SHIFT_U(y) (y >> 11)
	#define TEMPERING_SHIFT_S(y) (y << 7)
	#define TEMPERING_SHIFT_T(y) (y << 15)
	#define TEMPERING_SHIFT_L(y) (y >> 18)
	static unsigned int mt[N];
	static int mti = N + 1;

	Randomizer::Randomizer()
	{

	}
	Randomizer::Randomizer(const Randomizer& oRNG)
	{
		*this = (Randomizer&)oRNG;
	}
	Randomizer::~Randomizer()
	{

	}
	Randomizer& Randomizer::operator=(Randomizer& oRNG)
	{
		return *this;
	}
	double Randomizer::Random()
	{
		return MTRNG();
	}
	double Randomizer::Random(const double& dMin,const double& dMax)
	{
		if(dMax < dMin)
		{
			return 0.0;
		}
		return (dMin + Random()*(dMax - dMin));
	}
	int Randomizer::RandomInteger(const int& iMin,const int& iMax)
	{
		bool bAccepted = false;
		int iResult = 0;
		while(!bAccepted)
		{
			iResult = (int)floor(Random(iMin - 1,iMax + 1) + 0.5);
			if(iResult != (iMin - 1) && iResult != (iMax + 1))
			{
				bAccepted = true;
			}
		}
		return iResult;
	}
	int Randomizer::RandomSign()
	{
		if(RandomInteger(0,1) == 0)
		{
			return -1;
		}
		return 1;
	}
	double Randomizer::RandomNormal()
	{
		double x1, x2, w, y1, y2;
		do
		{
			x1 = 2.0 * Random() - 1.0;
			x2 = 2.0 * Random() - 1.0;
			w = x1 * x1 + x2 * x2;
		}
		while (w >= 1.0);

		w = sqrt((-2.0*log(w))/w);
		y1 = x1*w;
		y2 = x2*w;
		return y1;
	}
	double Randomizer::RandomNormal(const double& dMean,const double& dStandardDeviation)
	{
		return dMean + dStandardDeviation*RandomNormal();
	}
	double Randomizer::RandomExponential(const double& dMean)
	{
		return (-dMean*log(1.0 - Random()));
	}
	void Randomizer::Seed(unsigned int iSeed)
	{
		mt[0] = iSeed;
		for(mti = 1; mti < N ; mti++)
		{
			mt[mti] = (69069*mt[mti - 1]);
		}
	}
	vector<unsigned int> Randomizer::Shuffle(const unsigned int& iSize,const unsigned int& iPasses)
	{
		vector<unsigned int> viResult;
		unsigned int i = 0;
		for(i = 0 ; i < iSize ; i++)
		{
			viResult[i] = i + 1;
		}
		unsigned int j = 0;
		unsigned int iTempIndex = 0;
		unsigned int iTemp = 0;
		for(i = 0 ; i < iPasses ; i++)
		{
			for(j = 0 ; j < iSize ; j++)
			{
				iTempIndex = RandomInteger(1,iSize) - 1;
				iTemp = viResult[j];
				viResult[j] = viResult[iTempIndex];
				viResult[iTempIndex] = iTemp;
			}
		}
		return viResult;
	}
	double Randomizer::MTRNG()
	{
		unsigned int y = 0;
		static unsigned int mag01[2] = {0x0,MATRIX_A};
		if(mti >= N)
		{
			unsigned int kk = 0;
			if(mti == N + 1)
			{
				Seed((unsigned int)time(NULL));
			}
			for(kk = 0; kk < N - M ; kk++)
			{
				y = (mt[kk]&UPPER_MASK)|(mt[kk + 1]&LOWER_MASK);
				mt[kk] = mt[kk + M]^(y >> 1)^mag01[y & 0x1];
			}
			for(; kk < N - 1; kk++)
			{
				y = (mt[kk]&UPPER_MASK)|(mt[kk + 1]&LOWER_MASK);
				mt[kk] = mt[kk + (M - N)]^(y >> 1)^mag01[y & 0x1];
			}
			y = (mt[N - 1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
			mt[N - 1] = mt[M - 1]^(y >> 1)^mag01[y & 0x1];
			mti = 0;
		}

		y = mt[mti++];
		y = y^TEMPERING_SHIFT_U(y);
		y = y^TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
		y = y^TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
		y = y^TEMPERING_SHIFT_L(y);

		return ((double)y)/((unsigned int)0xffffffff);
	}
}

