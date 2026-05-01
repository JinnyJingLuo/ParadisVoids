// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef MATHSERVICES_H_
#define MATHSERVICES_H_

#include "vector"

using namespace std;

namespace EZ
{
	class MathServices
	{
	public:
		static void GenerateGaussPoints(vector<double>& vdLocations,vector<double>& vdWeights,const unsigned int& iCount);
		static int KroneckerDelta(const int& i,const int& j);
		static int PermutationSymbol(const int& i,const int& j,const int& k);
	private:
		
	protected:
		
	};
}

#endif


