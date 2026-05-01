// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "MathServices.h"
#include "Point.h"
#include "math.h"

namespace EZ
{
	void MathServices::GenerateGaussPoints(vector<double>& vdLocations,vector<double>& vdWeights,const unsigned int& iCount)
	{
		vdLocations.resize(iCount);
		vdWeights.resize(iCount);
		double dTemp = 0.0;
		// some of the known cases, just to save time
		if(iCount == 1)
		{
			vdLocations[0] = 0.0;
			vdWeights[0] = 2.0;
		}
		else if(iCount == 2)
		{
			dTemp = 1.0/sqrt(3.0);
			
			vdLocations[0] = -dTemp;
			vdLocations[1] = dTemp;
			
			vdWeights[0] = 1.0;
			vdWeights[1] = 1.0;
		}
		else if(iCount == 3)
		{
			dTemp = sqrt(3.0/5.0);
			vdLocations[0] = -dTemp;
			vdLocations[1] = 0.0;
			vdLocations[2] = dTemp;
			
			dTemp = 5.0/9.0;
			vdWeights[0] = dTemp;
			vdWeights[1] = 8.0/9.0;
			vdWeights[2] = dTemp;
		}
		else if(iCount == 4)
		{
			dTemp = sqrt((3.0 + 2.0*sqrt(6.0/5.0))/7.0);
			vdLocations[0] = -dTemp;
			vdLocations[3] = dTemp;
			
			dTemp = sqrt((3.0 - 2.0*sqrt(6.0/5.0))/7.0);
			vdLocations[1] = -dTemp;
			vdLocations[2] = dTemp;
			
			dTemp = (18.0 - sqrt(30.0))/36.0;
			vdWeights[0] = dTemp;
			vdWeights[3] = dTemp;
			
			dTemp = (18.0 + sqrt(30.0))/36.0;
			vdWeights[1] = dTemp;
			vdWeights[2] = dTemp;
		}
		else if(iCount == 5)
		{
			dTemp = 1.0/3.0*sqrt(5.0 + 2.0*sqrt(10.0/7.0));
			vdLocations[0] = -dTemp;
			vdLocations[4] = dTemp;
			
			dTemp = 1.0/3.0*sqrt(5.0 - 2.0*sqrt(10.0/7.0));
			vdLocations[1] = -dTemp;
			vdLocations[3] = dTemp;
			
			vdLocations[2] = 0.0;
			
			dTemp = 1.0/900.0*(322.0 - 13.0*sqrt(70.0));
			vdWeights[0] = dTemp;
			vdWeights[4] = dTemp;
			
			dTemp = 1.0/900.0*(322.0 + 13.0*sqrt(70.0));
			vdWeights[1] = dTemp;
			vdWeights[3] = dTemp;
			
			vdWeights[2] = 128.0/225.0;
		}
		else
		{
			unsigned int iM = (iCount + 1)/2;
			double dXM = 0.0;
			double dXL = 1.0;
			double dZ = 0.0;
			double dZ1 = 0.0;
			unsigned int i = 0;
			unsigned int j = 0;
			double dError = 0.0;
	
			double p1 = 0.0;
			double p2 = 0.0;
			double p3 = 0.0;
			double pp = 0.0;
	
			for(i = 1; i <= iM ; i++)
			{
				dError = 10.0;
				dZ=cos(PI*(i - 0.25)/(iCount + 0.5));
				while(dError > 3E-14)
				{
					p1=1.0;
					p2=0.0;
					for(j = 1; j <= iCount; j++)
					{
						p3 = p2;
						p2 = p1;
						p1 = ((2.0*j-1.0)*dZ*p2-(j-1.0)*p3)/j;
					}
					pp= iCount*(dZ*p1-p2)/(dZ*dZ-1.0);
					dZ1 = dZ;
					dZ = dZ1 - p1/pp;
					dError = fabs(dZ - dZ1);
				}
				vdLocations[i - 1] = -dZ;
				vdLocations[iCount - i] = dZ;
				vdWeights[i - 1] = 2.0/((1.0-dZ*dZ)*pp*pp);
				vdWeights[iCount - i] = vdWeights[i - 1];
			}
		}
	}
	int MathServices::KroneckerDelta(const int& i,const int& j)
	{
		if(i == j)
		{
			return 1;
		}
		return 0;
	}
	int MathServices::PermutationSymbol(const int& i,const int& j,const int& k)
	{
		int iResult = 0;
		if((i > 3 || i < 1) || (j > 3 || j < 1) || (k > 3 || k < 1))
		{
			iResult = 0;
		}
		if((i == j) || (i == k) || (j == k))
		{
			iResult = 0;
		}
		
		if(i == 1)
		{
			if(j == 2)
			{
				if(k == 3)
				{
					iResult = 1;
				}
			}
			else if(j == 3)
			{
				if(k == 2)
				{
					iResult = -1;
				}
			}
		}
		
		if(i == 2)
		{
			if(j == 1)
			{
				if(k == 3)
				{
					iResult = -1;
				}
			}
			else if(j == 3)
			{
				if(k == 1)
				{
					iResult = 1;
				}
			}
		}
		
		if(i == 3)
		{
			if(j == 1)
			{
				if(k == 2)
				{
					iResult = 1;
				}
			}
			else if(j == 2)
			{
				if(k == 1)
				{
					iResult = -1;
				}
			}
		}
		return iResult;
	}
}


