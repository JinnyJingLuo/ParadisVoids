#include "FourierTransformable.h"
#include "Point.h"
#include "math.h"

namespace EZ
{
	FourierTransformable::~FourierTransformable()
	{
		Reset();
	}
	FourierTransformable& FourierTransformable::operator=(const FourierTransformable& oTransformable)
	{
		return *this;
	}
	void FourierTransformable::Reset()
	{
	
	}
	void FourierTransformable::Initialize()
	{
	
	}
	unsigned int FourierTransformable::GetPaddingSize(const unsigned int& iSize)
	{
		double dPower = log(iSize)/log(2.0);
		double dTolerance = 1.0E-6;
		if(fabs(dPower - floor(dPower)) < dTolerance)
		{
			return iSize;
		}
		int iPower = (int)floor(dPower) + 1;
		unsigned int iNewSize = (unsigned int)floor(pow(2.0,iPower) + 0.5);
		return iNewSize;
	}
	unsigned int FourierTransformable::Pad(ComplexNumber*& poInput,const unsigned int& iSize)
	{
		unsigned int iNewSize = GetPaddingSize(iSize);
		// if the sizes are the same, don't do anything
		if(iSize == iNewSize)
		{
			return iSize;
		}
		// the size is not a power of two, pad zeros at the end of the array
		ComplexNumber* poPaddedInput = new ComplexNumber[iNewSize];
		unsigned int i = 0;
		for(i = 0 ; i < iNewSize ; i++)
		{
			if(i < iSize)		poPaddedInput[i] = poInput[i];
			else				poPaddedInput[i].Set(0.0,0.0);
		}
		// delete the original array and make its pointer point to the newly created padded array
		delete [] poInput;
		poInput = poPaddedInput;
		return iNewSize;
	}
	void FourierTransformable::Pad2D(ComplexNumber**& poInput,const unsigned int& iXSize,const unsigned int& iYSize,unsigned int& iNewXSize,unsigned int& iNewYSize)
	{
		iNewXSize = GetPaddingSize(iXSize);
		iNewYSize = GetPaddingSize(iYSize);
		// if all the sizes are the same, don't do anything
		if((iXSize == iNewXSize) && (iYSize == iNewYSize))
		{
			return;
		}
		// the size is not a power of two, pad zeros at the end of the array
		ComplexNumber** poPaddedInput = new ComplexNumber*[iNewXSize];
		unsigned int i = 0;
		for(i = 0 ; i < iNewXSize ; i++)
		{
			poPaddedInput[i] = new ComplexNumber[iNewYSize];
		}
		unsigned int j = 0;
		for(i = 0 ; i < iNewXSize ; i++)
		{
			if(i < iXSize)
			{
				for(j = 0 ; j < iNewYSize ; j++)
				{
					if(j < iYSize)		poPaddedInput[i][j] = poInput[i][j];
					else				poPaddedInput[i][j].Set(0.0,0.0);
				}
			}
			else
			{
				for(j = 0 ; j < iNewYSize ; j++)
				{
					poPaddedInput[i][j].Set(0.0,0.0);
				}
			}
		}
		// delete the original array
		for(i = 0 ; i < iXSize ; i++)
		{
			delete [] poInput[i];
		}
		delete [] poInput;
		// make its pointer point to the newly created padded array
		poInput = poPaddedInput;
	}
	ComplexNumber* FourierTransformable::DirectTransform1D(ComplexNumber* poInput,const unsigned int& iSize)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		ComplexNumber* poOutput = new ComplexNumber[iSize];
		ComplexNumber oFactor;
		double dFactor = -2.0*PI/(double)iSize;
		for(i = 0 ; i < iSize ; i++)
		{
			poOutput[i].Set(0.0,0.0);
			for(j = 0 ; j < iSize ; j++)
			{
				oFactor.SetAmplitudeAndPhase(1.0,dFactor*i*j);
				poOutput[i] = poOutput[i] + oFactor*poInput[j];
			}
			poOutput[i] = poOutput[i]/(double)iSize;
		}
		return poOutput;
	}
	ComplexNumber* FourierTransformable::DirectInvert1D(ComplexNumber* poInput,const unsigned int& iSize)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		ComplexNumber* poOutput = new ComplexNumber[iSize];
		ComplexNumber oFactor;
		double dFactor = 2.0*PI/(double)iSize;
		for(i = 0 ; i < iSize ; i++)
		{
			poOutput[i].Set(0.0,0.0);
			for(j = 0 ; j < iSize ; j++)
			{
				oFactor.SetAmplitudeAndPhase(1.0,dFactor*i*j);
				poOutput[i] = poOutput[i] + oFactor*poInput[j];
			}
		}
		return poOutput;
	}
	ComplexNumber* FourierTransformable::Transform1D(ComplexNumber* poInput,const unsigned int& iSize,const unsigned int& iStepSize)
	{
		// the input's size is guaranteed to be a power of two
		ComplexNumber* poOutput = new ComplexNumber[iSize];
		if(iSize == 1)
		{
			poOutput[0] = poInput[0];
		}
		else
		{
			ComplexNumber* poEven = Transform1D(poInput,iSize/2,2*iStepSize);
			ComplexNumber* poOdd = Transform1D(&(poInput[iStepSize]),iSize/2,2*iStepSize);
			unsigned int i = 0;
			ComplexNumber oTwiddleFactor;
			double dFactor = -2.0*PI/(double)iSize;
			unsigned int iJump = iSize/2;
			ComplexNumber oTemp;
			for(i = 0 ; i < iJump ; i++)
			{
				oTwiddleFactor.SetAmplitudeAndPhase(1.0,dFactor*(double)i);
				oTemp = oTwiddleFactor*poOdd[i];
				poOutput[i] = poEven[i] + oTemp;
				poOutput[i + iJump] = poEven[i] - oTemp;
			}
			delete [] poEven;
			delete [] poOdd;
		}
		// if this is the topmost function in the recursion heirarchy, divide all the output by N
		if(iStepSize == 1)
		{
			unsigned int i = 0;
			for(i = 0 ; i < iSize ; i++)
			{
				poOutput[i] = poOutput[i]/(double)iSize;
			}
		}
		return poOutput;
	}
	ComplexNumber* FourierTransformable::Invert1D(ComplexNumber* poInput,const unsigned int& iSize,const unsigned int& iStepSize)
	{
		// the input's size is guaranteed to be a power of two
		ComplexNumber* poOutput = new ComplexNumber[iSize];
		if(iSize == 1)
		{
			poOutput[0] = poInput[0];
		}
		else
		{
			ComplexNumber* poEven = Invert1D(poInput,iSize/2,2*iStepSize);
			ComplexNumber* poOdd = Invert1D(&(poInput[iStepSize]),iSize/2,2*iStepSize);
			unsigned int i = 0;
			ComplexNumber oTwiddleFactor;
			double dFactor = 2.0*PI/(double)iSize;
			unsigned int iJump = iSize/2;
			ComplexNumber oTemp;
			for(i = 0 ; i < iJump ; i++)
			{
				oTwiddleFactor.SetAmplitudeAndPhase(1.0,dFactor*(double)i);
				oTemp = oTwiddleFactor*poOdd[i];
				poOutput[i] = poEven[i] + oTemp;
				poOutput[i + iJump] = poEven[i] - oTemp;
			}
			delete [] poEven;
			delete [] poOdd;
		}
		return poOutput;
	}
	ComplexNumber** FourierTransformable::Transform2D(ComplexNumber** poInput,const unsigned int& iXSize,const unsigned int& iYSize)
	{
		// allocate arrays for the output and for the intermediate calculations
		ComplexNumber** poColumnFFT = new ComplexNumber*[iXSize];
		// get the 1D FFT for each column in the 2D signal (a column has a fixed x value and a size of iYSize)
		unsigned int i = 0;
		for(i = 0 ; i < iXSize ; i++)
		{
			poColumnFFT[i] = Transform1D(poInput[i],iYSize,1);
		}
		// transpose it
		ComplexNumber** poColumnFFTTranspose = new ComplexNumber*[iYSize];
		unsigned int j = 0;
		for(j = 0 ; j < iYSize ; j++)
		{
			poColumnFFTTranspose[j] = new ComplexNumber[iXSize];
			for(i = 0 ; i < iXSize ; i++)
			{
				poColumnFFTTranspose[j][i] = poColumnFFT[i][j];
			}
		}
		// delete the intermediate
		for(i = 0 ; i < iXSize ; i++)
		{
			delete [] poColumnFFT[i];
		}
		delete [] poColumnFFT;
		// do an FFT on the transposed ffts, these will be the output
		ComplexNumber** poOutputTranspose = new ComplexNumber*[iYSize];
		for(j = 0 ; j < iYSize ; j++)
		{
			poOutputTranspose[j] = Transform1D(poColumnFFTTranspose[j],iXSize,1);
			delete [] poColumnFFTTranspose[j];
		}
		delete [] poColumnFFTTranspose;
		// for organizational purposes, put it back in the same shape it came in, and 
		ComplexNumber** poOutput = new ComplexNumber*[iXSize];
		for(i = 0 ; i < iXSize ; i++)
		{
			poOutput[i] = new ComplexNumber[iYSize];
			for(j = 0 ; j < iYSize ; j++)
			{
				poOutput[i][j] = poOutputTranspose[j][i];
			}
		}
		// delete the output transpose
		for(j = 0 ; j < iYSize ; j++)
		{
			delete [] poOutputTranspose[j];
		}
		delete [] poOutputTranspose;
		// return the output
		return poOutput;
	}
	ComplexNumber** FourierTransformable::Invert2D(ComplexNumber** poInput,const unsigned int& iXSize,const unsigned int& iYSize)
	{
		// allocate arrays for the output and for the intermediate calculations
		ComplexNumber** poColumnIFFT = new ComplexNumber*[iXSize];
		// get the 1D FFT for each column in the 2D signal (a column has a fixed x value and a size of iYSize)
		unsigned int i = 0;
		for(i = 0 ; i < iXSize ; i++)
		{
			poColumnIFFT[i] = Invert1D(poInput[i],iYSize,1);
		}
		// transpose it
		ComplexNumber** poColumnIFFTTranspose = new ComplexNumber*[iYSize];
		unsigned int j = 0;
		for(j = 0 ; j < iYSize ; j++)
		{
			poColumnIFFTTranspose[j] = new ComplexNumber[iXSize];
			for(i = 0 ; i < iXSize ; i++)
			{
				poColumnIFFTTranspose[j][i] = poColumnIFFT[i][j];
			}
		}
		// delete the intermediate
		for(i = 0 ; i < iXSize ; i++)
		{
			delete [] poColumnIFFT[i];
		}
		delete [] poColumnIFFT;
		// do an FFT on the transposed ffts, these will be the output
		ComplexNumber** poOutputTranspose = new ComplexNumber*[iYSize];
		for(j = 0 ; j < iYSize ; j++)
		{
			poOutputTranspose[j] = Invert1D(poColumnIFFTTranspose[j],iXSize,1);
			delete [] poColumnIFFTTranspose[j];
		}
		delete [] poColumnIFFTTranspose;
		// for organizational purposes, put it back in the same shape it came in, and 
		ComplexNumber** poOutput = new ComplexNumber*[iXSize];
		for(i = 0 ; i < iXSize ; i++)
		{
			poOutput[i] = new ComplexNumber[iYSize];
			for(j = 0 ; j < iYSize ; j++)
			{
				poOutput[i][j] = poOutputTranspose[j][i];
			}
		}
		// delete the output transpose
		for(j = 0 ; j < iYSize ; j++)
		{
			delete [] poOutputTranspose[j];
		}
		delete [] poOutputTranspose;
		// return the output
		return poOutput;
	}
	void FourierTransformable::Free2DData(ComplexNumber**& poData,const unsigned int& iXSize,const unsigned int& iYSize)
	{
		if(poData == NULL)
		{
			return;
		}
		unsigned int i = 0;
		for(i = 0 ; i < iXSize ; i++)
		{
			delete [] poData[i];
		}
		delete [] poData;
		poData = NULL;
	}
	void FourierTransformable::Free2DData(double**& pdData,const unsigned int& iXSize,const unsigned int& iYSize)
	{
		if(pdData == NULL)
		{
			return;
		}
		unsigned int i = 0;
		for(i = 0 ; i < iXSize ; i++)
		{
			delete [] pdData[i];
		}
		delete [] pdData;
		pdData = NULL;
	}
}

