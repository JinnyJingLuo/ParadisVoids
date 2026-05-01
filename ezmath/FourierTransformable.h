#ifndef FOURIERTRANSFORMABLE_H_
#define FOURIERTRANSFORMABLE_H_

#include "ComplexNumber.h"

namespace EZ
{
	class FourierTransformable
	{
	public:
		virtual ~FourierTransformable();
		virtual FourierTransformable& operator=(const FourierTransformable& oTransformable);
		virtual void Reset();
		virtual unsigned int GetDimension() = 0;
		virtual void Read(const string& sFileName) = 0;
		virtual void WriteData(const string& sFileName) = 0;
		virtual void WriteTransform(const string& sFileName) = 0;
		virtual FourierTransformable* GetAutoCorrelation() = 0;
		virtual void Swap() = 0;
		virtual void Scale(const double& dFactor) = 0;

	private:
	
	protected:
		virtual void Initialize();
		static unsigned int GetPaddingSize(const unsigned int& iSize);
		static ComplexNumber* DirectTransform1D(ComplexNumber* poInput,const unsigned int& iSize);
		static ComplexNumber* DirectInvert1D(ComplexNumber* poInput,const unsigned int& iSize);
		static unsigned int Pad(ComplexNumber*& poInput,const unsigned int& iSize);
		static void Pad2D(ComplexNumber**& poInput,const unsigned int& iXSize,const unsigned int& iYSize,unsigned int& iNewXSize,unsigned int& iNewYSize);
		static ComplexNumber* Transform1D(ComplexNumber* poInput,const unsigned int& iSize,const unsigned int& iStepSize);
		static ComplexNumber* Invert1D(ComplexNumber* poInput,const unsigned int& iSize,const unsigned int& iStepSize);
		static ComplexNumber** Transform2D(ComplexNumber** poInput,const unsigned int& iXSize,const unsigned int& iYSize);
		static ComplexNumber** Invert2D(ComplexNumber** poInput,const unsigned int& iXSize,const unsigned int& iYSize);
		static void Free2DData(ComplexNumber**& poData,const unsigned int& iXSize,const unsigned int& iYSize);
		static void Free2DData(double**& poData,const unsigned int& iXSize,const unsigned int& iYSize);
	};
}

#endif

