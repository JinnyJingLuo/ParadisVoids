#ifndef FOURIERTRANSFORMABLE2D_H_
#define FOURIERTRANSFORMABLE2D_H_

#include "FourierTransformable.h"

namespace EZ
{
	class FourierTransformable2D : public FourierTransformable
	{
	public:
		FourierTransformable2D();
		FourierTransformable2D(const FourierTransformable2D& oTransformable);
		virtual ~FourierTransformable2D();
		virtual FourierTransformable2D& operator=(const FourierTransformable2D& oTransformable);
		virtual void Reset();
		
		unsigned int GetDimension();
		unsigned int GetXSize() const;
		unsigned int GetYSize() const;

		void Set(double** pdData,const unsigned int& iXSize,const unsigned int& iYSize);
		void Read(const string& sFileName);
		void WriteData(const string& sFileName);
		void WriteTransform(const string& sFileName);
		ComplexNumber** Shift();

		const ComplexNumber** GetData();
		const ComplexNumber** GetTransform();
		const ComplexNumber** Transform();
		const ComplexNumber** Invert();
		
		double* GetPowerSpectrum(unsigned int& iSpectrumSize);
		virtual FourierTransformable* GetAutoCorrelation();
		void Swap();
		void Scale(const double& dFactor);
		
	private:
	
	protected:
		virtual void Initialize();
		virtual void Free();
		unsigned int m_iXSize;
		unsigned int m_iYSize;
		ComplexNumber** m_poData;
		ComplexNumber** m_poFFT;
		bool m_bShifted;
	};
}


#endif


