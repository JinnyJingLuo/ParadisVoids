#include "FourierTransformable2D.h"
#include "math.h"
#include "Tools.h"

using namespace SupportSystem;

namespace EZ
{
	FourierTransformable2D::FourierTransformable2D()
	{
		Initialize();
	}
	FourierTransformable2D::FourierTransformable2D(const FourierTransformable2D& oTransformable)
	{
		*this = oTransformable;
	}
	FourierTransformable2D::~FourierTransformable2D()
	{
		Reset();
	}
	FourierTransformable2D& FourierTransformable2D::operator=(const FourierTransformable2D& oTransformable)
	{
		Reset();
		FourierTransformable::operator=(oTransformable);
		m_iXSize = oTransformable.m_iXSize;
		m_iYSize = oTransformable.m_iYSize;
		m_bShifted = oTransformable.m_bShifted;
		unsigned int i = 0;
		unsigned int j = 0;
		if(oTransformable.m_poData != NULL)
		{
			m_poData = new ComplexNumber*[m_iXSize];
			for(i = 0 ; i < m_iXSize ; i++)
			{
				m_poData[i] = new ComplexNumber[m_iYSize];
				for(j = 0 ; j < m_iYSize ; j++)
				{
					m_poData[i][j] = oTransformable.m_poData[i][j];
				}
			}
		}
		if(oTransformable.m_poFFT != NULL)
		{
			m_poFFT = new ComplexNumber*[m_iXSize];
			for(i = 0 ; i < m_iXSize ; i++)
			{
				m_poFFT[i] = new ComplexNumber[m_iYSize];
				for(j = 0 ; j < m_iYSize ; j++)
				{
					m_poFFT[i][j] = oTransformable.m_poFFT[i][j];
				}
			}
		}

		return *this;
	}
	void FourierTransformable2D::Reset()
	{
		FourierTransformable::Reset();
		Free();
		Initialize();
	}
	unsigned int FourierTransformable2D::GetDimension()
	{
		return 2;
	}
	unsigned int FourierTransformable2D::GetXSize() const
	{
		return m_iXSize;
	}
	unsigned int FourierTransformable2D::GetYSize() const
	{
		return m_iYSize;
	}
	void FourierTransformable2D::Set(double** pdData,const unsigned int& iXSize,const unsigned int& iYSize)
	{
		Reset();
		m_iXSize = iXSize;
		m_iYSize = iYSize;
		unsigned int i = 0;
		unsigned int j = 0;
		m_poData = new ComplexNumber*[m_iXSize];
		double dValue = 0.0;
		for(i = 0 ; i < m_iXSize ; i++)
		{
			m_poData[i] = new ComplexNumber[m_iYSize];
			for(j = 0 ; j < m_iYSize ; j++)
			{
				m_poData[i][j].Set(pdData[i][j],0.0);
			}
		}
	}
	void FourierTransformable2D::Read(const string& sFileName)
	{
		Reset();
		// file format
		// xres yres
		// x y z(x,y)
		FILE* fpFile = fopen(sFileName.c_str(),"r");
		string sRead = GetRealString(1024,fpFile);
		// get the x and y resolutions
		sscanf(sRead.c_str(),"%d\t%d\n",&m_iXSize,&m_iYSize);
		unsigned int i = 0;
		unsigned int j = 0;
		m_poData = new ComplexNumber*[m_iXSize];
		double dValue = 0.0;
		for(i = 0 ; i < m_iXSize ; i++)
		{
			m_poData[i] = new ComplexNumber[m_iYSize];
			for(j = 0 ; j < m_iYSize ; j++)
			{
				string sRead = GetRealString(1024,fpFile);
				sscanf(sRead.c_str(),"%*f\t%*f\t%lf\n",&dValue);
				m_poData[i][j].Set(dValue,0.0);
			}
		}
		fclose(fpFile);
	}
	void FourierTransformable2D::WriteData(const string& sFileName)
	{
		if(m_poData == NULL)
		{
			fprintf(stderr,"error: cannot write data for an empty transformable\n");
			return;
		}
		// file format
		// xres yres
		// x y z(x,y)
		FILE* fpFile = fopen(sFileName.c_str(),"w");
		// write the x and y resolutions
		fprintf(fpFile,"%d\t%d\n",m_iXSize,m_iYSize);
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iXSize ; i++)
		{
			for(j = 0 ; j < m_iYSize ; j++)
			{
				fprintf(fpFile,"%e\t%e\n",m_poData[i][j].GetReal(),m_poData[i][j].GetImaginary());
			}
		}
		fclose(fpFile);
	}
	void FourierTransformable2D::WriteTransform(const string& sFileName)
	{
		if(m_poFFT == NULL)
		{
			if(m_poData != NULL)
			{
				Transform();
			}
			else
			{
				fprintf(stderr,"error: cannot write transform for an empty transformable\n");
				return;
			}
		}
		// file format
		// xres yres
		// x y z(x,y)
		FILE* fpFile = fopen(sFileName.c_str(),"w");
		// get the x and y resolutions
		fprintf(fpFile,"%d\t%d\n",m_iXSize,m_iYSize);
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iXSize ; i++)
		{
			for(j = 0 ; j < m_iYSize ; j++)
			{
				fprintf(fpFile,"%e\t%e\n",m_poFFT[i][j].GetReal(),m_poFFT[i][j].GetImaginary());
			}
		}
		fclose(fpFile);
	}
	ComplexNumber** FourierTransformable2D::Shift()
	{
		if(m_poFFT == NULL)
		{
			if(m_poData != NULL)
			{
				Transform();
			}
			else
			{
				fprintf(stderr,"error: cannot shift an empty transformable\n");
				return NULL;
			}
		}
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int iXShiftIndex = 0;
		unsigned int iYShiftIndex = 0;
		unsigned int iXShift = m_iXSize/2;
		unsigned int iYShift = m_iYSize/2;
		ComplexNumber** poFFTShift = new ComplexNumber*[m_iXSize];
		for(i = 0 ; i < m_iXSize ; i++)
		{
			poFFTShift[i] = new ComplexNumber[m_iYSize];
		}
		for(i = 0 ; i < m_iXSize ; i++)
		{
			if(i < iXShift)			iXShiftIndex = i + iXShift;
			else					iXShiftIndex = i - iXShift;
			for(j = 0 ; j < m_iYSize ; j++)
			{
				if(j < iYShift)			iYShiftIndex = j + iYShift;
				else					iYShiftIndex = j - iYShift;
				poFFTShift[iXShiftIndex][iYShiftIndex] = m_poFFT[i][j];
			}
		}
		Free2DData(m_poFFT,m_iXSize,m_iYSize);
		m_poFFT = poFFTShift;
		m_bShifted = !m_bShifted;
		return m_poFFT;
	}
	const ComplexNumber** FourierTransformable2D::GetData()
	{
		return (const ComplexNumber**)m_poData;
	}
	const ComplexNumber** FourierTransformable2D::GetTransform()
	{
		return (const ComplexNumber**)m_poFFT;
	}
	const ComplexNumber** FourierTransformable2D::Transform()
	{
		if(m_poData == NULL)
		{
			fprintf(stderr,"error: cannot transform an empty transformable\n");
			return NULL;
		}
		Free2DData(m_poFFT,m_iXSize,m_iYSize);
		m_bShifted = false;
		unsigned int iNewXSize = m_iXSize;
		unsigned int iNewYSize = m_iYSize;
		Pad2D(m_poData,m_iXSize,m_iYSize,iNewXSize,iNewYSize);
		m_iXSize = iNewXSize;
		m_iYSize = iNewYSize;
		m_poFFT = Transform2D(m_poData,m_iXSize,m_iYSize);
		return (const ComplexNumber**)m_poFFT;
	}
	const ComplexNumber** FourierTransformable2D::Invert()
	{
		if(m_poFFT == NULL)
		{
			fprintf(stderr,"error: cannot invert an untransformed transformable\n");
			return NULL;
		}
		Free2DData(m_poData,m_iXSize,m_iYSize);
		bool bOriginallyShifted = m_bShifted;
		if(bOriginallyShifted)
		{
			Shift();
		}
		unsigned int iNewXSize = m_iXSize;
		unsigned int iNewYSize = m_iYSize;
		Pad2D(m_poFFT,m_iXSize,m_iYSize,iNewXSize,iNewYSize);
		m_iXSize = iNewXSize;
		m_iYSize = iNewYSize;
		m_poData = Invert2D(m_poFFT,m_iXSize,m_iYSize);
		if(bOriginallyShifted)
		{
			Shift();
		}
		return (const ComplexNumber**)m_poData;
	}
	double* FourierTransformable2D::GetPowerSpectrum(unsigned int& iSpectrumSize)
	{
		if(m_poFFT == NULL)
		{
			if(m_poData != NULL)
			{
				Transform();
			}
			else
			{
				fprintf(stderr,"error: cannot compute the power spectrum for an empty transformable\n");
				return NULL;
			}
		}
		// get the power array first
		unsigned int i = 0;
		unsigned int j = 0;
		double** pdPower = new double*[m_iXSize];
		for(i = 0 ; i < m_iXSize ; i++)
		{
			pdPower[i] = new double[m_iYSize];
			for(j = 0 ; j < m_iYSize ; j++)
			{
				pdPower[i][j] = m_poFFT[i][j].GetSquareAmplitude();
			}
		}

		if(m_bShifted)		iSpectrumSize = (unsigned int)floor(sqrt(m_iXSize*m_iXSize/4 + m_iYSize*m_iYSize/4) + 0.5);
		else				iSpectrumSize = (unsigned int)floor(sqrt(m_iXSize*m_iXSize + m_iYSize*m_iYSize) + 0.5);
		double* pdSpectrum = new double[iSpectrumSize];
		unsigned int* piCounts = new unsigned int[iSpectrumSize];
		for(i = 0 ; i < iSpectrumSize ; i++)
		{
			pdSpectrum[i] = 0.0;
			piCounts[i] = 0;
		}
		int iIndex = 0;
		unsigned int iXShift = m_iXSize/2;
		unsigned int iYShift = m_iYSize/2;
		for(i = 0 ; i < m_iXSize ; i++)
		{
			for(j = 0 ; j < m_iYSize ; j++)
			{	
				if(m_bShifted)	iIndex = (int)floor(sqrt((i - iXShift)*(i - iXShift) + (j - iYShift)*(j - iYShift)) + 0.5);
				else			iIndex = (int)floor(sqrt(i*i + j*j) + 0.5);
				if(iIndex < 0)					iIndex = 0;
				if(iIndex >= iSpectrumSize)		iIndex = iSpectrumSize - 1;
				pdSpectrum[iIndex] += pdPower[i][j];
				piCounts[iIndex]++;
			}
			delete [] pdPower[i];
		}
		delete [] pdPower;
		// get the average value for each index in the spectrum array
		for(i = 0 ; i < iSpectrumSize ; i++)
		{
			if(piCounts[i] == 0)
			{
				continue;
			}
			pdSpectrum[i] = pdSpectrum[i]/(double)piCounts[i];
		}
		delete [] piCounts;
		return pdSpectrum;
	}
	FourierTransformable* FourierTransformable2D::GetAutoCorrelation()
	{
		if(m_poFFT == NULL)
		{
			if(m_poData != NULL)
			{
				Transform();
			}
			else
			{
				fprintf(stderr,"error: cannot compute autocorrelation function for an empty transformable\n");
				return NULL;
			}
		}
		bool bOriginallyShifted = m_bShifted;
		if(bOriginallyShifted)
		{
			Shift();
		}
		// get the power array first
		unsigned int i = 0;
		unsigned int j = 0;
		double** pdPower = new double*[m_iXSize];
		double dSizeFactor = m_iXSize*m_iYSize;
		for(i = 0 ; i < m_iXSize ; i++)
		{
			pdPower[i] = new double[m_iYSize];
			for(j = 0 ; j < m_iYSize ; j++)
			{
				pdPower[i][j] = m_poFFT[i][j].GetSquareAmplitude()/dSizeFactor;
			}
		}
		FourierTransformable2D* poAutoCorrelation = new FourierTransformable2D;
		poAutoCorrelation->Set(pdPower,m_iXSize,m_iYSize);
		poAutoCorrelation->Swap();
		poAutoCorrelation->Invert();
		for(i = 0 ; i < m_iXSize ; i++)
		{
			delete [] pdPower[i];
		}
		delete [] pdPower;
		if(bOriginallyShifted)
		{
			Shift();
		}
		return poAutoCorrelation;
	}
	void FourierTransformable2D::Swap()
	{
		ComplexNumber** poTemp = m_poFFT;
		m_poFFT = m_poData;
		m_poData = poTemp;
	}
	void FourierTransformable2D::Scale(const double& dFactor)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		if(m_poData != NULL)
		{
			for(i = 0 ; i < m_iXSize ; i++)
			{
				for(j = 0 ; j < m_iYSize ; j++)
				{
					m_poData[i][j] = m_poData[i][j]*dFactor;
				}
			}
		}
		if(m_poFFT != NULL)
		{
			for(i = 0 ; i < m_iXSize ; i++)
			{
				for(j = 0 ; j < m_iYSize ; j++)
				{
					m_poFFT[i][j] = m_poFFT[i][j]*dFactor;
				}
			}
		}
	}
	void FourierTransformable2D::Initialize()
	{
		FourierTransformable::Initialize();
		m_iXSize = 0;
		m_iYSize = 0;
		m_poData = NULL;
		m_poFFT = NULL;
		m_bShifted = false;
	}
	void FourierTransformable2D::Free()
	{
		Free2DData(m_poData,m_iXSize,m_iYSize);
		Free2DData(m_poFFT,m_iXSize,m_iYSize);
	}
}

// namespace EZ
// {
// 	ComplexNumber* FFT::ForwardTransform1D(ComplexNumber*& poInput,unsigned int& iSize)
// 	{
// 		unsigned int iNewSize = Pad(poInput,iSize);
// 		iSize = iNewSize;
// 		ComplexNumber* poFFT = FFT::Transform1D(poInput,iSize,1);
// 		return poFFT;
// 	}
// 	ComplexNumber* FFT::InverseTransform1D(ComplexNumber*& poInput,unsigned int& iSize)
// 	{
// 		unsigned int iNewSize = Pad(poInput,iSize);
// 		iSize = iNewSize;
// 		ComplexNumber* poIFFT = FFT::Invert1D(poInput,iSize,1);
// 		return poIFFT;
// 	}
// }

