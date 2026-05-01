#include "SurfaceProfile.h"
#include "math.h"
#include "float.h"
#include "Vector.h"

SurfaceProfile::SurfaceProfile()
{
	Initialize();
}
SurfaceProfile::SurfaceProfile(const SurfaceProfile& oProfile)
{
	*this = oProfile;
}
SurfaceProfile::~SurfaceProfile()
{
	Reset();
}
SurfaceProfile& SurfaceProfile::operator=(const SurfaceProfile& oProfile)
{
	FourierTransformable2D::operator=(oProfile);
	m_dMean = oProfile.m_dMean;
	m_dVariance = oProfile.m_dVariance;
	m_dSkewness = oProfile.m_dSkewness;
	m_dKurtosis = oProfile.m_dKurtosis;
	m_dMinPit = oProfile.m_dMinPit;
	m_dMaxPeak = oProfile.m_dMaxPeak;
	m_dMeanHeight = oProfile.m_dMeanHeight;
	m_dHeight = oProfile.m_dHeight;
	m_dXMin = oProfile.m_dXMin;
	m_dXMax = oProfile.m_dXMax;
	m_dYMin = oProfile.m_dYMin;
	m_dYMax = oProfile.m_dYMax;
	return *this;
}
void SurfaceProfile::Reset()
{
	FourierTransformable2D::Reset();
}
FourierTransformable* SurfaceProfile::GetAutoCorrelation()
{
	FourierTransformable* poAutoCorrelation = FourierTransformable2D::GetAutoCorrelation();
	if(poAutoCorrelation == NULL)
	{
		return NULL;
	}
	poAutoCorrelation->Scale(1.0/m_dVariance);
	return poAutoCorrelation;
}
void SurfaceProfile::GetAmplitudeCharacteristics(const string& sFileName)
{
	if(m_poData == NULL)		return;
	unsigned int i = 0;
	unsigned int j = 0;
	
	m_dMeanHeight = 0.0;
	m_dMaxPeak = -DBL_MAX;
	m_dMinPit = DBL_MAX;
	
	double dX = 0.0;
	double dX2 = 0.0;
	double dX3 = 0.0;
	double dX4 = 0.0;
	double dValue = 0.0;
	
	for(i = 0 ; i < m_iXSize ; i++)
	{
		for(j = 0 ; j < m_iYSize ; j++)
		{
			dValue = m_poData[i][j].GetReal();
			dX += dValue;
			dX2 += dValue*dValue;
			dX3 += dValue*dValue*dValue;
			dX4 += dValue*dValue*dValue*dValue;
			m_dMeanHeight += fabs(dValue);
			if(dValue < m_dMinPit)
			{
				m_dMinPit = dValue;
			}
			if(dValue > m_dMaxPeak)
			{
				m_dMaxPeak = dValue;
			}
		}
	}
	double dSize = (double)m_iXSize*(double)m_iYSize;
	m_dMean = dX/dSize;
	m_dVariance = dX2/dSize - m_dMean*m_dMean;
	m_dSkewness = dX3/dSize - 3.0*m_dMean*m_dVariance - m_dMean*m_dMean*m_dMean;
	m_dKurtosis = dX4/dSize - 4.0*m_dMean*m_dSkewness - 6.0*m_dMean*m_dMean*m_dVariance*m_dVariance - m_dMean*m_dMean*m_dMean*m_dMean;
	m_dMeanHeight = m_dMeanHeight/dSize;
	m_dHeight = m_dMaxPeak - m_dMinPit;
	
	if(sFileName.empty())
	{
		printf("%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",m_dMean,m_dVariance,m_dSkewness,m_dKurtosis,m_dMinPit,m_dMaxPeak,m_dHeight,m_dMeanHeight);
	}
	else
	{
		FILE* fpFile = fopen(sFileName.c_str(),"a");
		fprintf(fpFile,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",m_dMean,m_dVariance,m_dSkewness,m_dKurtosis,m_dMinPit,m_dMaxPeak,m_dHeight,m_dMeanHeight);
		fclose(fpFile);
	}
}
void SurfaceProfile::GetSpatialCharacteristics(const string& sFileName)
{
	if(m_poData == NULL)		return;
	
	FourierTransformable* poAutoCorrelation = GetAutoCorrelation();
	((FourierTransformable2D*)poAutoCorrelation)->Shift();
	
	double dRMax = 0.0;
	double dRMin = sqrt(m_iXSize*m_iXSize + m_iYSize*m_iYSize);
	int i = 0;
	int j = 0;
	double dR = 0.0;
	const ComplexNumber** poData = ((FourierTransformable2D*)poAutoCorrelation)->GetData();
	// get the maximum and minimum autocorrelation values
	double dMin = DBL_MAX;
	double dMax = -DBL_MAX;
	double dTemp = 0.0;
	for(i = 0 ; i < m_iXSize ; i++)
	{
		for(j = 0 ; j < m_iYSize ; j++)
		{
			dTemp = poData[i][j].GetReal();
			if(dTemp < dMin)	dMin = dTemp;
			if(dTemp > dMax)	dMax = dTemp;
		}
	}
	// the target autocorrelation is at 20% of the autocorrelation range
	double dTargetAutoCorrelation = 0.2*(dMax - dMin) + dMin;
	for(i = 0 ; i < m_iXSize ; i++)
	{
		for(j = 0 ; j < m_iYSize ; j++)
		{
			if(poData[i][j].GetReal() < dTargetAutoCorrelation)
			{
				dR = sqrt(i*i + j*j);
				if(dR < dRMin)
				{
					dRMin = dR;
				}
				if(dR > dRMax)
				{
					dRMax = dR;
				}
			}
		}
	}
	// if the autocorrelation never reaches the target value, then the ratio is 1
	if(dRMax < dRMin)
	{
		dRMax = dRMin;
	}
	// get the angular spectral power density function
	const ComplexNumber** poFFT = ((FourierTransformable2D*)poAutoCorrelation)->GetTransform();
	double dAngle = 0.0;
	int iXCenter = m_iXSize/2;
	int iYCenter = m_iYSize/2;
	double dMaxPower = 0.0;
	double dMaxAngle = 0.0;
	for(i = 0 ; i < m_iXSize ; i++)
	{
		for(j = 0 ; j < m_iYSize ; j++)
		{
			// exclude the dc component
			if((i == 0) && (j == 0))
			{
				continue;
			}
			dAngle = atan2((double)(j - iYCenter),(double)(i - iXCenter));
			dTemp = poFFT[i][j].GetAmplitude();
			if(dTemp > dMaxPower)
			{
				dMaxPower = dTemp;
				dMaxAngle = dAngle;
			}
		}
	}

	delete poAutoCorrelation;
	
	if(sFileName.empty())
	{
		printf("%e\t\t%e\t\t%e\n",dRMax,dRMin/dRMax,dMaxAngle);
	}
	else
	{
		FILE* fpFile = fopen(sFileName.c_str(),"a");
		fprintf(fpFile,"%e\t\t%e\t\t%e\n",dRMax,dRMin/dRMax,dMaxAngle);
		fclose(fpFile);
	}
}
void SurfaceProfile::GetHybridCharacteristics(const string& sFileName)
{

}
void SurfaceProfile::SetXBounds(const double& dMin,const double& dMax)
{
	m_dXMin = dMin;
	m_dXMax = dMax;
}
void SurfaceProfile::SetYBounds(const double& dMin,const double& dMax)
{
	m_dYMin = dMin;
	m_dYMax = dMax;
}
double SurfaceProfile::GetValue(const double& dX,const double& dY) const
{
	double dTolerance = 1.0E-3;
	if((dX < (m_dXMin - dTolerance)) || (dX > (m_dXMax + dTolerance)))		return 0.0;
	if((dY < (m_dYMin - dTolerance)) || (dY > (m_dYMax + dTolerance)))		return 0.0;
	double dXStep = (m_dXMax - m_dXMin)/(double)(m_iXSize - 1);
	double dYStep = (m_dYMax - m_dYMin)/(double)(m_iXSize - 1);
	// get the space local point coordinates
	double dXi = (dX - m_dXMin)/(m_dXMax - m_dXMin);
	double dEta = (dY - m_dYMin)/(m_dYMax - m_dYMin);
	unsigned int iXMin = (unsigned int)floor(dXi*(m_iXSize - 1));
	unsigned int iXMax = (unsigned int)ceil(dXi*(m_iXSize - 1));
	unsigned int iYMin = (unsigned int)floor(dEta*(m_iYSize - 1));
	unsigned int iYMax = (unsigned int)ceil(dEta*(m_iYSize - 1));
	if(iXMin >= m_iXSize)			iXMin = m_iXSize - 1;
	if(iXMax >= m_iXSize)			iXMax = m_iXSize - 1;
	if(iYMin >= m_iYSize)			iYMin = m_iYSize - 1;
	if(iYMax >= m_iYSize)			iYMax = m_iYSize - 1;
	// get the cell local point coordinates
	dXi = (dX - dXStep*iXMin - m_dXMin)/dXStep;
	dEta = (dY - dYStep*iYMin - m_dYMin)/dYStep;
	// do a double interpolation
	double dLowMinXVal = m_poData[iXMin][iYMin].GetReal();
	double dLowMaxXVal = m_poData[iXMax][iYMin].GetReal();
	double dHighMinXVal = m_poData[iXMin][iYMax].GetReal();
	double dHighMaxXVal = m_poData[iXMax][iYMax].GetReal();
	double dMinYVal = (1.0 - dXi)*dLowMinXVal + dXi*dLowMaxXVal;
	double dMaxYVal = (1.0 - dXi)*dHighMinXVal + dXi*dHighMaxXVal;
	double dValue = (1.0 - dEta)*dMinYVal + dEta*dMaxYVal;
	return dValue;
}
double SurfaceProfile::GetHausdorffArea(const unsigned int& iXResolution,const unsigned int& iYResolution,double& dXScale,double& dYScale) const
{
	unsigned int i = 0;
	unsigned int j = 0;
	dXScale = (m_dXMax - m_dXMin)/(double)(iXResolution);
	dYScale = (m_dYMax - m_dYMin)/(double)(iYResolution);
	double dX = 0.0;
	double dY = 0.0;
	Vector oV1;
	Vector oV2;
	double dArea1 = 0.0;
	double dArea2 = 0.0;
	double dTotalArea = 0.0;
	double dZ1 = 0.0;
	double dZ2 = 0.0;
	double dZ3 = 0.0;
	double dZ4 = 0.0;
	for(i = 0 ; i < iXResolution ; i++)
	{
		dX = m_dXMin + i*dXScale;
		for(j = 0 ; j < iYResolution ; j++)
		{
			dY = m_dYMin + j*dYScale;
			//printf("%lf,%lf\n",dX,dY);fflush(NULL);
			// get the values at the 4 corners
			dZ1 = GetValue(dX,dY);
			dZ2 = GetValue(dX + dXScale,dY);
			dZ3 = GetValue(dX + dXScale,dY + dYScale);
			dZ4 = GetValue(dX,dY + dYScale);
			// get the area of the first triangle
			oV1.Set(dXScale,0.0,dZ2 - dZ1);
			oV2.Set(0.0,dYScale,dZ4 - dZ1);
			dArea1 = (oV1^oV2).Length();
			// get the area of the second triangle
			oV1.Set(0.0,-dYScale,dZ2 - dZ3);
			oV2.Set(-dXScale,0.0,dZ4 - dZ3);
			dArea2 = (oV1^oV2).Length();
			dTotalArea = dTotalArea + 0.5*(dArea1 + dArea2);
		}
	}
	return dTotalArea;
}
double SurfaceProfile::GetHausdorffDimension() const
{
	double dHausdorffArea = 0.0;
	double dXScale = 0.0;
	double dYScale = 0.0;
	double dScale = 0.0;
	double dSumX = 0.0;
	double dSumY = 0.0;
	double dSumX2 = 0.0;
	double dSumXY = 0.0;
	
	unsigned int iResolutionsCount = 10;
	unsigned int iaResolutions[10] = {25,50,100,500,600,700,800,900,1000,2000};
	unsigned int i = 0;
	unsigned int iResolution;
	for(i = 0 ; i < iResolutionsCount ; i++)
	{
		iResolution = iaResolutions[i];
		dHausdorffArea = GetHausdorffArea(iResolution,iResolution,dXScale,dYScale);
		dScale = 0.5*(dXScale + dYScale);
		dSumX += log(dScale);
		dSumY += log(dHausdorffArea);
		dSumX2 += log(dScale)*log(dScale);
		dSumXY += log(dScale)*log(dHausdorffArea);
	}
	// get a power law fit for the data
	double dBeta = (iResolutionsCount*dSumXY - dSumX*dSumY)/(iResolutionsCount*dSumX2 - dSumX*dSumX);
	double dAlpha = (dSumY - dBeta*dSumX)/iResolutionsCount;
	double dLambda = exp(dAlpha);
	double dN = dBeta;
	double dHausdorffDimension = -(dN - 2);
	return dHausdorffDimension;
}
void SurfaceProfile::Initialize()
{
	FourierTransformable2D::Initialize();
	m_dMean = 0.0;
	m_dVariance = 1.0;
	m_dSkewness = 0.0;
	m_dKurtosis = 0.0;
	m_dMinPit = 0.0;
	m_dMaxPeak = 0.0;
	m_dMeanHeight = 0.0;
	m_dHeight = 0.0;
	m_dXMin = 0.0;
	m_dXMax = 0.0;
	m_dYMin = 0.0;
	m_dYMax = 0.0;
}


