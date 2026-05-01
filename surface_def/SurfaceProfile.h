#ifndef SURFACEPROFILE_H_
#define SURFACEPROFILE_H_

#include "FourierTransformable2D.h"

using namespace EZ;

class SurfaceProfile : public FourierTransformable2D
{
public:
	SurfaceProfile();
	SurfaceProfile(const SurfaceProfile& oProfile);
	~SurfaceProfile();
	SurfaceProfile& operator=(const SurfaceProfile& oProfile);
	void Reset();
	virtual FourierTransformable* GetAutoCorrelation();
	void GetAmplitudeCharacteristics(const string& sFileName = "");
	void GetSpatialCharacteristics(const string& sFileName = "");
	void GetHybridCharacteristics(const string& sFileName = "");
	void SetXBounds(const double& dMin,const double& dMax);
	void SetYBounds(const double& dMin,const double& dMax);
	double GetValue(const double& dX,const double& dY) const;
	double GetHausdorffArea(const unsigned int& iXResolution,const unsigned int& iYResolution,double& dXScale,double& dYScale) const;
	double GetHausdorffDimension() const;
	
private:

protected:
	void Initialize();
	double m_dMean;
	double m_dVariance;
	double m_dSkewness;
	double m_dKurtosis;
	double m_dMinPit;
	double m_dMaxPeak;
	double m_dHeight;
	double m_dMeanHeight;
	double m_dXMin;
	double m_dXMax;
	double m_dYMin;
	double m_dYMax;
};


#endif


