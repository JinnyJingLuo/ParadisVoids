// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "AxisAlignedBoundingBox.h"
#include "Tools.h"
#include "Randomizer.h"
#include "math.h"

using namespace SupportSystem;

namespace GeometrySystem
{
	AxisAlignedBoundingBox::AxisAlignedBoundingBox()
	{
		Initialize();
	}
	AxisAlignedBoundingBox::AxisAlignedBoundingBox(const AxisAlignedBoundingBox& oBox)
	{
		*this = oBox;
	}
	AxisAlignedBoundingBox::~AxisAlignedBoundingBox()
	{
		
	}
	void AxisAlignedBoundingBox::Reset()
	{
		Initialize();
	}
	AxisAlignedBoundingBox& AxisAlignedBoundingBox::operator=(const AxisAlignedBoundingBox& oBox)
	{
		m_dXMin = oBox.m_dXMin;
		m_dXMax = oBox.m_dXMax;
		m_dYMin = oBox.m_dYMin;
		m_dYMax = oBox.m_dYMax;
		m_dZMin = oBox.m_dZMin;
		m_dZMax = oBox.m_dZMax;
		m_oSystem = oBox.m_oSystem;
		return *this;
	}
	bool AxisAlignedBoundingBox::IsPointInside(const Point& oPoint,const double& dToleranceFactor) const
	{
		Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
		double dX = oLocalPoint.GetX();
		double dY = oLocalPoint.GetY();
		double dZ = oLocalPoint.GetZ();
		double dActualTolerance = dToleranceFactor*GetMinimumDimension();
		if((dX >= (m_dXMin - dActualTolerance)) && (dX <= (m_dXMax + dActualTolerance)))
		{
			if((dY >= (m_dYMin - dActualTolerance)) && (dY <= (m_dYMax + dActualTolerance)))
			{
				if((dZ >= (m_dZMin - dActualTolerance)) && (dZ <= (m_dZMax + dActualTolerance)))
				{
					return true;
				}
			}
		}
		return false;
	}
	void AxisAlignedBoundingBox::Initialize()
	{
		m_dXMin = 0.0;
		m_dXMax = 0.0;
		m_dYMin = 0.0;
		m_dYMax = 0.0;
		m_dZMin = 0.0;
		m_dZMax = 0.0;
		m_oSystem.Reset();
	}
	double AxisAlignedBoundingBox::GetMinimumDimension() const
	{
		return min((m_dXMax - m_dXMin),min((m_dYMax - m_dYMin),(m_dZMax - m_dZMin)));
	}
	double AxisAlignedBoundingBox::GetMaximumDimension() const
	{
		return max((m_dXMax - m_dXMin),max((m_dYMax - m_dYMin),(m_dZMax - m_dZMin)));
	}
	double AxisAlignedBoundingBox::GetVolume() const
	{
		return ((m_dXMax - m_dXMin)*(m_dYMax - m_dYMin)*(m_dZMax - m_dZMin));
	}
	void AxisAlignedBoundingBox::SetXMin(const double& dValue)
	{
		m_dXMin = dValue;
	}
	void AxisAlignedBoundingBox::SetXMax(const double& dValue)
	{
		m_dXMax = dValue;
	}
	void AxisAlignedBoundingBox::SetYMin(const double& dValue)
	{
		m_dYMin = dValue;
	}
	void AxisAlignedBoundingBox::SetYMax(const double& dValue)
	{
		m_dYMax = dValue;
	}
	void AxisAlignedBoundingBox::SetZMin(const double& dValue)
	{
		m_dZMin = dValue;
	}
	void AxisAlignedBoundingBox::SetZMax(const double& dValue)
	{
		m_dZMax = dValue;
	}
	double AxisAlignedBoundingBox::GetXMin() const
	{
		return m_dXMin;
	}
	double AxisAlignedBoundingBox::GetXMax() const
	{
		return m_dXMax;
	}
	double AxisAlignedBoundingBox::GetYMin() const
	{
		return m_dYMin;
	}
	double AxisAlignedBoundingBox::GetYMax() const
	{
		return m_dYMax;
	}
	double AxisAlignedBoundingBox::GetZMin() const
	{
		return m_dZMin;
	}
	double AxisAlignedBoundingBox::GetZMax() const
	{
		return m_dZMax;
	}
	list<AxisAlignedBoundingBox*> AxisAlignedBoundingBox::UniformPartition(const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount) const
	{
		list<AxisAlignedBoundingBox*> lpoOldPartitions;
		list<AxisAlignedBoundingBox*> lpoNewPartitions;
		lpoOldPartitions.clear();
		lpoNewPartitions.clear();
		UniformXSlice(iXSlicesCount,&lpoNewPartitions);
		lpoOldPartitions = lpoNewPartitions;
		lpoNewPartitions.clear();
		list<AxisAlignedBoundingBox*>::iterator liPartitions;
		// Y direction partitioning
		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
		{
			(*liPartitions)->UniformYSlice(iYSlicesCount,&lpoNewPartitions);
			delete *liPartitions;
		}
		lpoOldPartitions.clear();
		lpoOldPartitions = lpoNewPartitions;
		lpoNewPartitions.clear();
		// Z direction partitioning
		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
		{
			(*liPartitions)->UniformZSlice(iZSlicesCount,&lpoNewPartitions);
			delete *liPartitions;
		}
		lpoOldPartitions.clear();
		return lpoNewPartitions;
	}
	list<AxisAlignedBoundingBox*> AxisAlignedBoundingBox::NonUniformPartition(const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount) const
	{
		list<AxisAlignedBoundingBox*> lpoOldPartitions;
		list<AxisAlignedBoundingBox*> lpoNewPartitions;
		lpoOldPartitions.clear();
		lpoNewPartitions.clear();
		NonUniformXSlice(iXSlicesCount,&lpoNewPartitions);
		lpoOldPartitions = lpoNewPartitions;
		lpoNewPartitions.clear();
		list<AxisAlignedBoundingBox*>::iterator liPartitions;
		// Y direction partitioning
		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
		{
			(*liPartitions)->NonUniformYSlice(iYSlicesCount,&lpoNewPartitions);
			delete *liPartitions;
		}
		lpoOldPartitions.clear();
		lpoOldPartitions = lpoNewPartitions;
		lpoNewPartitions.clear();
		// Z direction partitioning
		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
		{
			(*liPartitions)->NonUniformZSlice(iZSlicesCount,&lpoNewPartitions);
			delete *liPartitions;
		}
		lpoOldPartitions.clear();
		return lpoNewPartitions;
	}
	list<AxisAlignedBoundingBox*> AxisAlignedBoundingBox::PointBasedPartition(list<Point>* ploPoints,const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount) const
	{
		list<AxisAlignedBoundingBox*> lpoOldPartitions;
		list<AxisAlignedBoundingBox*> lpoNewPartitions;
		lpoOldPartitions.clear();
		lpoNewPartitions.clear();
		PointBasedXSlice(ploPoints,iXSlicesCount,&lpoNewPartitions);
		lpoOldPartitions = lpoNewPartitions;
		lpoNewPartitions.clear();
		list<AxisAlignedBoundingBox*>::iterator liPartitions;
		// Y direction partitioning
		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
		{
			(*liPartitions)->PointBasedYSlice(ploPoints,iYSlicesCount,&lpoNewPartitions);
			delete *liPartitions;
		}
		lpoOldPartitions.clear();
		lpoOldPartitions = lpoNewPartitions;
		lpoNewPartitions.clear();
		// Z direction partitioning
		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
		{
			(*liPartitions)->PointBasedZSlice(ploPoints,iZSlicesCount,&lpoNewPartitions);
			delete *liPartitions;
		}
		lpoOldPartitions.clear();
		return lpoNewPartitions;
	}
	void AxisAlignedBoundingBox::UniformXSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		double dTotalLength = m_dXMax - m_dXMin;
		double dPartitionLength = dTotalLength/(double)iSlicesCount;
		unsigned int i = 0;
		AxisAlignedBoundingBox* poBox = NULL;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new AxisAlignedBoundingBox;
			poBox->SetXMin(m_dXMin + i*dPartitionLength);
			poBox->SetXMax(m_dXMin + (i + 1)*dPartitionLength);
			poBox->SetYMin(m_dYMin);
			poBox->SetYMax(m_dYMax);
			poBox->SetZMin(m_dZMin);
			poBox->SetZMax(m_dZMax);
			poBox->SetSystem(m_oSystem);
			plpoBoxes->push_back(poBox);
		}
	}
	void AxisAlignedBoundingBox::UniformYSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		double dTotalLength = m_dYMax - m_dYMin;
		double dPartitionLength = dTotalLength/(double)iSlicesCount;
		unsigned int i = 0;
		AxisAlignedBoundingBox* poBox = NULL;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new AxisAlignedBoundingBox;
			poBox->SetXMin(m_dXMin);
			poBox->SetXMax(m_dXMax);
			poBox->SetYMin(m_dYMin + i*dPartitionLength);
			poBox->SetYMax(m_dYMin + (i + 1)*dPartitionLength);
			poBox->SetZMin(m_dZMin);
			poBox->SetZMax(m_dZMax);
			poBox->SetSystem(m_oSystem);
			plpoBoxes->push_back(poBox);
		}
	}
	void AxisAlignedBoundingBox::UniformZSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		double dTotalLength = m_dZMax - m_dZMin;
		double dPartitionLength = dTotalLength/(double)iSlicesCount;
		unsigned int i = 0;
		AxisAlignedBoundingBox* poBox = NULL;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new AxisAlignedBoundingBox;
			poBox->SetXMin(m_dXMin);
			poBox->SetXMax(m_dXMax);
			poBox->SetYMin(m_dYMin);
			poBox->SetYMax(m_dYMax);
			poBox->SetZMin(m_dZMin + i*dPartitionLength);
			poBox->SetZMax(m_dZMin + (i + 1)*dPartitionLength);
			poBox->SetSystem(m_oSystem);
			plpoBoxes->push_back(poBox);
		}
	}
	void AxisAlignedBoundingBox::NonUniformXSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		double dTotalLength = m_dXMax - m_dXMin;
		// get partition length mean and standard deviation
		double dMu = dTotalLength/(double)iSlicesCount;
		double dSigma = 0.3*dMu;
		unsigned int i = 0;
		AxisAlignedBoundingBox* poBox = NULL;
		double dLower = m_dXMin;
		double dUpper = 0.0;
		for(i = 0; i < iSlicesCount ; i++)
		{
			if(i == (iSlicesCount - 1))
			{
				dUpper = m_dXMax;
			}
			else
			{
				while(true)
				{
					dUpper = dLower + Randomizer::RandomNormal(dMu,dSigma);
					if((dUpper > dLower) && (dUpper < m_dXMax))
					{
						break;
					}
				}
			}
			poBox = new AxisAlignedBoundingBox;
			poBox->SetXMin(dLower);
			poBox->SetXMax(dUpper);
			poBox->SetYMin(m_dYMin);
			poBox->SetYMax(m_dYMax);
			poBox->SetZMin(m_dZMin);
			poBox->SetZMax(m_dZMax);
			poBox->SetSystem(m_oSystem);
			plpoBoxes->push_back(poBox);			
			dLower = dUpper;
		}
	}
	void AxisAlignedBoundingBox::NonUniformYSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		double dTotalLength = m_dYMax - m_dYMin;
		// get partition length mean and standard deviation
		double dMu = dTotalLength/(double)iSlicesCount;
		double dSigma = 0.3*dMu;
		unsigned int i = 0;
		AxisAlignedBoundingBox* poBox = NULL;
		double dLower = m_dYMin;
		double dUpper = 0.0;
		for(i = 0; i < iSlicesCount ; i++)
		{
			if(i == (iSlicesCount - 1))
			{
				dUpper = m_dYMax;
			}
			else
			{
				while(true)
				{
					dUpper = dLower + Randomizer::RandomNormal(dMu,dSigma);
					if((dUpper > dLower) && (dUpper < m_dYMax))
					{
						break;
					}
				}
			}
			poBox = new AxisAlignedBoundingBox;
			poBox->SetXMin(m_dXMin);
			poBox->SetXMax(m_dXMax);
			poBox->SetYMin(dLower);
			poBox->SetYMax(dUpper);
			poBox->SetZMin(m_dZMin);
			poBox->SetZMax(m_dZMax);
			poBox->SetSystem(m_oSystem);
			plpoBoxes->push_back(poBox);			
			dLower = dUpper;
		}
	}
	void AxisAlignedBoundingBox::NonUniformZSlice(const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		double dTotalLength = m_dZMax - m_dZMin;
		// get partition length mean and standard deviation
		double dMu = dTotalLength/(double)iSlicesCount;
		double dSigma = 0.3*dMu;
		unsigned int i = 0;
		AxisAlignedBoundingBox* poBox = NULL;
		double dLower = m_dZMin;
		double dUpper = 0.0;
		for(i = 0; i < iSlicesCount ; i++)
		{
			if(i == (iSlicesCount - 1))
			{
				dUpper = m_dZMax;
			}
			else
			{
				while(true)
				{
					dUpper = dLower + Randomizer::RandomNormal(dMu,dSigma);
					if((dUpper > dLower) && (dUpper < m_dZMax))
					{
						break;
					}
				}
			}
			poBox = new AxisAlignedBoundingBox;
			poBox->SetXMin(m_dXMin);
			poBox->SetXMax(m_dXMax);
			poBox->SetYMin(m_dYMin);
			poBox->SetYMax(m_dYMax);
			poBox->SetZMin(dLower);
			poBox->SetZMax(dUpper);
			poBox->SetSystem(m_oSystem);
			plpoBoxes->push_back(poBox);			
			dLower = dUpper;
		}
	}
	void AxisAlignedBoundingBox::PointBasedXSlice(list<Point>* ploPoints,const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		
		list<Point> loWorkingPoints;
		list<Point>::iterator liPoints;
		for(liPoints = ploPoints->begin() ; liPoints != ploPoints->end() ; liPoints++)
		{
			if(IsPointInside((*liPoints)))
			{
				loWorkingPoints.push_back((*liPoints));
			}
		}
		unsigned int iSize = (unsigned int)loWorkingPoints.size();
		
		unsigned int iPointsPerSlice = iSize/iSlicesCount;
		unsigned int iRemainingPoints = iSize%iSlicesCount;
		unsigned int i = 0;
		vector<unsigned int> viSliceLastPointIndex;
		viSliceLastPointIndex.resize(iSlicesCount);

		if(iRemainingPoints > 0)
		{
			viSliceLastPointIndex[0] = iPointsPerSlice + 1;
			iRemainingPoints = iRemainingPoints - 1;
		}
		else
		{
			viSliceLastPointIndex[0] = iPointsPerSlice;
		}
		for(i = 1 ; i < iSlicesCount - 1 ; i++)
		{
			if(iRemainingPoints > 0)
			{
				viSliceLastPointIndex[i] = viSliceLastPointIndex[i - 1] + iPointsPerSlice + 1;
				iRemainingPoints = iRemainingPoints - 1;
			}
			else
			{
				viSliceLastPointIndex[i] = viSliceLastPointIndex[i - 1] + iPointsPerSlice;
			}
		}
		viSliceLastPointIndex[iSlicesCount - 1] = iSize - 1;
		
		vector<double> vdValues;
		vdValues.resize(iSize);
		i = 0;
		for(liPoints = loWorkingPoints.begin() ; liPoints != loWorkingPoints.end() ; liPoints++)
		{
			vdValues[i] = (*liPoints).GetX();
			i = i + 1;
		}
		QuickSort(vdValues);
		AxisAlignedBoundingBox* poBox = NULL;
		double dLowerBound = m_dXMin;
		double dUpperBound = 0.0;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new AxisAlignedBoundingBox;
			if(i == iSlicesCount - 1)
			{
				dUpperBound = m_dXMax;
			}
			else
			{
				dUpperBound = 0.5*(vdValues[viSliceLastPointIndex[i]] + vdValues[viSliceLastPointIndex[i] + 1]);
			}
			poBox->SetXMin(dLowerBound);
			poBox->SetXMax(dUpperBound);
			poBox->SetYMin(m_dYMin);
			poBox->SetYMax(m_dYMax);
			poBox->SetZMin(m_dZMin);
			poBox->SetZMax(m_dZMax);
			plpoBoxes->push_back(poBox);
			dLowerBound = dUpperBound;
		}
	}
	void AxisAlignedBoundingBox::PointBasedYSlice(list<Point>* ploPoints,const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		
		list<Point> loWorkingPoints;
		list<Point>::iterator liPoints;
		for(liPoints = ploPoints->begin() ; liPoints != ploPoints->end() ; liPoints++)
		{
			if(IsPointInside((*liPoints)))
			{
				loWorkingPoints.push_back((*liPoints));
			}
		}
		unsigned int iSize = (unsigned int)loWorkingPoints.size();
		unsigned int iPointsPerSlice = iSize/iSlicesCount;
		unsigned int iRemainingPoints = iSize%iSlicesCount;
		unsigned int i = 0;
		vector<unsigned int> viSliceLastPointIndex;
		viSliceLastPointIndex.resize(iSlicesCount);
		
		if(iRemainingPoints > 0)
		{
			viSliceLastPointIndex[0] = iPointsPerSlice + 1;
			iRemainingPoints = iRemainingPoints - 1;
		}
		else
		{
			viSliceLastPointIndex[0] = iPointsPerSlice;
		}
		for(i = 1 ; i < iSlicesCount - 1 ; i++)
		{
			if(iRemainingPoints > 0)
			{
				viSliceLastPointIndex[i] = viSliceLastPointIndex[i - 1] + iPointsPerSlice + 1;
				iRemainingPoints = iRemainingPoints - 1;
			}
			else
			{
				viSliceLastPointIndex[i] = viSliceLastPointIndex[i - 1] + iPointsPerSlice;
			}
		}
		viSliceLastPointIndex[iSlicesCount - 1] = iSize - 1;
		
		vector<double> vdValues;
		vdValues.resize(iSize);
		i = 0;
		for(liPoints = loWorkingPoints.begin() ; liPoints != loWorkingPoints.end() ; liPoints++)
		{
			vdValues[i] = (*liPoints).GetY();
			i = i + 1;
		}
		QuickSort(vdValues);
		AxisAlignedBoundingBox* poBox = NULL;
		double dLowerBound = m_dYMin;
		double dUpperBound = 0.0;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new AxisAlignedBoundingBox;
			if(i == iSlicesCount - 1)
			{
				dUpperBound = m_dYMax;
			}
			else
			{
				dUpperBound = 0.5*(vdValues[viSliceLastPointIndex[i]] + vdValues[viSliceLastPointIndex[i] + 1]);
			}
			poBox->SetXMin(m_dXMin);
			poBox->SetXMax(m_dXMax);
			poBox->SetYMin(dLowerBound);
			poBox->SetYMax(dUpperBound);
			poBox->SetZMin(m_dZMin);
			poBox->SetZMax(m_dZMax);
			plpoBoxes->push_back(poBox);
			dLowerBound = dUpperBound;
		}
	}
	void AxisAlignedBoundingBox::PointBasedZSlice(list<Point>* ploPoints,const unsigned int& iSlicesCount,list<AxisAlignedBoundingBox*>* plpoBoxes) const
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		
		list<Point> loWorkingPoints;
		list<Point>::iterator liPoints;
		for(liPoints = ploPoints->begin() ; liPoints != ploPoints->end() ; liPoints++)
		{
			if(IsPointInside((*liPoints)))
			{
				loWorkingPoints.push_back((*liPoints));
			}
		}
		unsigned int iSize = (unsigned int)loWorkingPoints.size();
		unsigned int iPointsPerSlice = iSize/iSlicesCount;
		unsigned int iRemainingPoints = iSize%iSlicesCount;
		unsigned int i = 0;
		vector<unsigned int> viSliceLastPointIndex;
		viSliceLastPointIndex.resize(iSlicesCount);
		
		if(iRemainingPoints > 0)
		{
			viSliceLastPointIndex[0] = iPointsPerSlice + 1;
			iRemainingPoints = iRemainingPoints - 1;
		}
		else
		{
			viSliceLastPointIndex[0] = iPointsPerSlice;
		}
		for(i = 1 ; i < iSlicesCount - 1 ; i++)
		{
			if(iRemainingPoints > 0)
			{
				viSliceLastPointIndex[i] = viSliceLastPointIndex[i - 1] + iPointsPerSlice + 1;
				iRemainingPoints = iRemainingPoints - 1;
			}
			else
			{
				viSliceLastPointIndex[i] = viSliceLastPointIndex[i - 1] + iPointsPerSlice;
			}
		}
		viSliceLastPointIndex[iSlicesCount - 1] = iSize - 1;
		
		vector<double> vdValues;
		vdValues.resize(iSize);
		i = 0;
		for(liPoints = loWorkingPoints.begin() ; liPoints != loWorkingPoints.end() ; liPoints++)
		{
			vdValues[i] = (*liPoints).GetZ();
			i = i + 1;
		}
		QuickSort(vdValues);
		AxisAlignedBoundingBox* poBox = NULL;
		double dLowerBound = m_dZMin;
		double dUpperBound = 0.0;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new AxisAlignedBoundingBox;
			if(i == iSlicesCount - 1)
			{
				dUpperBound = m_dZMax;
			}
			else
			{
				dUpperBound = 0.5*(vdValues[viSliceLastPointIndex[i]] + vdValues[viSliceLastPointIndex[i] + 1]);
			}
			poBox->SetXMin(m_dXMin);
			poBox->SetXMax(m_dXMax);
			poBox->SetYMin(m_dYMin);
			poBox->SetYMax(m_dYMax);
			poBox->SetZMin(dLowerBound);
			poBox->SetZMax(dUpperBound);
			plpoBoxes->push_back(poBox);
			dLowerBound = dUpperBound;
		}
	}
	Point AxisAlignedBoundingBox::GenerateRandomPoint(const double& dMargin) const
	{
 		double dX = Randomizer::Random(m_dXMin + dMargin,m_dXMax - dMargin);
 		double dY = Randomizer::Random(m_dYMin + dMargin,m_dYMax - dMargin);
 		double dZ = Randomizer::Random(m_dZMin + dMargin,m_dZMax - dMargin);
 		return m_oSystem.GetInGlobalCoordinates(Point(dX,dY,dZ));
 	}
 	Point AxisAlignedBoundingBox::GenerateRandomExternalPoint(const double& dSpacingFactor) const
 	{
		double dRadius = dSpacingFactor*GetMaximumDimension();
		double dTheta = Randomizer::Random(0,2*PI);
		double dPhi = Randomizer::Random(0,PI);	
		double dX = 0.5*(m_dXMax + m_dXMin) + dRadius*sin(dPhi)*cos(dTheta);
		double dY = 0.5*(m_dYMax + m_dYMin) + dRadius*sin(dPhi)*sin(dTheta);
		double dZ = 0.5*(m_dZMax + m_dZMin) + dRadius*cos(dPhi);
		return m_oSystem.GetInGlobalCoordinates(Point(dX,dY,dZ));
 	}
 	void AxisAlignedBoundingBox::CenterAt(const Point& oPoint)
 	{
 		double dTolerance = 1.0E-10;
 		m_dXMin = oPoint.GetX() - dTolerance;
 		m_dXMax = oPoint.GetX() + dTolerance;
 		m_dYMin = oPoint.GetY() - dTolerance;
 		m_dYMax = oPoint.GetY() + dTolerance;
 		m_dZMin = oPoint.GetZ() - dTolerance;
 		m_dZMax = oPoint.GetZ() + dTolerance;
 	}
 	void AxisAlignedBoundingBox::ExpandToContain(const Point& oPoint)
 	{
 		double dX = oPoint.GetX();
 		double dY = oPoint.GetY();
 		double dZ = oPoint.GetZ();
 		if(dX > m_dXMax)
 		{
 			m_dXMax = dX;
 		}
 		if(dX < m_dXMin)
 		{
 			m_dXMin = dX;
 		}
 		if(dY > m_dYMax)
 		{
 			m_dYMax = dY;
 		}
 		if(dY < m_dYMin)
 		{
 			m_dYMin = dY;
 		}
 		if(dZ > m_dZMax)
 		{
 			m_dZMax = dZ;
 		}
 		if(dZ < m_dZMin)
 		{
 			m_dZMin = dZ;
 		}
 	}
 	int AxisAlignedBoundingBox::GetMaximumDimensionIndex() const
 	{
 	 	double dX = m_dXMax - m_dXMin;
 		double dY = m_dYMax - m_dYMin;
 		double dZ = m_dZMax - m_dZMin;
 		if(dX > dY)
 		{
 			if(dX > dZ)
 			{
 				return 1;
 			}
 			else
 			{
 				return 3;
 			}
 		}
 		else
 		{
 			if(dY > dZ)
 			{
 				return 2;
 			}
 		}
 		return 3;
 	}
 	Point AxisAlignedBoundingBox::GetCenter() const
 	{
 		return m_oSystem.GetInGlobalCoordinates(Point(0.5*(m_dXMin + m_dXMax),0.5*(m_dYMin + m_dYMax),0.5*(m_dZMin + m_dZMax)));
 	}
 	Geometry* AxisAlignedBoundingBox::Clone()
 	{
		return new AxisAlignedBoundingBox(*this);
	}
	string AxisAlignedBoundingBox::ToString() const
	{
		char cString[256];
		sprintf(cString,"%e -> %e,%e -> %e,%e -> %e",m_dXMin,m_dXMax,m_dYMin,m_dYMax,m_dZMin,m_dZMax);
		return string(cString);
	}
	Point AxisAlignedBoundingBox::GetRelativePosition(const Point& oPoint) const
	{
		Point oRelativePosition = oPoint;
		oRelativePosition.SetX(oRelativePosition.GetX()/(m_dXMax - m_dXMin));
		oRelativePosition.SetY(oRelativePosition.GetY()/(m_dYMax - m_dYMin));
		oRelativePosition.SetZ(oRelativePosition.GetZ()/(m_dZMax - m_dZMin));
		return oRelativePosition;
	}
	AxisAlignedBoundingBox AxisAlignedBoundingBox::GetSubVolume(const double& dVolumeFraction) const
	{
		double dLengthFactor = 0.5*pow(dVolumeFraction,(1.0/3.0));
		double dXDim = dLengthFactor*(m_dXMax - m_dXMin);
		double dYDim = dLengthFactor*(m_dYMax - m_dYMin);
		double dZDim = dLengthFactor*(m_dZMax - m_dZMin);
		double dXCenter = 0.5*(m_dXMax + m_dXMin);
		double dYCenter = 0.5*(m_dYMax + m_dYMin);
		double dZCenter = 0.5*(m_dZMax + m_dZMin);
		AxisAlignedBoundingBox oBox;
		oBox.SetXMin(dXCenter - dXDim);
		oBox.SetXMax(dXCenter + dXDim);
		oBox.SetYMin(dYCenter - dYDim);
		oBox.SetYMax(dYCenter + dYDim);
		oBox.SetZMin(dZCenter - dZDim);
		oBox.SetZMax(dZCenter + dZDim);
		return oBox;
	}
}


