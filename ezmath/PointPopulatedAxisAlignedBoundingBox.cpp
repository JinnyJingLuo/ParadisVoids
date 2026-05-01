// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "PointPopulatedAxisAlignedBoundingBox.h"
#include "Tools.h"

using namespace SupportSystem;

namespace GeometrySystem
{
	PointPopulatedAxisAlignedBoundingBox::PointPopulatedAxisAlignedBoundingBox()
	{
		Initialize();
	}
	PointPopulatedAxisAlignedBoundingBox::PointPopulatedAxisAlignedBoundingBox(const AxisAlignedBoundingBox& oBox)
	{
		AxisAlignedBoundingBox::operator=(oBox);
	}
	PointPopulatedAxisAlignedBoundingBox::PointPopulatedAxisAlignedBoundingBox(const PointPopulatedAxisAlignedBoundingBox& oBox)
	{
		*this = oBox;
	}
	PointPopulatedAxisAlignedBoundingBox::~PointPopulatedAxisAlignedBoundingBox()
	{
		m_lpoPoints.clear();
	}
	PointPopulatedAxisAlignedBoundingBox& PointPopulatedAxisAlignedBoundingBox::operator=(const PointPopulatedAxisAlignedBoundingBox& oBox)
	{
		AxisAlignedBoundingBox::operator=(oBox);
		m_lpoPoints = oBox.m_lpoPoints;
		return *this;
	}
	list<Point*>* PointPopulatedAxisAlignedBoundingBox::GetPoints()
	{
		return &m_lpoPoints;
	}
	void PointPopulatedAxisAlignedBoundingBox::AddPoint(Point* poPoint)
	{
		m_lpoPoints.push_back(poPoint);
	}
	void PointPopulatedAxisAlignedBoundingBox::SetPoints(list<Point*>* plpoPoints)
	{
		m_lpoPoints = *plpoPoints;
	}
	void PointPopulatedAxisAlignedBoundingBox::Initialize()
	{
		AxisAlignedBoundingBox::Initialize();
		m_lpoPoints.clear();
	}
	AxisAlignedBoundingBox PointPopulatedAxisAlignedBoundingBox::GetBox() const
	{
		AxisAlignedBoundingBox oBox;
		oBox = *this;
		return oBox;
	}
	list<AxisAlignedBoundingBox*> PointPopulatedAxisAlignedBoundingBox::PointBasedPartition(const unsigned int& iXSlicesCount,const unsigned int& iYSlicesCount,const unsigned int& iZSlicesCount)
	{
 		list<PointPopulatedAxisAlignedBoundingBox*> lpoOldPartitions;
 		list<PointPopulatedAxisAlignedBoundingBox*> lpoNewPartitions;
		
		lpoNewPartitions.clear();
 		PointBasedXSlice(iXSlicesCount,&lpoNewPartitions);
 		lpoOldPartitions.clear();
 		lpoOldPartitions = lpoNewPartitions;
  		lpoNewPartitions.clear();
 		
 		list<PointPopulatedAxisAlignedBoundingBox*>::iterator liPartitions;
 		
 		// Y direction partitioning
 		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
 		{
 			(*liPartitions)->PointBasedYSlice(iYSlicesCount,&lpoNewPartitions); 			
 			delete *liPartitions;
 		}
 		lpoOldPartitions.clear();
 		lpoOldPartitions = lpoNewPartitions;
 		lpoNewPartitions.clear();
 		
 		// Z direction partitioning
 		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
 		{
 			(*liPartitions)->PointBasedZSlice(iZSlicesCount,&lpoNewPartitions);
 			delete *liPartitions;
 		}
 		lpoOldPartitions.clear();
 		lpoOldPartitions = lpoNewPartitions;
 		lpoNewPartitions.clear();
 		
 		list<AxisAlignedBoundingBox*> lpoPartitions;
 		AxisAlignedBoundingBox* poPartition = NULL;
 		for(liPartitions = lpoOldPartitions.begin() ; liPartitions != lpoOldPartitions.end() ; liPartitions++)
 		{
 			poPartition = new AxisAlignedBoundingBox;
 			*poPartition = *(*liPartitions);
 			lpoPartitions.push_back(poPartition);
 			delete *liPartitions;
 		}
 		
 		return lpoPartitions;
	}
	void PointPopulatedAxisAlignedBoundingBox::PointBasedXSlice(const unsigned int& iSlicesCount,list<PointPopulatedAxisAlignedBoundingBox*>* plpoBoxes)
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		
		unsigned int iSize = (unsigned int)m_lpoPoints.size();
		if(iSize == 0)
		{
			double dTotalLength = m_dXMax - m_dXMin;
			double dPartitionLength = dTotalLength/(double)iSlicesCount;
			unsigned int i = 0;
			PointPopulatedAxisAlignedBoundingBox* poBox = NULL;
			for(i = 0; i < iSlicesCount ; i++)
			{
				poBox = new PointPopulatedAxisAlignedBoundingBox;
				poBox->SetXMin(m_dXMin + i*dPartitionLength);
				poBox->SetXMax(m_dXMin + (i + 1)*dPartitionLength);
				poBox->SetYMin(m_dYMin);
				poBox->SetYMax(m_dYMax);
				poBox->SetZMin(m_dZMin);
				poBox->SetZMax(m_dZMax);
				plpoBoxes->push_back(poBox);
			}
		}
		unsigned int iPointsPerSlice = iSize/iSlicesCount;
		unsigned int iRemainingPoints = iSize%iSlicesCount;
		unsigned int i = 0;
		vector<unsigned int> viSliceLastPointIndex;
		viSliceLastPointIndex.resize(iSlicesCount);
		
		if(iRemainingPoints > 0)
		{
			viSliceLastPointIndex[0] = iPointsPerSlice;
			iRemainingPoints = iRemainingPoints - 1;
		}
		else
		{
			viSliceLastPointIndex[0] = iPointsPerSlice - 1;
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
		list<Point*>::iterator liPoints;
		i = 0;
		
		list<bool> lbIsPointTaken;
		for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ; liPoints++)
		{
			vdValues[i] = (*liPoints)->GetX();
			lbIsPointTaken.push_back(false);
			i = i + 1;
		}
		QuickSort(vdValues);
		
		double dValue = 0.0;
		PointPopulatedAxisAlignedBoundingBox* poBox = NULL;
		double dLowerBound = m_dXMin;
		double dUpperBound = 0.0;
		
		list<bool>::iterator lbTakenPoints;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new PointPopulatedAxisAlignedBoundingBox;
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
			
			lbTakenPoints = lbIsPointTaken.begin();
			for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ; liPoints++)
			{
				if(!(*lbTakenPoints))
				{
					dValue = (*liPoints)->GetX();
					if(dValue >= dLowerBound && dValue <= dUpperBound)
					{
						poBox->AddPoint(*liPoints);
						(*lbTakenPoints) = true;
					}
				}
				lbTakenPoints++;
			}
			dLowerBound = dUpperBound;
			plpoBoxes->push_back(poBox);
		}
	}
	void PointPopulatedAxisAlignedBoundingBox::PointBasedYSlice(const unsigned int& iSlicesCount,list<PointPopulatedAxisAlignedBoundingBox*>* plpoBoxes)
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		
		unsigned int iSize = (unsigned int)m_lpoPoints.size();
		if(iSize == 0)
		{
			double dTotalLength = m_dYMax - m_dYMin;
			double dPartitionLength = dTotalLength/(double)iSlicesCount;
			unsigned int i = 0;
			PointPopulatedAxisAlignedBoundingBox* poBox = NULL;
			for(i = 0; i < iSlicesCount ; i++)
			{
				poBox = new PointPopulatedAxisAlignedBoundingBox;
				poBox->SetXMin(m_dXMin);
				poBox->SetXMax(m_dXMax);
				poBox->SetYMin(m_dYMin + i*dPartitionLength);
				poBox->SetYMax(m_dYMin + (i + 1)*dPartitionLength);
				poBox->SetZMin(m_dZMin);
				poBox->SetZMax(m_dZMax);
				plpoBoxes->push_back(poBox);
			}
		}
		unsigned int iPointsPerSlice = iSize/iSlicesCount;
		unsigned int iRemainingPoints = iSize%iSlicesCount;
		unsigned int i = 0;
		vector<unsigned int> viSliceLastPointIndex;
		viSliceLastPointIndex.resize(iSlicesCount);
		
		if(iRemainingPoints > 0)
		{
			viSliceLastPointIndex[0] = iPointsPerSlice;
			iRemainingPoints = iRemainingPoints - 1;
		}
		else
		{
			viSliceLastPointIndex[0] = iPointsPerSlice - 1;
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
		list<Point*>::iterator liPoints;
		i = 0;
		
		list<bool> lbIsPointTaken;
		for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ; liPoints++)
		{
			vdValues[i] = (*liPoints)->GetY();
			lbIsPointTaken.push_back(false);
			i = i + 1;
		}
		QuickSort(vdValues);
		
		double dValue = 0.0;
		PointPopulatedAxisAlignedBoundingBox* poBox = NULL;
		double dLowerBound = m_dYMin;
		double dUpperBound = 0.0;
		
		list<bool>::iterator lbTakenPoints;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new PointPopulatedAxisAlignedBoundingBox;
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
			
			lbTakenPoints = lbIsPointTaken.begin();
			for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ; liPoints++)
			{
				if(!(*lbTakenPoints))
				{
					dValue = (*liPoints)->GetY();
					if(dValue >= dLowerBound && dValue <= dUpperBound)
					{
						poBox->AddPoint(*liPoints);
						(*lbTakenPoints) = true;
					}
				}
				lbTakenPoints++;
			}
			dLowerBound = dUpperBound;
			plpoBoxes->push_back(poBox);
		}
	}
	void PointPopulatedAxisAlignedBoundingBox::PointBasedZSlice(const unsigned int& iSlicesCount,list<PointPopulatedAxisAlignedBoundingBox*>* plpoBoxes)
	{
		if(iSlicesCount == 0)
		{
			return;
		}
		
		unsigned int iSize = (unsigned int)m_lpoPoints.size();
		if(iSize == 0)
		{
			double dTotalLength = m_dZMax - m_dZMin;
			double dPartitionLength = dTotalLength/(double)iSlicesCount;
			unsigned int i = 0;
			PointPopulatedAxisAlignedBoundingBox* poBox = NULL;
			for(i = 0; i < iSlicesCount ; i++)
			{
				poBox = new PointPopulatedAxisAlignedBoundingBox;
				poBox->SetXMin(m_dXMin);
				poBox->SetXMax(m_dXMax);
				poBox->SetYMin(m_dYMin);
				poBox->SetYMax(m_dYMax);
				poBox->SetZMin(m_dZMin + i*dPartitionLength);
				poBox->SetZMax(m_dZMin + (i + 1)*dPartitionLength);
				plpoBoxes->push_back(poBox);
			}
		}
		unsigned int iPointsPerSlice = iSize/iSlicesCount;
		unsigned int iRemainingPoints = iSize%iSlicesCount;
		unsigned int i = 0;
		vector<unsigned int> viSliceLastPointIndex;
		viSliceLastPointIndex.resize(iSlicesCount);
		
		if(iRemainingPoints > 0)
		{
			viSliceLastPointIndex[0] = iPointsPerSlice;
			iRemainingPoints = iRemainingPoints - 1;
		}
		else
		{
			viSliceLastPointIndex[0] = iPointsPerSlice - 1;
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
		list<Point*>::iterator liPoints;
		i = 0;

		list<bool> lbIsPointTaken;
		for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ; liPoints++)
		{
			vdValues[i] = (*liPoints)->GetZ();
			lbIsPointTaken.push_back(false);
			i = i + 1;
		}
		QuickSort(vdValues);
		
		double dValue = 0.0;
		PointPopulatedAxisAlignedBoundingBox* poBox = NULL;
		double dLowerBound = m_dZMin;
		double dUpperBound = 0.0;
		
		list<bool>::iterator lbTakenPoints;
		for(i = 0; i < iSlicesCount ; i++)
		{
			poBox = new PointPopulatedAxisAlignedBoundingBox;
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

			lbTakenPoints = lbIsPointTaken.begin();
			for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ; liPoints++)
			{
				if(!(*lbTakenPoints))
				{
					dValue = (*liPoints)->GetZ();
					if(dValue >= dLowerBound && dValue <= dUpperBound)
					{
						poBox->AddPoint(*liPoints);
						(*lbTakenPoints) = true;
					}
				}
				lbTakenPoints++;
			}
			dLowerBound = dUpperBound;
			plpoBoxes->push_back(poBox);
		}
	}
}

