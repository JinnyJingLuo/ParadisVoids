#include "Curve.h"


namespace GeometrySystem
{
	Curve::Curve()
	{
		Initialize();
	}
	Curve::Curve(const Curve& oCurve)
	{
		*this = oCurve;
	}
	Curve::~Curve()
	{
		Reset();
	}
	Curve& Curve::operator=(const Curve& oCurve)
	{
		Reset();
		list<Segment*>::const_iterator liSegments;
		for(liSegments = oCurve.m_lpoSegments.begin() ; liSegments != oCurve.m_lpoSegments.end() ; liSegments++)
		{
			AddSegment(*(*liSegments)->GetPoint1(),*(*liSegments)->GetPoint1());
		}
		return *this;
	}
	void Curve::Reset()
	{
		list<Segment*>::iterator liSegments;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			if((*liSegments) != NULL)			delete (*liSegments);
		}
		m_lpoSegments.clear();
	}
	void Curve::AddSegment(const Point& oStartPoint,const Point& oEndPoint)
	{
		GenericNode* poStartPoint = new GenericNode(oStartPoint);
		GenericNode* poEndPoint = new GenericNode(oEndPoint);
		Segment* poSegment = new Segment(poStartPoint,poEndPoint,true);
		m_lpoSegments.push_back(poSegment);
	}
	void Curve::FrontPushSegment(const Point& oStartPoint,const Point& oEndPoint)
	{
		GenericNode* poStartPoint = new GenericNode(oStartPoint);
		GenericNode* poEndPoint = new GenericNode(oEndPoint);
		Segment* poSegment = new Segment(poStartPoint,poEndPoint,true);
		m_lpoSegments.push_front(poSegment);
	}
	void Curve::WriteVTK(const string sFileName) const
	{
		FILE* fpFile = fopen(sFileName.c_str(),"w");

		fprintf(fpFile,"# vtk DataFile Version 1.0\n");
		fprintf(fpFile,"curve\n");
		fprintf(fpFile,"ASCII\n");
		fprintf(fpFile,"DATASET POLYDATA\n");
		
		unsigned int iEdgesCount = (unsigned int)m_lpoSegments.size();
		unsigned int iNodesCount = 2*iEdgesCount;
		fprintf(fpFile,"POINTS %d float\n",iNodesCount);
		list<Segment*>::const_iterator liSegments;
		GenericNode* poPoint = NULL;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			poPoint = (*liSegments)->GetPoint1();
			fprintf(fpFile,"%E\t\t%E\t\t%E\n",poPoint->GetX(),poPoint->GetY(),poPoint->GetZ());
			poPoint = (*liSegments)->GetPoint2();
			fprintf(fpFile,"%E\t\t%E\t\t%E\n",poPoint->GetX(),poPoint->GetY(),poPoint->GetZ());
		}		
		
		fprintf(fpFile,"LINES %d %d\n",iEdgesCount,3*iEdgesCount);
		unsigned int iCounter = 0;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			fprintf(fpFile,"2\t\t%d\t\t%d\n",iCounter++,iCounter++);
		}

		fprintf(fpFile,"POINT_DATA %d\n",iNodesCount);
		fprintf(fpFile,"scalars NodeType integer\n");
		fprintf(fpFile,"LOOKUP_TABLE default\n");
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			fprintf(fpFile,"1\n");
			fprintf(fpFile,"1\n");
		}
		
		fprintf(fpFile,"CELL_DATA %d\n",iEdgesCount);
		fprintf(fpFile,"SCALARS SlipSystem integer\n");
		fprintf(fpFile,"LOOKUP_TABLE default\n");
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			fprintf(fpFile,"1\n");
		}
		
		fclose(fpFile);
	}
	void Curve::CyclicSort()
	{
		if(m_lpoSegments.size() < 2)			return;
		list<Segment*>::iterator liSegments = m_lpoSegments.begin();
		list<Segment*> lpoSortedSegments;
		Segment* poActiveSegment = (*liSegments);
		lpoSortedSegments.push_back(poActiveSegment);
		m_lpoSegments.erase(liSegments);
		liSegments = m_lpoSegments.begin();
		bool bEndStartTouching = false;
		bool bEndEndTouching = false;
		Segment* poCurrentSegment = NULL;
		double dTolerance = 1.0E-6;
		GenericNode* poActiveEndPoint = poActiveSegment->GetPoint2();
		while(m_lpoSegments.size() > 1)
		{
			poCurrentSegment = (*liSegments);
			bEndStartTouching = (poActiveEndPoint->Distance(*(poCurrentSegment->GetPoint1()))) < dTolerance;
			bEndEndTouching = (poActiveEndPoint->Distance(*(poCurrentSegment->GetPoint2()))) < dTolerance;
			if(bEndEndTouching)
			{
				poCurrentSegment->Flip();
				bEndStartTouching = true;
			}
			if(bEndStartTouching)
			{
				poActiveSegment = poCurrentSegment;
				lpoSortedSegments.push_back(poActiveSegment);
				liSegments = m_lpoSegments.erase(liSegments);
				poActiveEndPoint = poActiveSegment->GetPoint2();
			}
			else
			{
				liSegments++;
			}
			if(liSegments == m_lpoSegments.end())	liSegments = m_lpoSegments.begin();
		}
		// handle the last segment
		poActiveSegment = m_lpoSegments.front();
		lpoSortedSegments.push_back(poActiveSegment);
		m_lpoSegments.clear();
		poActiveEndPoint = poActiveSegment->GetPoint2();
				
		// see if the last segment needs flipping
		if(poActiveEndPoint->Distance(*(lpoSortedSegments.front()->GetPoint1())) > dTolerance)		poActiveSegment->Flip();
		// add it to the sorted segments list
		
		liSegments = lpoSortedSegments.begin();
		while(liSegments != lpoSortedSegments.end())
		{
			m_lpoSegments.push_back((*liSegments));
			liSegments++;
		}
		lpoSortedSegments.clear();
	}
	void Curve::Localize(CartesianOrthogonalCoordinateSystem* poSystem)
	{
		list<Segment*>::iterator liSegments;
		Point oTempPoint;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			oTempPoint = poSystem->GetInLocalCoordinates(*(*liSegments)->GetPoint1());
			(*liSegments)->GetPoint1()->Set(oTempPoint);
			oTempPoint = poSystem->GetInLocalCoordinates(*(*liSegments)->GetPoint2());
			(*liSegments)->GetPoint2()->Set(oTempPoint);
		}
	}
	void Curve::Globalize(CartesianOrthogonalCoordinateSystem* poSystem)
	{
		list<Segment*>::iterator liSegments;
		Point oTempPoint;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			oTempPoint = poSystem->GetInGlobalCoordinates(*(*liSegments)->GetPoint1());
			(*liSegments)->GetPoint1()->Set(oTempPoint);
			oTempPoint = poSystem->GetInGlobalCoordinates(*(*liSegments)->GetPoint2());
			(*liSegments)->GetPoint2()->Set(oTempPoint);
		}
	}
	void Curve::Expand(const Vector& oPlaneNormal,const double& dDistance)
	{
		if(m_lpoSegments.empty())				return;
		list<Segment*>::iterator liSegments = m_lpoSegments.begin();
		Vector oShiftVector;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			oShiftVector = (*liSegments)->GetDirection()^oPlaneNormal;
			oShiftVector.Normalize();
			oShiftVector = oShiftVector*dDistance;
			(*liSegments)->GetPoint1()->Shift(oShiftVector.GetX(),oShiftVector.GetY(),oShiftVector.GetZ());
			(*liSegments)->GetPoint2()->Shift(oShiftVector.GetX(),oShiftVector.GetY(),oShiftVector.GetZ());
		}
		liSegments = m_lpoSegments.begin();
		Segment* poPreviousSegment = NULL;
		Segment* poThisSegment = (*liSegments);
		GenericNode* poNewStartPoint = NULL;
		GenericNode* poNewEndPoint = NULL;
		Segment* poNewSegment = NULL;
		while(true)
		{
			poPreviousSegment = poThisSegment;
			liSegments++;
			if(liSegments == m_lpoSegments.end())		break;
			poThisSegment = (*liSegments);
			poNewStartPoint = new GenericNode(*(poPreviousSegment->GetPoint2()));
			poNewEndPoint = new GenericNode(*(poThisSegment->GetPoint1()));
			poNewSegment = new Segment(poNewStartPoint,poNewEndPoint,true);
			liSegments = m_lpoSegments.insert(liSegments,poNewSegment);
			liSegments++;
		}
		// insert the last segment
		liSegments = m_lpoSegments.begin();
		poThisSegment = (*liSegments);
		poNewStartPoint = new GenericNode(*(poPreviousSegment->GetPoint2()));
		poNewEndPoint = new GenericNode(*(poThisSegment->GetPoint1()));
		poNewSegment = new Segment(poNewStartPoint,poNewEndPoint,true);
		liSegments = m_lpoSegments.insert(liSegments,poNewSegment);
	}
	list<Point> Curve::GetIntersectionPoints(const Segment& oSegment) const
	{
		list<Segment*>::const_iterator liSegments;
		Point oIntersectionPoint;
		list<Point> loIntersectionPoints;
		loIntersectionPoints.clear();
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			if((*liSegments)->GetIntersectionPoint(oSegment,oIntersectionPoint))
			{
				loIntersectionPoints.push_back(oIntersectionPoint);
			}
		}
		return loIntersectionPoints;
	}
	list<Point> Curve::GetIntersectionPoints(const Curve& oCurve) const
	{
		list<Segment*>::const_iterator liSegments;
		list<Point> loIntersectionPoints;
		list<Point> loTempPoints;
		list<Point>::iterator liPoints;
		loIntersectionPoints.clear();
		for(liSegments = oCurve.m_lpoSegments.begin() ; liSegments != oCurve.m_lpoSegments.end() ; liSegments++)
		{
			loTempPoints = GetIntersectionPoints(*(*liSegments));
			for(liPoints = loTempPoints.begin() ; liPoints != loTempPoints.end() ; liPoints++)
			{
				loIntersectionPoints.push_back((*liPoints));
			}
		}
		return loIntersectionPoints;
	}
	void Curve::Split(const Point& oPoint)
	{
		list<Segment*>::iterator liSegments;
		GenericNode* poNewStartNode = NULL;
		GenericNode* poNewEndNode = NULL;
		Segment* poNewSegment = NULL;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			if((*liSegments)->IsPointOnSegment(oPoint))
			{
				// segment found, create the new segment nodes and move the original segment's end point
				poNewStartNode = new GenericNode(oPoint);
				poNewEndNode = new GenericNode(*(*liSegments)->GetPoint2());
				(*liSegments)->GetPoint2()->Set(oPoint);
				// create the new segment and insert it AFTER the current segment
				poNewSegment = new Segment(poNewStartNode,poNewEndNode,true);
				liSegments++;
				m_lpoSegments.insert(liSegments,poNewSegment);
				break;
			}
		}
	}
	Curve* Curve::ExtractSubCurve(const Point& oStartPoint,const Point& oEndPoint)
	{
		list<Segment*>::iterator liSegments;
		bool bStartPointFound = false;
		double dTolerance = 1.0E-8;
		list<Segment*> lpoSubCurveSegments;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			if(!bStartPointFound)
			{
				if((*liSegments)->GetPoint1()->GetDistanceSquared(oStartPoint) < dTolerance)
				{
					bStartPointFound = true;
					lpoSubCurveSegments.push_back((*liSegments));
				}
			}
			else
			{
				if((*liSegments)->GetPoint1()->GetDistanceSquared(oEndPoint) < dTolerance)
				{
					break;
				}
				else
				{
					lpoSubCurveSegments.push_back((*liSegments));
				}
			}
		}
		// create the subcurve, no need to sort
		Curve* poSubCurve = new Curve;
		for(liSegments = lpoSubCurveSegments.begin() ; liSegments != lpoSubCurveSegments.end() ; liSegments++)
		{
			poSubCurve->AddSegment(*((*liSegments)->GetPoint1()),*((*liSegments)->GetPoint2()));
		}
		lpoSubCurveSegments.clear();
		return poSubCurve;
	}
	GenericNode* Curve::GetStartPoint() const
	{
		return m_lpoSegments.front()->GetPoint1();
	}
	GenericNode* Curve::GetEndPoint() const
	{
		return m_lpoSegments.back()->GetPoint2();
	}
	list<GenericNode*> Curve::GetPoints() const
	{
		list<GenericNode*> lpoPoints;
		list<Segment*>::const_iterator liSegments;
		for(liSegments = m_lpoSegments.begin() ; liSegments != m_lpoSegments.end() ; liSegments++)
		{
			lpoPoints.push_back((*liSegments)->GetPoint1());
			lpoPoints.push_back((*liSegments)->GetPoint2());
		}
		return lpoPoints;
	}
	void Curve::Flip()
	{
		list<Segment*>::iterator liSegments = m_lpoSegments.begin();
		while(liSegments != m_lpoSegments.end())
		{
			(*liSegments)->Flip();
			m_lpoSegments.push_front((*liSegments));
			liSegments = m_lpoSegments.erase(liSegments);
		}
	}
	void Curve::Initialize()
	{
		m_lpoSegments.clear();
	}
}


