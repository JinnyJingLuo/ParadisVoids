// Ahmed M. Hussein

#include "SurfaceServer.h"
#include "Randomizer.h"
#include "float.h"
#include "cmath"

using namespace std;
using namespace EZ;

SurfaceServer* SurfaceServer::m_poInstance = NULL;
SurfaceServer* SurfaceServer::GetInstance()
{
	if(m_poInstance == NULL)
	{
		m_poInstance = new SurfaceServer;
	}
	return m_poInstance;
}
SurfaceServer::~SurfaceServer()
{
	Reset();
}
void SurfaceServer::Reset()
{
	list<TriPatch*>::iterator liTriangles;
	for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
	{
		if((*liTriangles) != NULL)
		{
			delete (*liTriangles);
		}
	}
	m_lpoTriangles.clear();
	
	list<GenericNode*>::iterator liPoints;
	for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ; liPoints++)
	{
		if((*liPoints) != NULL)
		{
			delete (*liPoints);
		}
	}
	m_lpoPoints.clear();
}
void SurfaceServer::SetDataStructure(MainDataStructure* poDataStructure)
{
	Reset();
	m_poDataStructure = poDataStructure;
	AxisAlignedBoundingBox* poBox = m_poDataStructure->GetSimulationBox();
	m_dXMin = poBox->GetXMin();
	m_dXMax = poBox->GetXMax();
	m_dYMin = poBox->GetYMin();
	m_dYMax = poBox->GetYMax();
	m_dZMin = poBox->GetZMin();
	m_dZMax = poBox->GetZMax();
	double dToleranceFactor = 1.0E-6;
	double dMinimumDimension = m_dXMax - m_dXMin;
	if(dMinimumDimension > (m_dYMax - m_dYMin))
	{
		dMinimumDimension = m_dYMax - m_dZMin;
	}
	if(dMinimumDimension > (m_dZMax - m_dYMin))
	{
		dMinimumDimension = m_dZMax - m_dZMin;
	}
	double dTolerance = dToleranceFactor*dMinimumDimension;
	m_dXMin = m_dXMin - dTolerance;
	m_dXMax = m_dXMax + dTolerance;
	m_dYMin = m_dYMin - dTolerance;
	m_dYMax = m_dYMax + dTolerance;
	m_dZMin = m_dZMin - dTolerance;
	m_dZMax = m_dZMax + dTolerance;
	GenerateTriangulations();
}
void SurfaceServer::CheckNodes()
{
	DislocationNode* poNode = NULL;
	Point oNodePoint;
	Point oNearestSurfacePoint;
	TriPatch* poTriangle = NULL;
	Vector oSurfaceNormal;
	unsigned int iConstraintType;
	Vector oConstraintVector;
	double dSurfaceAnnihilationDistance = 1.0;
	bool bIsSurfaceNode = false;
	list<DislocationNode*>* plpoNodes = m_poDataStructure->GetLocalNodes();
	list<DislocationNode*>::iterator liNodes;
	for(liNodes = plpoNodes->begin() ; liNodes != plpoNodes->end() ; liNodes++)
	{
		poNode = (*liNodes);
		if(poNode->GetCategory() == PinnedNode)
		{
			continue;
		}
		oNodePoint = *poNode;
		// a node is converted to a surface node if it is outside the box or the minimum
		// distance to the surface is less than the surface annihilation distance
		bIsSurfaceNode = false;
		//bIsSurfaceNode = !IsPointInside(oNodePoint);
		// faster than the above function for cuboidal volumes
		bIsSurfaceNode = !IsPointInBox(oNodePoint);
		bIsSurfaceNode = bIsSurfaceNode || (GetLeastSurfaceDistance(oNodePoint) < dSurfaceAnnihilationDistance);
		if(bIsSurfaceNode)
		{
			poNode->GetDynamicConstraint(iConstraintType,oConstraintVector);
			poTriangle = NULL;
			if(iConstraintType == 0)
			{
				// node is free to move in the space
				poTriangle = GetNearestTriangle(oNodePoint,oNearestSurfacePoint);
			}
			else if(iConstraintType == 1)
			{
				// node moves in a plane
				Plane oPlane(oConstraintVector,oNodePoint);
				poTriangle = GetNearestTriangleOnPlane(oPlane,oNodePoint,oNearestSurfacePoint);
			}
			else if(iConstraintType == 2)
			{
				// node moves on a line
				Line oLine(oConstraintVector,oNodePoint);
				poTriangle = GetNearestTriangleOnLine(oLine,oNodePoint,oNearestSurfacePoint);
			}
			else if(iConstraintType == 3)
			{
				// node cannot move at all
				oNearestSurfacePoint = oNodePoint;
			}
			if(poTriangle == NULL)
			{
				// in case the triangle was not obtained at all, just
				// get the nearest surface point. this is not the best 
				// solution in this case. 
				poTriangle = GetNearestTriangle(oNodePoint,oNearestSurfacePoint);
			}
			poNode->Set(oNearestSurfacePoint.GetX(),oNearestSurfacePoint.GetY(),oNearestSurfacePoint.GetZ());
			poNode->SetCategory(SurfaceNode);
			poNode->SetSurfaceNormal(poTriangle->GetNormal());
		}
	}
}
void SurfaceServer::RemoveSurfaceArms()
{	
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	ClearOpList(poHome);
// 	unsigned int iArmsCount = 0;
// 	list<Node_t*> lpoNeighbours;
// 	list<Node_t*>::iterator liNeighbours;
// 	unsigned int i = 0;
// 	unsigned int j = 0;
// 	unsigned int iSize = (unsigned int)poHome->newNodeKeyPtr;
// 	Node_t* poNode = NULL;
// 	Node_t* poNeighbour = NULL;
// 	for(i = 0 ; i < iSize ; i++)
// 	{
// 		poNode = poHome->nodeKeys[i];
// 		if(poNode == NULL)
// 		{
// 			continue;
// 		}
// 		if(poNode->constraint == SURFACE_NODE)
// 		{
// 			iArmsCount = poNode->numNbrs;
// 			lpoNeighbours.clear();
// 			for(j = 0 ; j < iArmsCount ; j++)
// 			{
// 				poNeighbour = GetNeighborNode(poHome,poNode,j);
// 				if(poNeighbour == NULL)
// 				{
// 					continue;
// 				}
// 				if((poNeighbour->constraint == SURFACE_NODE) || (poNeighbour->constraint == VIRTUAL_NODE))
// 				{
// 					lpoNeighbours.push_back(poNeighbour);
// 				}
// 			}
// 			for(liNeighbours = lpoNeighbours.begin() ; liNeighbours != lpoNeighbours.end() ; liNeighbours++)
// 			{
// 				poNeighbour = *liNeighbours;
// 				ChangeArmBurg(poHome,poNode,&poNeighbour->myTag,0.0,0.0,0.0,0.0,0.0,0.0,1,DEL_SEG_NONE);
// 				ChangeArmBurg(poHome,poNeighbour,&poNode->myTag,0.0,0.0,0.0,0.0,0.0,0.0,1,DEL_SEG_NONE);
// 			}
// 		}
// 	}
// 	// get rid of any nodes that have no arms
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	CommSendRemesh(poHome);
// 	FixRemesh(poHome);
// 	RemoveOrphanedNodes(poHome);
// 	CommSendRemesh(poHome);
// 	FixRemesh(poHome);
// 	MPI_Barrier(MPI_COMM_WORLD);
}
SurfaceServer::SurfaceServer()
{
	Initialize();
}
void SurfaceServer::Initialize()
{
	m_lpoTriangles.clear();
	m_lpoPoints.clear();
	m_dXMin = 0.0;
	m_dXMax = 0.0;
	m_dYMin = 0.0;
	m_dYMax = 0.0;
	m_dZMin = 0.0;
	m_dZMax = 0.0;
}
TriPatch* SurfaceServer::GetNearestTriangle(const Point& oPoint,Point& oNearestPoint) const
{
	list<TriPatch*>::const_iterator liTriangles;
	double dTemp = 0.0;
	double dMinDistance = DBL_MAX;
	TriPatch* poTriangle = NULL;
	Point oTempPoint;
	for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
	{
		oTempPoint = (*liTriangles)->GetNearestPoint(oPoint,dTemp);
		if(dTemp < dMinDistance)
		{
			dMinDistance = dTemp;
			poTriangle = *liTriangles;
			oNearestPoint = oTempPoint;
		}
	}
	return poTriangle;
}
TriPatch* SurfaceServer::GetNearestTriangleOnPlane(const Plane& oPlane,const Point& oPoint,Point& oNearestPoint) const
{
	list<TriPatch*>::const_iterator liTriangles;
	double dTemp = 0.0;
	double dMinDistance = DBL_MAX;
	TriPatch* poTriangle = NULL;
	Point oTempPoint;
	for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
	{
		if((*liTriangles)->GetNearestPointOnPlane(oPoint,oPlane,oTempPoint,dTemp))
		{
			if(dTemp < dMinDistance)
			{
				dMinDistance = dTemp;
				poTriangle = *liTriangles;
				oNearestPoint = oTempPoint;
			}
		}
	}
	return poTriangle;
}
TriPatch* SurfaceServer::GetNearestTriangleOnLine(const Line& oLine,const Point& oPoint,Point& oNearestPoint) const
{
	list<TriPatch*>::const_iterator liTriangles;
	double dTemp = 0.0;
	double dMinDistance = DBL_MAX;
	TriPatch* poTriangle = NULL;
	Point oTempPoint;
	for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
	{
		if((*liTriangles)->GetNearestPointOnLine(oPoint,oLine,oTempPoint,dTemp))
		{
			if(dTemp < dMinDistance)
			{
				dMinDistance = dTemp;
				poTriangle = *liTriangles;
				oNearestPoint = oTempPoint;
			}
		}
	}
	return poTriangle;
}
bool SurfaceServer::IsPointInBox(const Point& oPoint)
{
	double dX = oPoint.GetX();
	if((dX >= m_dXMin) && (dX <= m_dXMax))
	{
		double dY = oPoint.GetY();
		if((dY >= m_dYMin) && (dY <= m_dYMax))
		{
			double dZ = oPoint.GetZ();
			if((dZ >= m_dZMin) && (dZ <= m_dZMax))
			{
				return true;
			}
		}
	}
	return false;
}
bool SurfaceServer::IsPointInside(const Point& oPoint)
{
	unsigned int iTestPointsCount = 5;
	double dAmplificationFactor = 1.5;
	unsigned int i = 0;
	unsigned int j = 0;
	Point oTestPoint;
	double dMaxDimension = max((m_dXMax - m_dXMin),max((m_dYMax - m_dYMin),(m_dZMax - m_dZMin)));
	double dXMean = 0.5*(m_dXMax + m_dXMin);
	double dYMean = 0.5*(m_dYMax + m_dYMin);
	double dZMean = 0.5*(m_dZMax + m_dZMin);
	Point oCenter(dXMean,dYMean,dZMean);
	double dRadius = dAmplificationFactor*dMaxDimension;
	double dTheta = 0.0;
	double dPhi = 0.0;
	list<TriPatch*>::const_iterator liTriangles;
	double dIntersectionParameter = 0.0;
	Point oIntersectionPoint;
	unsigned int iCount = 0;
	unsigned int iIsPointInsideCount = 0;
	for(i = 0 ; i < iTestPointsCount ; i++)
	{
		dPhi = Randomizer::Random(0,PI);
		dTheta = Randomizer::Random(0,2*PI);
		oTestPoint.Set(dRadius*sin(dPhi)*cos(dTheta),dRadius*sin(dPhi)*sin(dTheta),dRadius*cos(dPhi));
		oTestPoint = oTestPoint + oCenter;
		iCount = 0;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			if((*liTriangles)->GetLineIntersection(oPoint,oTestPoint,oIntersectionPoint,dIntersectionParameter))
			{
				iCount = iCount + 1;
			}
		}
		if(iCount%2 != 0)
		{
			iIsPointInsideCount = iIsPointInsideCount + 1;
		}
	}
	if(2*iIsPointInsideCount > iTestPointsCount)
	{
		return true;
	}
	return false;
}
bool SurfaceServer::IsPointOnSurface(const Point& oPoint) const
{
	list<TriPatch*>::const_iterator liTriangles;
	for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
	{
		if((*liTriangles)->IsPointInTriangle(oPoint))
		{
			return true;
		}
	}
	return false;
}
double SurfaceServer::GetLeastSurfaceDistance(const Point& oPoint) const
{
	list<TriPatch*>::const_iterator liTriangles;
	double dMin = DBL_MAX;
	double dTemp = 0.0;
	for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
	{
		dTemp = (*liTriangles)->GetNormalDistance(oPoint);
		if(dTemp < dMin)
		{
			dMin = dTemp;
		}
	}
	return dMin;
}
void SurfaceServer::GenerateTriangulations()
{
	GenericNode* poPointNNN = new GenericNode(m_dXMin,m_dYMin,m_dZMin);
	GenericNode* poPointPNN = new GenericNode(m_dXMax,m_dYMin,m_dZMin);
	GenericNode* poPointPPN = new GenericNode(m_dXMax,m_dYMax,m_dZMin);
	GenericNode* poPointNPN = new GenericNode(m_dXMin,m_dYMax,m_dZMin);
	GenericNode* poPointNNP = new GenericNode(m_dXMin,m_dYMin,m_dZMax);
	GenericNode* poPointPNP = new GenericNode(m_dXMax,m_dYMin,m_dZMax);
	GenericNode* poPointPPP = new GenericNode(m_dXMax,m_dYMax,m_dZMax);
	GenericNode* poPointNPP = new GenericNode(m_dXMin,m_dYMax,m_dZMax);
	m_lpoPoints.push_back(poPointNNN);
	m_lpoPoints.push_back(poPointPNN);
	m_lpoPoints.push_back(poPointPPN);
	m_lpoPoints.push_back(poPointNPN);
	m_lpoPoints.push_back(poPointNNP);
	m_lpoPoints.push_back(poPointPNP);
	m_lpoPoints.push_back(poPointPPP);
	m_lpoPoints.push_back(poPointNPP);
	
	m_lpoTriangles.push_back(new TriPatch(poPointNNN,poPointPNN,poPointNNP));
	m_lpoTriangles.push_back(new TriPatch(poPointPNN,poPointPNP,poPointNNP));
	
	m_lpoTriangles.push_back(new TriPatch(poPointNPN,poPointNPP,poPointPPP));
	m_lpoTriangles.push_back(new TriPatch(poPointNPN,poPointPPP,poPointPPN));
	
	m_lpoTriangles.push_back(new TriPatch(poPointNPN,poPointNNN,poPointNNP));
	m_lpoTriangles.push_back(new TriPatch(poPointNPN,poPointNNP,poPointNPP));
	
	m_lpoTriangles.push_back(new TriPatch(poPointPNN,poPointPPN,poPointPPP));
	m_lpoTriangles.push_back(new TriPatch(poPointPNN,poPointPPP,poPointPNP));
	
	m_lpoTriangles.push_back(new TriPatch(poPointNNN,poPointPPN,poPointPNN));
	m_lpoTriangles.push_back(new TriPatch(poPointNNN,poPointNPN,poPointPPN));
	
	m_lpoTriangles.push_back(new TriPatch(poPointPNP,poPointNPP,poPointNNP));
	m_lpoTriangles.push_back(new TriPatch(poPointPNP,poPointPPP,poPointNPP));
}



