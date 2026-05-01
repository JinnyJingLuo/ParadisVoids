// Ahmed M. Hussein
#include "ParadisTriangulatedSurface.h"

// ParadisSurface::ParadisSurface()
// {
// 	Initialize();
// }
// ParadisSurface::ParadisSurface(const ParadisSurface& oSurface)
// {
// 	*this = oSurface;
// }
// ParadisSurface::~ParadisSurface()
// {
// 	Reset();
// }
// ParadisSurface& ParadisSurface::operator=(const ParadisSurface& oSurface)
// {
// 	Reset();
// 	m_dXMin = oSurface.m_dXMin;
// 	m_dXMax = oSurface.m_dXMax;
// 	m_dYMin = oSurface.m_dYMin;
// 	m_dYMax = oSurface.m_dYMax;
// 	m_dZMin = oSurface.m_dZMin;
// 	m_dZMax = oSurface.m_dZMax;
// 	m_dVolume = oSurface.m_dVolume;
// 	m_dTolerance = oSurface.m_dTolerance;
// 	GenerateTriangulations();
// 	return *this;
// }
// void ParadisSurface::Reset()
// {
// 	ClearTriangulations();
// 	Initialize();
// }
// void ParadisSurface::Set(Home_t* poHome)
// {
// 	Reset();
// 	m_dXMin = poHome->param->minCoordinates[0];
// 	m_dXMax = poHome->param->maxCoordinates[0];
// 	m_dYMin = poHome->param->minCoordinates[1];
// 	m_dYMax = poHome->param->maxCoordinates[1];
// 	m_dZMin = poHome->param->minCoordinates[2];
// 	m_dZMax = poHome->param->maxCoordinates[2];
// 	double dToleranceFactor = 1.0E-6;
// 	double dMinimumDimension = m_dXMax - m_dXMin;
// 	if(dMinimumDimension > (m_dYMax - m_dYMin))
// 	{
// 		dMinimumDimension = m_dYMax - m_dZMin;
// 	}
// 	if(dMinimumDimension > (m_dZMax - m_dYMin))
// 	{
// 		dMinimumDimension = m_dZMax - m_dZMin;
// 	}
// 	m_dTolerance = dToleranceFactor*dMinimumDimension;
// 	m_dXMin = m_dXMin - m_dTolerance;
// 	m_dXMax = m_dXMax + m_dTolerance;
// 	m_dYMin = m_dYMin - m_dTolerance;
// 	m_dYMax = m_dYMax + m_dTolerance;
// 	m_dZMin = m_dZMin - m_dTolerance;
// 	m_dZMax = m_dZMax + m_dTolerance;
// 	m_dVolume = (m_dXMax - m_dXMin)*(m_dYMax - m_dYMin)*(m_dZMax - m_dZMin);
// 	GenerateTriangulations();
// 	poHome->param->simVol = m_dVolume;
// 	poHome->param->burgVolFactor = 1.0
// /(poHome->param->burgMag*poHome->param->burgMag*poHome->param->simVol);
// }
// void ParadisSurface::CheckNodes(Home_t* poHome)
// {
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	unsigned int iSize = (unsigned int)poHome->newNodeKeyPtr;
// 	unsigned int i = 0;
// 	Node_t* poNode = NULL;
// 	Point oNodePoint;
// 	Point oNearestSurfacePoint;
// 	TriPatch* poTriangle = NULL;
// 	Vector oSurfaceNormal;
// 	double daPos[3] = {0.0,0.0,0.0};
// 	ClearOpList(poHome);
// 	unsigned int iConstraintType;
// 	Vector oConstraintVector;
// 	double dSurfaceAnnihilationDistance = 1.0;
// 	bool bIsSurfaceNode = false;
// 	for(i = 0 ; i < iSize ; i++)
// 	{
// 		poNode = poHome->nodeKeys[i];
// 		if(poNode == NULL)
// 		{
// 			continue;
// 		}
// 		if(poNode->constraint == PINNED_NODE)
// 		{
// 			continue;
// 		}
// 		oNodePoint.Set(poNode->x,poNode->y,poNode->z);
// 		// a node is converted to a surface node if it is outside the
// box or the minimum
// 		// distance to the surface is less than the surface annihilation
// distance 		bIsSurfaceNode = false;
// 		//bIsSurfaceNode = !IsPointInside(oNodePoint);
// 		// faster than the above function for cuboidal volumes
// 		bIsSurfaceNode = !IsPointInBox(oNodePoint);
// 		bIsSurfaceNode = bIsSurfaceNode ||
// (GetLeastSurfaceDistance(oNodePoint) < dSurfaceAnnihilationDistance);
// 		if(bIsSurfaceNode)
// 		{
// 			GetNodeDynamicConstraint(poNode,iConstraintType,oConstraintVector);
// 			poTriangle = NULL;
// 			if(iConstraintType == 0)
// 			{
// 				// node is free to move in the space
// 				poTriangle =
// GetNearestTriangle(oNodePoint,oNearestSurfacePoint);
// 			}
// 			else if(iConstraintType == 1)
// 			{
// 				// node moves in a plane
// 				Plane oPlane(oConstraintVector,oNodePoint);
// 				poTriangle =
// GetNearestTriangleOnPlane(oPlane,oNodePoint,oNearestSurfacePoint);
// 			}
// 			else if(iConstraintType == 2)
// 			{
// 				// node moves on a line
// 				Line oLine(oConstraintVector,oNodePoint);
// 				poTriangle =
// GetNearestTriangleOnLine(oLine,oNodePoint,oNearestSurfacePoint);
// 			}
// 			else if(iConstraintType == 3)
// 			{
// 				// node cannot move at all
// 				oNearestSurfacePoint = oNodePoint;
// 			}
// 			if(poTriangle == NULL)
// 			{
// 				// in case the triangle was not obtained at all,
// just
// 				// get the nearest surface point. this is not
// the best
// 				// solution in this case.
// 				poTriangle =
// GetNearestTriangle(oNodePoint,oNearestSurfacePoint);
// 			}
// 			daPos[0] = oNearestSurfacePoint.GetX();
// 			daPos[1] = oNearestSurfacePoint.GetY();
// 			daPos[2] = oNearestSurfacePoint.GetZ();
//
// 			RepositionNode(poHome,daPos,&(poNode->myTag),1);
//
// 			poNode->constraint = SURFACE_NODE;
// 			oSurfaceNormal = poTriangle->GetNormal();
// 			poNode->dNx = oSurfaceNormal.GetX();
// 			poNode->dNy = oSurfaceNormal.GetY();
// 			poNode->dNz = oSurfaceNormal.GetZ();
// 		}
// 	}
//
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	CommSendRemesh(poHome);
// 	FixRemesh(poHome);
// 	MPI_Barrier(MPI_COMM_WORLD);
// }
// void ParadisSurface::RemoveSurfaceArms(Home_t* poHome)
// {
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
// 		if(poNode->constraint >=SURFACE_NODE)
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
// 				if(poNeighbour->constraint >=SURFACE_NODE)
// 				{
// 					lpoNeighbours.push_back(poNeighbour);
// 				}
// 			}
// 			for(liNeighbours = lpoNeighbours.begin() ; liNeighbours
// != lpoNeighbours.end() ; liNeighbours++)
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
// }
// void ParadisSurface::Initialize()
// {
// 	m_lpoTriangles.clear();
// 	m_lpoPoints.clear();
// 	m_dXMin = 0.0;
// 	m_dXMax = 0.0;
// 	m_dYMin = 0.0;
// 	m_dYMax = 0.0;
// 	m_dZMin = 0.0;
// 	m_dZMax = 0.0;
// 	m_dVolume = 0.0;
// 	m_dTolerance = 0.0;
// }
// void ParadisSurface::GetNodeDynamicConstraint(Node_t* poNode,unsigned int&
// iConstraintType,Vector& oConstraintVector)
// {
// 	oConstraintVector.Set(0.0,0.0,0.0);
// 	iConstraintType = 0;
// 	double dTolerance = 1.0E-6;
// 	unsigned int i = 0;
// 	// get the node dynamic constraint
// 	if(poNode->numNbrs == 0)
// 	{
// 		oConstraintVector.Set(0.0,0.0,0.0);
// 		iConstraintType = 0;
// 	}
// 	else if(poNode->numNbrs == 1)
// 	{
// 		oConstraintVector.Set(poNode->nx[0],poNode->ny[0],poNode->nz[0]);
// 		iConstraintType = 1;
// 	}
// 	else if(poNode->numNbrs == 2)
// 	{
// 		Vector oNormal1(poNode->nx[0],poNode->ny[0],poNode->nz[0]);
// 		Vector oNormal2(poNode->nx[1],poNode->ny[1],poNode->nz[1]);
// 		oNormal1.Normalize();
// 		oNormal2.Normalize();
// 		Vector oCrossProduct = oNormal1^oNormal2;
// 		if(oCrossProduct.Length() < dTolerance)
// 		{
// 			iConstraintType = 1;
// 			oConstraintVector = oNormal1;
// 		}
// 		else
// 		{
// 			iConstraintType = 2;
// 			oConstraintVector = oCrossProduct;
// 		}
// 	}
// 	else
// 	{
// 		list<Vector> loNormals;
// 		list<Vector>::iterator liNormals;
// 		loNormals.push_back(Vector(poNode->nx[0],poNode->ny[0],poNode->nz[0]));
// 		loNormals.front().Normalize();
// 		Vector oTempNormal;
// 		Vector oCrossProduct;
// 		bool bAddNormal = false;
// 		for(i = 1 ; i < poNode->numNbrs ; i++)
// 		{
// 			oTempNormal.Set(poNode->nx[i],poNode->ny[i],poNode->nz[i]);
// 			oTempNormal.Normalize();
// 			bAddNormal = false;
// 			for(liNormals = loNormals.begin() ; liNormals !=
// loNormals.end() ; liNormals++)
// 			{
// 				oCrossProduct = oTempNormal^(*liNormals);
// 				if(oCrossProduct.Length() > dTolerance)
// 				{
// 					bAddNormal = true;
// 					break;
// 				}
// 			}
// 			if(bAddNormal)
// 			{
// 				loNormals.push_back(oTempNormal);
// 			}
// 		}
// 		unsigned int iDynamicConstraintsCount = (unsigned
// int)loNormals.size(); 		if(iDynamicConstraintsCount == 1)
// 		{
// 			iConstraintType = 1;
// 			oConstraintVector.Set(poNode->nx[0],poNode->ny[0],poNode->nz[0]);
// 		}
// 		else if(iDynamicConstraintsCount == 2)
// 		{
// 			iConstraintType = 2;
// 			oConstraintVector = loNormals.front()^loNormals.back();
// 		}
// 		else
// 		{
// 			iConstraintType = 3;
// 			oConstraintVector.Set(0.0,0.0,0.0);
// 		}
// 	}
// 	oConstraintVector.Normalize();
// }
// TriPatch* ParadisSurface::GetNearestTriangle(const Point& oPoint,Point&
// oNearestPoint) const
// {
// 	list<TriPatch*>::const_iterator liTriangles;
// 	double dTemp = 0.0;
// 	double dMinDistance = DBL_MAX;
// 	TriPatch* poTriangle = NULL;
// 	Point oTempPoint;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		oTempPoint = (*liTriangles)->GetNearestPoint(oPoint,dTemp);
// 		if(dTemp < dMinDistance)
// 		{
// 			dMinDistance = dTemp;
// 			poTriangle = *liTriangles;
// 			oNearestPoint = oTempPoint;
// 		}
// 	}
// 	return poTriangle;
// }
// TriPatch* ParadisSurface::GetNearestTriangleOnPlane(const Plane& oPlane,const
// Point& oPoint,Point& oNearestPoint) const
// {
// 	list<TriPatch*>::const_iterator liTriangles;
// 	double dTemp = 0.0;
// 	double dMinDistance = DBL_MAX;
// 	TriPatch* poTriangle = NULL;
// 	Point oTempPoint;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		if((*liTriangles)->GetNearestPointOnPlane(oPoint,oPlane,oTempPoint,dTemp))
// 		{
// 			if(dTemp < dMinDistance)
// 			{
// 				dMinDistance = dTemp;
// 				poTriangle = *liTriangles;
// 				oNearestPoint = oTempPoint;
// 			}
// 		}
// 	}
// 	return poTriangle;
// }
// TriPatch* ParadisSurface::GetNearestTriangleOnLine(const Line& oLine,const
// Point& oPoint,Point& oNearestPoint) const
// {
// 	list<TriPatch*>::const_iterator liTriangles;
// 	double dTemp = 0.0;
// 	double dMinDistance = DBL_MAX;
// 	TriPatch* poTriangle = NULL;
// 	Point oTempPoint;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		if((*liTriangles)->GetNearestPointOnLine(oPoint,oLine,oTempPoint,dTemp))
// 		{
// 			if(dTemp < dMinDistance)
// 			{
// 				dMinDistance = dTemp;
// 				poTriangle = *liTriangles;
// 				oNearestPoint = oTempPoint;
// 			}
// 		}
// 	}
// 	return poTriangle;
// }
// bool ParadisSurface::IsPointInBox(const Point& oPoint)
// {
// 	double dX = oPoint.GetX();
// 	if((dX >= m_dXMin) && (dX <= m_dXMax))
// 	{
// 		double dY = oPoint.GetY();
// 		if((dY >= m_dYMin) && (dY <= m_dYMax))
// 		{
// 			double dZ = oPoint.GetZ();
// 			if((dZ >= m_dZMin) && (dZ <= m_dZMax))
// 			{
// 				return true;
// 			}
// 		}
// 	}
// 	return false;
// }
// bool ParadisSurface::IsPointInside(const Point& oPoint)
// {
// 	unsigned int iTestPointsCount = 5;
// 	double dAmplificationFactor = 1.5;
// 	unsigned int i = 0;
// 	unsigned int j = 0;
// 	Point oTestPoint;
// 	double dMaxDimension = max((m_dXMax - m_dXMin),max((m_dYMax -
// m_dYMin),(m_dZMax - m_dZMin))); 	double dXMean = 0.5*(m_dXMax + m_dXMin);
// 	double dYMean = 0.5*(m_dYMax + m_dYMin);
// 	double dZMean = 0.5*(m_dZMax + m_dZMin);
// 	Point oCenter(dXMean,dYMean,dZMean);
// 	double dRadius = dAmplificationFactor*dMaxDimension;
// 	double dTheta = 0.0;
// 	double dPhi = 0.0;
// 	list<TriPatch*>::const_iterator liTriangles;
// 	double dIntersectionParameter = 0.0;
// 	Point oIntersectionPoint;
// 	unsigned int iCount = 0;
// 	unsigned int iIsPointInsideCount = 0;
// 	for(i = 0 ; i < iTestPointsCount ; i++)
// 	{
// 		dPhi = Randomizer::Random(0,PI);
// 		dTheta = Randomizer::Random(0,2*PI);
// 		oTestPoint.Set(dRadius*sin(dPhi)*cos(dTheta),dRadius*sin(dPhi)*sin(dTheta),dRadius*cos(dPhi));
// 		oTestPoint = oTestPoint + oCenter;
// 		iCount = 0;
// 		for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 		{
// 			if((*liTriangles)->GetLineIntersection(oPoint,oTestPoint,oIntersectionPoint,dIntersectionParameter))
// 			{
// 				iCount = iCount + 1;
// 			}
// 		}
// 		if(iCount%2 != 0)
// 		{
// 			iIsPointInsideCount = iIsPointInsideCount + 1;
// 		}
// 	}
// 	if(2*iIsPointInsideCount > iTestPointsCount)
// 	{
// 		return true;
// 	}
// 	return false;
// }
// bool ParadisSurface::IsPointOnSurface(const Point& oPoint) const
// {
// 	list<TriPatch*>::const_iterator liTriangles;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		if((*liTriangles)->IsPointInTriangle(oPoint))
// 		{
// 			return true;
// 		}
// 	}
// 	return false;
// }
// double ParadisSurface::GetLeastSurfaceDistance(const Point& oPoint) const
// {
// 	list<TriPatch*>::const_iterator liTriangles;
// 	double dMin = DBL_MAX;
// 	double dTemp = 0.0;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		dTemp = (*liTriangles)->GetNormalDistance(oPoint);
// 		if(dTemp < dMin)
// 		{
// 			dMin = dTemp;
// 		}
// 	}
// 	return dMin;
// }
// void ParadisSurface::GenerateTriangulations()
// {
// 	GenericNode* poPointNNN = new GenericNode(m_dXMin,m_dYMin,m_dZMin);
// 	GenericNode* poPointPNN = new GenericNode(m_dXMax,m_dYMin,m_dZMin);
// 	GenericNode* poPointPPN = new GenericNode(m_dXMax,m_dYMax,m_dZMin);
// 	GenericNode* poPointNPN = new GenericNode(m_dXMin,m_dYMax,m_dZMin);
// 	GenericNode* poPointNNP = new GenericNode(m_dXMin,m_dYMin,m_dZMax);
// 	GenericNode* poPointPNP = new GenericNode(m_dXMax,m_dYMin,m_dZMax);
// 	GenericNode* poPointPPP = new GenericNode(m_dXMax,m_dYMax,m_dZMax);
// 	GenericNode* poPointNPP = new GenericNode(m_dXMin,m_dYMax,m_dZMax);
// 	m_lpoPoints.push_back(poPointNNN);
// 	m_lpoPoints.push_back(poPointPNN);
// 	m_lpoPoints.push_back(poPointPPN);
// 	m_lpoPoints.push_back(poPointNPN);
// 	m_lpoPoints.push_back(poPointNNP);
// 	m_lpoPoints.push_back(poPointPNP);
// 	m_lpoPoints.push_back(poPointPPP);
// 	m_lpoPoints.push_back(poPointNPP);
//
// 	m_lpoTriangles.push_back(new
// TriPatch(poPointNNN,poPointPNN,poPointNNP));
// m_lpoTriangles.push_back(new TriPatch(poPointPNN,poPointPNP,poPointNNP));
//
// 	m_lpoTriangles.push_back(new
// TriPatch(poPointNPN,poPointNPP,poPointPPP));
// m_lpoTriangles.push_back(new TriPatch(poPointNPN,poPointPPP,poPointPPN));
//
// 	m_lpoTriangles.push_back(new
// TriPatch(poPointNPN,poPointNNN,poPointNNP));
// m_lpoTriangles.push_back(new TriPatch(poPointNPN,poPointNNP,poPointNPP));
//
// 	m_lpoTriangles.push_back(new
// TriPatch(poPointPNN,poPointPPN,poPointPPP));
// m_lpoTriangles.push_back(new TriPatch(poPointPNN,poPointPPP,poPointPNP));
//
// 	m_lpoTriangles.push_back(new
// TriPatch(poPointNNN,poPointPPN,poPointPNN));
// m_lpoTriangles.push_back(new TriPatch(poPointNNN,poPointNPN,poPointPPN));
//
// 	m_lpoTriangles.push_back(new
// TriPatch(poPointPNP,poPointNPP,poPointNNP));
// m_lpoTriangles.push_back(new TriPatch(poPointPNP,poPointPPP,poPointNPP));
// }
// void ParadisSurface::ClearTriangulations()
// {
// 	list<TriPatch*>::iterator liTriangles;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		if((*liTriangles) != NULL)
// 		{
// 			delete (*liTriangles);
// 		}
// 	}
// 	m_lpoTriangles.clear();
//
// 	list<GenericNode*>::iterator liPoints;
// 	for(liPoints = m_lpoPoints.begin() ; liPoints != m_lpoPoints.end() ;
// liPoints++)
// 	{
// 		if((*liPoints) != NULL)
// 		{
// 			delete (*liPoints);
// 		}
// 	}
// 	m_lpoPoints.clear();
// }
