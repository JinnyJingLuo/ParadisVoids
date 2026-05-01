// Ahmed M. Hussein
#include "ParadisGrainSurfaceServer.h"
#include "list"
#include "float.h"
#include "Randomizer.h"
#include "mpi.h"
#include "Comm.h"
#include "QueueOps.h"
#include "map"
#include "Dodecahedron.h"

using namespace GeometrySystem;

// ParadisGrainSurfaceServer::ParadisGrainSurfaceServer()
// {
// 	Initialize();
// }
// ParadisGrainSurfaceServer::ParadisGrainSurfaceServer(const
// ParadisGrainSurfaceServer& oSurface)
// {
// 	*this = oSurface;
// }
// ParadisGrainSurfaceServer::~ParadisGrainSurfaceServer()
// {
// 	Reset();
// }
// ParadisGrainSurfaceServer& ParadisGrainSurfaceServer::operator=(const
// ParadisGrainSurfaceServer& oSurface)
// {
// 	Reset();
// 	ParadisGrainSurfaceServer::operator=(oSurface);
// 	m_dGrainDiameter = oSurface.m_dGrainDiameter;
// 	m_bIsHalfGrain = oSurface.m_bIsHalfGrain;
// 	// override triangulation
// 	ClearTriangulations();
// 	Dodecahedron oGrain;
// 	oGrain.Set(m_dGrainDiameter,m_bIsHalfGrain);
// 	m_lpoPoints = *oGrain.GetTriangulationPoints();
// 	m_lpoTriangles = *oGrain.GetTriangles();
// 	// the grain object has to release the triangulations so that the server
// can use them after
// 	// destroying the object when it goes out of scope
// 	oGrain.ReleaseTriangulations();
// 	return *this;
// }
// void ParadisGrainSurfaceServer::Reset()
// {
// 	ParadisSurface::Reset();
// 	Initialize();
// }
// void ParadisGrainSurfaceServer::Set(Home_t* poHome)
// {
// 	ParadisSurface::Set(poHome);
// 	if(poHome->param->GeometryType == 2)
// 	{
// 		m_bIsHalfGrain = false;
// 	}
// 	else if(poHome->param->GeometryType == 3)
// 	{
// 		m_bIsHalfGrain = true;
// 	}
// 	else
// 	{
// 		printf("error: unknown grain geometry type
// %d\n",poHome->param->GeometryType); 		return;
// 	}
// 	m_dGrainDiameter = poHome->param->GrainDiameter;
// 	// override triangulation
// 	ClearTriangulations();
// 	Dodecahedron oGrain;
// 	oGrain.Set(m_dGrainDiameter,m_bIsHalfGrain);
// 	m_lpoPoints = *oGrain.GetTriangulationPoints();
// 	m_lpoTriangles = *oGrain.GetTriangles();
// 	m_dVolume = oGrain.GetVolume();
//
// 	list<TriPatch*>::iterator liTriangles;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		if((*liTriangles)->IsConstrained())
// 		{
// 			printf("rigid\n");
// 		}
// 		else
// 		{
// 			printf("free\n");
// 		}
// 	}
//
// 	// the grain object has to release the triangulations so that the server
// can use them after
// 	// destroying the object when it goes out of scope
// 	oGrain.ReleaseTriangulations();
// 	// update paradis parameters
// 	poHome->param->simVol = m_dVolume;
// 	poHome->param->burgVolFactor = 1.0
// /(poHome->param->burgMag*poHome->param->burgMag*poHome->param->simVol);
// }
// void ParadisGrainSurfaceServer::CheckNodes(Home_t* poHome)
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
// distance 		bIsSurfaceNode = false; 		bIsSurfaceNode =
// !IsPointInside(oNodePoint); 		if(bIsSurfaceNode)
// 		{
// 			printf("out !!!!\n");
// 		}
// 		else
// 		{
// 			printf("in !!!!\n");
// 		}
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
// 			oSurfaceNormal = poTriangle->GetNormal();
// 			if(poTriangle->IsConstrained())
// 			{
// 				// for rigid surfaces, the node goes 1 unit into
// the
// volume 				daPos[0] = daPos[0] -
// oSurfaceNormal.GetX();
// daPos[1] = daPos[1] - oSurfaceNormal.GetY();
// daPos[2] = daPos[2] - oSurfaceNormal.GetZ();
// 			}
// 			else
// 			{
// 				// for free surfaces, the node stays on the
// surface
// and
// becomes a surface node 				poNode->constraint =
// SURFACE_NODE;
// poNode->dNx =
// oSurfaceNormal.GetX(); 				poNode->dNy =
// oSurfaceNormal.GetY(); 				poNode->dNz =
// oSurfaceNormal.GetZ();
// 			}
// 			RepositionNode(poHome,daPos,&(poNode->myTag),1);
// 		}
// 	}
//
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	CommSendRemesh(poHome);
// 	FixRemesh(poHome);
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	fflush(NULL);
// }
// void ParadisGrainSurfaceServer::Initialize()
// {
// 	ParadisSurface::Initialize();
// 	m_dGrainDiameter = 0.0;
// 	m_bIsHalfGrain = false;
// }
// bool ParadisGrainSurfaceServer::IsPointInside(const Point& oPoint)
// {
// 	list<TriPatch*>::const_iterator liTriangles;
// 	Vector oRelativePosition;
// 	double dTemp = 0.0;
// 	for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// 	{
// 		oRelativePosition.SetByPoints(*(*liTriangles)->GetPoint1(),oPoint);
// 		dTemp = oRelativePosition*(*liTriangles)->GetNormal();
// 		if(dTemp > m_dTolerance)
// 		{
// 			return false;
// 		}
// 	}
// 	return true;
//
// // 	unsigned int iTestPointsCount = 5;
// // 	double dAmplificationFactor = 1.5;
// // 	unsigned int i = 0;
// // 	unsigned int j = 0;
// // 	Point oTestPoint;
// // 	double dMaxDimension = max((m_dXMax - m_dXMin),max((m_dYMax -
// m_dYMin),(m_dZMax - m_dZMin)));
// // 	double dXMean = 0.5*(m_dXMax + m_dXMin);
// // 	double dYMean = 0.5*(m_dYMax + m_dYMin);
// // 	double dZMean = 0.5*(m_dZMax + m_dZMin);
// // 	Point oCenter(dXMean,dYMean,dZMean);
// // 	double dRadius = dAmplificationFactor*dMaxDimension;
// // 	double dTheta = 0.0;
// // 	double dPhi = 0.0;
// // 	list<TriPatch*>::const_iterator liTriangles;
// // 	double dIntersectionParameter = 0.0;
// // 	Point oIntersectionPoint;
// // 	unsigned int iCount = 0;
// // 	unsigned int iIsPointInsideCount = 0;
// // 	for(i = 0 ; i < iTestPointsCount ; i++)
// // 	{
// // 		dPhi = Randomizer::Random(0,PI);
// // 		dTheta = Randomizer::Random(0,2*PI);
// //
// oTestPoint.Set(dRadius*sin(dPhi)*cos(dTheta),dRadius*sin(dPhi)*sin(dTheta),dRadius*cos(dPhi));
// // 		oTestPoint = oTestPoint + oCenter;
// // 		iCount = 0;
// // 		for(liTriangles = m_lpoTriangles.begin() ; liTriangles !=
// m_lpoTriangles.end() ; liTriangles++)
// // 		{
// //
// if((*liTriangles)->GetLineIntersection(oPoint,oTestPoint,oIntersectionPoint,dIntersectionParameter))
// // 			{
// // 				iCount = iCount + 1;
// // 			}
// // 		}
// // 		printf("count %d\n",iCount);
// // 		if(iCount%2 != 0)
// // 		{
// // 			iIsPointInsideCount = iIsPointInsideCount + 1;
// // 		}
// // 	}
// // 	if(2*iIsPointInsideCount > iTestPointsCount)
// // 	{
// // 		return true;
// // 	}
// // 	return false;
// }
//
//
//
