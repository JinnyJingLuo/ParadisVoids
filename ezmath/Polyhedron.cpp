#include "Polyhedron.h"
#include "vector"
#include "float.h"
#include "Randomizer.h"
#include "stdlib.h"
#include "math.h"

using namespace std;
using namespace EZ;

namespace GeometrySystem
{
	Polyhedron::Polyhedron()
	{
		Initialize();
	}
	Polyhedron::Polyhedron(const Polyhedron& oPolyhedron)
	{
		*this = oPolyhedron;
	}
	Polyhedron::~Polyhedron()
	{
		Reset();
	}
	Polyhedron& Polyhedron::operator=(const Polyhedron& oPolyhedron)
	{
		Reset();
		// create a new copy of everything the polyhedron has
		vector< GenericNode* > vpoNodes;
		vector< Edge* > vpoEdges;
		vector< PolyhedronFace* > vpoFaces;
		
		// copy the nodes
		list<GenericNode*> loNodes = oPolyhedron.m_lpoNodes.GetList();
		list<GenericNode*>::const_iterator liNodes;
		vpoNodes.resize(loNodes.size());
		GenericNode* poOriginalNode = NULL;
		GenericNode* poNewNode = NULL;
		for(liNodes = loNodes.begin() ; liNodes != loNodes.end() ; liNodes++)
		{
			poOriginalNode = (*liNodes);
			poNewNode = new GenericNode(*poOriginalNode);
			poNewNode->SetID(poOriginalNode->GetID());
			vpoNodes[poOriginalNode->GetID() - 1] = poNewNode;
			AddNode(poNewNode);
		}

		// copy the edges and connect them to their nodes
		Edge* poOriginalEdge = NULL;
		Edge* poNewEdge = NULL;
		list< Edge* >::const_iterator liEdges;
		vpoEdges.resize(oPolyhedron.m_lpoEdges.size());
		for(liEdges = oPolyhedron.m_lpoEdges.begin() ; liEdges != oPolyhedron.m_lpoEdges.end() ; liEdges++)
		{
			poOriginalEdge = (*liEdges);
			poNewEdge = new Edge(*poOriginalEdge);
			vpoEdges[poOriginalEdge->GetID() - 1] = poNewEdge;
			poNewEdge->SetEndPoints(vpoNodes[poOriginalEdge->GetStartPoint()->GetID() - 1],vpoNodes[poOriginalEdge->GetEndPoint()->GetID() - 1]);
			AddEdge(poNewEdge);
		}

		// copy the faces and connect them to their nodes and edges
		PolyhedronFace* poOriginalFace = NULL;
		PolyhedronFace* poNewFace = NULL;
		list< PolyhedronFace* >::const_iterator liFaces;
		vpoFaces.resize(oPolyhedron.m_lpoFaces.size());
		for(liFaces = oPolyhedron.m_lpoFaces.begin() ; liFaces != oPolyhedron.m_lpoFaces.end() ; liFaces++)
		{
			poOriginalFace = (*liFaces);
			poNewFace = new PolyhedronFace(*poOriginalFace);
			vpoFaces[poOriginalFace->GetID() - 1] = poNewFace;
			poNewFace->SetPoints(vpoNodes[poOriginalFace->GetPoint1()->GetID() - 1],vpoNodes[poOriginalFace->GetPoint2()->GetID() - 1],vpoNodes[poOriginalFace->GetPoint3()->GetID() - 1]);
			poNewFace->SetEdges(vpoEdges[poOriginalFace->GetEdge1()->GetID() - 1],vpoEdges[poOriginalFace->GetEdge2()->GetID() - 1],vpoEdges[poOriginalFace->GetEdge3()->GetID() - 1]);
			AddFace(poNewFace);
		}

		// finally connect the edges to the faces
		for(liEdges = oPolyhedron.m_lpoEdges.begin() ; liEdges != oPolyhedron.m_lpoEdges.end() ; liEdges++)
		{
			poOriginalEdge = (*liEdges);
			poNewEdge = vpoEdges[poOriginalEdge->GetID() - 1];
			poNewEdge->SetRightFace(vpoFaces[poOriginalEdge->GetRightFace()->GetID() - 1]);
			poNewEdge->SetLeftFace(vpoFaces[poOriginalEdge->GetLeftFace()->GetID() - 1]);
		}

		// clear all the vectors and return
		vpoFaces.clear();
		vpoEdges.clear();
		vpoNodes.clear();
		m_bIsConvex = oPolyhedron.m_bIsConvex;
		m_iID = oPolyhedron.m_iID;
		return *this;
	}
	void Polyhedron::Reset()
	{
		list< PolyhedronFace* >::iterator liFaces;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces) != NULL)
			{
				delete (*liFaces);
				(*liFaces) = NULL;
			}
		}
		list< Edge* >::iterator liEdges;
		for(liEdges = m_lpoEdges.begin() ; liEdges != m_lpoEdges.end() ; liEdges++)
		{
			if((*liEdges) != NULL)
			{
				delete (*liEdges);
				(*liEdges) = NULL;
			}
		}
		unsigned int i = 0;
		unsigned int iSize = m_lpoNodes.GetSize();
		GenericNode* poNode = NULL;
		for(i = 0 ; i < iSize ; i++)
		{
			poNode = m_lpoNodes.GetCurrentItem();
			if(poNode != NULL)
			{
				delete poNode;
				poNode = NULL;
			}
			m_lpoNodes.IncrementIterator();
		}
		Initialize();
	}
	void Polyhedron::AddNode(GenericNode* poNode)
	{
		poNode->SetID(++m_iLastNodeID);
		m_lpoNodes.Append(poNode);
		// update the bounding box
		double dX = poNode->GetX();
		double dY = poNode->GetY();
		double dZ = poNode->GetZ();
		if(dX > m_oBox.GetXMax())
		{
			m_oBox.SetXMax(dX);
		}
		else if(dX < m_oBox.GetXMin())
		{
			m_oBox.SetXMin(dX);
		}

		if(dY > m_oBox.GetYMax())
		{
			m_oBox.SetYMax(dY);
		}
		else if(dY < m_oBox.GetYMin())
		{
			m_oBox.SetYMin(dY);
		}

		if(dZ > m_oBox.GetZMax())
		{
			m_oBox.SetZMax(dZ);
		}
		else if(dZ < m_oBox.GetZMin())
		{
			m_oBox.SetZMin(dZ);
		}
	}
	void Polyhedron::AddEdge(Edge* poEdge)
	{
		poEdge->SetID(++m_iLastEdgeID);
		m_lpoEdges.push_back(poEdge);
	}
	void Polyhedron::AddFace(PolyhedronFace* poFace)
	{
		poFace->SetID(++m_iLastFaceID);
		m_lpoFaces.push_back(poFace);
	}
	CircularLinkedList< GenericNode* >* Polyhedron::GetNodes()
	{
		return &m_lpoNodes;
	}
	list< Edge* >* Polyhedron::GetEdges()
	{
		return &m_lpoEdges;
	}
	list< PolyhedronFace* >* Polyhedron::GetFaces()
	{
		return &m_lpoFaces;
	}
	string Polyhedron::ToString()
	{
		string sString = "";

		sString = sString + "nodes : \n";
		//m_lpoNodes.ResetIterator();
		GenericNode* poNode = NULL;
		unsigned int iSize = m_lpoNodes.GetSize();
		unsigned int i = 0;
		for(i = 0 ; i < iSize ; i++)
		{
			poNode = m_lpoNodes.GetCurrentItem();
			sString = sString + poNode->ToString() + "\n";
			m_lpoNodes.IncrementIterator();
		}
		// return the iterator to its original state
		m_lpoNodes.IncrementIterator();

		sString = sString + "edges : \n";
		list< Edge* >::iterator liEdges;
		for(liEdges = m_lpoEdges.begin() ; liEdges != m_lpoEdges.end() ; liEdges++)
		{
			sString = sString + (*liEdges)->ToString() + "\n";
		}

		sString = sString + "faces : \n";
		list< PolyhedronFace* >::iterator liFaces;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			sString = sString + (*liFaces)->ToString() + "\n";
		}

		return sString;
	}
	bool Polyhedron::CreateAsHullFromPoints(list< Point >* ploPoints)
	{
		Reset();
		m_bHullConstructionError = false;
		// first, create the list of points that might end up in the polyhedron
		list< Point >::iterator liPoints;
		GenericNode* poNode = NULL;
		for(liPoints = ploPoints->begin() ; liPoints != ploPoints->end() ; liPoints++)
		{
			poNode = new GenericNode((*liPoints));
			poNode->SetCategory(FRESH_NODE);
			AddNode(poNode);
		}
		
		if(!CreateDihedron())
		{
			printf("error: couldn't create initial dihedron\n");
			return false;
		}
		if(!ConstructHull())
		{
			printf("error: couldn't construct hull\n");
			return false;
		}
		CleanUpNodes();
		UpdateIDs();
		m_bIsConvex = true;
		return true;
	}
	void Polyhedron::WriteNodes(FILE* fpFile)
	{
		unsigned int iNodesCount = m_lpoNodes.GetSize();
		fprintf(fpFile,"%d\n",iNodesCount);
		GenericNode* poNode = NULL;
		unsigned int i = 0;
		m_lpoNodes.ResetIterator();
		for(i = 0 ; i < iNodesCount ; i++)
		{
			poNode = m_lpoNodes.GetCurrentItem();
			fprintf(fpFile,"%E\t\t%E\t\t%E\n",poNode->GetX(),poNode->GetY(),poNode->GetZ());
			m_lpoNodes.IncrementIterator();
		}
	}
	void Polyhedron::WriteParaview(const string& sFileName)
	{
		FILE* fpFile = fopen(sFileName.c_str(),"w");

		fprintf(fpFile,"# vtk DataFile Version 1.0\n");
		fprintf(fpFile,"polyhedron\n");
		fprintf(fpFile,"ASCII\n");
		fprintf(fpFile,"DATASET POLYDATA\n");

		unsigned int iNodesCount = m_lpoNodes.GetSize();
		fprintf(fpFile,"POINTS %d float\n",iNodesCount);
		GenericNode* poNode = NULL;
		unsigned int i = 0;
		m_lpoNodes.ResetIterator();
		for(i = 0 ; i < iNodesCount ; i++)
		{
			poNode = m_lpoNodes.GetCurrentItem();
			fprintf(fpFile,"%E\t\t%E\t\t%E\n",poNode->GetX(),poNode->GetY(),poNode->GetZ());
			m_lpoNodes.IncrementIterator();
		}

//		unsigned int iEdgesCount = (unsigned int)m_lpoEdges.size();
//		fprintf(fpFile,"LINES %d %d\n",iEdgesCount,3*iEdgesCount);
//		list< Edge* >::iterator liEdges;
//		for(liEdges = m_lpoEdges.begin() ; liEdges != m_lpoEdges.end() ; liEdges++)
//		{
//			fprintf(fpFile,"2\t\t%d\t\t%d\n",(*liEdges)->GetStartPoint()->GetID() - 1,(*liEdges)->GetEndPoint()->GetID() - 1);
//		}

		unsigned int iFacesCount = (unsigned int)m_lpoFaces.size();
		fprintf(fpFile,"POLYGONS %d %d\n",iFacesCount,4*iFacesCount);
		list< PolyhedronFace* >::iterator liFaces;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			fprintf(fpFile,"3\t\t%d\t\t%d\t\t%d\n",(*liFaces)->GetPoint1()->GetID() - 1,(*liFaces)->GetPoint2()->GetID() - 1,(*liFaces)->GetPoint3()->GetID() - 1);
		}

		fprintf(fpFile,"POINT_DATA %d\n",iNodesCount);
		fprintf(fpFile,"scalars NodeType integer\n");
		fprintf(fpFile,"LOOKUP_TABLE default\n");
		for(i = 0 ; i < iNodesCount ; i++)
		{
			fprintf(fpFile,"%d\n",1);
		}

		//fprintf(fpFile,"CELL_DATA %d\n",iEdgesCount);
		fprintf(fpFile,"CELL_DATA %d\n",iFacesCount);
		fprintf(fpFile,"SCALARS SlipSystem integer\n");
		fprintf(fpFile,"LOOKUP_TABLE default\n");
		//for(liEdges = m_lpoEdges.begin() ; liEdges != m_lpoEdges.end() ; liEdges++)
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			fprintf(fpFile,"%d\n",1);
		}

		fclose(fpFile);
	}
	bool Polyhedron::IsPointInside(const Point& oPoint) const
	{
		if(!m_oBox.IsPointInside(oPoint))
		{
			return false;
		}
		// non-convex polyhedra are not supported yet
		if(!m_bIsConvex)
		{
			return false;
		}
		// if any of the faces is visible from the given point, then it is outside
		list< PolyhedronFace* >::const_iterator liFaces;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces)->IsVisibleFromPoint(oPoint))
			{
				return false;
			}
		}
		return true;
	}
	void Polyhedron::SetID(const unsigned int& iID)
	{
		m_iID = iID;
	}
	unsigned int Polyhedron::GetID() const
	{
		return m_iID;
	}
	double Polyhedron::GetDistance(const Point& oPoint) const
	{
		list< PolyhedronFace* >::const_iterator liFaces;
		double dMinDistance = DBL_MAX;
		double dDistance = 0.0;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			(*liFaces)->GetNearestPoint(oPoint,dDistance);
			if(dDistance < dMinDistance)
			{
				dMinDistance = dDistance;
			}
		}
		return dMinDistance;
	}
	double Polyhedron::GetSignedDistance(const Point& oPoint) const
	{
		double dMinDistance = GetDistance(oPoint);
		if(!IsPointInside(oPoint))			dMinDistance = -dMinDistance;
		return dMinDistance;
	}
	Point Polyhedron::GetNearestPoint(const Point& oPoint,double& dDistance) const
	{
		list< PolyhedronFace* >::const_iterator liFaces;
		Point oNearestPoint;
		Point oTempPoint;
		double dMinDistance = DBL_MAX;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			oTempPoint = (*liFaces)->GetNearestPoint(oPoint,dDistance);
			if(dDistance < dMinDistance)
			{
				dMinDistance = dDistance;
				oNearestPoint = oTempPoint;
			}
		}
		return oNearestPoint;
	}
	bool Polyhedron::GetNearestPointOnPlane(const Point& oPoint,const Plane& oPlane,Point& oNearestPoint,double& dDistance) const
	{
		list< PolyhedronFace* >::const_iterator liFaces;
		Point oTempPoint;
		double dMinDistance = DBL_MAX;
		bool bFound = false;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces)->GetNearestPointOnPlane(oPoint,oPlane,oTempPoint,dDistance))
			{
				if(dDistance < dMinDistance)
				{
					dMinDistance = dDistance;
					oNearestPoint = oTempPoint;
					bFound = true;
				}
			}
		}
		return bFound;
	}
	bool Polyhedron::GetNearestPointOnLine(const Point& oPoint,const Line& oLine,Point& oNearestPoint,double& dDistance) const
	{
		list< PolyhedronFace* >::const_iterator liFaces;
		Point oTempPoint;
		double dMinDistance = DBL_MAX;
		bool bFound = false;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces)->GetNearestPointOnLine(oPoint,oLine,oTempPoint,dDistance))
			{
				if(dDistance < dMinDistance)
				{
					dMinDistance = dDistance;
					oNearestPoint = oTempPoint;
					bFound = true;
				}
			}
		}
		return bFound;
	}
	Curve* Polyhedron::GetIntersectionCurve(const Plane& oPlane) const
	{
		list< PolyhedronFace* >::const_iterator liFaces;
		list< PolyhedronFace* > lpoIntersectingFaces;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces)->DoesPlaneCut(oPlane))		lpoIntersectingFaces.push_back((*liFaces));
		}
		list<Segment*> lpoSegments;
		Point oPoint1;
		Point oPoint2;
		Curve* poCurve = new Curve;
		for(liFaces = lpoIntersectingFaces.begin() ; liFaces != lpoIntersectingFaces.end() ; liFaces++)
		{
			if((*liFaces)->GetPlaneIntersection(oPlane,oPoint1,oPoint2))
			{
				poCurve->AddSegment(oPoint1,oPoint2);
			}
		}
		lpoIntersectingFaces.clear();
		// now sort the segments so that they form a closed curve
		poCurve->CyclicSort();
		return poCurve;
	}
	bool Polyhedron::IsIntersecting(const Plane& oPlane) const
	{
		list< PolyhedronFace* >::const_iterator liFaces;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces)->DoesPlaneCut(oPlane))		return true;
		}
		return false;
	}
	int Polyhedron::ClassifyPlane(const Plane& oPlane) const
	{
		// returns 0 if the precipitate intersects the plane, 1 if it is completely above
		// the plane and -1 if it is completely below it
		list<GenericNode*> loNodes = m_lpoNodes.GetList();
		list<GenericNode*>::const_iterator liNodes;
		bool bPointsAbove = false;
		bool bPointsBelow = false;
		int iResult = 0;
		for(liNodes = loNodes.begin() ; liNodes != loNodes.end() ; liNodes++)
		{
			iResult = oPlane.ClassifyPoint(*(*liNodes));
			if(iResult > 0)			bPointsAbove = true;
			if(iResult < 0)			bPointsBelow = true;
		}
		if(bPointsAbove && !bPointsBelow)		return 1;
		if(!bPointsAbove && bPointsBelow)		return -1;
		return 0;
	}
	void Polyhedron::Initialize()
	{
		m_lpoNodes.Reset();
		m_lpoEdges.clear();
		m_lpoFaces.clear();
		m_iLastNodeID = 0;
		m_iLastEdgeID = 0;
		m_iLastFaceID = 0;
		m_bIsConvex = false;
		m_bHullConstructionError = false;
		m_iID = 0;
		m_oBox.Reset();
	}
	bool Polyhedron::CreateDihedron()
	{
		// look for 3 non collinear points
		int iSize = m_lpoNodes.GetSize();
		GenericNode* poPreviousPoint = NULL;
		GenericNode* poThisPoint = NULL;
		GenericNode* poNextPoint = NULL;
		unsigned int i = 0;
		double dTolerance = 1.0E-6;
		bool bFound = false;
		for(i = 0 ; i < iSize ; i++)
		{
			poPreviousPoint = m_lpoNodes.GetPreviousItem();
			poThisPoint = m_lpoNodes.GetCurrentItem();
			poNextPoint = m_lpoNodes.GetNextItem();
			if((Vector(*poPreviousPoint,*poThisPoint)^Vector(*poThisPoint,*poNextPoint)).Length() > dTolerance)
			{
				bFound = true;
				break;
			}
			m_lpoNodes.IncrementIterator();
		}
		if(!bFound)
		{
			printf("error: all input points are collinear\n");
			return false;
		}

		PolyhedronFace* poFace1 = new PolyhedronFace;
		PolyhedronFace* poFace2 = new PolyhedronFace;
		Edge* poEdge1 = new Edge;
		Edge* poEdge2 = new Edge;
		Edge* poEdge3 = new Edge;
		// set the edges
		poEdge1->SetEndPoints(poPreviousPoint,poThisPoint);
		poEdge2->SetEndPoints(poThisPoint,poNextPoint);
		poEdge3->SetEndPoints(poNextPoint,poPreviousPoint);
		poEdge1->SetLeftFace(poFace1);
		poEdge1->SetRightFace(poFace2);
		poEdge2->SetLeftFace(poFace1);
		poEdge2->SetRightFace(poFace2);
		poEdge3->SetLeftFace(poFace1);
		poEdge3->SetRightFace(poFace2);
		poEdge1->SetCategory(INVISIBLE_EDGE);
		poEdge2->SetCategory(INVISIBLE_EDGE);
		poEdge3->SetCategory(INVISIBLE_EDGE);
		// add them to the hull
		AddEdge(poEdge1);
		AddEdge(poEdge2);
		AddEdge(poEdge3);
		// set the faces
		poFace1->SetPoints(poPreviousPoint,poThisPoint,poNextPoint);
		poFace2->SetPoints(poPreviousPoint,poNextPoint,poThisPoint);
		poFace1->SetEdges(poEdge1,poEdge2,poEdge3);
		poFace2->SetEdges(poEdge1,poEdge3,poEdge2);
		poFace1->SetCategory(INVISIBLE_FACE);
		poFace2->SetCategory(INVISIBLE_FACE);
		// add them to the hull
		AddFace(poFace1);
		AddFace(poFace2);
		// mark nodes as processed
		poPreviousPoint->SetCategory(PROCESSED_NODE);
		poThisPoint->SetCategory(PROCESSED_NODE);
		poNextPoint->SetCategory(PROCESSED_NODE);
		// finally, move the iterator of the points list to the first point that would give a nonzero volume with the
		// three points just obtained
		bFound = false;
		for(i = 0 ; i < iSize ; i++)
		{
			m_lpoNodes.IncrementIterator();
			poThisPoint = m_lpoNodes.GetCurrentItem();
			if(!poFace1->IsPointOnPlane(*poThisPoint))
			{
				bFound = true;
				break;
			}
		}
		if(!bFound)
		{
			printf("error: all input points are coplanar\n");
			return false;
		}
		return true;
	}
	bool Polyhedron::ConstructHull()
	{
		unsigned int iSize = m_lpoNodes.GetSize();
		unsigned int i = 0;
		GenericNode* poNode = NULL;
		// this part is very sensitive, if the iterator was slightly modified between the calls to 
		// the CreateDihedron() and this functions, everything will be messed up
		for(i = 0 ; i < iSize ; i++)
		{
			poNode = m_lpoNodes.GetCurrentItem();
			if(poNode->GetCategory() == PROCESSED_NODE)
			{
				m_lpoNodes.IncrementIterator();
				continue;
			}
			ProcessHullPoint(poNode);
			if(m_bHullConstructionError)
			{
				return false;
			}
			poNode->SetCategory(PROCESSED_NODE);
			CleanUp();
			m_lpoNodes.IncrementIterator();
		}
		return true;
	}
	bool Polyhedron::ProcessHullPoint(GenericNode* poPoint)
	{
		// 1. compile a list of the visible faces for this node, if there are no visible faces, then the
		// point is an internal node. a visible face is a face which the point is in its positive half space
		// all the current faces are invisible initially due to the clean up process that follows the call to
		// this function
		bool bIsExternal = false;
		list< PolyhedronFace* >::iterator liFaces;
		// check all the faces for visibility
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces)->IsVisibleFromPoint(*poPoint))
			{
				bIsExternal = true;
				(*liFaces)->SetCategory(VISIBLE_FACE);
			}
		}

		// if the point is internal, remove it
		if(!bIsExternal)
		{
			return false;
		}

		// detect the internal and the boundary edges, notice that some edges are neither internal nor boundary
		list< Edge* >::iterator liEdges;
		int iRightFaceCategory = 0;
		int iLeftFaceCategory = 0;
		// create an array to hold the nodes that were used in the cone faces creation, the C type allocations
		// will be used for efficiency
		Edge** ppoNodeNewEdgeIncidence = (Edge**)calloc(m_iLastNodeID,sizeof(Edge*));
		PolyhedronFace* poNewFace = NULL;
		Edge* poNewEdge1 = NULL;
		Edge* poNewEdge2 = NULL;
		GenericNode* poEdgeStartNode = NULL;
		GenericNode* poEdgeEndNode = NULL;
		PolyhedronFace* poHiddenFace = NULL;
		for(liEdges = m_lpoEdges.begin() ; liEdges != m_lpoEdges.end() ; liEdges++)
		{
			Face* poRightFace = (*liEdges)->GetRightFace();
			Face* poLeftFace = (*liEdges)->GetLeftFace();
			if((poRightFace == NULL) || (poLeftFace == NULL))
			{
				m_bHullConstructionError = true;
				free(ppoNodeNewEdgeIncidence);
				return false;
			}
			iRightFaceCategory = poRightFace->GetCategory();
			iLeftFaceCategory = poLeftFace->GetCategory();
			if((iRightFaceCategory == VISIBLE_FACE) && (iLeftFaceCategory == VISIBLE_FACE))
			{
				// both faces are visible, edge is internal
				(*liEdges)->SetCategory(INTERNAL_EDGE);
			}
			else if((iRightFaceCategory == VISIBLE_FACE) || (iLeftFaceCategory == VISIBLE_FACE))
			{
				if(iRightFaceCategory == VISIBLE_FACE)
				{
					poHiddenFace = (PolyhedronFace*)(*liEdges)->GetRightFace();
				}
				else
				{
					poHiddenFace = (PolyhedronFace*)(*liEdges)->GetLeftFace();
				}
				// only one side is visible, the edge is a boundary edge
				(*liEdges)->SetCategory(BOUNDARY_EDGE);
				poEdgeStartNode = (*liEdges)->GetStartPoint();
				poEdgeEndNode = (*liEdges)->GetEndPoint();
				// create a new face between the point in question and this edge
				// 1. see if the first node in the boundary edge has a new incident edge, if so, use it
				if(ppoNodeNewEdgeIncidence[poEdgeStartNode->GetID() - 1] != NULL)
				{
					poNewEdge1 = ppoNodeNewEdgeIncidence[poEdgeStartNode->GetID() - 1];
				}
				else
				{
					poNewEdge1 = new Edge;
					poNewEdge1->SetEndPoints(poPoint,poEdgeStartNode);
					poNewEdge1->SetCategory(INVISIBLE_EDGE);
					AddEdge(poNewEdge1);
					ppoNodeNewEdgeIncidence[poEdgeStartNode->GetID() - 1] = poNewEdge1;
				}
				// 2. see if the first node in the boundary edge has a new incident edge, if so, use it
				if(ppoNodeNewEdgeIncidence[poEdgeEndNode->GetID() - 1] != NULL)
				{
					poNewEdge2 = ppoNodeNewEdgeIncidence[poEdgeEndNode->GetID() - 1];
				}
				else
				{
					poNewEdge2 = new Edge;
					poNewEdge2->SetEndPoints(poPoint,poEdgeEndNode);
					poNewEdge2->SetCategory(INVISIBLE_EDGE);
					AddEdge(poNewEdge2);
					ppoNodeNewEdgeIncidence[poEdgeEndNode->GetID() - 1] = poNewEdge2;
				}

				poNewFace = new PolyhedronFace;
				if(poHiddenFace->IsEdgeAligned((*liEdges)))
				{
					poNewFace->SetPoints(poPoint,poEdgeStartNode,poEdgeEndNode);
				}
				else
				{
					poNewFace->SetPoints(poPoint,poEdgeEndNode,poEdgeStartNode);
				}

				poNewFace->SetEdges(poNewEdge1,(*liEdges),poNewEdge2);
				poNewFace->SetCategory(INVISIBLE_FACE);
				AddFace(poNewFace);

				// set the edges, no need to check which faces are set, everything is taken care of because
				// if the edge was previously created, the correct faces had been already set at the time
				if(poHiddenFace->IsEdgeAligned((*liEdges)))
				{
					poNewEdge1->SetLeftFace(poNewFace);
					poNewEdge2->SetRightFace(poNewFace);
				}
				else
				{
					poNewEdge1->SetRightFace(poNewFace);
					poNewEdge2->SetLeftFace(poNewFace);
				}

				// finally, update the right or left faces of the boundary edge and mark it as invisible
				(*liEdges)->ReplaceFace(poHiddenFace,poNewFace);
				(*liEdges)->SetCategory(INVISIBLE_EDGE);
			}
			else
			{
				// the edge is invisible
				(*liEdges)->SetCategory(INVISIBLE_EDGE);
			}
		}
		free(ppoNodeNewEdgeIncidence);
		return true;
	}
	void Polyhedron::CleanUp()
	{
		list< PolyhedronFace* >::iterator liFaces = m_lpoFaces.begin();
		while(liFaces != m_lpoFaces.end())
		{
			if((*liFaces)->GetCategory() == VISIBLE_FACE)
			{
				delete (*liFaces);
				liFaces = m_lpoFaces.erase(liFaces);
			}
			else
			{
				liFaces++;
			}
		}

		list< Edge* >::iterator liEdges = m_lpoEdges.begin();
		while(liEdges != m_lpoEdges.end())
		{
			if((*liEdges)->GetCategory() == INTERNAL_EDGE)
			{
				delete (*liEdges);
				liEdges = m_lpoEdges.erase(liEdges);
			}
			else
			{
				liEdges++;
			}
		}
	}
	void Polyhedron::CleanUpNodes()
	{
		list< PolyhedronFace* >::iterator liFaces;
		char* piNodeUsage = (char*)calloc(m_iLastNodeID,sizeof(char));
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			piNodeUsage[(*liFaces)->GetPoint1()->GetID() - 1] = 1;
			piNodeUsage[(*liFaces)->GetPoint2()->GetID() - 1] = 1;
			piNodeUsage[(*liFaces)->GetPoint3()->GetID() - 1] = 1;
		}

		unsigned int i = 0;
		unsigned int iSize = m_lpoNodes.GetSize();
		GenericNode* poNode = NULL;
		for(i = 0 ; i < iSize ; i++)
		{
			poNode = m_lpoNodes.GetCurrentItem();
			if(piNodeUsage[poNode->GetID() - 1] == 0)
			{
				delete poNode;
				m_lpoNodes.DropItem();
				m_lpoNodes.DecrementIterator();
			}
			m_lpoNodes.IncrementIterator();
		}
		free(piNodeUsage);
	}
	void Polyhedron::UpdateIDs()
	{
		list< PolyhedronFace* >::iterator liFaces;
		m_iLastFaceID = 0;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			(*liFaces)->SetID(++m_iLastFaceID);
		}
		list< Edge* >::iterator liEdges;
		m_iLastEdgeID = 0;
		for(liEdges = m_lpoEdges.begin() ; liEdges != m_lpoEdges.end() ; liEdges++)
		{
			(*liEdges)->SetID(++m_iLastEdgeID);
		}

		unsigned int i = 0;
		unsigned int iSize = m_lpoNodes.GetSize();
		m_iLastNodeID = 0;
		m_lpoNodes.ResetIterator();
		for(i = 0 ; i < iSize ; i++)
		{
			m_lpoNodes.GetCurrentItem()->SetID(++m_iLastNodeID);
			m_lpoNodes.IncrementIterator();
		}
	}
	bool Polyhedron::DoesSegmentPierce(const Point& oStartPoint,const Point& oEndPoint) const
	{
		list< PolyhedronFace* >::const_iterator liFaces;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			if((*liFaces)->DoesSegmentIntersect(oStartPoint,oEndPoint))
			{
				return true;
			}
		}
		return false;
	}
	const AxisAlignedBoundingBox* Polyhedron::GetBox() const
	{
		return &m_oBox;
	}
	bool Polyhedron::IsIntersecting(Polyhedron& oPolyhedron)
	{
		// this check is not robust because high aspect ratio polyhedra can intersect without
		// having points inside each other, fix that
		// loop over this polyhedron nodes and check whether any of them is inside the other polyhedron
		unsigned int i = 0;
		unsigned int iSize = m_lpoNodes.GetSize();
		GenericNode* poNode = NULL;
		for(i = 0 ; i < iSize ; i++)
		{
			poNode = m_lpoNodes.GetCurrentItem();
			if(oPolyhedron.IsPointInside(*poNode))
			{
				return true;
			}
			m_lpoNodes.IncrementIterator();
		}
		// none of this polyhedron points is inside the other polyhedron, check the other case
		iSize = oPolyhedron.m_lpoNodes.GetSize();
		for(i = 0 ; i < iSize ; i++)
		{
			poNode = oPolyhedron.m_lpoNodes.GetCurrentItem();
			if(IsPointInside(*poNode))
			{
				return true;
			}
			oPolyhedron.m_lpoNodes.IncrementIterator();
		}
		// none of the other polyhedron's points in inside this one, then they are not intersecting
		return false;
	}
	bool Polyhedron::Contains(Polyhedron& oPolyhedron)
	{
		// all of the points of the other polyhedron should lie inside this polyhedron
		unsigned int iSize = oPolyhedron.m_lpoNodes.GetSize();
		unsigned int i = 0;
		GenericNode* poNode = NULL;
		for(i = 0 ; i < iSize ; i++)
		{
			poNode = oPolyhedron.m_lpoNodes.GetCurrentItem();
			if(!IsPointInside(*poNode))
			{
				return false;
			}
			oPolyhedron.m_lpoNodes.IncrementIterator();
		}
		return true;
	}
	Polyhedron* Polyhedron::PlaneCut(const Plane& oPlane,const bool& bPerturb,const double& dPerturbationDistance)
	{
		// this function currently does not support non-convex polyhedra
		Polyhedron* poCutPolyhedron = NULL;
		if(!m_bIsConvex)				return poCutPolyhedron;
		// the cut polyhedron contains all the original polyhedron points which are on the negative 
		// side of the plane in addition to the points on the intersection curve of the plane and
		// the original polyhedron, start by building a list of these points
		vector<Point> voNewPoints;
		unsigned int iSize = m_lpoNodes.GetSize();
		// just an estimate for the final size
		voNewPoints.reserve(iSize);
		
		
		// 1. the points below the cutting plane
		unsigned int i = 0;
		for(i = 0 ; i < iSize ; i++)
		{
			if(!oPlane.IsPointOnOrAbove(*m_lpoNodes.GetCurrentItem()))
			{
				voNewPoints.push_back(Point(*m_lpoNodes.GetCurrentItem()));
			}
			m_lpoNodes.IncrementIterator();
		}
		
		// 2. the points on the intersection curve of that plane with the cutting plane
		Curve* poIntersectionCurve = GetIntersectionCurve(oPlane);
		list<GenericNode*> lpoIntersectionPoints = poIntersectionCurve->GetPoints();
		list<GenericNode*>::iterator liPoints;
		Point oPerturbedPoint;
		Vector oNormal = oPlane.GetNormal();
		double dShiftDistance = 0.0;
		for(liPoints = lpoIntersectionPoints.begin() ; liPoints != lpoIntersectionPoints.end() ; liPoints++)
		{
			oPerturbedPoint = Point(*(*liPoints));
			if(bPerturb)
			{
				dShiftDistance = fabs(Randomizer::RandomNormal(dPerturbationDistance,0.3*dPerturbationDistance));
				oPerturbedPoint = oPerturbedPoint + oNormal*dShiftDistance;
			}
			voNewPoints.push_back(oPerturbedPoint);
		}
		delete poIntersectionCurve;
		
		// for best results, shuffle the points list before creating new polyhedron
		iSize = (unsigned int)voNewPoints.size();
		
		Point oTempPoint;
		unsigned int iRandomIndex = 0;
		for(i = 0 ; i < iSize ; i++)
		{
			iRandomIndex = Randomizer::RandomInteger(0,iSize - 1);
			oTempPoint = voNewPoints[iRandomIndex];
			voNewPoints[iRandomIndex] = voNewPoints[i];
			voNewPoints[i] = oTempPoint;
		}
		
		// pack the nodes in a list
		list<Point> loNewPoints;
		for(i = 0 ; i < iSize ; i++)
		{
			loNewPoints.push_back(voNewPoints[i]);
		}
		voNewPoints.clear();

		// build the new polyhedron
		poCutPolyhedron = new Polyhedron;
		poCutPolyhedron->CreateAsHullFromPoints(&loNewPoints);
		loNewPoints.clear();
		return poCutPolyhedron;
	}
	Point Polyhedron::GenerateInternalPoint() const
	{
		Point oPoint;
		while(true)
		{
			oPoint = m_oBox.GenerateRandomPoint();
			if(IsPointInside(oPoint))		break;
		}
		return oPoint;
	}
	double Polyhedron::GetVolume() const
	{
		// compute the volume using surface integration
		list< PolyhedronFace* >::const_iterator liFaces;
		double dVolume = 0.0;
		for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
		{
			dVolume += (*liFaces)->GetArea()*(*liFaces)->GetCentroid().GetX()*(*liFaces)->GetNormal().GetX();
		}
		return dVolume;
	}
}

