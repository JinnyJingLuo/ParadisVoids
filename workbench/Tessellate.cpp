#include "stdio.h"
#include "AxisAlignedBoundingBox.h"
#include "GenericNode.h"
#include "Polyhedron.h"
#include "PlanarPolygon.h"
#include "math.h"

using namespace GeometrySystem;


/*
void WriteVoronoiParaview(const string& sFileName,vector< GenericNode* >* pvpoVertices,vector< Edge* >* pvpoEdges)
{
	FILE* fpFile = fopen(sFileName.c_str(),"w");

	fprintf(fpFile,"# vtk DataFile Version 1.0\n");
	fprintf(fpFile,"Voronoi tessellation\n");
	fprintf(fpFile,"ASCII\n");
	fprintf(fpFile,"DATASET POLYDATA\n");

	unsigned int iNodesCount = (unsigned int)pvpoVertices->size();
	fprintf(fpFile,"POINTS %d float\n",iNodesCount);
	GenericNode* poNode = NULL;
	unsigned int i = 0;
	for(i = 0 ; i < iNodesCount ; i++)
	{
		poNode = pvpoVertices->at(i);
		fprintf(fpFile,"%E\t\t%E\t\t%E\n",poNode->GetX(),poNode->GetY(),poNode->GetZ());
	}

	unsigned int iEdgesCount = (unsigned int)pvpoEdges->size();
	fprintf(fpFile,"LINES %d %d\n",iEdgesCount,3*iEdgesCount);
	Edge* poEdge = NULL;
	for(i = 0 ; i < iEdgesCount ; i++)
	{
		poEdge = pvpoEdges->at(i);
		fprintf(fpFile,"2\t\t%d\t\t%d\n",poEdge->GetStartPoint()->GetID() - 1,poEdge->GetEndPoint()->GetID() - 1);
	}

// 	unsigned int iFacesCount = (unsigned int)m_lpoFaces.size();
// 	fprintf(fpFile,"POLYGONS %d %d\n",iFacesCount,4*iFacesCount);
// 	list< PolyhedronFace* >::iterator liFaces;
// 	for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
// 	{
// 		fprintf(fpFile,"3\t\t%d\t\t%d\t\t%d\n",(*liFaces)->GetPoint1()->GetID() - 1,(*liFaces)->GetPoint2()->GetID() - 1,(*liFaces)->GetPoint3()->GetID() - 1);
// 	}

	fprintf(fpFile,"POINT_DATA %d\n",iNodesCount);
	fprintf(fpFile,"scalars NodeType integer\n");
	fprintf(fpFile,"LOOKUP_TABLE default\n");
	for(i = 0 ; i < iNodesCount ; i++)
	{
		fprintf(fpFile,"%d\n",1);
	}

	fprintf(fpFile,"CELL_DATA %d\n",iEdgesCount);
// 	fprintf(fpFile,"CELL_DATA %d\n",iFacesCount);
	fprintf(fpFile,"SCALARS SlipSystem integer\n");
	fprintf(fpFile,"LOOKUP_TABLE default\n");
	for(i = 0 ; i < iEdgesCount ; i++)
// 	for(liFaces = m_lpoFaces.begin() ; liFaces != m_lpoFaces.end() ; liFaces++)
	{
		fprintf(fpFile,"%d\n",1);
	}

	fclose(fpFile);
}
	
void Tessalate(const list< GenericNode* >& lpoNodes)
{
	Polyhedron oPolyhedron;
	list< Point > loPoints;
	list< GenericNode* >::const_iterator liNodes;
	Point oTemp;
	double dX = 0.0;
	double dY = 0.0;
	double dZ = 0.0;
	for(liNodes = lpoNodes.begin() ; liNodes != lpoNodes.end() ; liNodes++)
	{
		oTemp = *(*liNodes);
		dX = oTemp.GetX();
		dY = oTemp.GetY();
		dZ = dX*dX + dY*dY;
		oTemp.Set(dX,dY,dZ);
		loPoints.push_back(oTemp);
	}
	printf("before creation\n");
	if(!oPolyhedron.CreateAsHullFromPoints(&loPoints))
	{
		printf("error: convex hull generation failed\n");
	}
	printf("after creation\n");
	oPolyhedron.WriteParaview("polyhedron.vtk");
	oPolyhedron.CutTop(Vector(0.0,0.0,1.0));
	oPolyhedron.WriteParaview("base.vtk");
	oPolyhedron.Project(Plane(Vector(0.0,0.0,1.0),Point(0.0,0.0,0.0)));
	oPolyhedron.WriteParaview("tessalation.vtk");
	
	// now the polyhedron is just a Delaunay triangulation, get its dual Voronoi diagram
	
	// 1. loop over all the vertices and create their corresponding Voronoi polygons that will 
	// be set later with Voronoi edges
	printf("creating cells\n");
	CircularLinkedList< GenericNode* >* plpoNodes = oPolyhedron.GetNodes();
	vector<PlanarPolygon*> vpoVoronoiCells;
	unsigned int iSize = plpoNodes->GetSize();
	unsigned int i = 0;
	vpoVoronoiCells.resize(iSize);
	unsigned int iID = 0;
	PlanarPolygon* poCell = NULL;
	for(i = 0 ; i < iSize ; i++)
	{
		iID = plpoNodes->GetCurrentItem()->GetID();
		vpoVoronoiCells[iID - 1] = new PlanarPolygon;
		vpoVoronoiCells[iID - 1]->SetID(iID);
		plpoNodes->IncrementIterator();
	}
	
	// 2. loop over all the triangles and get their circumcenters
	printf("creating vertices\n");
	list< PolyhedronFace* >* plpoFaces = oPolyhedron.GetFaces();
	list< PolyhedronFace* >::iterator liFaces;
	vector< GenericNode* > vpoVoronoiVertices;
	vpoVoronoiVertices.resize(plpoFaces->size());
	for(liFaces = plpoFaces->begin() ; liFaces != plpoFaces->end() ; liFaces++)
	{
		iID = (*liFaces)->GetID();
		vpoVoronoiVertices[iID - 1] = new GenericNode((*liFaces)->GetCircumcenter());
		vpoVoronoiVertices[iID - 1]->SetID(iID);
	}
	
	// 3. get the Delaunay edges and form the Voronoi edges from them
	printf("creating edges\n");
	list< Edge* >* plpoEdges = oPolyhedron.GetEdges();
	list< Edge* >::iterator liEdges;
	vector< Edge* > vpoVoronoiEdges;
	vpoVoronoiEdges.reserve(plpoEdges->size());
	PolyhedronFace* poLeftFace = NULL;
	PolyhedronFace* poRightFace = NULL;
	GenericNode* poStartNode = NULL;
	GenericNode* poEndNode = NULL;
	unsigned int iStartID = 0;
	unsigned int iEndID = 0;
	PlanarPolygon* poLeftCell = NULL;
	PlanarPolygon* poRightCell = NULL;
	iID = 0;
	Edge* poEdge = NULL;
	Vector oTestVector;
	for(liEdges = plpoEdges->begin() ; liEdges != plpoEdges->end() ; liEdges++)
	{
		poLeftFace = (PolyhedronFace*)(*liEdges)->GetLeftFace();
		poRightFace = (PolyhedronFace*)(*liEdges)->GetRightFace();
		if((poLeftFace == NULL) || (poRightFace == NULL))
		{
			continue;
		}
		poStartNode = (*liEdges)->GetStartPoint();
		poEndNode = (*liEdges)->GetStartPoint();
		poEdge = new Edge;
		poEdge->SetEndPoints(vpoVoronoiVertices[poLeftFace->GetID() - 1],vpoVoronoiVertices[poRightFace->GetID() - 1]);
		poEdge->SetID(++iID);
		// now we need to set the edge cells, assume that the cell of the start point is the left 
		// cell, check, and switch the cells if the check fails
		poLeftCell = vpoVoronoiCells[poStartNode->GetID() - 1];
		poRightCell = vpoVoronoiCells[poEndNode->GetID() - 1];
		// note that all of the faces have the same normal at this point
		oTestVector.SetByPoints(*vpoVoronoiVertices[poLeftFace->GetID() - 1],*poStartNode);
		if(oTestVector*(poLeftFace->GetNormal()^poEdge->GetVector()) < 0.0)
		{
			poLeftCell = vpoVoronoiCells[poEndNode->GetID() - 1];
			poRightCell = vpoVoronoiCells[poStartNode->GetID() - 1];
		}
		// now add the edge to its cells
		poLeftCell->AddEdge(poEdge);
		poRightCell->AddEdge(poEdge);
		vpoVoronoiEdges.push_back(poEdge);
	}
	
	// at this point, we have already set all the edges in their cells, finalize them
	// the cells that fail to finalize are open (infinite cells) and will be discarded
	printf("finalizing cells\n");
	iSize = (unsigned int)vpoVoronoiCells.size();
	for(i = 0 ; i < iSize ; i++)
	{
		if(!vpoVoronoiCells[i]->FinalizeEdges())
		{
			// the cell is open, remove it
			printf("cell %d failed\n",vpoVoronoiCells[i]->GetID());
			vpoVoronoiCells[i]->Detach();
			
			delete vpoVoronoiCells[i];
			vpoVoronoiCells[i] = NULL;
		}
	}
	printf("writing\n");
	WriteVoronoiParaview("voronoi.vtk",&vpoVoronoiVertices,&vpoVoronoiEdges);
}
*/

int main(int argc,char** argv)
{
	printf("tessellator running\n");
	AxisAlignedBoundingBox oBox;
	oBox.SetXMin(-1000.0);
	oBox.SetXMax(1000.0);
	oBox.SetYMin(-1000.0);
	oBox.SetYMax(1000.0);
	oBox.SetZMin(-1000.0);
	oBox.SetZMax(1000.0);
	
	unsigned int iSize = 100;
	unsigned int i = 0;
	list< GenericNode* > lpoNodes;

	for(i = 0 ; i < iSize ; i++)
	{
		lpoNodes.push_back(new GenericNode(oBox.GenerateRandomPoint(1.0)));
	}
	//Tessalate(lpoNodes);
	list< GenericNode* >::iterator liNodes;
	for(liNodes = lpoNodes.begin() ; liNodes != lpoNodes.end() ; liNodes++)
	{
		delete (*liNodes);
	}
	return 0;
}


