// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#include "Boundary.h"
#include "Tools.h"
#include "float.h"
#include "stdio.h"

using namespace SupportSystem;


namespace GeometrySystem
{
	Boundary::Boundary()
	{
		Initialize();
	}
	Boundary::~Boundary()
	{
		Reset();
	}
	Boundary::Boundary(const Boundary& oBoundary)
	{
		*this = oBoundary;
	}
	Boundary& Boundary::operator=(const Boundary& oBoundary)
	{
		Reset();
		// make an entirely new copy of the input structure
		list<GenericNode*>* plpoPoints = &(((Boundary&)oBoundary).m_lpoTriangulationPoints);
		list<TriPatch*>* plpoTriangles = &(((Boundary&)oBoundary).m_lpoTriangles);
		Set(plpoPoints,plpoTriangles);
		return *this;
	}
	void Boundary::Reset()
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
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			if((*liPoints) != NULL)
			{
				delete (*liPoints);
			}
		}
		m_lpoTriangulationPoints.clear();
	}
	void Boundary::Set(const string& sFileName)
	{
		FILE* fpFile = fopen(sFileName.c_str(),"r");
		unsigned int iSize = 0;
		string sRead;
		sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%d\n",&iSize);
		unsigned int i = 0;
		vector<GenericNode*> vpoPoints;
		vpoPoints.resize(iSize);
		double dX = 0.0;
		double dY = 0.0;
		double dZ = 0.0;
		// read the points
		for(i = 0 ; i < iSize ; i++)
		{
			sRead = GetRealString(500,fpFile);
			sscanf(sRead.c_str(),"%lf\t\t%lf\t\t%lf\n",&dX,&dY,&dZ);
			vpoPoints[i] = new GenericNode(dX,dY,dZ);
			vpoPoints[i]->SetID(i + 1);
		}
		// read the triangles
		unsigned int iPoint1Index = 0;
		unsigned int iPoint2Index = 0;
		unsigned int iPoint3Index = 0;
		sRead = GetRealString(500,fpFile);
		sscanf(sRead.c_str(),"%d\n",&iSize);
		for(i = 0 ; i < iSize ; i++)
		{
			sRead = GetRealString(500,fpFile);
			sscanf(sRead.c_str(),"%d\t\t%d\t\t%d\n",&iPoint1Index,&iPoint2Index,&iPoint3Index);
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[iPoint1Index - 1],vpoPoints[iPoint2Index - 1],vpoPoints[iPoint3Index - 1]));
		}
		fclose(fpFile);
		iSize = vpoPoints.size();
		for(i = 0 ; i < iSize ; i++)
		{
			m_lpoTriangulationPoints.push_back(vpoPoints[i]);
		}
		vpoPoints.clear();
	}
	void Boundary::Set(list<GenericNode*>* plpoPoints,list<TriPatch*>* plpoTriangles)
	{
		list<GenericNode*>::iterator liPoints;
		unsigned int i = 0;
		unsigned int iSize = plpoPoints->size();
		vector<GenericNode*> vpoPoints;
		vpoPoints.resize(iSize);
		vector<unsigned int> viOriginalIndices;
		viOriginalIndices.resize(iSize);
		for(liPoints = plpoPoints->begin() ; liPoints != plpoPoints->end() ; liPoints++)
		{
			vpoPoints[i] = new GenericNode((*(*liPoints)));
			vpoPoints[i]->SetID(i + 1);
			viOriginalIndices[i] = (*liPoints)->GetID();
			(*liPoints)->SetID(i + 1);
			i = i + 1;
		}
	
		unsigned int iPoint1Index = 0;
		unsigned int iPoint2Index = 0;
		unsigned int iPoint3Index = 0;
		list<TriPatch*>::iterator liTriangles;
		for(liTriangles = plpoTriangles->begin() ; liTriangles != plpoTriangles->end() ; liTriangles++)
		{
			iPoint1Index = (*liTriangles)->GetPoint1()->GetID() - 1;
			iPoint2Index = (*liTriangles)->GetPoint2()->GetID() - 1;
			iPoint3Index = (*liTriangles)->GetPoint3()->GetID() - 1;
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[iPoint1Index],vpoPoints[iPoint2Index],vpoPoints[iPoint3Index]));
		}
		iSize = vpoPoints.size();
		liPoints = plpoPoints->begin();
		for(i = 0 ; i < iSize ; i++)
		{
			vpoPoints[i]->SetID(viOriginalIndices[i]);
			(*liPoints)->SetID(viOriginalIndices[i]);
			liPoints++;
			m_lpoTriangulationPoints.push_back(vpoPoints[i]);
		}
		vpoPoints.clear();
		viOriginalIndices.clear();
	}
	void Boundary::Initialize()
	{
		m_lpoTriangles.clear();
		m_lpoTriangulationPoints.clear();
	}
	list<TriPatch*>* Boundary::GetTriangles()
	{
		return &m_lpoTriangles;
	}
	list<GenericNode*>* Boundary::GetTriangulationPoints()
	{
		return &m_lpoTriangulationPoints;
	}
	void Boundary::Refine(const unsigned int& iRefinementsCount)
	{
		unsigned int i = 0;
		for(i = 0; i < iRefinementsCount ; i++)
		{
			RefineTriangulation();
		}
	}
	void Boundary::WriteTriangles(FILE* fpFile)
	{
		unsigned int iSize = m_lpoTriangulationPoints.size();
		char cWrite[500];
		sprintf(cWrite,"%d\n",iSize);
		fputs(cWrite,fpFile);
		list<GenericNode*>::iterator liPoints;
		unsigned int i = 0;
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			sprintf(cWrite,"%d\t\t%E\t\t%E\t\t%E\n",i,(*liPoints)->GetX(),(*liPoints)->GetY(),(*liPoints)->GetZ());
			fputs(cWrite,fpFile);
			i = i + 1;
		}
		
		iSize = m_lpoTriangles.size();
		sprintf(cWrite,"%d\n",iSize);
		fputs(cWrite,fpFile);
		list<TriPatch*>::iterator liTriangles;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			sprintf(cWrite,"%d\t\t%d\t\t%d\n",(*liTriangles)->GetPoint1()->GetID() - 1,(*liTriangles)->GetPoint2()->GetID() - 1,(*liTriangles)->GetPoint3()->GetID() - 1);
			fputs(cWrite,fpFile);
		}
	}
	void Boundary::WriteTriangles(const string& sFileName)
	{
		FILE* fpFile = fopen(sFileName.c_str(),"w");
		WriteTriangles(fpFile);
		fclose(fpFile);
	}
	void Boundary::WriteParaviewTriangles(const string& sFilePath,const unsigned int& iPointsCategory,const unsigned int& iTrianglesCategory)
	{
		printf("writing\n");
		FILE* fpFile = fopen(sFilePath.c_str(),"w");
		list<GenericNode*>::iterator liNodes;
		
		fprintf(fpFile,"# vtk DataFile Version 1.0\n");
		fprintf(fpFile,"boundary triangulation\n");
		fprintf(fpFile,"ASCII\n");
		fprintf(fpFile,"DATASET POLYDATA\n");
		
		unsigned int iNodesCount = (unsigned int)m_lpoTriangulationPoints.size();
		fprintf(fpFile,"POINTS %d float\n",iNodesCount);
		GenericNode* poNode = NULL;
		for(liNodes = m_lpoTriangulationPoints.begin() ; liNodes != m_lpoTriangulationPoints.end() ; liNodes++)
		{
			poNode = (*liNodes);
			fprintf(fpFile,"%E\t\t%E\t\t%E\n",poNode->GetX(),poNode->GetY(),poNode->GetZ());
		}
		
		unsigned int iTrianglesCount = (unsigned int)m_lpoTriangles.size();
		fprintf(fpFile,"POLYGONS %d %d\n",iTrianglesCount,4*iTrianglesCount);
		list<TriPatch*>::iterator liTriangles;
		TriPatch* poTriangles = NULL;
		unsigned int iID1 = 0;
		unsigned int iID2 = 0;
		unsigned int iID3 = 0;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			poTriangles = (*liTriangles);
			iID1 = poTriangles->GetPoint1()->GetID() - 1;
			iID2 = poTriangles->GetPoint2()->GetID() - 1;
			iID3 = poTriangles->GetPoint3()->GetID() - 1;
			fprintf(fpFile,"3\t\t%d\t\t%d\t\t%d\n",iID1,iID2,iID3);
		}

// 		fprintf(fpFile,"POINT_DATA %d\n",iNodesCount);
// 		fprintf(fpFile,"scalars NodeType integer\n");
// 		fprintf(fpFile,"LOOKUP_TABLE default\n");
// 		for(liNodes = plpoGraphNodes->begin() ; liNodes != plpoGraphNodes->end() ; liNodes++)
// 		{
// 			poNode = (*liNodes)->GetDataPointer();
// 			fprintf(fpFile,"%d\n",(*liNodes)->GetDataPointer()->GetCategory());
// 		}
// 		
// 		fprintf(fpFile,"CELL_DATA %d\n",iEdgesCount);
// 		fprintf(fpFile,"SCALARS SlipSystem integer\n");
// 		fprintf(fpFile,"LOOKUP_TABLE default\n");
// 		for(liEdges = plpoGraphEdges->begin() ; liEdges != plpoGraphEdges->end() ; liEdges++)
// 		{
// 			fprintf(fpFile,"%d\n",IndentifyFCCSlipSystem((*liEdges)->GetDataPointer()->GetSlipPlaneNormal(),(*liEdges)->GetDataPointer()->GetBurgersVector()));
// 		}
		
		fclose(fpFile);
	}
	list<Point*> Boundary::GenerateSurfacePoints(const unsigned int& iResolution)
	{
		list<Point*> lpoPoints;
		lpoPoints.clear();
		list<TriPatch*>::iterator liTriangles;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			(*liTriangles)->GeneratePoints(iResolution,&lpoPoints);
		}
		return lpoPoints;
	}
	void Boundary::RefineTriangulation()
	{
		list<TriPatch*>::iterator liTriangles;
		double dMinSideLength = DBL_MAX;
		double dTemp = 0.0;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			dTemp = (*liTriangles)->GetMinSideLength();
			if(dTemp < dMinSideLength)
			{
				dMinSideLength = dTemp;
			}
		}
		double dTolerance = 1E-6*dMinSideLength;
		list<GenericNode*> lpoNewPoints;
		lpoNewPoints.clear();
		list<TriPatch*> lpoNewTriangles;
		lpoNewTriangles.clear();
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			(*liTriangles)->Refine(lpoNewPoints,lpoNewTriangles,dTolerance);
			delete *liTriangles;
		}
		m_lpoTriangles.clear();
		list<GenericNode*>::iterator liPoints;
		// add the new points to the old ones
		for(liPoints = lpoNewPoints.begin() ; liPoints != lpoNewPoints.end() ; liPoints++)
		{
			m_lpoTriangulationPoints.push_back(*liPoints);
		}
		lpoNewPoints.clear();
		// add the new triangles to the now empty list
		for(liTriangles = lpoNewTriangles.begin() ; liTriangles != lpoNewTriangles.end() ; liTriangles++)
		{
			m_lpoTriangles.push_back(*liTriangles);
		}
			
		unsigned int iIndex = 0;
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			iIndex = iIndex + 1;
			(*liPoints)->SetID(iIndex);
		}
	}
}



