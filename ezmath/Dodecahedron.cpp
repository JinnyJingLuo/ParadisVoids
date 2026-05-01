// Paradis Processor Project
// Ahmed M. Hussein (mailto : am.hussin@gmail.com)
// June 2012

#include "Dodecahedron.h"
#include "math.h"

using namespace EZ;

namespace GeometrySystem
{
	Dodecahedron::Dodecahedron()
	{
		Initialize();
	}
	Dodecahedron::Dodecahedron(const Dodecahedron& oDodecahedron)
	{
		*this = oDodecahedron;
	}
	Dodecahedron::~Dodecahedron()
	{
		Reset();
	}
	Dodecahedron& Dodecahedron::operator=(const Dodecahedron& oDodecahedron)
	{
		Set(oDodecahedron.m_dDiameter,oDodecahedron.m_bIsHalfGrain);
		m_oSystem = oDodecahedron.m_oSystem;
		return *this;
	}
	void Dodecahedron::Reset()
	{
		list<TriPatch*>::iterator liTriangles;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			if((*liTriangles) != NULL)
			{
				delete (*liTriangles);
			}
		}
		
		list<GenericNode*>::iterator liPoints;
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			if((*liPoints) != NULL)
			{
				delete (*liPoints);
			}
		}
		Initialize();
	}
	void Dodecahedron::Set(const double& dDiameter,const bool& bHalfGrain)
	{
		Reset();
		m_dDiameter = dDiameter;
		m_bIsHalfGrain = bHalfGrain;
		GenerateTriangulations();
	}
	double Dodecahedron::GetDiameter() const
	{
		return m_dDiameter;
	}
	list<GenericNode*>* Dodecahedron::GetTriangulationPoints()
	{
		return &m_lpoTriangulationPoints;
	}
	list<TriPatch*>* Dodecahedron::GetTriangles()
	{
		return &m_lpoTriangles;
	}
	bool Dodecahedron::IsPointInside(const Point& oPoint,const double& dTolerance) const
	{
		Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
		list<TriPatch*>::const_iterator liTriangles;
		Vector oRelativePosition;
		double dTemp = 0.0;
		double dThreshold = -1.0E-6*m_dDiameter;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			oRelativePosition.SetByPoints(*(*liTriangles)->GetPoint1(),oLocalPoint);
			dTemp = oRelativePosition*(*liTriangles)->GetNormal();
			if(dTemp > dThreshold)
			{
				return false;
			}
		}
		return true;
	}
	double Dodecahedron::GetVolume() const
	{
		list<TriPatch*>::const_iterator liTriangles;
		double dTemp = 0.0;
		double dVolume = 0.0;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			dTemp = (*liTriangles)->GetPoint1()->GetX() + (*liTriangles)->GetPoint2()->GetX() + (*liTriangles)->GetPoint3()->GetX();
			dTemp = dTemp/3.0;
			dTemp = dTemp*(*liTriangles)->GetNormal().GetX();
			dTemp = dTemp*(*liTriangles)->GetArea();
			dVolume = dVolume + dTemp;
		}
		return dVolume;
	}
	void Dodecahedron::WriteParaviewFile(const string& sFileName) const
	{
		FILE* fpFile = fopen(sFileName.c_str(),"w");
		
		fprintf(fpFile,"# vtk DataFile Version 2.0\n");
		fprintf(fpFile,"# Dodecahedron\n");
		fprintf(fpFile,"ASCII\n");
		fprintf(fpFile,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(fpFile,"POINTS 20 float\n");
		
		list<GenericNode*>::const_iterator liPoints;
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			fprintf(fpFile,"%25.20f\t\t%25.20f\t\t%25.20f\n",(*liPoints)->GetX(),(*liPoints)->GetY(),(*liPoints)->GetZ());
		}
				
		list<TriPatch*>::const_iterator liTriangles;
		fprintf(fpFile,"CELLS %d %d\n",(unsigned int)m_lpoTriangles.size(),4*(unsigned int)m_lpoTriangles.size());
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
 			fprintf(fpFile,"3\t\t%d\t\t%d\t\t%d\n",(*liTriangles)->GetPoint1()->GetID() - 1,(*liTriangles)->GetPoint2()->GetID() - 1,(*liTriangles)->GetPoint3()->GetID() - 1);
		}
		fprintf(fpFile,"CELL_TYPES %d\n",(unsigned int)m_lpoTriangles.size());
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			fprintf(fpFile,"5\n");
		}
		// write some dummy data
		
		fprintf(fpFile,"POINT_DATA 20\n");
		fprintf(fpFile,"scalars VertexType integer\n");
		fprintf(fpFile,"LOOKUP_TABLE default\n");

		unsigned int iCount = 0;
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			if(iCount < 8)
			{
				fprintf(fpFile,"1\n");
			}
			else if(iCount < 12)
			{
				fprintf(fpFile,"2\n");
			}
			else if(iCount < 16)
			{
				fprintf(fpFile,"3\n");
			}
			else
			{
				fprintf(fpFile,"4\n");
			}
			iCount = iCount + 1;
		}
		
		fclose(fpFile);
	}
	void Dodecahedron::GenerateTriangulations()
	{
		double dScale = m_dDiameter/2.0/sqrt(3.0);
		double dOne = dScale;
		double dPhi = 1.61803398874989484820458683436563811772*dScale;
		double dOneOverPhi = dScale*dScale/dPhi;
		
		// generate the points
		vector<GenericNode*> vpoPoints;
		vpoPoints.reserve(20);
		vpoPoints.push_back(new GenericNode(-dOne,-dOne,-dOne));
		vpoPoints.push_back(new GenericNode(-dOne,dOne,-dOne));
		vpoPoints.push_back(new GenericNode(dOne,dOne,-dOne));
		vpoPoints.push_back(new GenericNode(dOne,-dOne,-dOne));
		vpoPoints.push_back(new GenericNode(-dOne,-dOne,dOne));
		vpoPoints.push_back(new GenericNode(-dOne,dOne,dOne));
		vpoPoints.push_back(new GenericNode(dOne,dOne,dOne));
		vpoPoints.push_back(new GenericNode(dOne,-dOne,dOne));

		vpoPoints.push_back(new GenericNode(0.0,-dOneOverPhi,-dPhi));
		vpoPoints.push_back(new GenericNode(0.0,-dOneOverPhi,dPhi));
		vpoPoints.push_back(new GenericNode(0.0,dOneOverPhi,-dPhi));
		vpoPoints.push_back(new GenericNode(0.0,dOneOverPhi,dPhi));
		
		vpoPoints.push_back(new GenericNode(-dPhi,0.0,-dOneOverPhi));
		vpoPoints.push_back(new GenericNode(-dPhi,0.0,dOneOverPhi));
		vpoPoints.push_back(new GenericNode(dPhi,0.0,-dOneOverPhi));
		vpoPoints.push_back(new GenericNode(dPhi,0.0,dOneOverPhi));
		
		vpoPoints.push_back(new GenericNode(-dOneOverPhi,-dPhi,0.0));
		vpoPoints.push_back(new GenericNode(-dOneOverPhi,dPhi,0.0));
		vpoPoints.push_back(new GenericNode(dOneOverPhi,-dPhi,0.0));
		vpoPoints.push_back(new GenericNode(dOneOverPhi,dPhi,0.0));
		unsigned int i = 0;
		unsigned int iSize = (unsigned int)vpoPoints.size();
		for(i = 0 ; i < iSize ; i++)
		{
			vpoPoints[i]->SetID(i + 1);
			m_lpoTriangulationPoints.push_back(vpoPoints[i]);
		}

		// generate the triangles
		
		// bottom
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[11],vpoPoints[6],vpoPoints[19],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[11],vpoPoints[19],vpoPoints[17],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[11],vpoPoints[17],vpoPoints[5],true));
		
		// bottom connection
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[11],vpoPoints[9],vpoPoints[6],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[6],vpoPoints[9],vpoPoints[15],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[6],vpoPoints[15],vpoPoints[19],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[19],vpoPoints[15],vpoPoints[2],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[19],vpoPoints[2],vpoPoints[17],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[17],vpoPoints[2],vpoPoints[1],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[17],vpoPoints[1],vpoPoints[5],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[5],vpoPoints[1],vpoPoints[13],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[5],vpoPoints[13],vpoPoints[11],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[11],vpoPoints[13],vpoPoints[9],true));
		
		// mid section
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[9],vpoPoints[4],vpoPoints[7],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[9],vpoPoints[7],vpoPoints[15],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[15],vpoPoints[7],vpoPoints[14],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[15],vpoPoints[14],vpoPoints[2],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[2],vpoPoints[14],vpoPoints[10],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[2],vpoPoints[10],vpoPoints[1],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[1],vpoPoints[10],vpoPoints[12],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[1],vpoPoints[12],vpoPoints[13],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[13],vpoPoints[12],vpoPoints[4],true));
		m_lpoTriangles.push_back(new TriPatch(vpoPoints[13],vpoPoints[4],vpoPoints[9],true));
		
		// continue based on whether the grain is full or not
		if(!m_bIsHalfGrain)
		{
			//top connection
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[7],vpoPoints[18],vpoPoints[14],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[14],vpoPoints[18],vpoPoints[3],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[14],vpoPoints[3],vpoPoints[10],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[10],vpoPoints[3],vpoPoints[8],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[10],vpoPoints[8],vpoPoints[12],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[12],vpoPoints[8],vpoPoints[0],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[12],vpoPoints[0],vpoPoints[4],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[4],vpoPoints[0],vpoPoints[16],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[4],vpoPoints[16],vpoPoints[7],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[7],vpoPoints[16],vpoPoints[18],true));
	
			//top
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[0],vpoPoints[18],vpoPoints[16],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[0],vpoPoints[3],vpoPoints[18],true));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[0],vpoPoints[8],vpoPoints[3],true));
		}
		else
		{
			//truncation face
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[4],vpoPoints[14],vpoPoints[7],false));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[4],vpoPoints[10],vpoPoints[14],false));
			m_lpoTriangles.push_back(new TriPatch(vpoPoints[4],vpoPoints[12],vpoPoints[10],false));
		}

		vpoPoints.clear();
	}
	void Dodecahedron::WriteTriangulation(const string& sFileName) const
	{
		FILE* fpFile = fopen(sFileName.c_str(),"w");
		fprintf(fpFile,"* points count\n");
		fprintf(fpFile,"20\n");
		list<GenericNode*>::const_iterator liPoints;
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			fprintf(fpFile,"%25.20f\t\t%25.20f\t\t%25.20f\n",(*liPoints)->GetX(),(*liPoints)->GetY(),(*liPoints)->GetZ());
		}
		
		fprintf(fpFile,"* traingles count\n");
		fprintf(fpFile,"36\n");
		
		// all the triangles should be constrained unless their normal is in the positive x direction
		list<TriPatch*>::const_iterator liTriangles;
		unsigned int iConstrained = 0;
		double dTolerance = 1.0E-6;
		for(liTriangles = m_lpoTriangles.begin() ; liTriangles != m_lpoTriangles.end() ; liTriangles++)
		{
			iConstrained = 1;
			if(fabs((*liTriangles)->GetNormal().GetX() - 1.0) < dTolerance)
			{
				iConstrained = 0;
			}
			fprintf(fpFile,"%d\t\t%d\t\t%d\t\t%d\n",(*liTriangles)->GetPoint1()->GetID(),(*liTriangles)->GetPoint2()->GetID(),(*liTriangles)->GetPoint3()->GetID(),iConstrained);
		}
		fprintf(fpFile,"\n\n");
		
		fclose(fpFile);
	}
	void Dodecahedron::SetHalfGrain(const bool& bHalfGrain)
	{
		m_bIsHalfGrain = bHalfGrain;
	}
	void Dodecahedron::SetOrigin(const Point& oOrigin)
	{
		double dX = oOrigin.GetX();
		double dY = oOrigin.GetY();
		double dZ = oOrigin.GetZ();
		list<GenericNode*>::iterator liPoints;
		for(liPoints = m_lpoTriangulationPoints.begin() ; liPoints != m_lpoTriangulationPoints.end() ; liPoints++)
		{
			(*liPoints)->Shift(dX,dY,dZ);
		}
	}
	void Dodecahedron::ReleaseTriangulations()
	{
		m_lpoTriangulationPoints.clear();
		m_lpoTriangles.clear();
	}
	Geometry* Dodecahedron::Clone()
	{
		return new Dodecahedron(*this);
	}
	void Dodecahedron::Initialize()
	{
		m_dDiameter = 0.0;
		m_lpoTriangulationPoints.clear();
		m_lpoTriangles.clear();
		m_bIsHalfGrain = false;
		m_oSystem.Reset();
	}
}


