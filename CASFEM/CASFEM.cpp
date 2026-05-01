#include "iostream"
#include "vector"
#include "CASFEM.h"
#include "cmath"
#include "Cylinder.h"
#include "Block.h"
#include "ctime"
#include "FEMMesh.h"
#include "Tools.h"
#include "Randomizer.h"
#include "FEMSolidSolver.h"
#include "FEMPotentialSolver.h"

using namespace std;
using namespace FEMSystem;
using namespace SupportSystem;


// test case generation function
void GenerateCylinderTestCase()
{
	Cylinder* poCylinder = new Cylinder;
	double dRadius = 3.0;
	double dLength = 10.0;
	unsigned int iRadialResolution = 6;
	unsigned int iCircumferentialResolution = 8;
	unsigned int iAxialResolution = 12;
	poCylinder->SetRadius(dRadius);
	poCylinder->SetLength(dLength);
	poCylinder->SetResolution(iRadialResolution,iCircumferentialResolution,iAxialResolution);
	vector<FEMNode*> vpoNodes;
	vector<FEMElement*> vpoElements;
	FEMMesh::GenerateMeshFromCylinder(poCylinder,vpoNodes,vpoElements,SolidFEMNode,SolidFEMElement);
	unsigned int iNodesCount = (unsigned int)vpoNodes.size();
	unsigned int iElementsCount = (unsigned int)vpoElements.size();
	unsigned int i = 0;
	unsigned int iTemp = 0;

	FILE* fpInput = fopen("cyl_TC.txt","w");
	char cTempString[1024];
	double dTolerance = 1.0E-6;
	// file organization
	// 1. header
	// 2. loads
	// 3. nodes (type, position, surface, initial conditions, boundary conditions)
	// 4. elements (connectivity, boundary conditions)
	
	// header section
	{
		fprintf(fpInput,"* output base file name\n");
		fprintf(fpInput,"fem_out\n");
		fprintf(fpInput,"* problem type, physics, number of steps and total time (ignored for (pseudo)static problems)\n");
		fprintf(fpInput,"2,2,1000000,50.0,10\n");
	}

	// loads section
	{
		fprintf(fpInput,"* Loads : \n");
		unsigned int iLoadsCount = 6;
		unsigned int iConstantLoadType = 1;
		fprintf(fpInput,"%d\n",iLoadsCount);
		// first load (constant with zero)
		double dConstantLoadValue = 0.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
	
		// second load (constant with value)
		dConstantLoadValue = 1000000.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
	
		// third load (x torsion)
		dConstantLoadValue = 100000.0;
		fprintf(fpInput,"3\n");
		fprintf(fpInput,"%e\n",dConstantLoadValue);
	
		// fourth load (y torsion)
		dConstantLoadValue = 100000.0;
		fprintf(fpInput,"4\n");
		fprintf(fpInput,"%e\n",dConstantLoadValue);
		
		// fifth load (nodal torsion)
		dConstantLoadValue = 5000.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
		
		// sixth load (time linear displacement)
		fprintf(fpInput,"%d\n",5);
		fprintf(fpInput,"%e\t\t%e\n",0.0,1.0);
	}
	
	// nodes section
	{
		fprintf(fpInput,"* Nodes : \n");
		fprintf(fpInput,"%d\n",iNodesCount);
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iTemp = 0;
			if(vpoNodes[i]->IsOnSurface())
			{
				iTemp = 1;
			}
			if(fabs(vpoNodes[i]->GetZ() + dLength/2.0) < dTolerance)
			{
				sprintf(cTempString,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1,1,1,1,1,1");
			}
			else if(fabs(vpoNodes[i]->GetZ() - dLength/2.0) < dTolerance)
			{
				sprintf(cTempString,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,1,0,1,1,6");
			}
			else
			{
				sprintf(cTempString,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,1,0,1,0,1");
			}
			
			fprintf(fpInput,"2\n");
			// x, y, z, ux, uy, uz, vx, vy, vz, ax, ay, az, cx, lx, cy, ly, cz, lz, surf
			fprintf(fpInput,"%e,%e,%e,%s,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),cTempString,iTemp);
		}
	}

	// elements section
	{
		unsigned int j = 0;
		vector<FEMNode*>* pvpoNodes = NULL;
		fprintf(fpInput,"* Elements : \n");
		fprintf(fpInput,"%d\n",iElementsCount);
		unsigned int iFacesCount = 6;
		Point oCenter;
		for(i = 0 ; i < iElementsCount ; i++)
		{
			pvpoNodes = vpoElements[i]->GetNodes();
			iTemp = (unsigned int)pvpoNodes->size();
			fprintf(fpInput,"2,2\n");
			fprintf(fpInput,"%d",pvpoNodes->at(0)->GetID());
			for(j = 1 ; j < iTemp ; j++)
			{
				fprintf(fpInput,",%d",pvpoNodes->at(j)->GetID());
			}
			fprintf(fpInput,"\n");
			
			sprintf(cTempString,"");
			for(j = 1 ; j <= iFacesCount ; j++)
			{
				oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
				if(fabs(oCenter.GetZ() - dLength/2.0) < dTolerance)
				{
					sprintf(cTempString,"%s1,1,1,",cTempString);
				}
				else
				{
					sprintf(cTempString,"%s1,1,1,",cTempString);
				}
	// 			if(fabs(oCenter.GetX()*oCenter.GetX() + oCenter.GetY()*oCenter.GetY() - dRadius*dRadius) < dTolerance)
	// 			{
	// 				sprintf(cTempString,"%s1,1,1,",cTempString);	// surface faces
	// 			}
	// 			else
	// 			{
	// 				sprintf(cTempString,"%s1,1,1,",cTempString);	// internal faces
	// 			}
			}
			sprintf(cTempString,"%s1,1,1",cTempString);			// body forces
			fprintf(fpInput,"%s\n",cTempString);
			vpoElements[i]->DeleteGeometry();
			delete vpoElements[i];
		}
	}

	for(i = 0 ; i < iNodesCount ; i++)
	{
		delete vpoNodes[i];
	}
	fclose(fpInput);
}

void GenerateBlockTestCase()
{
	Block* poBlock = new Block;
	double dXLength = 1.0;
	double dYLength = 1.0;
	double dZLength = 3.0;
	unsigned int iXResolution = 5;
	unsigned int iYResolution = 5;
	unsigned int iZResolution = 15;
	poBlock->SetXLength(dXLength);
	poBlock->SetYLength(dYLength);
	poBlock->SetZLength(dZLength);
	poBlock->SetResolution(iXResolution,iYResolution,iZResolution);
	vector<FEMNode*> vpoNodes;
	vector<FEMElement*> vpoElements;
	FEMMesh::GenerateMeshFromBlock(poBlock,vpoNodes,vpoElements,SolidFEMNode,SolidFEMElement);
	unsigned int iNodesCount = (unsigned int)vpoNodes.size();
	unsigned int iElementsCount = (unsigned int)vpoElements.size();
	unsigned int i = 0;
	unsigned int iTemp = 0;
	FILE* fpInput = fopen("block_const.txt","w");
	char cTempString[1024];
	double dTolerance = 1.0E-6;
	// file organization
	// 1. header
	// 2. loads
	// 3. nodes (type, position, surface, initial conditions, boundary conditions)
	// 4. elements (connectivity, boundary conditions)
	
	// header section
	{
		fprintf(fpInput,"* output base file name\n");
		fprintf(fpInput,"fem_out\n");
		fprintf(fpInput,"* problem type, physics, current time, time step, target time, this output count\n");
		fprintf(fpInput,"1,2,0.0,1.0E-3,1.0,0\n");
	}
	
	// loads section
	{
		fprintf(fpInput,"* Loads : \n");
		unsigned int iLoadsCount = 3;
		unsigned int iConstantLoadType = 1;
		fprintf(fpInput,"%d\n",iLoadsCount);

		// first load (constant with zero)
		double dConstantLoadValue = 0.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
	
		// second load (constant with value)
		dConstantLoadValue = 0.5;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
		
		// third load (time sinusoidal displacement)
		fprintf(fpInput,"%d\n",2);
		fprintf(fpInput,"%e\t\t%e\t\t%e\n",0.5,628.3185307,0.0);
	}
	
	// materials section
	{
		fprintf(fpInput,"* Materials : \n");
		unsigned int iMaterailsCount = 1;
		unsigned int iMaterialType = 1;					// J2 plastic material
		fprintf(fpInput,"%d\n",iMaterailsCount);
		// first material
		fprintf(fpInput,"%d\n",iMaterialType);
		// properties for nickel E nu rho k alpha c To sig_y K H
		fprintf(fpInput,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",2.0E11,0.31,8912.0,90.9,1.3E-6,440.0,298.0,1.5E8,0.0,0.0);
	}
	
	// nodes section
	{		
		fprintf(fpInput,"* Nodes : \n");
		fprintf(fpInput,"%d\n",iNodesCount);
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iTemp = 0;
			if(vpoNodes[i]->IsOnSurface())
			{
				iTemp = 1;
			}
			fprintf(fpInput,"2\n");
			if(fabs(vpoNodes[i]->GetZ() + dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,1,1,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else if(fabs(vpoNodes[i]->GetZ() - dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,2,1,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else
			{
				fprintf(fpInput,"%e,%e,%e,0,1,0,1,0,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			// initial displacements, velocities and accelerations
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
			// forces and stresses
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
		}
	}
	
	// elements section
	{
		unsigned int j = 0;
		vector<FEMNode*>* pvpoNodes = NULL;
		fprintf(fpInput,"* Elements : \n");
		fprintf(fpInput,"%d\n",iElementsCount);
		unsigned int iFacesCount = 6;
		Point oCenter;
		for(i = 0 ; i < iElementsCount ; i++)
		{
			pvpoNodes = vpoElements[i]->GetNodes();
			iTemp = (unsigned int)pvpoNodes->size();
			// element geometry type, element physics type, element material type
			fprintf(fpInput,"2,2,1\n");
			fprintf(fpInput,"%d",pvpoNodes->at(0)->GetID());
			for(j = 1 ; j < iTemp ; j++)
			{
				fprintf(fpInput,",%d",pvpoNodes->at(j)->GetID());
			}
			fprintf(fpInput,"\n");
			
			sprintf(cTempString,"");
			for(j = 1 ; j <= iFacesCount ; j++)
			{
				oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
				if(fabs(oCenter.GetZ() - dZLength/2.0) < dTolerance)
				{
					sprintf(cTempString,"%s1,1,1,",cTempString);
				}
				else
				{
					sprintf(cTempString,"%s1,1,1,",cTempString);
				}
			}
			sprintf(cTempString,"%s1,1,1",cTempString);			// body forces
			fprintf(fpInput,"%s\n",cTempString);
			vpoElements[i]->DeleteGeometry();
			delete vpoElements[i];
		}
	}

	for(i = 0 ; i < iNodesCount ; i++)
	{
		delete vpoNodes[i];
	}
	fclose(fpInput);
}

void GenerateCylinderThermalTestCase()
{
	Cylinder* poCylinder = new Cylinder;
	double dRadius = 3.0;
	double dLength = 10.0;
	unsigned int iRadialResolution = 8;
	unsigned int iCircumferentialResolution = 16;
	unsigned int iAxialResolution = 10;
	poCylinder->SetRadius(dRadius);
	poCylinder->SetLength(dLength);
	poCylinder->SetResolution(iRadialResolution,iCircumferentialResolution,iAxialResolution);
	vector<FEMNode*> vpoNodes;
	vector<FEMElement*> vpoElements;
	FEMMesh::GenerateMeshFromCylinder(poCylinder,vpoNodes,vpoElements,PotentialFEMNode,PotentialFEMElement);
	unsigned int iNodesCount = (unsigned int)vpoNodes.size();
	unsigned int iElementsCount = (unsigned int)vpoElements.size();
	unsigned int i = 0;
	unsigned int iTemp = 0;

	FILE* fpInput = fopen("cyl_LIN_Therm.txt","w");
	
	// file organization
	// 1. header
	// 2. loads
	// 3. nodes (type, position, surface, initial conditions, boundary conditions)
	// 4. elements (connectivity, boundary conditions)
	
	// header section
	fprintf(fpInput,"* output base file name\n");
	fprintf(fpInput,"fem_out\n");
	fprintf(fpInput,"* problem type, physics, number of steps and total time (ignored for (pseudo)static problems)\n");
	fprintf(fpInput,"1,1,1,0.0\n");
	fprintf(fpInput,"* Gauss points counts : face - body\n");
	fprintf(fpInput,"3,3\n");

	// loads section
	fprintf(fpInput,"* Loads : \n");
	unsigned int iLoadsCount = 3;
	unsigned int iConstantLoadType = 1;
	fprintf(fpInput,"%d\n",iLoadsCount);
	// first load (constant with zero)
	double dConstantLoadValue = 0.0;
	fprintf(fpInput,"%d\n",iConstantLoadType);
	fprintf(fpInput,"%e\n",dConstantLoadValue);

	// second load (constant with value)
	dConstantLoadValue = 0.0;
	fprintf(fpInput,"%d\n",iConstantLoadType);
	fprintf(fpInput,"%e\n",dConstantLoadValue);
	
	// third load (constant with value)
	dConstantLoadValue = 100.0;
	fprintf(fpInput,"%d\n",iConstantLoadType);
	fprintf(fpInput,"%e\n",dConstantLoadValue);
	
	// nodes section
	fprintf(fpInput,"* Nodes : \n");
	fprintf(fpInput,"%d\n",iNodesCount);
	double dTolerance = 1.0E-6;
	char cTempString[1024];
	for(i = 0 ; i < iNodesCount ; i++)
	{
		iTemp = 0;
		if(vpoNodes[i]->IsOnSurface())
		{
			iTemp = 1;
		}
		if(fabs(vpoNodes[i]->GetZ() + dLength/2.0) < dTolerance)
		{
			sprintf(cTempString,"1,2");
		}
		else if(fabs(vpoNodes[i]->GetZ() - dLength/2.0) < dTolerance)
		{
			sprintf(cTempString,"1,3");
		}
		else
		{
			sprintf(cTempString,"0,1");
		}
		
		fprintf(fpInput,"1\n");
		// x, y, z, cx, lx, cy, ly, cz, lz, surf
		fprintf(fpInput,"%e,%e,%e,%s,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),cTempString,iTemp);
	}

	// elements section
	unsigned int j = 0;
	vector<FEMNode*>* pvpoNodes = NULL;
	fprintf(fpInput,"* Elements : \n");
	fprintf(fpInput,"%d\n",iElementsCount);
	unsigned int iFacesCount = 6;
	Point oCenter;
	for(i = 0 ; i < iElementsCount ; i++)
	{
		pvpoNodes = vpoElements[i]->GetNodes();
		iTemp = (unsigned int)pvpoNodes->size();
		fprintf(fpInput,"2,1\n");
		fprintf(fpInput,"%d",pvpoNodes->at(0)->GetID());
		for(j = 1 ; j < iTemp ; j++)
		{
			fprintf(fpInput,",%d",pvpoNodes->at(j)->GetID());
		}
		fprintf(fpInput,"\n");
		
		sprintf(cTempString,"");
		for(j = 1 ; j <= iFacesCount ; j++)
		{
			oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
			if(fabs(oCenter.GetZ() - dLength/2.0) < dTolerance)
			{
				sprintf(cTempString,"%s1,1,1,",cTempString);
			}
			else
			{
				sprintf(cTempString,"%s1,1,1,",cTempString);
			}
// 			if(fabs(oCenter.GetX()*oCenter.GetX() + oCenter.GetY()*oCenter.GetY() - dRadius*dRadius) < dTolerance)
// 			{
// 				sprintf(cTempString,"%s1,1,1,",cTempString);	// surface faces
// 			}
// 			else
// 			{
// 				sprintf(cTempString,"%s1,1,1,",cTempString);	// internal faces
// 			}
		}
		sprintf(cTempString,"%s1,1,1",cTempString);			// body forces
		fprintf(fpInput,"%s\n",cTempString);
		vpoElements[i]->DeleteGeometry();
		delete vpoElements[i];
	}

	for(i = 0 ; i < iNodesCount ; i++)
	{
		delete vpoNodes[i];
	}
	fclose(fpInput);
}

void GenerateTMBlockTestCase()
{
	Block* poBlock = new Block;
	double dXLength = 0.01;
	double dYLength = 0.01;
	double dZLength = 0.03;
	unsigned int iXResolution = 5;
	unsigned int iYResolution = 5;
	unsigned int iZResolution = 15;
	poBlock->SetXLength(dXLength);
	poBlock->SetYLength(dYLength);
	poBlock->SetZLength(dZLength);
	poBlock->SetResolution(iXResolution,iYResolution,iZResolution);
	vector<FEMNode*> vpoNodes;
	vector<FEMElement*> vpoElements;
	FEMMesh::GenerateMeshFromBlock(poBlock,vpoNodes,vpoElements,ThermoMechanicalFEMNode,ThermoMechanicalFEMElement);
	unsigned int iNodesCount = (unsigned int)vpoNodes.size();
	unsigned int iElementsCount = (unsigned int)vpoElements.size();
	unsigned int i = 0;
	unsigned int iTemp = 0;
	FILE* fpInput = fopen("block_const.txt","w");
	char cTempString[1024];
	double dTolerance = 1.0E-6;
	// file organization
	// 1. header
	// 2. loads
	// 3. nodes (type, position, surface, initial conditions, boundary conditions)
	// 4. elements (connectivity, boundary conditions)
	
	// header section
	{
		fprintf(fpInput,"* output base file name\n");
		fprintf(fpInput,"fem_out\n");
		fprintf(fpInput,"* problem type, physics, current time, time step, target time, this output count\n");
		fprintf(fpInput,"3,3,0.0,1.0E-3,1.0,0\n");
	}
	
	// loads section
	{
		fprintf(fpInput,"* Loads : \n");
		unsigned int iLoadsCount = 4;
		unsigned int iConstantLoadType = 1;
		fprintf(fpInput,"%d\n",iLoadsCount);
		// first load (constant with zero)
		double dConstantLoadValue = 0.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
		
		// second load (time sinusoidal displacement)
		fprintf(fpInput,"%d\n",2);
		fprintf(fpInput,"%e\t\t%e\t\t%e\n",0.005,628.3185307,0.0);
		
		// third load (themral linear load)
		fprintf(fpInput,"%d\n",5);
		fprintf(fpInput,"%e\t\t%e\n",0.0,25.0);
		
		// fourth load (constant themral load)
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",500.0);
	}
	
	// materials section
	{
		fprintf(fpInput,"* Materials : \n");
		unsigned int iMaterailsCount = 1;
		unsigned int iLinearIsotropicMaterialType = 1;
		fprintf(fpInput,"%d\n",iMaterailsCount);
		// first material
		fprintf(fpInput,"%d\n",iLinearIsotropicMaterialType);
		// properties for nickel E nu rho k alpha c To
		fprintf(fpInput,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",2.0E11,0.31,8912.0,90.9,1.3E-6,440.0,298.0);
	}
	
	// nodes section
	{		
		fprintf(fpInput,"* Nodes : \n");
		fprintf(fpInput,"%d\n",iNodesCount);
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iTemp = 0;
			if(vpoNodes[i]->IsOnSurface())
			{
				iTemp = 1;
			}
			fprintf(fpInput,"3\n");
			if(fabs(vpoNodes[i]->GetZ() + dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,1,1,1,1,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else if(fabs(vpoNodes[i]->GetZ() - dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,2,1,1,1,1,1,4,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else
			{
				fprintf(fpInput,"%e,%e,%e,0,1,0,1,0,1,0,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			// initial displacements, velocities, accelerations, temperatures, heating rates, ...
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
			// forces, fluxes and stresses
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
		}
	}
	
	// elements section
	{
		unsigned int j = 0;
		vector<FEMNode*>* pvpoNodes = NULL;
		fprintf(fpInput,"* Elements : \n");
		fprintf(fpInput,"%d\n",iElementsCount);
		unsigned int iFacesCount = 6;
		Point oCenter;
		for(i = 0 ; i < iElementsCount ; i++)
		{
			pvpoNodes = vpoElements[i]->GetNodes();
			iTemp = (unsigned int)pvpoNodes->size();
			// element geometry type, element physics type, element material type
			fprintf(fpInput,"2,3,1\n");
			fprintf(fpInput,"%d",pvpoNodes->at(0)->GetID());
			for(j = 1 ; j < iTemp ; j++)
			{
				fprintf(fpInput,",%d",pvpoNodes->at(j)->GetID());
			}
			fprintf(fpInput,"\n");
			
			sprintf(cTempString,"");
			for(j = 1 ; j <= iFacesCount ; j++)
			{
				oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
				if(fabs(oCenter.GetZ() - dZLength/2.0) < dTolerance)
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
				else
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
			}
			sprintf(cTempString,"%s1,1,1,1",cTempString);			// body forces
			fprintf(fpInput,"%s\n",cTempString);
			vpoElements[i]->DeleteGeometry();
			delete vpoElements[i];
		}
	}
	
	for(i = 0 ; i < iNodesCount ; i++)
	{
		delete vpoNodes[i];
	}
	fclose(fpInput);
}

void GenerateTMBlockTestCaseInnerSource()
{
	Block* poBlock = new Block;
	double dXLength = 1.0;
	double dYLength = 1.0;
	double dZLength = 3.0;
	unsigned int iXResolution = 5;
	unsigned int iYResolution = 5;
	unsigned int iZResolution = 15;
	poBlock->SetXLength(dXLength);
	poBlock->SetYLength(dYLength);
	poBlock->SetZLength(dZLength);
	poBlock->SetResolution(iXResolution,iYResolution,iZResolution);
	vector<FEMNode*> vpoNodes;
	vector<FEMElement*> vpoElements;
	FEMMesh::GenerateMeshFromBlock(poBlock,vpoNodes,vpoElements,ThermoMechanicalFEMNode,ThermoMechanicalFEMElement);
	unsigned int iNodesCount = (unsigned int)vpoNodes.size();
	unsigned int iElementsCount = (unsigned int)vpoElements.size();
	unsigned int i = 0;
	unsigned int iTemp = 0;
	FILE* fpInput = fopen("block_const.txt","w");
	char cTempString[1024];
	double dTolerance = 1.0E-6;
	// file organization
	// 1. header
	// 2. loads
	// 3. nodes (type, position, surface, initial conditions, boundary conditions)
	// 4. elements (connectivity, boundary conditions)
	
	// header section
	{
		fprintf(fpInput,"* output base file name\n");
		fprintf(fpInput,"fem_out\n");
		fprintf(fpInput,"* problem type, physics, number of steps and total time (ignored for (pseudo)static problems)\n");
		fprintf(fpInput,"3,3,1000000,50.0,1\n");
	}
	
	// loads section
	{
		fprintf(fpInput,"* Loads : \n");
		unsigned int iLoadsCount = 2;
		unsigned int iConstantLoadType = 1;
		fprintf(fpInput,"%d\n",iLoadsCount);
		// first load (constant with zero)
		double dConstantLoadValue = 0.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
		
		// third load (time linear temperature)
		fprintf(fpInput,"%d\n",5);
		fprintf(fpInput,"%e\t\t%e\n",0.0,2.0);
	}
	
	// materials section
	{
		fprintf(fpInput,"* Materials : \n");
		unsigned int iMaterailsCount = 1;
		unsigned int iLinearIsotropicMaterialType = 1;
		fprintf(fpInput,"%d\n",iMaterailsCount);
		// first material
		double dConstantLoadValue = 0.0;
		fprintf(fpInput,"%d\n",iLinearIsotropicMaterialType);
		// properties fo nickel E nu rho k c To
		fprintf(fpInput,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",2.0E11,0.31,8912.0,90.9,440.0,298.0);
	}
	
	// nodes section
	{		
		fprintf(fpInput,"* Nodes : \n");
		fprintf(fpInput,"%d\n",iNodesCount);
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iTemp = 0;
			if(vpoNodes[i]->IsOnSurface())
			{
				iTemp = 1;
			}
			if(fabs(vpoNodes[i]->GetZ() + dZLength/2.0) < dTolerance)
			{
				sprintf(cTempString,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1,1,1,1,1,1,1,1");
			}
			else if(fabs(vpoNodes[i]->GetZ() - dZLength/2.0) < dTolerance)
			{
				sprintf(cTempString,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1,1,1,1,1,1,0,1");
			}
			else
			{
				sprintf(cTempString,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,1,0,1,0,1,0,1");
			}
			fprintf(fpInput,"3\n");
			// x, y, z, ux, uy, uz, vx, vy, vz, ax, ay, az, cx, lx, cy, ly, cz, lz, surf
			fprintf(fpInput,"%e,%e,%e,%s,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),cTempString,iTemp);
		}
	}
	
	// elements section
	{
		unsigned int j = 0;
		vector<FEMNode*>* pvpoNodes = NULL;
		fprintf(fpInput,"* Elements : \n");
		fprintf(fpInput,"%d\n",iElementsCount);
		unsigned int iFacesCount = 6;
		Point oCenter;
		for(i = 0 ; i < iElementsCount ; i++)
		{
			pvpoNodes = vpoElements[i]->GetNodes();
			iTemp = (unsigned int)pvpoNodes->size();
			// element geometry type, element physics type, element material type
			fprintf(fpInput,"2,3,1\n");
			fprintf(fpInput,"%d",pvpoNodes->at(0)->GetID());
			for(j = 1 ; j < iTemp ; j++)
			{
				fprintf(fpInput,",%d",pvpoNodes->at(j)->GetID());
			}
			fprintf(fpInput,"\n");
			
			sprintf(cTempString,"");
			for(j = 1 ; j <= iFacesCount ; j++)
			{
				oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
				if(fabs(oCenter.GetZ() - dZLength/2.0) < dTolerance)
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
				else
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
			}
			if(vpoElements[i]->GetGeometry()->GetCenterPoint().Distance(Point(0.0,0.0,0.0)) < 1.0)
			{
				sprintf(cTempString,"%s1,1,1,2",cTempString);			// body forces
			}
			else
			{
				sprintf(cTempString,"%s1,1,1,1",cTempString);			// body forces
			}
			fprintf(fpInput,"%s\n",cTempString);
			vpoElements[i]->DeleteGeometry();
			delete vpoElements[i];
		}
	}

	for(i = 0 ; i < iNodesCount ; i++)
	{
		delete vpoNodes[i];
	}
	fclose(fpInput);
}

void GenerateBlockAndSphere()
{

	Block* poBlock = new Block;
	double dXLength = 1.25;
	double dYLength = 1.25;
	double dZLength = 1.25;
	unsigned int iXResolution = 10;
	unsigned int iYResolution = 10;
	unsigned int iZResolution = 10;
	poBlock->SetXLength(dXLength);
	poBlock->SetYLength(dYLength);
	poBlock->SetZLength(dZLength);
	poBlock->SetResolution(iXResolution,iYResolution,iZResolution);
	vector<FEMNode*> vpoNodes;
	vector<FEMElement*> vpoElements;
	FEMMesh::GenerateMeshFromBlock(poBlock,vpoNodes,vpoElements,ThermoMechanicalFEMNode,ThermoMechanicalFEMElement);
	unsigned int iNodesCount = (unsigned int)vpoNodes.size();
	unsigned int iElementsCount = (unsigned int)vpoElements.size();
	unsigned int i = 0;
	unsigned int iTemp = 0;


	FILE* fpInput = fopen("spherical_precipitate.txt","w");
	char cTempString[1024];
	double dTolerance = 1.0E-6;
	// file organization
	// 1. header
	// 2. loads
	// 3. nodes (type, position, surface, initial conditions, boundary conditions)
		// 4. elements (connectivity, boundary conditions)
	
	// header section
	{
		fprintf(fpInput,"* output base file name\n");
		fprintf(fpInput,"fem_out\n");
		fprintf(fpInput,"* problem type, physics, current time, time step, target time, this output count\n");
		fprintf(fpInput,"1,3,0.0,1.0E-3,1.0,0\n");
	}
	
	// loads section
	{
		fprintf(fpInput,"* Loads : \n");
		unsigned int iLoadsCount = 3;
		unsigned int iConstantLoadType = 1;
		fprintf(fpInput,"%d\n",iLoadsCount);

		// first load (constant with zero)
		double dConstantLoadValue = 0.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
	
		// second load (constant with value)
		dConstantLoadValue = 500.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
		
		// third load (time sinusoidal displacement)
		fprintf(fpInput,"%d\n",2);
		fprintf(fpInput,"%e\t\t%e\t\t%e\n",0.5,628.3185307,0.0);
	}
	
	// materials section
	{
		fprintf(fpInput,"* Materials : \n");
		unsigned int iMaterailsCount = 1;
		unsigned int iMaterialType = 1;					// LEI material
		fprintf(fpInput,"%d\n",iMaterailsCount);
		// first material
		fprintf(fpInput,"%d\n",iMaterialType);
		// properties for nickel E nu rho k alpha c To sig_y K H
		fprintf(fpInput,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",2.0E11,0.31,8912.0,0.0,1.3E-6,440.0,298.0);
	}
	
	// nodes section
	{		
		fprintf(fpInput,"* Nodes : \n");
		fprintf(fpInput,"%d\n",iNodesCount);

		for(i = 0 ; i < iNodesCount ; i++)
		{
			iTemp = 0;
			if(vpoNodes[i]->IsOnSurface())
			{
				iTemp = 1;
			}
			fprintf(fpInput,"3\n");
			if(fabs(vpoNodes[i]->GetZ() + dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,1,1,1,1,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else if(fabs(vpoNodes[i]->GetZ() - dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,1,1,1,1,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else
			{
				// nodes inside the central sphere of radius 1.5 are hotter
				double dX = vpoNodes[i]->GetX();
				double dY = vpoNodes[i]->GetY();
				double dZ = vpoNodes[i]->GetZ();
				double dR = sqrt(dX*dX + dY*dY + dZ*dZ);
				if(dR > -0.5)
				{
					fprintf(fpInput,"%e,%e,%e,0,1,0,1,0,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
				}
				else
				{
					fprintf(fpInput,"%e,%e,%e,0,1,0,1,0,1,1,2,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
				}
			}
			// initial displacements, velocities and accelerations
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
			// forces and stresses
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
		}
	}
	
	// elements section
	{
		unsigned int j = 0;
		vector<FEMNode*>* pvpoNodes = NULL;
		fprintf(fpInput,"* Elements : \n");
		fprintf(fpInput,"%d\n",iElementsCount);
		unsigned int iFacesCount = 6;
		Point oCenter;
		for(i = 0 ; i < iElementsCount ; i++)
		{
			pvpoNodes = vpoElements[i]->GetNodes();
			iTemp = (unsigned int)pvpoNodes->size();
			// element geometry type, element physics type, element material type
			fprintf(fpInput,"2,3,1\n");
			fprintf(fpInput,"%d",pvpoNodes->at(0)->GetID());
			for(j = 1 ; j < iTemp ; j++)
			{
				fprintf(fpInput,",%d",pvpoNodes->at(j)->GetID());
			}
			fprintf(fpInput,"\n");
			
			sprintf(cTempString,"");
			for(j = 1 ; j <= iFacesCount ; j++)
			{
				oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
				if(fabs(oCenter.GetZ() - dZLength/2.0) < dTolerance)
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
				else
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
			}
			sprintf(cTempString,"%s1,1,1,1",cTempString);			// body forces
			fprintf(fpInput,"%s\n",cTempString);
			vpoElements[i]->DeleteGeometry();
			delete vpoElements[i];
		}
	}

	for(i = 0 ; i < iNodesCount ; i++)
	{
		delete vpoNodes[i];
	}
	fclose(fpInput);
}


void GenerateInnerHeatBlockTestCase()
{
	Block* poBlock = new Block;
	double dXLength = 1.0;
	double dYLength = 1.0;
	double dZLength = 3.0;
	unsigned int iXResolution = 10;
	unsigned int iYResolution = 10;
	unsigned int iZResolution = 30;
	poBlock->SetXLength(dXLength);
	poBlock->SetYLength(dYLength);
	poBlock->SetZLength(dZLength);
	poBlock->SetResolution(iXResolution,iYResolution,iZResolution);
	vector<FEMNode*> vpoNodes;
	vector<FEMElement*> vpoElements;
	FEMMesh::GenerateMeshFromBlock(poBlock,vpoNodes,vpoElements,ThermoMechanicalFEMNode,ThermoMechanicalFEMElement);
	unsigned int iNodesCount = (unsigned int)vpoNodes.size();
	unsigned int iElementsCount = (unsigned int)vpoElements.size();
	unsigned int i = 0;
	unsigned int iTemp = 0;
	FILE* fpInput = fopen("inner_heat_block.txt","w");
	char cTempString[1024];
	double dTolerance = 1.0E-6;
	// file organization
	// 1. header
	// 2. loads
	// 3. nodes (type, position, surface, initial conditions, boundary conditions)
	// 4. elements (connectivity, boundary conditions)
	
	// header section
	{
		fprintf(fpInput,"* output base file name\n");
		fprintf(fpInput,"fem_out\n");
		fprintf(fpInput,"* problem type, physics, current time, time step, target time, this output count\n");
		fprintf(fpInput,"1,3,0.0,1.0E-3,1.0,0\n");
	}
	
	// loads section
	{
		fprintf(fpInput,"* Loads : \n");
		unsigned int iLoadsCount = 3;
		unsigned int iConstantLoadType = 1;
		fprintf(fpInput,"%d\n",iLoadsCount);

		// first load (constant with zero)
		double dConstantLoadValue = 0.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
	
		// second load (constant with value)
		dConstantLoadValue = 500.0;
		fprintf(fpInput,"%d\n",iConstantLoadType);
		fprintf(fpInput,"%e\n",dConstantLoadValue);
		
		// third load (time sinusoidal displacement)
		fprintf(fpInput,"%d\n",2);
		fprintf(fpInput,"%e\t\t%e\t\t%e\n",0.5,628.3185307,0.0);
	}
	
	// materials section
	{
		fprintf(fpInput,"* Materials : \n");
		unsigned int iMaterailsCount = 1;
		unsigned int iMaterialType = 1;					// LEI material
		fprintf(fpInput,"%d\n",iMaterailsCount);
		// first material
		fprintf(fpInput,"%d\n",iMaterialType);
		// properties for nickel E nu rho k alpha c To sig_y K H
		fprintf(fpInput,"%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\t\t%e\n",2.0E11,0.31,8912.0,90.9,1.3E-6,440.0,298.0);
	}
	
	// nodes section
	{		
		fprintf(fpInput,"* Nodes : \n");
		fprintf(fpInput,"%d\n",iNodesCount);
		for(i = 0 ; i < iNodesCount ; i++)
		{
			iTemp = 0;
			if(vpoNodes[i]->IsOnSurface())
			{
				iTemp = 1;
			}
			fprintf(fpInput,"3\n");
			if(fabs(vpoNodes[i]->GetZ() + dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,1,1,1,1,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else if(fabs(vpoNodes[i]->GetZ() - dZLength/2.0) < dTolerance)
			{
				fprintf(fpInput,"%e,%e,%e,1,1,1,1,1,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
			}
			else
			{
				if(vpoNodes[i]->GetZ() > 0.5)
				{
					fprintf(fpInput,"%e,%e,%e,0,1,0,1,0,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
				}
				else if(vpoNodes[i]->GetZ() < -0.5)
				{
					fprintf(fpInput,"%e,%e,%e,0,1,0,1,0,1,1,1,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
				}
				else
				{
					fprintf(fpInput,"%e,%e,%e,0,1,0,1,0,1,1,2,%d\n",vpoNodes[i]->GetX(),vpoNodes[i]->GetY(),vpoNodes[i]->GetZ(),iTemp);
				}
			}
			// initial displacements, velocities and accelerations
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
			// forces and stresses
			fprintf(fpInput,"0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0\n");
		}
	}
	
	// elements section
	{
		unsigned int j = 0;
		vector<FEMNode*>* pvpoNodes = NULL;
		fprintf(fpInput,"* Elements : \n");
		fprintf(fpInput,"%d\n",iElementsCount);
		unsigned int iFacesCount = 6;
		Point oCenter;
		for(i = 0 ; i < iElementsCount ; i++)
		{
			pvpoNodes = vpoElements[i]->GetNodes();
			iTemp = (unsigned int)pvpoNodes->size();
			// element geometry type, element physics type, element material type
			fprintf(fpInput,"2,3,1\n");
			fprintf(fpInput,"%d",pvpoNodes->at(0)->GetID());
			for(j = 1 ; j < iTemp ; j++)
			{
				fprintf(fpInput,",%d",pvpoNodes->at(j)->GetID());
			}
			fprintf(fpInput,"\n");
			
			sprintf(cTempString,"");
			for(j = 1 ; j <= iFacesCount ; j++)
			{
				oCenter = vpoElements[i]->GetGeometry()->GetFaceCenter(j);
				if(fabs(oCenter.GetZ() - dZLength/2.0) < dTolerance)
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
				else
				{
					sprintf(cTempString,"%s1,1,1,1,",cTempString);
				}
			}
			sprintf(cTempString,"%s1,1,1,1",cTempString);			// body forces
			fprintf(fpInput,"%s\n",cTempString);
			vpoElements[i]->DeleteGeometry();
			delete vpoElements[i];
		}
	}

	for(i = 0 ; i < iNodesCount ; i++)
	{
		delete vpoNodes[i];
	}
	fclose(fpInput);
}



CASFEM* CASFEM::CASFEMInstance = NULL;
CASFEM* CASFEM::CreateInstance()
{
	if(CASFEMInstance == NULL)
	{
		CASFEMInstance = new CASFEM();
	}
	return CASFEMInstance;
}
CASFEM::CASFEM()
{
 	m_poData = NULL;
 	m_poFEMSolver = NULL;
 	char cFileName[500];
 	sprintf(cFileName,"CASFEM_log.txt");
 	m_fpLog = fopen(cFileName,"w");
 	WriteInLog("CASFEM Object Created");
}
CASFEM::~CASFEM()
{
	Reset();
}
void CASFEM::Reset()
{
	if(m_poData != NULL)
	{
		m_poData->Reset();
		delete m_poData;
	}
	if(m_poFEMSolver != NULL)
	{
		m_poFEMSolver->Reset();
		delete m_poFEMSolver;
	}

	WriteInLog("Destroying CASFEM Object");
	fclose(m_fpLog);
	m_fpLog = NULL;
}
bool CASFEM::Initialize(const string& sFileName)
{
	m_poData = MainDataStructure::CreateInstance();
	m_poData->SetInputFileName(sFileName);
	// read FEM problem
	PrintOnScreen("Generating FEM Problem");
 	m_poData->ReadInput();
 	PrintOnScreen("Done generating FEM problem");
 	WriteInLog("FEM Problem Generated");
	m_poFEMSolver = FEMSolver::CreateSolverByPhysicsIndex(m_poData->GetProblemPhysics(),(const ProblemType&)(m_poData->GetProblemType()));
	if(m_poFEMSolver == NULL)
	{
		printf("error: unknown problem physics %d\n",m_poData->GetProblemPhysics());
		return false;
	}
	m_poFEMSolver->SetDataStructure(m_poData);
	return true;
}
MainDataStructure* CASFEM::GetMainDataStructure() const
{
	return m_poData;
}
void CASFEM::WriteInLog(const string& sLogString)
{
	fputs(sLogString.c_str(),m_fpLog);
	fputs("\n",m_fpLog);
	fflush(m_fpLog);
}
FEMSolver* CASFEM::GetFEMSolver() const
{
	return m_poFEMSolver;
}
void CASFEM::Run()
{
	WriteInLog("CASFEM Running");
	char cWrite[500];
	double dTime = m_poData->GetCurrentTime();
	double dTargetTime = m_poData->GetTargetTime();
	m_poFEMSolver->InitializeMatrices(m_poData->GetTimeStep());
	// correct the time step if required
	double dTimeStep = m_poFEMSolver->GetTimeStep();
	unsigned int iOutputCount = 0;
	while(dTime <= dTargetTime)
	{
		dTime = dTime + dTimeStep;
		m_poData->UpdateTime(dTime);
		m_poData->IncrementCurrentOutputCount();
		printf("time = %lf\n",dTime);
		m_poData->ApplyLoads(dTime);
		m_poFEMSolver->Solve(dTime);
		iOutputCount = m_poData->GetCurrentOutputCount();
		m_poFEMSolver->WriteFEMSolution(iOutputCount);
		m_poFEMSolver->WriteFEMSolutionToParaview(iOutputCount);
	}
}



void Fire()
{
	unsigned int iRows = 200;
	unsigned int iColumns = 200;
	SparseMatrix A(iRows);
	Matrix B(iRows,1);
	//Matrix oOrig = Matrix::GenerateRandomSymmetricPositiveDefiniteMatrix(iRows);
	
	A.Set(1,1,1.0);
	A.Set(1,2,1.0);
	A.Set(1,3,1.0);
	
	A.Set(2,1,2.0);
	A.Set(2,2,3.0);
	A.Set(2,3,4.0);
	
	A.Set(3,1,1.0);
	A.Set(3,2,7.0);
	A.Set(3,3,5.0);
	
	B.Set(1,1,3.0);
	B.Set(2,1,9.0);
	B.Set(3,1,13.0);
	
	unsigned int i = 0;
	unsigned int j = 0;
	double dTemp = 0.0;
	double dDensity = 0.8;
	double dMin = -1000.0;
	double dMax = 1000.0;
	double dDiagonalFactor = 1.0;
	
	for(i = 1 ; i <= iRows ; i++)
	{
		for(j = 1 ; j <= iColumns ; j++)
		{
			if(i == j)
			{
				dTemp = Randomizer::Random(dMin,dMax);
				A.Set(i,j,dDiagonalFactor*fabs(dTemp));
			}
			else
			{
				if(Randomizer::Random() <= dDensity)
				{
					A.Set(i,j,Randomizer::Random(dMin,dMax));
				}
			}
			//A.Set(i,j,oOrig.Get(i,j));
		}
		B.Set(i,1,Randomizer::Random(dMin,dMax));
	}
	
	A.InitializeGMRESMatrices();
	Matrix X = A.SolveGMRES(B);
	printf("%E\n",(A*X - B).GetNorm());
	
	for(i = 1 ; i <= iRows ; i++)
	{
		for(j = 1 ; j <= iColumns ; j++)
		{
			//printf("%25.20f\t",A.Get(i,j));
		}
		//printf("\t\t%25.20f\n",X.Get(i,1));
	}
	
// 	Matrix C = B*B.GetTranspose();
// 	C.SubtractFromIdentity();
//	printf("%E\n",(oQ*oR - A).GetNorm());
}


int main(int argc,char** argv)
{
//	GenerateCylinderTestCase();
        GenerateBlockTestCase();
 	if(argc != 2)
 	{
 		printf("arguments error\n");
 		printf("usage : casfem fem_input_file\n");
 		return 1;
	}
	string sFEMInputFile = argv[1];
 	clock_t tTime = clock();
 	
  	printf("CASFEM running\n");
  	CASFEM* poCASFEM = CASFEM::CreateInstance();
	if(!poCASFEM->Initialize(sFEMInputFile))
	{
		return 1;
	}
	poCASFEM->Run();
 	delete poCASFEM;
 
 	tTime = clock() - tTime;
 	printf("CASFEM has ended successfully in %E seconds.\n",(double)tTime/(double)CLOCKS_PER_SEC);
 	return 0;
}
 


