// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "FEMMesh.h"
#include "Tools.h"
#include "cmath"

namespace FEMSystem
{
	FEMMesh::FEMMesh()
	{
		Initialize();
 	}
	FEMMesh::FEMMesh(const FEMMesh& oMesh)
	{
 		*this = oMesh;
 	}
	FEMMesh::~FEMMesh()
	{
		Reset();
 	}
	FEMMesh& FEMMesh::operator=(const FEMMesh& oMesh)
	{
 		return *this;
 	}
 	void FEMMesh::Initialize()
 	{
 	
 	}
 	void FEMMesh::Reset()
 	{
 	
 	}
	void FEMMesh::GenerateMeshFromFile(FILE* fpFile,vector<FEMNode*>* pvpoNodes,vector<FEMElement*>* pvpoElements,vector<FEMLoad*>* pvpoLoads,vector<FEMMaterial*>* pvpoMaterials)
	{
 		unsigned int iNodesCount = 0;
 		string sRead = SupportSystem::GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%d\n",&iNodesCount);
 		unsigned int i = 0;
 		pvpoNodes->resize(iNodesCount);
 		double dX = 0.0;
 		double dY = 0.0;
 		double dZ = 0.0;
 		unsigned int iTemp = 0;
 		FEMNode* poNode = NULL;
 		for(i = 0; i < iNodesCount ; i++)
 		{
 			sRead = SupportSystem::GetRealString(500,fpFile);
 			sscanf(sRead.c_str(),"%d\n",&iTemp);
 			poNode = FEMNode::CreateNodeByTypeIndex(iTemp);
 			poNode->ReadNode(fpFile,pvpoLoads);
 			poNode->SetID(i + 1);
			pvpoNodes->at(i) = poNode;
 		}
 		unsigned int iElementsCount = 0;
 		sRead = SupportSystem::GetRealString(500,fpFile);
 		sscanf(sRead.c_str(),"%d\n",&iElementsCount);
 		pvpoElements->resize(iElementsCount);
		unsigned int iElementTypeIndex = 0;
		unsigned int iElementMaterial = 0;
		FEMElementGeometry* poGeometry = NULL;
		FEMElement* poElement = NULL;
 		for(i = 0; i < iElementsCount ; i++)
 		{
			sRead = SupportSystem::GetRealString(500,fpFile);
			// read the geometry
 			sscanf(sRead.c_str(),"%d,%d,%d\n",&iTemp,&iElementTypeIndex,&iElementMaterial);
			poGeometry = FEMElementGeometry::CreateElementGeometryByTypeIndex(iTemp);
			poGeometry->Read(fpFile,pvpoNodes);
			// create the element
			poElement = FEMElement::CreateElementByTypeIndex(iElementTypeIndex);
			poElement->SetGeometry(poGeometry);
			// set element material
			poElement->SetMaterial(pvpoMaterials->at(iElementMaterial - 1));
			// read the loads
			poElement->ReadLoads(fpFile,pvpoLoads);
			pvpoElements->at(i) = poElement;
			pvpoElements->at(i)->InitializeGaussPoints();
 		}
 		for(i = 0; i < iElementsCount ; i++)
 		{
 			pvpoElements->at(i)->GetGeometry()->UpdateAxisAlignedBoundingBox();
 		}
 	}
	void FEMMesh::GenerateMeshFromCylinder(Cylinder* poCylinder,vector<FEMNode*>& vpoNodes,vector<FEMElement*>& vpoElements,FEMNodeType eNodeType,FEMElementType eElementType)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		unsigned int iRadialResolution = poCylinder->GetRadialResolution();
		unsigned int iCircumferentialResolution = poCylinder->GetCircumferentialResolution();
		unsigned int iAxialResolution = poCylinder->GetAxialResolution();
		double dRadius = poCylinder->GetRadius();
		double dLength = poCylinder->GetLength();
		unsigned int iWorkingRLines = 2*iRadialResolution;
		unsigned int iWorkingThetaLines = 4*iCircumferentialResolution;
		unsigned int iWorkingZLines = 2*iAxialResolution - 1;

		unsigned int iMainPlaneSize = iCircumferentialResolution*(3*iWorkingRLines - 1) + 1;
		unsigned int iLayerSize = (4*iWorkingRLines - 1)*iCircumferentialResolution + 2;
		vpoNodes.clear();
		vpoElements.clear();
		vpoNodes.reserve(iLayerSize*(iAxialResolution - 1) + iMainPlaneSize);
		double dX = 0.0;
		double dY = 0.0;
		double dZ = 0.0;
		double dRIncrement = dRadius/(double)iWorkingRLines;
		double dThetaIncrement = 2*PI/(double)iWorkingThetaLines;
		double dZIncrement = dLength/(double)(iWorkingZLines - 1);
		double dTheta = 0.0;
		double dR = 0.0;
		FEMNode* poNode = NULL;
		for(k = 0; k < iWorkingZLines ; k++)
		{
			dZ = k*dZIncrement - 0.5*dLength;
			poNode = FEMNode::CreateNodeByType(eNodeType);
			poNode->Set(0.0,0.0,dZ);
			vpoNodes.push_back(poNode);
			if(k%2 == 0)
			{
				for(i = 1; i <= iWorkingRLines ; i++)
				{
					dR = i*dRIncrement;
					if(i == 1)
					{
						for(j = 0; j < iWorkingThetaLines ; j = j + 4)
						{
							dTheta = j*dThetaIncrement;
							dX = dR*cos(dTheta);
							dY = dR*sin(dTheta);
							poNode = FEMNode::CreateNodeByType(eNodeType);
							poNode->Set(dX,dY,dZ);
							vpoNodes.push_back(poNode);
						}
						continue;
					}
					if(i%2 != 0)
					{
						for(j = 0; j < iWorkingThetaLines ; j = j + 2)
						{
							dTheta = j*dThetaIncrement;
							dX = dR*cos(dTheta);
							dY = dR*sin(dTheta);
							poNode = FEMNode::CreateNodeByType(eNodeType);
							poNode->Set(dX,dY,dZ);
							vpoNodes.push_back(poNode);
						}
						continue;
					}
					for(j = 0; j < iWorkingThetaLines ; j = j++)
					{
						dTheta = j*dThetaIncrement;
						dX = dR*cos(dTheta);
						dY = dR*sin(dTheta);
						poNode = FEMNode::CreateNodeByType(eNodeType);
						poNode->Set(dX,dY,dZ);
						vpoNodes.push_back(poNode);
					}
				}
			}
			else
			{
				for(i = 2; i <= iWorkingRLines ; i = i + 2)
				{
					dR = i*dRIncrement;
					for(j = 0; j < iWorkingThetaLines ; j = j + 2)
					{
						dTheta = j*dThetaIncrement;
						dX = dR*cos(dTheta);
						dY = dR*sin(dTheta);
						poNode = FEMNode::CreateNodeByType(eNodeType);
						poNode->Set(dX,dY,dZ);
						vpoNodes.push_back(poNode);
					}
				}
			}
		}
		
		vpoElements.reserve(((iRadialResolution - 1)*2*iCircumferentialResolution + iCircumferentialResolution)*(iAxialResolution - 1));
		vector<FEMNode*> vpoTempNodes;
		unsigned int iNodesPerElementCount = 20;
		vpoTempNodes.resize(iNodesPerElementCount);
		FEMElementGeometry* poGeometry = NULL;
		// generate the cylinder core elements
		for(k = 0; k < iAxialResolution - 1 ; k++)
		{
			for(j = 0; j < iCircumferentialResolution ; j++)
			{
				if(j == iCircumferentialResolution - 1)
				{
					vpoTempNodes[0] = vpoNodes[k*iLayerSize];
					vpoTempNodes[1] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j];
					vpoTempNodes[2] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j + 2];
					vpoTempNodes[3] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1];
					vpoTempNodes[4] = vpoNodes[(k + 1)*iLayerSize];
					vpoTempNodes[5] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j];
					vpoTempNodes[6] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j + 2];
					vpoTempNodes[7] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1];
					vpoTempNodes[8] = vpoNodes[k*iLayerSize + 1 + j];
					vpoTempNodes[9] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j + 1];
					vpoTempNodes[10] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j + 3];
					vpoTempNodes[11] = vpoNodes[k*iLayerSize + 1];
					vpoTempNodes[12] = vpoNodes[(k + 1)*iLayerSize + 1 + j];
					vpoTempNodes[13] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j + 1];
					vpoTempNodes[14] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j + 3];
					vpoTempNodes[15] = vpoNodes[(k + 1)*iLayerSize + 1];
					vpoTempNodes[16] = vpoNodes[k*iLayerSize + iMainPlaneSize];
					vpoTempNodes[17] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*j];
					vpoTempNodes[18] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*j + 1];
					vpoTempNodes[19] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1];
				}
				else
				{
					vpoTempNodes[0] = vpoNodes[k*iLayerSize];
					vpoTempNodes[1] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j];
					vpoTempNodes[2] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j + 2];
					vpoTempNodes[3] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j + 4];
					vpoTempNodes[4] = vpoNodes[(k + 1)*iLayerSize];
					vpoTempNodes[5] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j];
					vpoTempNodes[6] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j + 2];
					vpoTempNodes[7] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j + 4];
					vpoTempNodes[8] = vpoNodes[k*iLayerSize + 1 + j];
					vpoTempNodes[9] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j + 1];
					vpoTempNodes[10] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 4*j + 3];
					vpoTempNodes[11] = vpoNodes[k*iLayerSize + 1 + j + 1];
					vpoTempNodes[12] = vpoNodes[(k + 1)*iLayerSize + 1 + j];
					vpoTempNodes[13] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j + 1];
					vpoTempNodes[14] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 4*j + 3];
					vpoTempNodes[15] = vpoNodes[(k + 1)*iLayerSize + 1 + j + 1];
					vpoTempNodes[16] = vpoNodes[k*iLayerSize + iMainPlaneSize];
					vpoTempNodes[17] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*j];
					vpoTempNodes[18] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*j + 1];
					vpoTempNodes[19] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*j + 2];
				}
				poGeometry = FEMElementGeometry::CreateElementGeometryByType(HexahedralFEMElement);
				poGeometry->Set(vpoTempNodes);
				vpoElements.push_back(FEMElement::CreateElementByType(eElementType));
				vpoElements.back()->SetGeometry(poGeometry);
			}
		}

		// generate the rest of the elements
		for(k = 0; k < iAxialResolution - 1 ; k++)
		{
			for(i = 1; i < iRadialResolution ; i++)
			{
				for(j = 0; j < iCircumferentialResolution ; j++)
				{
					// generate the first brick element
					vpoTempNodes[0] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j];
					vpoTempNodes[1] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j];
					vpoTempNodes[2] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2];
					vpoTempNodes[3] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2];
					vpoTempNodes[4] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j];
					vpoTempNodes[5] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j];
					vpoTempNodes[6] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2];
					vpoTempNodes[7] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2];
					vpoTempNodes[8] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j];
					vpoTempNodes[9] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 1];
					vpoTempNodes[10] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1];
					vpoTempNodes[11] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 1];
					vpoTempNodes[12] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j];
					vpoTempNodes[13] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 1];
					vpoTempNodes[14] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1];
					vpoTempNodes[15] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 1];
					vpoTempNodes[16] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j];
					vpoTempNodes[17] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 2*iCircumferentialResolution];
					vpoTempNodes[18] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 2*iCircumferentialResolution + 1];
					vpoTempNodes[19] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 1];
					poGeometry = FEMElementGeometry::CreateElementGeometryByType(HexahedralFEMElement);
					poGeometry->Set(vpoTempNodes);
					vpoElements.push_back(FEMElement::CreateElementByType(eElementType));
					vpoElements.back()->SetGeometry(poGeometry);		
					// generate the second brick element
					if(j == iCircumferentialResolution - 1)
					{
						vpoTempNodes[0] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[1] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[2] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution ];
						vpoTempNodes[3] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution];
						vpoTempNodes[4] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[5] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[6] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution];
						vpoTempNodes[7] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution];
						vpoTempNodes[8] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1];
						vpoTempNodes[9] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[10] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution];
						vpoTempNodes[11] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[12] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1];
						vpoTempNodes[13] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[14] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution];
						vpoTempNodes[15] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[16] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 1];
						vpoTempNodes[17] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 2*iCircumferentialResolution + 1];
						vpoTempNodes[18] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*iCircumferentialResolution];
						vpoTempNodes[19] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution];
					}
					else
					{
						vpoTempNodes[0] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[1] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[2] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2 + 2];
						vpoTempNodes[3] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2 + 2];
						vpoTempNodes[4] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[5] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2];
						vpoTempNodes[6] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 2 + 2];
						vpoTempNodes[7] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 2 + 2];
						vpoTempNodes[8] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1];
						vpoTempNodes[9] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[10] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1 + 1];
						vpoTempNodes[11] = vpoNodes[k*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[12] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1];
						vpoTempNodes[13] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*i*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[14] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*iCircumferentialResolution + 2*j + 1 + 1];
						vpoTempNodes[15] = vpoNodes[(k + 1)*iLayerSize + iCircumferentialResolution + 1 + 6*(i - 1)*iCircumferentialResolution + 4*j + 1 + 2];
						vpoTempNodes[16] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 1];
						vpoTempNodes[17] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 2*iCircumferentialResolution + 1];
						vpoTempNodes[18] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 2*iCircumferentialResolution + 1 + 1];
						vpoTempNodes[19] = vpoNodes[k*iLayerSize + iMainPlaneSize + 1 + 2*(i - 1)*iCircumferentialResolution + 2*j + 1 + 1];
					}

					poGeometry = FEMElementGeometry::CreateElementGeometryByType(HexahedralFEMElement);
					poGeometry->Set(vpoTempNodes);
					vpoElements.push_back(FEMElement::CreateElementByType(eElementType));
					vpoElements.back()->SetGeometry(poGeometry);
				}
			}
		}

//		// check surface nodes and apply coordinate system transformations
		unsigned int iSize = vpoNodes.size();
		double dTolerance = 1E-6*min(dRadius,dLength);
		for(i = 0; i < iSize ; i++)
		{
			dX = vpoNodes[i]->GetX();
			dY = vpoNodes[i]->GetY();
			dZ = vpoNodes[i]->GetZ();
			if(fabs(dZ - 0.5*dLength) < dTolerance || fabs(dZ + 0.5*dLength) < dTolerance)
			{
				vpoNodes[i]->SetOnSurface();
			}
			else if(fabs(dX*dX + dY*dY - dRadius*dRadius) < dTolerance)
			{
				vpoNodes[i]->SetOnSurface();
			}
			vpoNodes[i]->SetID(i + 1);
			vpoNodes[i]->SetX(vpoNodes[i]->GetX());
			vpoNodes[i]->SetY(vpoNodes[i]->GetY());
			vpoNodes[i]->SetZ(vpoNodes[i]->GetZ());
		}
		iSize = vpoElements.size();
		for(i = 0; i < iSize ; i++)
		{
			vpoElements[i]->GetGeometry()->UpdateAxisAlignedBoundingBox();
		}
	}
	void FEMMesh::GenerateMeshFromBlock(Block* poBlock,vector<FEMNode*>& vpoNodes,vector<FEMElement*>& vpoElements,FEMNodeType eNodeType,FEMElementType eElementType)
	{
 		vpoNodes.clear();
 		vpoElements.clear();
 		unsigned int i = 0;
 		unsigned int j = 0;
 		unsigned int k = 0;
 		unsigned int iXResolution = poBlock->GetXResolution();
		unsigned int iYResolution = poBlock->GetYResolution();
		unsigned int iZResolution = poBlock->GetZResolution();
		double dXLength = poBlock->GetXLength();
		double dYLength = poBlock->GetYLength();
		double dZLength = poBlock->GetZLength();
		
 		unsigned int iWorkingXLines = 2*iXResolution - 1;
 		unsigned int iWorkingYLines = 2*iYResolution - 1;
 		unsigned int iWorkingZLines = 2*iZResolution - 1;
 		double dXIncrement = dXLength/(double)(iWorkingXLines - 1);
		double dYIncrement = dYLength/(double)(iWorkingYLines - 1);
 		double dZIncrement = dZLength/(double)(iWorkingZLines - 1);
 		double dX = 0.0;
 		double dY = 0.0;
 		double dZ = 0.0;
 		unsigned int iIIncrement = 0;
 		unsigned int iJIncrement = 0;
 		double dHalfXLength = 0.5*dXLength;
 		double dHalfYLength = 0.5*dYLength;
 		double dHalfZLength = 0.5*dZLength;
 
 		unsigned int iMainPlaneSize = iXResolution*iWorkingYLines + (iXResolution - 1)*iYResolution;
 		unsigned int iLayerSize = iMainPlaneSize + iXResolution*iYResolution;
 		unsigned int iTotalBlockSize = iLayerSize*(iZResolution - 1) + iMainPlaneSize;
 		vpoNodes.reserve(iTotalBlockSize);

 		// generate the block nodes
 		FEMNode* poNode = NULL;
 		for(k = 0; k < iWorkingZLines ; k++)
 		{
 			dZ = k*dZIncrement - dHalfZLength;
 			iIIncrement = k%2 + 1;
 			for(i = 0; i < iWorkingXLines ; i = i + iIIncrement)
 			{
 				dX = i*dXIncrement - dHalfXLength;
 				if(k%2 == 0)
 				{
 					iJIncrement = i%2 + 1;
 				}
 				else
 				{
 					iJIncrement = iIIncrement;
 				}
 				for(j = 0; j < iWorkingYLines ; j = j + iJIncrement)
 				{
 					dY = j*dYIncrement - dHalfYLength;
					poNode = FEMNode::CreateNodeByType(eNodeType);
					poNode->Set(dX,dY,dZ);
 					vpoNodes.push_back(poNode);

 				}
 			}
 		}

 		unsigned int iTotalBlockElements = (iXResolution - 1)*(iYResolution - 1)*(iZResolution - 1);

				
 		vpoElements.reserve(iTotalBlockElements);
 		vector<FEMNode*> vpoTempNodes;
 		vpoTempNodes.resize(20);
		FEMElement* poElement = NULL;
		unsigned int iStripSize = iWorkingYLines + iYResolution;
		FEMElementGeometry* poGeometry = NULL;
 		// generate the block elements
 		for(k = 0; k < iZResolution - 1 ; k++)
 		{
 			for(i = 0; i < iXResolution - 1 ; i++)
 			{
 				for(j = 0; j < iYResolution - 1 ; j++)
 				{
 					vpoTempNodes[0] = vpoNodes[k*iLayerSize + i*iStripSize + 2*j];
 					vpoTempNodes[1] = vpoNodes[k*iLayerSize + (i + 1)*iStripSize + 2*j];
 					vpoTempNodes[2] = vpoNodes[k*iLayerSize + (i + 1)*iStripSize + 2*j + 2];
 					vpoTempNodes[3] = vpoNodes[k*iLayerSize + i*iStripSize + 2*j + 2];
 					vpoTempNodes[4] = vpoNodes[(k + 1)*iLayerSize + i*iStripSize + 2*j];
 					vpoTempNodes[5] = vpoNodes[(k + 1)*iLayerSize + (i + 1)*iStripSize + 2*j];
 					vpoTempNodes[6] = vpoNodes[(k + 1)*iLayerSize + (i + 1)*iStripSize + 2*j + 2];
 					vpoTempNodes[7] = vpoNodes[(k + 1)*iLayerSize + i*iStripSize + 2*j + 2];
 					vpoTempNodes[8] = vpoNodes[k*iLayerSize + i*iStripSize + iWorkingYLines + j];
 					vpoTempNodes[9] = vpoNodes[k*iLayerSize + (i + 1)*iStripSize + 2*j + 1];
 					vpoTempNodes[10] = vpoNodes[k*iLayerSize + i*iStripSize + iWorkingYLines + j + 1];
 					vpoTempNodes[11] = vpoNodes[k*iLayerSize + i*iStripSize + 2*j + 1];
 					vpoTempNodes[12] = vpoNodes[(k + 1)*iLayerSize + i*iStripSize + iWorkingYLines + j];
 					vpoTempNodes[13] = vpoNodes[(k + 1)*iLayerSize + (i + 1)*iStripSize + 2*j + 1];
 					vpoTempNodes[14] = vpoNodes[(k + 1)*iLayerSize + i*iStripSize + iWorkingYLines + j + 1];
 					vpoTempNodes[15] = vpoNodes[(k + 1)*iLayerSize + i*iStripSize + 2*j + 1];
 					vpoTempNodes[16] = vpoNodes[k*iLayerSize + iMainPlaneSize + i*iYResolution + j];
 					vpoTempNodes[17] = vpoNodes[k*iLayerSize + iMainPlaneSize + (i + 1)*iYResolution + j];
 					vpoTempNodes[18] = vpoNodes[k*iLayerSize + iMainPlaneSize + (i + 1)*iYResolution + j + 1];
 					vpoTempNodes[19] = vpoNodes[k*iLayerSize + iMainPlaneSize + i*iYResolution + j + 1];
 					poGeometry = FEMElementGeometry::CreateElementGeometryByType(HexahedralFEMElement);
					poGeometry->Set(vpoTempNodes);
					vpoElements.push_back(FEMElement::CreateElementByType(eElementType));
					vpoElements.back()->SetGeometry(poGeometry);
 				}
 			}
 		}

		// check surface nodes and apply coordinate system transformations
		unsigned int iSize = vpoNodes.size();
		double dTolerance = 1E-6*min(min(dXLength,dYLength),dZLength);
		for(i = 0; i < iSize ; i++)
		{
			dX = vpoNodes[i]->GetX();
			dY = vpoNodes[i]->GetY();
			dZ = vpoNodes[i]->GetZ();
			if(fabs(dX - dHalfXLength) < dTolerance || fabs(dX + dHalfXLength) < dTolerance)
			{
				vpoNodes[i]->SetOnSurface();
			}
			
			if(fabs(dY - dHalfYLength) < dTolerance || fabs(dY + dHalfYLength) < dTolerance)
			{
				vpoNodes[i]->SetOnSurface();
			}
			
			if(fabs(dZ - dHalfZLength) < dTolerance || fabs(dZ + dHalfZLength) < dTolerance)
			{
				vpoNodes[i]->SetOnSurface();
			}
			vpoNodes[i]->SetID(i + 1);
			vpoNodes[i]->SetX(vpoNodes[i]->GetX());
			vpoNodes[i]->SetY(vpoNodes[i]->GetY());
			vpoNodes[i]->SetZ(vpoNodes[i]->GetZ());
		}
		iSize = vpoElements.size();
		for(i = 0; i < iSize ; i++)
		{
			vpoElements[i]->GetGeometry()->UpdateAxisAlignedBoundingBox();
		}

	}

}




 
