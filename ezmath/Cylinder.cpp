// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Cylinder.h"
#include "math.h"
#include "Tools.h"


namespace GeometrySystem
{
	Cylinder::Cylinder()
	{
		Reset();
	}
	Cylinder::~Cylinder()
	{
		// do nothing
	}
	Cylinder::Cylinder(const Cylinder& oCylinder)
	{
		*this = oCylinder;
	}
	Cylinder& Cylinder::operator=(const Cylinder& oCylinder)
	{
		m_oSystem = oCylinder.m_oSystem;
		m_dRadius = oCylinder.m_dRadius;
		m_dLength = oCylinder.m_dLength;
		m_iRadialResolution = oCylinder.m_iRadialResolution;
		m_iCircumferentialResolution = oCylinder.m_iCircumferentialResolution;
		m_iAxialResolution = oCylinder.m_iAxialResolution;
		return *this;
	}
	void Cylinder::Reset()
	{
		m_dRadius = 0.0;
		m_dLength = 0.0;
		m_iRadialResolution = 0;
		m_iCircumferentialResolution = 0;
		m_iAxialResolution = 0;
		m_oSystem.Reset();
	}
//	// void Cylinder::GenerateMesh(vector<Node*>& vpoNodes,vector<TetrahedralElement*>& vpoElements)
//	// {
//		// unsigned int i = 0;
//		// unsigned int j = 0;
//		// unsigned int k = 0;
//		// unsigned int iWorkingRLines = 2*m_iRadialResolution;
//		// unsigned int iWorkingThetaLines = 2*m_iCircumferentialResolution;
//		// unsigned int iWorkingZLines = 2*m_iAxialResolution - 1;

//		// unsigned int iPlaneSize = (iWorkingRLines - 0.5)*iWorkingThetaLines + 1;
//		// vpoNodes.clear();
//		// vpoElements.clear();
//		// vpoNodes.reserve(iPlaneSize*iWorkingZLines);
//		// double dX = 0.0;
//		// double dY = 0.0;
//		// double dZ = 0.0;
//		// double dRIncrement = m_dRadius/(double)iWorkingRLines;
//		// double dThetaIncrement = 2*PI/(double)iWorkingThetaLines;
//		// double dZIncrement = m_dLength/(double)(iWorkingZLines - 1);
//		// double dTheta = 0.0;
//		// double dR = 0.0;
//		// double dTemp = 0.0;
//		// Node* poTempNode1 = NULL;
//		// Node* poTempNode2 = NULL;
//		// Node* poTempNode3 = NULL;

//		// for(k = 0; k < iWorkingZLines ; k++)
//		// {
//			// dZ = k*dZIncrement - 0.5*m_dLength;
//			// vpoNodes.push_back(new Node(0.0,0.0,dZ));
//			// for(i = 1; i <= iWorkingRLines ; i++)
//			// {
//				// dR = i*dRIncrement;
//				// if(i == 1)
//				// {
//					// for(j = 0; j < iWorkingThetaLines ; j = j + 2)
//					// {
//						// dTheta = j*dThetaIncrement;
//						// dX = dR*cos(dTheta);
//						// dY = dR*sin(dTheta);
//						// vpoNodes.push_back(new Node(dX,dY,dZ));
//					// }
//					// continue;
//				// }
//				// for(j = 0; j < iWorkingThetaLines ; j = j + 2)
//				// {
//					// dTheta = j*dThetaIncrement;
//					// dX = dR*cos(dTheta);
//					// dY = dR*sin(dTheta);
//					// vpoNodes.push_back(new Node(dX,dY,dZ));

//					// dTemp = dTheta + 2*dThetaIncrement;
//					// dX = 0.5*dR*(cos(dTheta) + cos(dTemp));
//					// dY = 0.5*dR*(sin(dTheta) + sin(dTemp));
//					// vpoNodes.push_back(new Node(dX,dY,dZ));
//				// }
//			// }
//			// // adjust nodal positions to generate flat surface elements
//			// for(i = 3; i <= iWorkingRLines ; i = i + 2)
//			// {
//				// for(j = 1; j < iWorkingThetaLines ; j = j + 2)
//				// {
//					// poTempNode1 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i - 3)*iWorkingThetaLines + j - 1];
//					// if(j == iWorkingThetaLines - 1)
//					// {
//						// poTempNode2 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i - 1)*iWorkingThetaLines];
//					// }
//					// else
//					// {
//						// poTempNode2 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i - 1)*iWorkingThetaLines + j + 1];
//					// }
//					// poTempNode3 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i - 2)*iWorkingThetaLines + j];
//					// poTempNode3->SetX(0.5*(poTempNode1->GetX() + poTempNode2->GetX()));
//					// poTempNode3->SetY(0.5*(poTempNode1->GetY() + poTempNode2->GetY()));
//				// }
//			// }
//		// }

//		// vpoElements.reserve(3*(m_iAxialResolution - 1)*(2*m_iRadialResolution - 1)*m_iCircumferentialResolution);

//		// Node* poNode1 = NULL;
//		// Node* poNode2 = NULL;
//		// Node* poNode3 = NULL;
//		// Node* poNode4 = NULL;
//		// Node* poNode5 = NULL;
//		// Node* poNode6 = NULL;
//		// Node* poNode7 = NULL;
//		// Node* poNode8 = NULL;
//		// Node* poNode9 = NULL;
//		// Node* poNode10 = NULL;
//		// Node* poNode11 = NULL;
//		// Node* poNode12 = NULL;
//		// Node* poNode13 = NULL;
//		// Node* poNode14 = NULL;
//		// Node* poNode15 = NULL;
//		// Node* poNode16 = NULL;
//		// Node* poNode17 = NULL;
//		// Node* poNode18 = NULL;
//		// Node* poNode19 = NULL;
//		// Node* poNode20 = NULL;
//		// Node* poNode21 = NULL;
//		// Node* poNode22 = NULL;
//		// Node* poNode23 = NULL;
//		// Node* poNode24 = NULL;
//		// Node* poNode25 = NULL;
//		// Node* poNode26 = NULL;
//		// Node* poNode27 = NULL;
//		// vector<Node*> vpoTempNodes;
//		// vpoTempNodes.resize(NodesPerElementCount);

//		// // generate cylinder core elements
//		// for(k = 0; k < iWorkingZLines - 1 ; k = k + 2)
//		// {
//			// for(j = 0; j < iWorkingThetaLines ; j = j + 2)
//			// {
//				// poNode1 = vpoNodes[k*iPlaneSize];
//				// poNode2 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + j + 1];
//				// poNode3 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + j + 3];
//				// poNode4 = vpoNodes[(k + 2)*iPlaneSize];
//				// poNode5 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + j + 1];
//				// poNode6 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + j + 3];

//				// poNode7 = vpoNodes[(k + 1)*iPlaneSize];
//				// poNode8 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + j + 1];
//				// poNode9 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + j + 3];

//				// poNode10 = vpoNodes[k*iPlaneSize + 0.5*j + 1];
//				// poNode11 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + j + 2];
//				// poNode12 = vpoNodes[k*iPlaneSize + 0.5*j + 2];

//				// poNode13 = vpoNodes[(k + 2)*iPlaneSize + 0.5*j + 1];
//				// poNode14 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + j + 2];
//				// poNode15 = vpoNodes[(k + 2)*iPlaneSize + 0.5*j + 2];

//				// poNode16 = vpoNodes[(k + 1)*iPlaneSize + 0.5*j + 1];
//				// poNode17 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + j + 2];
//				// poNode18 = vpoNodes[(k + 1)*iPlaneSize + 0.5*j + 2];

//				// if(j == iWorkingThetaLines - 2)
//				// {
//					// poNode3 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1];
//					// poNode6 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1];
//					// poNode9 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1];
//					// poNode12 = vpoNodes[k*iPlaneSize + 1];
//					// poNode15 = vpoNodes[(k + 2)*iPlaneSize + 1];
//					// poNode18 = vpoNodes[(k + 1)*iPlaneSize + 1];
//				// }
//				// vpoTempNodes[0] = poNode1;
//				// vpoTempNodes[1] = poNode2;
//				// vpoTempNodes[2] = poNode3;
//				// vpoTempNodes[3] = poNode5;
//				// vpoTempNodes[4] = poNode10;
//				// vpoTempNodes[5] = poNode11;
//				// vpoTempNodes[6] = poNode12;
//				// vpoTempNodes[7] = poNode16;
//				// vpoTempNodes[8] = poNode8;
//				// vpoTempNodes[9] = poNode17;
//				// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));

//				// vpoTempNodes[0] = poNode1;
//				// vpoTempNodes[1] = poNode4;
//				// vpoTempNodes[2] = poNode5;
//				// vpoTempNodes[3] = poNode6;
//				// vpoTempNodes[4] = poNode7;
//				// vpoTempNodes[5] = poNode13;
//				// vpoTempNodes[6] = poNode16;
//				// vpoTempNodes[7] = poNode18;
//				// vpoTempNodes[8] = poNode15;
//				// vpoTempNodes[9] = poNode14;
//				// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));

//				// vpoTempNodes[0] = poNode1;
//				// vpoTempNodes[1] = poNode3;
//				// vpoTempNodes[2] = poNode5;
//				// vpoTempNodes[3] = poNode6;
//				// vpoTempNodes[4] = poNode12;
//				// vpoTempNodes[5] = poNode17;
//				// vpoTempNodes[6] = poNode16;
//				// vpoTempNodes[7] = poNode18;
//				// vpoTempNodes[8] = poNode9;
//				// vpoTempNodes[9] = poNode14;
//				// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));
//			// }
//		// }
//		// // generate the rest of the cylinder elements
//		// for(k = 0; k < iWorkingZLines - 1 ; k = k + 2)
//		// {
//			// for(i = 0; i < iWorkingRLines - 2 ; i = i + 2)
//			// {
//				// for(j = 0; j < iWorkingThetaLines ; j = j + 2)
//				// {
//					// poNode1 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j];
//					// poNode2 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j];
//					// //poNode3 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 2];
//					// //poNode4 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 2];
//					// poNode5 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j];
//					// poNode6 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j];
//					// //poNode7 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 2];
//					// //poNode8 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 2];

//					// poNode9 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j];
//					// poNode10 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 1];
//					// //poNode11 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 2];
//					// poNode12 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 1];

//					// poNode13 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j];
//					// poNode14 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 1];
//					// //poNode15 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 2];
//					// poNode16 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 1];

//					// poNode17 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j];
//					// poNode18 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j];
//					// //poNode19 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 2];
//					// //poNode20 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 2];

//					// poNode21 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j];
//					// poNode22 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 1];
//					// //poNode23 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 2];
//					// poNode24 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 1];
//					// poNode25 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 1];
//					// poNode26 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 1];
//					// poNode27 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 1];

//					// if(j == iWorkingThetaLines - 2)
//					// {
//						// poNode3 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines];
//						// poNode4 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines];
//						// poNode7 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines];
//						// poNode8 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines];
//						// poNode11 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines];
//						// poNode15 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines];
//						// poNode19 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines];
//						// poNode20 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines];
//						// poNode23 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines];
//					// }
//					// else
//					// {
//						// poNode3 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 2];
//						// poNode4 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 2];
//						// poNode7 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 2];
//						// poNode8 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 2];
//						// poNode11 = vpoNodes[k*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 2];
//						// poNode15 = vpoNodes[(k + 2)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 2];
//						// poNode19 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 2)*iWorkingThetaLines + j + 2];
//						// poNode20 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + i*iWorkingThetaLines + j + 2];
//						// poNode23 = vpoNodes[(k + 1)*iPlaneSize + 0.5*iWorkingThetaLines + 1 + (i + 1)*iWorkingThetaLines + j + 2];
//					// }

//					// vpoTempNodes[0] = poNode1;
//					// vpoTempNodes[1] = poNode2;
//					// vpoTempNodes[2] = poNode3;
//					// vpoTempNodes[3] = poNode6;
//					// vpoTempNodes[4] = poNode9;
//					// vpoTempNodes[5] = poNode10;
//					// vpoTempNodes[6] = poNode25;
//					// vpoTempNodes[7] = poNode21;
//					// vpoTempNodes[8] = poNode18;
//					// vpoTempNodes[9] = poNode22;
//					// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));

//					// vpoTempNodes[0] = poNode1;
//					// vpoTempNodes[1] = poNode5;
//					// vpoTempNodes[2] = poNode6;
//					// vpoTempNodes[3] = poNode7;
//					// vpoTempNodes[4] = poNode17;
//					// vpoTempNodes[5] = poNode13;
//					// vpoTempNodes[6] = poNode21;
//					// vpoTempNodes[7] = poNode27;
//					// vpoTempNodes[8] = poNode26;
//					// vpoTempNodes[9] = poNode14;
//					// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));

//					// vpoTempNodes[0] = poNode1;
//					// vpoTempNodes[1] = poNode3;
//					// vpoTempNodes[2] = poNode6;
//					// vpoTempNodes[3] = poNode7;
//					// vpoTempNodes[4] = poNode25;
//					// vpoTempNodes[5] = poNode22;
//					// vpoTempNodes[6] = poNode21;
//					// vpoTempNodes[7] = poNode27;
//					// vpoTempNodes[8] = poNode19;
//					// vpoTempNodes[9] = poNode14;
//					// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));

//					// vpoTempNodes[0] = poNode1;
//					// vpoTempNodes[1] = poNode3;
//					// vpoTempNodes[2] = poNode4;
//					// vpoTempNodes[3] = poNode8;
//					// vpoTempNodes[4] = poNode25;
//					// vpoTempNodes[5] = poNode11;
//					// vpoTempNodes[6] = poNode12;
//					// vpoTempNodes[7] = poNode24;
//					// vpoTempNodes[8] = poNode23;
//					// vpoTempNodes[9] = poNode20;
//					// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));

//					// vpoTempNodes[0] = poNode1;
//					// vpoTempNodes[1] = poNode5;
//					// vpoTempNodes[2] = poNode8;
//					// vpoTempNodes[3] = poNode7;
//					// vpoTempNodes[4] = poNode17;
//					// vpoTempNodes[5] = poNode16;
//					// vpoTempNodes[6] = poNode24;
//					// vpoTempNodes[7] = poNode27;
//					// vpoTempNodes[8] = poNode26;
//					// vpoTempNodes[9] = poNode15;
//					// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));

//					// vpoTempNodes[0] = poNode1;
//					// vpoTempNodes[1] = poNode3;
//					// vpoTempNodes[2] = poNode7;
//					// vpoTempNodes[3] = poNode8;
//					// vpoTempNodes[4] = poNode25;
//					// vpoTempNodes[5] = poNode19;
//					// vpoTempNodes[6] = poNode27;
//					// vpoTempNodes[7] = poNode24;
//					// vpoTempNodes[8] = poNode23;
//					// vpoTempNodes[9] = poNode15;
//					// vpoElements.push_back(new FEMSystem::TetrahedralElement(vpoTempNodes));
//				// }
//			// }
//		// }
//		// // check surface nodes and apply coordinate system transformations
//		// unsigned int iSize = vpoNodes.size();
//		// double dTolerance = 1E-6*min(m_dRadius,m_dLength);
//		// Point oP;
//		// for(i = 0; i < iSize ; i++)
//		// {
//			// dX = vpoNodes[i]->GetX();
//			// dY = vpoNodes[i]->GetY();
//			// dZ = vpoNodes[i]->GetZ();
//			// if(fabs(dZ - 0.5*m_dLength) < dTolerance || fabs(dZ + 0.5*m_dLength) < dTolerance)
//			// {
//				// vpoNodes[i]->SetOnSurface();
//			// }
//			// else if(fabs(dX*dX + dY*dY - m_dRadius*m_dRadius) < dTolerance)
//			// {
//				// vpoNodes[i]->SetOnSurface();
//			// }
//			// vpoNodes[i]->SetID(i + 1);
//			// oP = m_oSystem.GetInGlobalCoordinates(*vpoNodes[i]);
//			// vpoNodes[i]->SetX(oP.GetX());
//			// vpoNodes[i]->SetY(oP.GetY());
//			// vpoNodes[i]->SetZ(oP.GetZ());
//		// }
//		// iSize = vpoElements.size();
//		// for(i = 0; i < iSize ; i++)
//		// {
//			// vpoElements[i]->SetID(i + 1);
//		// }
//	// }
	bool Cylinder::IsOnSurface(const Point& oPoint) const
	{
		Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
		double dX = oLocalPoint.GetX();
		double dY = oLocalPoint.GetY();
		double dZ = oLocalPoint.GetZ();
		double dR = dX*dX + dY*dY;
		double dToleranceFactor = 1.0E-6;
		double dZTolerance = dToleranceFactor*m_dLength;
		double dRTolerance = dToleranceFactor*m_dRadius;
		if(fabs(dZ - 0.5*m_dLength) < dZTolerance || fabs(dZ + 0.5*m_dLength) < dZTolerance)
		{
			if(dR <= m_dRadius*m_dRadius)
			{
				return true;
			}
			return false;
		}
		
		if(fabs(dR - m_dRadius*m_dRadius) < dRTolerance)
		{
			if(fabs(dZ) < 0.5*m_dLength)
			{
				return true;
			}
			return false;
		}
		return false;
	}
	double Cylinder::GetRadius() const
	{
		return m_dRadius;
	}
	double Cylinder::GetLength() const
	{
		return m_dLength;
	}
	void Cylinder::SetRadius(const double& dRadius)
	{
		m_dRadius = dRadius;
	}
	void Cylinder::SetLength(const double& dLength)
	{
		m_dLength = dLength;
	}
	void Cylinder::SetResolution(const unsigned int& iRadialResolution,const unsigned int& iCircumferentialResolution,const unsigned int& iAxialResolution)
	{
		m_iRadialResolution = iRadialResolution;
		m_iCircumferentialResolution = iCircumferentialResolution;
		m_iAxialResolution = iAxialResolution;
	}
	unsigned int Cylinder::GetRadialResolution() const
	{
		return m_iRadialResolution;
	}
	unsigned int Cylinder::GetCircumferentialResolution() const
	{
		return m_iCircumferentialResolution;
	}
	unsigned int Cylinder::GetAxialResolution() const
	{
		return m_iAxialResolution;
	}
	bool Cylinder::IsOnLowerFace(const Point& oPoint) const
	{
		if(IsOnSurface(oPoint))
		{
			Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
			double dZTolerance = m_dLength/(double)m_iAxialResolution/100.0;
			if(fabs(oLocalPoint.GetZ() + 0.5*m_dLength) < dZTolerance)
			{
				return true;
			}
		}
		return false;
	}
	bool Cylinder::IsOnUpperFace(const Point& oPoint) const
	{
		if(IsOnSurface(oPoint))
		{
			Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
			double dZTolerance = m_dLength/(double)m_iAxialResolution/100.0;
			if(fabs(oLocalPoint.GetZ() - 0.5*m_dLength) < dZTolerance)
			{
				return true;
			}
		}
		return false;
	}
	bool Cylinder::IsOnLateralFace(const Point& oPoint) const
	{
		if(IsOnSurface(oPoint))
		{
			Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
			double dX = oLocalPoint.GetX();
			double dY = oLocalPoint.GetY();
			double dR = dX*dX + dY*dY;
			dR = sqrt(dR);
			double dRTolerance = m_dRadius/(double)m_iRadialResolution/5.0;
			if(fabs(dR - m_dRadius) < dRTolerance)
			{
				return true;
			}
		}
		return false;
	}
	bool Cylinder::IsCenterOfLowerFace(const Point& oPoint) const
	{
		if(IsOnLowerFace(oPoint))
		{
			Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
			double dTolerance = 1E-6*min(m_dRadius,m_dLength);
			double dX = oLocalPoint.GetX();
			double dY = oLocalPoint.GetY();
			if((dX*dX + dY*dY) < dTolerance*dTolerance)
			{
				return true;
			}
		}
		return false;
	}
	bool Cylinder::IsCenterOfUpperFace(const Point& oPoint) const
	{
		if(IsOnUpperFace(oPoint))
		{
			Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
			double dTolerance = 1E-6*min(m_dRadius,m_dLength);
			double dX = oLocalPoint.GetX();
			double dY = oLocalPoint.GetY();
			if((dX*dX + dY*dY) < dTolerance*dTolerance)
			{
				return true;
			}
		}
		return false;
	}
	Geometry* Cylinder::Clone()
	{
		return new Cylinder(*this);
	}
	double Cylinder::GetVolume() const
	{
		double dBaseArea = PI*m_dRadius*m_dRadius;
		return dBaseArea*m_dLength;
	}
	bool Cylinder::IsPointInside(const Point& oPoint,const double& dTolerance) const
	{
		Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
		double dX = oLocalPoint.GetX();
		double dY = oLocalPoint.GetY();
		double dRadiusSquared = dX*dX + dY*dY;
		if(dRadiusSquared > m_dRadius*m_dRadius)
		{
			return false;
		}
		double dZ = oLocalPoint.GetZ();
		if(fabs(dZ) > 0.5*m_dLength)
		{
			return false;
		}
		return true;
	}
}



