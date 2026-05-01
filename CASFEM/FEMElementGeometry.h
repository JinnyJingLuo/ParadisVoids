// Ahmed M. Hussein

#ifndef FEMELEMENTGEOMETRY_H_
#define FEMELEMENTGEOMETRY_H_

#include "Matrix.h"
#include "FEMNode.h"
#include "AxisAlignedBoundingBox.h"
#include "GenericNode.h"


using namespace GeometrySystem;

namespace FEMSystem
{
	enum FEMElementGeometryType
	{
		NullFEMElementGeometry = 0,
		TetrahedralFEMElement = 1,
		HexahedralFEMElement = 2
	};
	
	class FEMElementGeometry
	{
	public:
		virtual ~FEMElementGeometry();
		virtual FEMElementGeometry& operator=(const FEMElementGeometry& oElementGeometry);
		virtual void Reset();
		vector<FEMNode*>* GetNodes();
		FEMNode* GetNode(const unsigned int& iNodeIndex) const;
		bool IsInAxisAlignedBoundingBox(Point* poPoint) const;
		virtual FEMElementGeometryType GetType() const = 0;
		virtual unsigned int GetNodesCount() const = 0;
		virtual unsigned int GetFacesCount() const = 0;
		
		void Set(const vector<FEMNode*>& vpoNodes);
		virtual double GetVolume() const = 0;
		virtual Point GetPoint(const Matrix& oNaturalCoordinates) const = 0;
		virtual vector<Point> GetPoints(const Matrix& oNaturalCoordinates) const = 0;
		virtual bool IsPointInside(Point* poPoint,vector<double>& vdNaturalCoordinates) const = 0;
		virtual double GetDistanceToElementCenterPoint(Point* poPoint) const = 0;
		virtual vector<double> GetNearestNodeNaturalCoordinates(Point* poPoint) const = 0;
		virtual bool IsFaceOnSurface(const unsigned int& iFaceID) const = 0;
		virtual Point GetFaceCenter(const unsigned int& iFaceID) const = 0;
		virtual vector<double> GetNaturalCoordinates(Point* poPoint,const unsigned int& iMaxIterationsCount = 750) const = 0;
		virtual bool IsSurfaceElement() const = 0;
		virtual bool IsValid() const = 0;
		virtual Point GetCenterPoint() const = 0;
		void UpdateAxisAlignedBoundingBox();
		virtual vector<unsigned int> GetSurfaceFacesIndices() const = 0;
		virtual Vector GetFaceOuterNormal(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const = 0;
		virtual Point GetPointOnFace(const unsigned int& iFaceID,const Matrix& oNaturalCoordinates) const = 0;
		virtual void Read(FILE* fpFile,const vector<FEMNode*>* pvpoNodes) = 0;
		virtual void Write(FILE* fpFile) const = 0;
		virtual void GenerateSurfacePatches(vector< vector<unsigned int> >& vviNodesIndices,vector<GenericNode*>& vpoMidPoints) const = 0; 
		static FEMElementGeometry* CreateElementGeometryByType(FEMElementGeometryType eType);
		static FEMElementGeometry* CreateElementGeometryByTypeIndex(const unsigned int& iIndex);
		virtual Matrix GetShapeFunctions(const double& dXi,const double& dEta,const double& dZeta) const = 0;
		virtual Matrix GetShapeFunctionsXYZDerivatives(const double& dXi,const double& dEta,const double& dZeta,double& dJacobian) const = 0;
		virtual double GetJacobianMatrixDeterminant(const double& dXi,const double& dEta,const double& dZeta) const = 0;
		virtual Matrix GetShapeFunctionsExtrapolationsMatrix() const = 0;
		AxisAlignedBoundingBox* GetBox();
		virtual bool IsMidSideNode(const unsigned int& iIndex) const = 0;
		virtual double GetMassLumpingFactor() const = 0;
		
		virtual vector< vector<double> > GetBodyGaussPointsCoordinates() const = 0;
		virtual vector<double> GetBodyGaussPointsWeights() const = 0;
		virtual vector< vector<double> > GetFaceGaussPointsCoordinates(const unsigned int& iFaceIndex) const = 0;
		virtual vector<double> GetFaceGaussPointsWeights() const = 0;
		virtual Matrix GetFaceNaturalCoordinates(const unsigned int& iFaceID,const double& dXi,const double& dEta) const = 0;
		virtual FEMElementGeometry* Clone() const = 0;

	private:
	
	protected:
		virtual void Initialize();
		vector<FEMNode*> m_vpoNodes;
		AxisAlignedBoundingBox m_oBoundingBox;
	};
}

#endif



