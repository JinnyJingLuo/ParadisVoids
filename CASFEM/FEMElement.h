// Ahmed M. Hussein

#ifndef FEMELEMENT_H_
#define FEMELEMENT_H_

#include "FEMElementGeometry.h"
#include "FEMMaterial.h"
#include "FEMGaussPoint.h"

using namespace GeometrySystem;

namespace FEMSystem
{
	enum FEMElementType
	{
		NullFEMElement = 0,
		PotentialFEMElement = 1,
		SolidFEMElement = 2,
		ThermoMechanicalFEMElement = 3
	};

	class FEMElement
	{
	public:
		virtual ~FEMElement();
		virtual void Reset();
		virtual FEMElement& operator=(const FEMElement& oElement);
		vector<FEMNode*>* GetNodes();
		FEMNode* GetNode(const unsigned int& iNodeIndex) const;
		unsigned int GetNodesCount() const;
		double GetVolume() const;
		virtual FEMElementType GetType() const = 0;
		virtual unsigned int GetDegreesOfFreedomCount() const = 0;
		virtual void ApplyLoads(const double& dTime) = 0;
		virtual void ReadLoads(FILE* fpFile,vector<FEMLoad*>* pvpoLoads) = 0;
		virtual void Write(FILE* fpFile) const = 0;
		FEMElementGeometry* GetGeometry() const;
		void SetGeometry(FEMElementGeometry* poElementGeometry);
		static FEMElement* CreateElementByType(FEMElementType eType);
		static FEMElement* CreateElementByTypeIndex(const unsigned int& iIndex);
		void DeleteGeometry();
		virtual void SetFacesLoadsSize() = 0;
		virtual vector<unsigned int> GetDOFIndices() const = 0;
		void SetMaterial(FEMMaterial* poMaterial);
		FEMMaterial* GetMaterial() const;
		void InitializeGaussPoints();
		void UpdateGaussPoints();
		
	private:
	
	protected:
		virtual void Initialize();
		FEMElementGeometry* m_poElementGeometry;
		FEMMaterial* m_poMaterial;
		vector< FEMGaussPoint > m_voGaussPoints;
	};
}


#endif




