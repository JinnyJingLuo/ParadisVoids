// Ahmed M. Hussein

#ifndef DISLOCATIONNODE_H_
#define DISLOCATIONNODE_H_

#include "CategorizedGenericNode.h"
#include "Vector.h"
#include "list"
#include "DislocationNodeTag.h"

using namespace EZ;
using namespace std;

namespace DislocationSystem
{
	enum NodeCategory
	{
		UnconstrainedNode = 0,
		PinnedNode = 7,
		SurfaceNode = 8
	};
	class DislocationSegment;
	class DislocationNode : public CategorizedGenericNode
	{
	public:
	    DislocationNode();
        DislocationNode(const double& dX,const double& dY,const double& dZ);
        DislocationNode(const DislocationNode& oNode);
        DislocationNode(const CategorizedGenericNode& oNode);
        DislocationNode(const GenericNode& oNode);
        DislocationNode(const Point& oNode);
        virtual ~DislocationNode();
        DislocationNode& operator=(const DislocationNode& oNode);
        DislocationNode& operator=(const CategorizedGenericNode& oNode);
        DislocationNode& operator=(const GenericNode& oNode);
        DislocationNode& operator=(const Point& oNode);
        virtual void Reset();
        void SetDomainID(const unsigned int& iID);
        unsigned int GetDomainID() const;
        void SetSurfaceNormal(const Vector& oNormal);
        Vector GetSurfaceNormal() const;
        void SetForce(const Vector& oForce);
        Vector GetForce() const;
        void SetVelocity(const Vector& oVelocity);
        Vector GetVelocity() const;
		list<DislocationSegment*>* GetArms();
		void AddArm(DislocationSegment* poArm);
		void RemoveArmByTag(DislocationNodeTag* poTag);
		void RemoveArmByTag(const unsigned int& iDomainID,const unsigned int& iNodeID);
		void RemoveArmByEndNode(DislocationNode* poNode);
		DislocationSegment* GetArmByTag(const unsigned int& iDomainID,const unsigned int& iNodeID);
		void Write(FILE* fpFile) const;
		void Pack(list<double>* pldData) const;
		void Unpack(list<double>* pldData);
		void PackVelocity(list<double>* pldData) const;
		static DislocationNodeTag UnpackVelocity(list<double>* pldData,Vector& oVelocity);
		void ResetForces();
		void AddForce(const Vector& oForceIncrement);
		bool IsNodeSuperior(DislocationNode* poNode) const;
		void GetDynamicConstraint(unsigned int& iConstraintType,Vector& oConstraintVector) const;
		bool IsAllowedToCollide() const;
		void EnableCollision();
		void DisableCollision();
		unsigned int GetNeighboursCount() const;
		bool IsSameNode(DislocationNode* poNode) const;
		bool IsSameAsNode(DislocationNodeTag* poTag) const;
		bool IsMasterNode() const;
		
	private:
	
	protected:
		void Initialize();
		unsigned int m_iDomainID;
		Vector m_oSurfaceNormal;
		Vector m_oForce;
		Vector m_oVelocity;
		bool m_bAllowedToCollide;
		list<DislocationSegment*> m_lpoArms;
	};
}


#endif



