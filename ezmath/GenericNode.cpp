// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "GenericNode.h"
#include "string"
#include "Tools.h"

using namespace std;
using namespace SupportSystem;

namespace EZ
{
	GenericNode::GenericNode()
	{
		Initialize();
	}
	GenericNode::GenericNode(const double& dX,const double& dY,const double& dZ)
	{
		Initialize();
		Set(dX,dY,dZ);
	}
	GenericNode::GenericNode(const GenericNode& oNode)
	{
		*this = oNode;
	}
	GenericNode::GenericNode(const Point& oNode)
	{
		*this = oNode;
	}
	GenericNode& GenericNode::operator=(const Point& oPoint)
	{
		Initialize();
		Point::operator =(oPoint);
		return *this;
	}
	GenericNode::~GenericNode()
	{
		Reset();
	}
	GenericNode& GenericNode::operator=(const GenericNode& oNode)
	{
		Point::operator =(oNode);
		GeometricComponent::operator =(oNode);
		return *this;
	}
	void GenericNode::Reset()
	{
		Point::Reset();
		GeometricComponent::Reset();

	}
	void GenericNode::Initialize()
	{
		m_dX = 0.0;
		m_dY = 0.0;
		m_dZ = 0.0;
		GeometricComponent::Initialize();
	}
}

