#include "Face.h"


namespace GeometrySystem
{
	Face::Face()
	{
		Initialize();
	}
	Face::Face(const Face& oFace)
	{
		*this = oFace;
	}
	Face::~Face()
	{
		Reset();
	}
	Face& Face::operator=(const Face& oFace)
	{
		GeometricComponent::operator=(oFace);
		return *this;
	}
	void Face::Reset()
	{
		GeometricComponent::Reset();
	}
	void Face::Initialize()
	{
		GeometricComponent::Initialize();
	}
}


