// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Block.h"
#include "math.h"
#include "Tools.h"


namespace GeometrySystem
{
	Block::Block()
	{
		Reset();
	}
	Block::~Block()
	{
		// do nothing
	}
	Block::Block(const Block& oBlock)
	{
		*this = oBlock;
	}
	Block& Block::operator=(const Block& oBlock)
	{
		m_oSystem = oBlock.m_oSystem;
		m_dXLength = oBlock.m_dXLength;
		m_dYLength = oBlock.m_dYLength;
		m_dZLength = oBlock.m_dZLength;
		m_iXResolution = oBlock.m_iXResolution;
		m_iYResolution = oBlock.m_iYResolution;
		m_iZResolution = oBlock.m_iZResolution;
		return *this;
	}
	void Block::Reset()
	{
		m_dXLength = 0.0;
		m_dYLength = 0.0;
		m_dZLength = 0.0;
		m_iXResolution = 0;
		m_iYResolution = 0;
		m_iZResolution = 0;
		m_oSystem.Reset();
	}
	double Block::GetXLength() const
	{
		return m_dXLength;
	}
	double Block::GetYLength() const
	{
		return m_dYLength;
	}
	double Block::GetZLength() const
	{
		return m_dZLength;
	}
	void Block::SetXLength(const double& dLength)
	{
		m_dXLength = dLength;
	}
	void Block::SetYLength(const double& dLength)
	{
		m_dYLength = dLength;
	}
	void Block::SetZLength(const double& dLength)
	{
		m_dZLength = dLength;
	}
	void Block::SetResolution(const unsigned int& iXResolution,const unsigned int& iYResolution,const unsigned int& iZResolution)
	{
		m_iXResolution = iXResolution;
		m_iYResolution = iYResolution;
		m_iZResolution = iZResolution;
	}
	unsigned int Block::GetXResolution() const
	{
		return m_iXResolution;
	}
	unsigned int Block::GetYResolution() const
	{
		return m_iYResolution;
	}
	unsigned int Block::GetZResolution() const
	{
		return m_iZResolution;
	}
	Geometry* Block::Clone()
	{
		return new Block(*this);
	}
	double Block::GetVolume() const
	{
		return m_dXLength*m_dYLength*m_dZLength;
	}
	bool Block::IsPointInside(const Point& oPoint,const double& dTolerance) const
	{
		Point oLocalPoint = m_oSystem.GetInLocalCoordinates(oPoint);
		if(fabs(oLocalPoint.GetX()) > 0.5*m_dXLength)
		{
			return false;
		}
		if(fabs(oLocalPoint.GetY()) > 0.5*m_dYLength)
		{
			return false;
		}
		if(fabs(oLocalPoint.GetZ()) > 0.5*m_dZLength)
		{
			return false;
		}
		return true;
	}
}



