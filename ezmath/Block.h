// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef BLOCK_H_
#define BLOCK_H_

#include "Geometry.h"

using namespace EZ;

namespace GeometrySystem
{
	class Block : public Geometry
	{
	public:
		Block();
		~Block();
		Block(const Block& oBlock);
		Block& operator=(const Block& oBlock);
		void Reset();
		
		double GetXLength() const;
		double GetYLength() const;
		double GetZLength() const;
		void SetXLength(const double& dLength);
		void SetYLength(const double& dLength);
		void SetZLength(const double& dLength);
		void SetResolution(const unsigned int& iXResolution,const unsigned int& iYResolution,const unsigned int& iZResolution);
		unsigned int GetXResolution() const;
		unsigned int GetYResolution() const;
		unsigned int GetZResolution() const;
		
		virtual Geometry* Clone();
		double GetVolume() const;
		bool IsPointInside(const Point& oPoint,const double& dTolerance = 1.0E-6) const;
	private:

	protected:
		double m_dXLength;
		double m_dYLength;
		double m_dZLength;
		unsigned int m_iXResolution;
		unsigned int m_iYResolution;
		unsigned int m_iZResolution;
	};
}

#endif



