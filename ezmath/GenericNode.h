// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef GENERICNODE_H
#define	GENERICNODE_H

#include "Point.h"
#include "GeometricComponent.h"
#include "iostream"

using namespace std;
using namespace GeometrySystem;

namespace EZ
{
    class GenericNode : public Point, public GeometricComponent
    {
    public:
        GenericNode();
        GenericNode(const double& dX,const double& dY,const double& dZ);
        GenericNode(const GenericNode& oNode);
        GenericNode(const Point& oNode);
        virtual ~GenericNode();
        GenericNode& operator=(const GenericNode& oNode);
        GenericNode& operator=(const Point& oPoint);
        virtual void Reset();
    private:
        
    protected:
		virtual void Initialize();
    };
}

#endif


