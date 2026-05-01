#ifndef VOXELARRAY_H_
#define VOXELARRAY_H_

#include "AxisAlignedBoundingBox.h"
#include "math.h"

using namespace GeometrySystem;

template <typename Type> class VoxelArray
{
public:
	VoxelArray()
	{
		Initialize();
	}
	VoxelArray(const VoxelArray& oArray)
	{
		*this = oArray;
	}
	~VoxelArray()
	{
		Reset();
	}
	VoxelArray& operator=(const VoxelArray& oArray)
	{
		Reset();
		m_dXMin = oArray.m_dXMin;
		m_dYMin = oArray.m_dYMin;
		m_dZMin = oArray.m_dZMin;
		m_dDeltaX = oArray.m_dDeltaX;
		m_dDeltaY = oArray.m_dDeltaY;
		m_dDeltaZ = oArray.m_dDeltaZ;
		SetResolution(oArray.m_iXResolution,oArray.m_iYResolution,oArray.m_iZResolution);
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		for(i = 0 ; i < m_iXResolution ; i++)
		{
			for(j = 0 ; j < m_iYResolution ; j++)
			{
				for(k = 0 ; k < m_iZResolution ; k++)
				{
					m_poArray[i][j][k] = oArray.m_poArray[i][j][k];
				}
			}
		}
		return *this;
	}
	void Reset()
	{
		DeleteArray();
		Initialize();
	}
	void SetBox(const AxisAlignedBoundingBox& oBox)
	{
		m_dXMin = oBox.GetXMin();
		m_dYMin = oBox.GetYMin();
		m_dZMin = oBox.GetZMin();
		m_dDeltaX = (oBox.GetXMax() - m_dXMin)/(double)(m_iXResolution);
		m_dDeltaY = (oBox.GetYMax() - m_dYMin)/(double)(m_iYResolution);
		m_dDeltaZ = (oBox.GetZMax() - m_dZMin)/(double)(m_iZResolution);
	}
	void SetBox(const double& dXDim,const double& dYDim,const double& dZDim)
	{
		SetBox(-0.5*dXDim,0.5*dXDim,-0.5*dYDim,0.5*dYDim,-0.5*dZDim,0.5*dZDim);
	}
	void SetBox(const double& dXMin,const double& dXMax,const double& dYMin,const double& dYMax,const double& dZMin,const double& dZMax)
	{
		AxisAlignedBoundingBox oBox;
		oBox.SetXMin(dXMin);
		oBox.SetXMax(dXMax);
		oBox.SetYMin(dYMin);
		oBox.SetYMax(dYMax);
		oBox.SetZMin(dZMin);
		oBox.SetZMax(dZMax);
		SetBox(oBox);
	}
	void SetResolution(const unsigned int& iXResolution,const unsigned int& iYResolution,const unsigned int& iZResolution)
	{
		DeleteArray();
		m_iXResolution = iXResolution;
		m_iYResolution = iYResolution;
		m_iZResolution = iZResolution;
		// allocate array memory
		unsigned int i = 0;
		unsigned int j = 0;
		m_poArray = new Type**[m_iXResolution];
		for(i = 0 ; i < m_iXResolution ; i++)
		{
			m_poArray[i] = new Type*[m_iYResolution];
			for(j = 0 ; j < m_iYResolution ; j++)
			{
				m_poArray[i][j] = new Type[m_iZResolution];
			}
		}
	}
	AxisAlignedBoundingBox GetVoxel(const unsigned int& i,const unsigned int& j,const unsigned int& k) const
	{
		AxisAlignedBoundingBox oVoxel;
		oVoxel.SetXMin(m_dXMin + i*m_dDeltaX);
		oVoxel.SetXMax(m_dXMin + (i + 1)*m_dDeltaX);
		oVoxel.SetYMin(m_dYMin + j*m_dDeltaY);
		oVoxel.SetYMax(m_dYMin + (j + 1)*m_dDeltaY);
		oVoxel.SetZMin(m_dZMin + k*m_dDeltaZ);
		oVoxel.SetZMax(m_dZMin + (k + 1)*m_dDeltaZ);
		return oVoxel;
	}
	void GetContainingVoxelIndices(const Point& oPoint,int& i,int& j,int& k)
	{
		i = (int)floor((oPoint.GetX() - m_dXMin)/m_dDeltaX);
		j = (int)floor((oPoint.GetY() - m_dYMin)/m_dDeltaY);
		k = (int)floor((oPoint.GetZ() - m_dZMin)/m_dDeltaZ);
		if(i >= m_iXResolution)			i = m_iXResolution - 1;
		if(j >= m_iYResolution)			j = m_iYResolution - 1;
		if(k >= m_iZResolution)			k = m_iZResolution - 1;
		if(i < 0)			i = 0;
		if(j < 0)			j = 0;
		if(k < 0)			k = 0;
	}
	Type GetVoxelData(const unsigned int& i,const unsigned int& j,const unsigned int& k) const
	{
		return m_poArray[i][j][k];
	}
	void SetVoxelData(const unsigned int& i,const unsigned int& j,const unsigned int& k,const Type& oData) const
	{
		m_poArray[i][j][k] = oData;
	}
	void WriteVTKGeometry(FILE* fpFile)
	{		
		WritePointGeometry(fpFile);
		WriteVoxelGeometry(fpFile);
	}
	Point GetVoxelCenter(const unsigned int& i,const unsigned int& j,const unsigned int& k) const
	{
		Point oCenter;
		oCenter.SetX(m_dXMin + (i + 0.5)*m_dDeltaX);
		oCenter.SetY(m_dYMin + (j + 0.5)*m_dDeltaY);
		oCenter.SetZ(m_dZMin + (k + 0.5)*m_dDeltaZ);
		return oCenter;
	}
	unsigned int GetXResolution() const
	{
		return m_iXResolution;
	}
	unsigned int GetYResolution() const
	{
		return m_iYResolution;
	}
	unsigned int GetZResolution() const
	{
		return m_iZResolution;
	}
	AxisAlignedBoundingBox GetBox() const
	{
		AxisAlignedBoundingBox oBox;
		oBox.SetXMin(m_dXMin);
		oBox.SetXMax(m_dXMin + m_iXResolution*m_dDeltaX);
		oBox.SetYMin(m_dYMin);
		oBox.SetYMax(m_dYMin + m_iYResolution*m_dDeltaY);
		oBox.SetZMin(m_dZMin);
		oBox.SetZMax(m_dZMin + m_iZResolution*m_dDeltaZ);
		return oBox;
	}
	
private:

protected:
	void Initialize()
	{
		m_poArray = NULL;
		m_iXResolution = 0;
		m_iYResolution = 0;
		m_iZResolution = 0;
		m_dXMin = 0.0;
		m_dYMin = 0.0;
		m_dZMin = 0.0;
		m_dDeltaX = 0.0;
		m_dDeltaY = 0.0;
		m_dDeltaZ = 0.0;
	}
	void DeleteArray()
	{
		unsigned int i = 0;
		unsigned int j = 0;
		if(m_poArray != NULL)
		{
			for(i = 0 ; i < m_iXResolution ; i++)
			{
				if(m_poArray[i] != NULL)
				{
					for(j = 0 ; j < m_iYResolution ; j++)
					{
						if(m_poArray[i][j] != NULL)
						{
							delete [] m_poArray[i][j];
						}
					}
					delete [] m_poArray[i];
				}
			}
			delete [] m_poArray;
		}
		m_poArray = NULL;
	}
	void WritePointGeometry(FILE* fpFile)
	{
		double dX = 0.0;
		double dY = 0.0;
		double dZ = 0.0;
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		fputs("      <Points>\n",fpFile);
		fputs("        <DataArray type=\"Float32\" NumberOfComponents=\" 3\" format=\"ascii\">\n",fpFile);
		for(i = 0 ; i <= m_iXResolution ; i++)
		{
			dX = m_dXMin + i*m_dDeltaX;
			for(j = 0 ; j <= m_iYResolution ; j++)
			{
				dY = m_dYMin + j*m_dDeltaY;
				for(k = 0 ; k <= m_iZResolution ; k++)
				{
					dZ = m_dZMin + k*m_dDeltaZ;
					fprintf(fpFile,"%e\t\t%e\t\t%e\n",dX,dY,dZ);
				}
			}
		}
		fputs("        </DataArray>\n",fpFile);
		fputs("      </Points>\n",fpFile);
	}
	void WriteVoxelGeometry(FILE* fpFile)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		fputs("      <Cells>\n",fpFile);
		{
			fputs("        <DataArray type=\"Float32\" Name=\"connectivity\" format=\"ascii\">\n",fpFile);		
			unsigned int iLayerSize = (m_iYResolution + 1)*(m_iZResolution + 1);
			unsigned int iRowSize = m_iZResolution + 1;
			unsigned int iNode1Index = 0;
			unsigned int iNode2Index = 0;
			unsigned int iNode3Index = 0;
			unsigned int iNode4Index = 0;
			unsigned int iNode5Index = 0;
			unsigned int iNode6Index = 0;
			unsigned int iNode7Index = 0;
			unsigned int iNode8Index = 0;
			for(i = 0 ; i < m_iXResolution ; i++)
			{
				for(j = 0 ; j < m_iYResolution ; j++)
				{
					for(k = 0 ; k < m_iZResolution ; k++)
					{
						iNode1Index = i*iLayerSize + j*iRowSize + k;
						iNode2Index = (i + 1)*iLayerSize + j*iRowSize + k;
						iNode3Index = i*iLayerSize + (j + 1)*iRowSize + k;
						iNode4Index = (i + 1)*iLayerSize + (j + 1)*iRowSize + k;
						iNode5Index = i*iLayerSize + j*iRowSize + k + 1;
						iNode6Index = (i + 1)*iLayerSize + j*iRowSize + k + 1;
						iNode7Index = i*iLayerSize + (j + 1)*iRowSize + k + 1;
						iNode8Index = (i + 1)*iLayerSize + (j + 1)*iRowSize + k + 1;
						fprintf(fpFile,"%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\t\t%d\n",iNode1Index,iNode2Index,iNode3Index,iNode4Index,iNode5Index,iNode6Index,iNode7Index,iNode8Index);
					}
				}
			}
			fputs("        </DataArray>\n",fpFile);
		}
		{
			fputs("        <DataArray type=\"Float32\" Name=\"offsets\" format=\"ascii\">\n",fpFile);
			unsigned int iCount = 0;
			unsigned int iOffset = 0;
			for(i = 0 ; i < m_iXResolution ; i++)
			{
				for(j = 0 ; j < m_iYResolution ; j++)
				{
					for(k = 0 ; k < m_iZResolution ; k++)
					{
						iOffset = 8*(iCount + 1);
						fprintf(fpFile,"%d\n",iOffset);
						iCount = iCount + 1;
					}
				}
			}
			fputs("        </DataArray>\n",fpFile);
		}
		{
			fputs("        <DataArray type=\"Float32\" Name=\"types\" format=\"ascii\">\n",fpFile);
			unsigned int iCellType = 11;
			for(i = 0 ; i < m_iXResolution ; i++)
			{
				for(j = 0 ; j < m_iYResolution ; j++)
				{
					for(k = 0 ; k < m_iZResolution ; k++)
					{
						fprintf(fpFile,"%d\n",iCellType);
					}
				}
			}
			fputs("        </DataArray>\n",fpFile);
		}
		fputs("      </Cells>\n",fpFile);
	}
	Type*** m_poArray;
	unsigned int m_iXResolution;
	unsigned int m_iYResolution;
	unsigned int m_iZResolution;
	double m_dXMin;
	double m_dYMin;
	double m_dZMin;
	double m_dDeltaX;
	double m_dDeltaY;
	double m_dDeltaZ;
};

#endif

