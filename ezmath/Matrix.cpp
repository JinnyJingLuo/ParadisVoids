// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "Matrix.h"
#include "math.h"
#include "float.h"
#include "Randomizer.h"

namespace EZ
{
	Matrix::Matrix()
	{
		m_iRowsCount = 0;
		m_iColumnsCount = 0;
		m_ppdData = NULL;
		Reset();
	}
	Matrix::Matrix(const Matrix& oMatrix)
	{
		m_iRowsCount = 0;
		m_iColumnsCount = 0;
		m_ppdData = NULL;
		*this = oMatrix;
	}
	Matrix::Matrix(const unsigned int& iRowsCount,const unsigned int& iColumnsCount)
	{
		m_iRowsCount = 0;
		m_iColumnsCount = 0;
		m_ppdData = NULL;
		SetSize(iRowsCount,iColumnsCount);
	}
	void Matrix::Reset()
	{
		unsigned int i = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			if(m_ppdData[i] != NULL)
			{
				delete[] m_ppdData[i];
				m_ppdData[i] = NULL;
			}
		}
		if(m_ppdData != NULL)
		{
			delete[] m_ppdData;
			m_ppdData = NULL;
		}
		m_iRowsCount = 0;
		m_iColumnsCount = 0;
	}
	Matrix::~Matrix()
	{
		Reset();
	}
	void Matrix::Set(const unsigned int& iRow,const unsigned int& iColumn,const double& dVal)
	{
		m_ppdData[iRow - 1][iColumn - 1] = dVal;
	}
	double Matrix::Get(const unsigned int& iRow,const unsigned int& iColumn) const
	{
		return m_ppdData[iRow - 1][iColumn - 1];
	}
	void Matrix::SetSize(const unsigned int& iRowsCount,const unsigned int& iColumnsCount)
	{
		Reset();
		m_ppdData = new double*[iRowsCount];
		double* pdTemp;
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0; i < iRowsCount ; i++)
		{
			pdTemp = new double[iColumnsCount];
			for(j = 0; j < iColumnsCount ; j++)
			{
				pdTemp[j] = 0.0;
			}
			m_ppdData[i] = pdTemp;
		}
		m_iRowsCount = iRowsCount;
		m_iColumnsCount = iColumnsCount;
	}
	Matrix& Matrix::operator=(const Matrix& oMatrix)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		if((m_iRowsCount != oMatrix.m_iRowsCount) || (m_iColumnsCount != oMatrix.m_iColumnsCount))
		{
			SetSize(oMatrix.m_iRowsCount,oMatrix.m_iColumnsCount);
		}
		for(i = 0; i < m_iRowsCount ; i++)
		{
			for(j = 0; j < m_iColumnsCount ; j++)
			{
				m_ppdData[i][j] = oMatrix.m_ppdData[i][j];
			}
		}
		return *this;
	}
	Matrix Matrix::operator*(const Matrix& oMatrix) const
	{
		Matrix oResult;
		unsigned int iRows = GetRowsCount();
		unsigned int iColumns = oMatrix.GetColumnsCount();
		unsigned int iCommon = GetColumnsCount();
		oResult.SetSize(iRows,iColumns);
		if(iCommon != oMatrix.GetRowsCount())
		{
			printf("error: matrix-matrix multiplication size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		double dSum = 0.0;
		for(i = 0 ; i < iRows ; i++)
		{
			for(j = 0 ; j < iColumns ; j++)
			{
				dSum = 0.0;
				for(k = 0; k < iCommon ; k++)
				{
					dSum = dSum + m_ppdData[i][k]*oMatrix.m_ppdData[k][j];
				}
				oResult.m_ppdData[i][j] = dSum;
			}
		}
		return oResult;
	}
	Vector Matrix::operator*(const Vector& o3DVector) const
	{
		Vector oResult(0.0,0.0,0.0);
		if(m_iColumnsCount != 3 || m_iColumnsCount != 3)
		{
			printf("error: matrix-3d vector multiplication size mismatch\n");
			return oResult;
		}
		double dX = o3DVector.GetX();
		double dY = o3DVector.GetY();
		double dZ = o3DVector.GetZ();
		oResult.SetX(Get(1,1)*dX + Get(1,2)*dY + Get(1,3)*dZ);
		oResult.SetY(Get(2,1)*dX + Get(2,2)*dY + Get(2,3)*dZ);
		oResult.SetZ(Get(3,1)*dX + Get(3,2)*dY + Get(3,3)*dZ);
		return oResult;
	}
	Matrix Matrix::operator*(const double& dFactor) const
	{
		unsigned int i = 0; 
		unsigned int j = 0;
		Matrix oResult(m_iRowsCount,m_iColumnsCount);
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = 0 ; j < m_iColumnsCount ; j++)
			{
				oResult.m_ppdData[i][j] = m_ppdData[i][j]*dFactor;
			}
		}
		return oResult;
	}
	unsigned int Matrix::GetRowsCount() const
	{
		return m_iRowsCount;
	}
	unsigned int Matrix::GetColumnsCount() const
	{
		return m_iColumnsCount;
	}
	Matrix Matrix::Solve(const Matrix& oRHS) const
	{
		////////////////////////////////////////////
		// the solution using Gauss elimination with partial
		// pivoting. very robust but very slow
		////////////////////////////////////////////
		Matrix oTemp;
		unsigned int iRHSCount = oRHS.GetColumnsCount();
		Matrix oResult(m_iRowsCount,iRHSCount);
		if(m_iRowsCount != oRHS.GetRowsCount())
		{
			printf("error: matrix solution size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		vector<double> vdScales;
		vdScales.resize(m_iRowsCount);
		double dTemp = 0.0;
		// get the maximum for each row and place it in the scales vector
		for(i = 0; i < m_iRowsCount ; i++)
		{
			dTemp = DBL_MIN;
			for(j = 0; j < m_iColumnsCount ; j++)
			{
				if(fabs(m_ppdData[i][j]) > dTemp)
				{
					dTemp = fabs(m_ppdData[i][j]);
				}
			}
			vdScales[i] = dTemp;
		}
		// set the temp matrix with the row wise scaled version
		// of the original system matrix
		oTemp.SetSize(m_iRowsCount,m_iColumnsCount);
		Matrix oTempRHS;
		oTempRHS.SetSize(oRHS.GetRowsCount(),oRHS.GetColumnsCount());
		vector<unsigned int> viIndices;		// used for virtual row swapping
		viIndices.resize(m_iRowsCount);
		for(i = 0; i < m_iRowsCount ; i++)
		{
			for(j = 0; j < m_iColumnsCount ; j++)
			{
				oTemp.Set(i + 1,j + 1,(m_ppdData[i][j])/vdScales[i]);
			}
			for(j = 0; j < iRHSCount ; j++)
			{
				oTempRHS.Set(i + 1,j + 1,oRHS.Get(i + 1,j + 1)/vdScales[i]);
			}
			viIndices[i] = i;
		}
		unsigned int iPivotIndex = 0;
		unsigned int k = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			dTemp = fabs(oTemp.Get(viIndices[i] + 1,i + 1));
			iPivotIndex = i;
			for(j = i + 1; j < m_iRowsCount ; j++)
			{
				if(fabs(oTemp.Get(viIndices[j] + 1,i + 1)) > dTemp)
				{
					dTemp = fabs(oTemp.Get(viIndices[j] + 1,i + 1));
					iPivotIndex = j;
				}
			}
			if(dTemp < 1E-50)
			{
				return oResult;
			}
			j = viIndices[i];
			viIndices[i] = viIndices[iPivotIndex];
			viIndices[iPivotIndex] = j;

			for(j = i + 1; j < m_iRowsCount; j++)
			{
				dTemp = -oTemp.Get(viIndices[j] + 1,i + 1)/oTemp.Get(viIndices[i] + 1,i + 1);
				for(k = i + 1; k < m_iColumnsCount; k++)
				{
					oTemp.Set(viIndices[j] + 1,k + 1,oTemp.Get(viIndices[j] + 1,k + 1) + dTemp*oTemp.Get(viIndices[i] + 1,k + 1));
				}
				for(k = 0; k < iRHSCount; k++)
				{
					oTempRHS.Set(viIndices[j] + 1,k + 1,oTempRHS.Get(viIndices[j] + 1,k + 1) + dTemp*oTempRHS.Get(viIndices[i] + 1,k + 1));
				}
			}
		}
		// back substitution to get the final solution
		for(int i = m_iRowsCount - 1 ; i >= 0; i--)
		{
			for(k = 0; k < iRHSCount ; k++)
			{	
				dTemp = 0.0;
				for(j = i + 1; j < m_iColumnsCount; j++)
				{
					dTemp = dTemp + oTemp.Get(viIndices[i] + 1,j + 1)*oTempRHS.Get(viIndices[j] + 1,k + 1);
				}
				oTempRHS.Set(viIndices[i] + 1,k + 1,(oTempRHS.Get(viIndices[i] + 1,k + 1) - dTemp)/oTemp.Get(viIndices[i] + 1,i + 1));
			}
		}

		for(i = 0; i < m_iColumnsCount; i++)
		{
			for(j = 0; j < iRHSCount ; j++)
			{
				oResult.Set(i + 1,j + 1,oTempRHS.Get(viIndices[i] + 1,j + 1));
			}
		}
		return oResult;
	}
	void Matrix::FillMatrixFromVector(vector<double>* poVector)
	{
		unsigned int i = 0;
		unsigned int iSize = poVector->size();
		SetSize(iSize,1);
		for(i = 0; i < iSize ; i++)
		{
			Set(i + 1,1,poVector->at(i));
		}
	}
	vector<double> Matrix::Vectorize()
	{
		vector<double> vdResult;
		vdResult.resize(m_iRowsCount);
		unsigned int i = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			vdResult[i] = Get(i+1,1);
		}
		return vdResult;
	}
	void Matrix::SetColumn(const unsigned int& iIndex,Matrix oColumn)
	{
		if(m_iRowsCount != oColumn.GetRowsCount() || iIndex > m_iColumnsCount)
		{
			printf("error: matrix column setting size/index mismatch\n");
			return;
		}

		unsigned int i = 0;
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			Set(i,iIndex,oColumn.Get(i,1));
		}
	}
	void Matrix::SetRow(const unsigned int& iIndex,Matrix oRow)
	{
		if(m_iColumnsCount != oRow.GetColumnsCount() || iIndex > m_iRowsCount)
		{
			printf("error: matrix row setting size/index mismatch\n");
			return;
		}

		unsigned int i = 0;
		for(i = 1; i <= m_iColumnsCount ; i++)
		{
			Set(iIndex,i,oRow.Get(1,i));
		}
	}
	void Matrix::PlaceSubMatrix(const unsigned int& iEntryRowIndex,const unsigned int& iEntryColumnIndex,Matrix* poMatrix)
	{
		unsigned int iHeight = poMatrix->GetRowsCount();
		unsigned int iWidth = poMatrix->GetColumnsCount();

		if(iEntryRowIndex + iHeight > m_iRowsCount || iEntryColumnIndex + iWidth > m_iColumnsCount)
		{
			printf("error: matrix-submatrix placement size/index mismatch\n");
			return;
		}

		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0; i < iHeight ; i++)
		{
			for(j = 0; j < iWidth ; j++)
			{
				Set(iEntryRowIndex + i + 1,iEntryColumnIndex + j + 1,poMatrix->Get(i + 1,j + 1));
			}
		}
	}
	Matrix Matrix::GetRow(const unsigned int& iIndex) const
	{
		Matrix oResult(1,m_iColumnsCount);
		if(iIndex > m_iRowsCount)
		{
			printf("error: matrix row getting index out of bounds\n");
			return oResult;
		}
		unsigned int i = 0;
		for(i = 1 ; i <= m_iColumnsCount ; i++)
		{
			oResult.Set(1,i,Get(iIndex,i));
		}
		return oResult;
	}
	Matrix Matrix::GetColumn(const unsigned int& iIndex) const
	{
		Matrix oResult(m_iRowsCount,1);
		if(iIndex > m_iColumnsCount)
		{
			printf("error: matrix column getting index out of bounds\n");
			return oResult;
		}
		unsigned int i = 0;
		for(i = 1 ; i <= m_iRowsCount ; i++)
		{
			oResult.Set(i,1,Get(i,iIndex));
		}
		return oResult;
	}
	Matrix Matrix::GetInverse() const
	{
		unsigned int i = 0;
		Matrix oRHS;
		oRHS.SetSize(m_iRowsCount,m_iColumnsCount);
		for(i = 0; i < m_iRowsCount ; i++)
		{
			oRHS.Set(i + 1,i + 1,1.0);
		}
		Matrix oResult = Solve(oRHS);
		return oResult;
	}
	double Matrix::GetAbsoluteMaximum() const
	{
		unsigned int i = 0;
		unsigned int j = 0;
		double dMax = 0.0;
		double dTemp = 0.0;
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			for(j = 1; j <= m_iColumnsCount ; j++)
			{
				dTemp = fabs(Get(i,j));
				if(dTemp > dMax)
				{
					dMax = dTemp;
				}
			}
		}
		return dMax;
	}
	void Matrix::Filter(const double& dTolerance)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		double dMax = GetAbsoluteMaximum();
		for(i = 0; i < m_iRowsCount ; i++)
		{
			for(j = 0; j < m_iColumnsCount ; j++)
			{
				if(fabs(Get(i + 1,j + 1))/dMax < dTolerance)
				{
					Set(i + 1,j + 1,0.0);
				}
			}
		}
	}
	Matrix Matrix::operator+(const Matrix& oMatrix) const
	{
		Matrix oResult;
		oResult.SetSize(m_iRowsCount,m_iColumnsCount);
		if(m_iRowsCount != oMatrix.GetRowsCount() || m_iColumnsCount != oMatrix.GetColumnsCount())
		{
			printf("error: matrix-matrix addition size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = 0 ; j < m_iColumnsCount ; j++)
			{
				oResult.m_ppdData[i][j] = m_ppdData[i][j] + oMatrix.m_ppdData[i][j];
			}
		}
		return oResult;
	}
	Matrix Matrix::operator-(const Matrix& oMatrix) const
	{
		Matrix oResult;
		oResult.SetSize(m_iRowsCount,m_iColumnsCount);
		if(m_iRowsCount != oMatrix.GetRowsCount() || m_iColumnsCount != oMatrix.GetColumnsCount())
		{
			printf("error: matrix-matrix subtraction size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = 0 ; j < m_iColumnsCount ; j++)
			{
				oResult.m_ppdData[i][j] = m_ppdData[i][j] - oMatrix.m_ppdData[i][j];
			}
		}
		return oResult;
	}
	void Matrix::SubtractFromIdentity()
	{
		if(!IsSquare())
		{
			return;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			for(j = 1; j <= m_iColumnsCount ; j++)
			{
				Set(i,j,-Get(i,j));
			}
		}
		AddToIdentity();
	}
	void Matrix::AddToIdentity()
	{
		unsigned int i = 0;
		if(m_iRowsCount != m_iColumnsCount)
		{
			return;
		}
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			Set(i,i,1.0 + Get(i,i));
		}
	}
	void Matrix::AddToEntry(const unsigned int& iRow,const unsigned int& iColumn,const double& dVal)
	{
		Set(iRow,iColumn,Get(iRow,iColumn) + dVal);
	}
	Matrix Matrix::GetTranspose() const
	{
		Matrix oTranspose(m_iColumnsCount,m_iRowsCount);
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			for(j = 1; j <= m_iColumnsCount ; j++)
			{
				oTranspose.Set(j,i,Get(i,j));
			}
		}
		return oTranspose;
	}
	void Matrix::ResetToZeros()
	{
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			for(j = 1; j <= m_iColumnsCount ; j++)
			{
				Set(i,j,0.0);
			}
		}
	}
	double Matrix::GetNorm() const
	{
		unsigned int i = 0;
		unsigned int j = 0;
		double dSum = 0.0;
		double dTemp = 0.0;
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			for(j = 1; j <= m_iColumnsCount ; j++)
			{
				dTemp = Get(i,j);
				dSum = dSum + dTemp*dTemp;
			}
		}
		return sqrt(dSum);
	}
	unsigned int Matrix::GetSSPDSize(const double& dTolerance) const
	{
		unsigned int i = 0;
		unsigned int j = 0;
		double dThreshold = dTolerance*GetAbsoluteMaximum();
		unsigned int iCount = 0;
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			for(j = i; j <= m_iColumnsCount ; j++)
			{
				if(fabs(Get(i,j)) > dThreshold)
				{
					iCount = iCount + 1;
				}
			}
		}
		return iCount;
	}
	void Matrix::StoreSSPD(double*& pdValues,unsigned int*& piRowsStartingLocations,unsigned int*& piColumnsIndices,const double& dTolerance) const
	{
		unsigned int iSSPDSize = GetSSPDSize(dTolerance);
		unsigned int i = 0;
		unsigned int j = 0;
		double dTemp = 0.0;
		unsigned int iCount = 0;
		pdValues = new double[iSSPDSize];
		piColumnsIndices = new unsigned int[iSSPDSize];
		piRowsStartingLocations = new unsigned int[m_iRowsCount + 1];
		for(i = 1; i <= m_iRowsCount ; i++)
		{
			piRowsStartingLocations[i - 1] = iCount;
			for(j = i; j <= m_iColumnsCount ; j++)
			{
				dTemp = Get(i,j);
				if(fabs(dTemp) > dTolerance)
				{
					pdValues[iCount] = dTemp;
					piColumnsIndices[iCount] = j;
					iCount = iCount + 1;
				}
			}
		}
		piRowsStartingLocations[m_iRowsCount] = iSSPDSize;
		printf("actual nnz size is : %d\n",iSSPDSize);
	}
	Matrix Matrix::Invert3x3Matrix(const Matrix& oMatrix,double& dDeterminant)
	{
		Matrix oInverse(3,3);
		if(oMatrix.GetRowsCount() != 3 || oMatrix.GetColumnsCount() != 3)
		{
			printf("error: 3x3 matrix inversion size mismatch\n");
			return oInverse;
		}
		double dDet11 = oMatrix.Get(2,2)*oMatrix.Get(3,3) - oMatrix.Get(2,3)*oMatrix.Get(3,2);
		double dDet12 = oMatrix.Get(2,3)*oMatrix.Get(3,1) - oMatrix.Get(2,1)*oMatrix.Get(3,3);
		double dDet13 = oMatrix.Get(2,1)*oMatrix.Get(3,2) - oMatrix.Get(3,1)*oMatrix.Get(2,2);

		double dDet21 = oMatrix.Get(3,2)*oMatrix.Get(1,3) - oMatrix.Get(1,2)*oMatrix.Get(3,3);
		double dDet22 = oMatrix.Get(3,3)*oMatrix.Get(1,1) - oMatrix.Get(3,1)*oMatrix.Get(1,3);
		double dDet23 = oMatrix.Get(3,1)*oMatrix.Get(1,2) - oMatrix.Get(3,2)*oMatrix.Get(1,1);

		double dDet31 = oMatrix.Get(1,2)*oMatrix.Get(2,3) - oMatrix.Get(1,3)*oMatrix.Get(2,2);
		double dDet32 = oMatrix.Get(1,3)*oMatrix.Get(2,1) - oMatrix.Get(2,3)*oMatrix.Get(1,1);
		double dDet33 = oMatrix.Get(1,1)*oMatrix.Get(2,2) - oMatrix.Get(1,2)*oMatrix.Get(2,1);

		dDeterminant = oMatrix.Get(1,1)*dDet11 + oMatrix.Get(1,2)*dDet12 + oMatrix.Get(1,3)*dDet13;

		oInverse.Set(1,1,dDet11);
		oInverse.Set(1,2,dDet12);
		oInverse.Set(1,3,dDet13);

		oInverse.Set(2,1,dDet21);
		oInverse.Set(2,2,dDet22);
		oInverse.Set(2,3,dDet23);
		
		oInverse.Set(3,1,dDet31);
		oInverse.Set(3,2,dDet32);
		oInverse.Set(3,3,dDet33);
		oInverse = oInverse.GetTranspose();
		oInverse = oInverse*(1.0/dDeterminant);

		return oInverse;
	}
	Matrix Matrix::Invert2x2Matrix(const Matrix& oMatrix,double& dDeterminant)
	{
		Matrix oInverse(2,2);
		if(oMatrix.GetRowsCount() != 2 || oMatrix.GetColumnsCount() != 2)
		{
			printf("error: 2x2 matrix inversion size mismatch\n");
			return oInverse;
		}
		
		dDeterminant = oMatrix.Get(1,1)*oMatrix.Get(2,2) - oMatrix.Get(1,2)*oMatrix.Get(2,1);

		oInverse.Set(1,1,oMatrix.Get(2,2));
		oInverse.Set(1,2,-oMatrix.Get(1,2));

		oInverse.Set(2,1,-oMatrix.Get(2,1));
		oInverse.Set(2,2,oMatrix.Get(1,1));

		oInverse = oInverse*(1.0/dDeterminant);
		return oInverse;
	}
	Matrix Matrix::Solve3x3System(const Matrix& oSystem,const Matrix& oRHS)
	{
		if(oRHS.GetRowsCount() != 3)
		{
			printf("error: 3x3 matrix solution size mismatch\n");
			return Matrix(oRHS.GetRowsCount(),oRHS.GetColumnsCount());
		}
		double dDummyDeterminant = 0.0;
		Matrix oInverse = Invert3x3Matrix(oSystem,dDummyDeterminant);
		Matrix oResult = oInverse*oRHS;
		return oResult;
	}
	Matrix Matrix::Solve2x2System(const Matrix& oSystem,const Matrix& oRHS)
	{
		if(oRHS.GetRowsCount() != 2)
		{
			printf("error: 2x2 matrix solution size mismatch\n");
			return Matrix(oRHS.GetRowsCount(),oRHS.GetColumnsCount());
		}
		double dDummyDeterminant = 0.0;
		Matrix oInverse = Invert2x2Matrix(oSystem,dDummyDeterminant);
		Matrix oResult = oInverse*oRHS;
		return oResult;
	}
	void Matrix::NormalizeColumns()
	{
		unsigned int i = 0;
		unsigned int j = 0;
		double dNorm = 0.0;
		for(j = 0; j < m_iColumnsCount ; j++)
		{
			dNorm = 0.0;
			for(i = 0; i < m_iRowsCount ; i++)
			{
				dNorm = dNorm + m_ppdData[i][j]*m_ppdData[i][j];
			}
			dNorm = sqrt(dNorm);
			if(dNorm < 1E-10)
			{
				continue;
			}
			for(i = 0; i < m_iRowsCount ; i++)
			{
				m_ppdData[i][j] = m_ppdData[i][j]/dNorm;
			}
		}
	}
	void Matrix::Randomize(const double& dMin,const double& dMax)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			for(j = 0; j < m_iColumnsCount ; j++)
			{
				m_ppdData[i][j] = Randomizer::Random(dMin,dMax);
			}
		}
	}
	void Matrix::RandomizeNormal(const double& dMean,const double& dStandardDeviation)
	{
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			for(j = 0; j < m_iColumnsCount ; j++)
			{
				m_ppdData[i][j] = Randomizer::RandomNormal(dMean,dStandardDeviation);
			}
		}
	}
	bool Matrix::IsSquare() const
	{
		if(m_iRowsCount == m_iColumnsCount)
		{
			return true;
		}
		return false;
	}
	Matrix Matrix::GenerateGramSchmidtOrthogonalMatrix() const
	{
		Matrix oOrthogonalMatrix(m_iRowsCount,m_iColumnsCount);
		Matrix oTempColumn(m_iRowsCount,1);
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		double dProjection = 0.0;
		for(j = 1 ; j <= m_iColumnsCount ; j++)
		{
			for(i = 1 ; i <= m_iRowsCount ; i++)
			{
				oTempColumn.Set(i,1,Get(i,j));
			}
			for(k = 1 ; k <= j - 1 ; k++)
			{
				// calculate the projection length
				dProjection = 0.0;
				for(i = 1; i <= m_iRowsCount ; i++)
				{
					dProjection = dProjection + oTempColumn.Get(i,1)*oOrthogonalMatrix.Get(i,k);
				}
				// update the new vector
				for(i = 1; i <= m_iRowsCount ; i++)
				{
					oTempColumn.AddToEntry(i,1,-dProjection*oOrthogonalMatrix.Get(i,k));
				}
			}
			oTempColumn.NormalizeColumns();
			for(i = 1; i <= m_iRowsCount ; i++)
			{
				oOrthogonalMatrix.Set(i,j,oTempColumn.Get(i,1));
			}
		}
		return oOrthogonalMatrix;
	}
	Matrix Matrix::GenerateRandomOrthogonalMatrix(const unsigned int& iSize)
	{
		Matrix oRandomMatrix(iSize,iSize);
		Matrix oOrthogonal;
		Matrix oTestMatrix;
		double dError = 0.0;
		double dTolerance = 1.0E-6;
		while(true)
		{
			oRandomMatrix.Randomize(-1000,1000);
			oOrthogonal = oRandomMatrix.GenerateGramSchmidtOrthogonalMatrix();
			oTestMatrix = oOrthogonal.GetTranspose()*oOrthogonal;
			oTestMatrix.Filter();
			oTestMatrix.SubtractFromIdentity();
			dError = oTestMatrix.GetNorm();
			if(dError < dTolerance)
			{
				break;
			}
		}
		return oOrthogonal;
	}
	Matrix Matrix::GenerateRandomSymmetricPositiveDefiniteMatrix(const unsigned int& iSize)
	{
		Matrix oResult;
		Matrix oTemp(iSize,iSize);
		unsigned int i = 0;
		unsigned int j = 0;
		Matrix oRotation = GenerateRandomOrthogonalMatrix(iSize);
		// an optimized version for the product of the diagonal
		// matrix with eigenvalues = 1,2,3,... and the transpose
		// of the rotation matrix
		for(i = 1; i <= iSize ; i++)
		{
			for(j = 1; j <= iSize ; j++)
			{
				oTemp.Set(i,j,i*oRotation.Get(j,i));
			}
		}
		oResult = oRotation*oTemp;
		return oResult;
	}
	Matrix Matrix::MultiplyDiagonalMatrixAndVector(const Matrix& oDiagonal,const Matrix& oVector)
	{
		unsigned int iRowsCount = oVector.GetRowsCount();
		Matrix oResult(iRowsCount,1);
		// the diagonal matrix should ALWAYS be given as a row vector
		if(oDiagonal.GetColumnsCount() != iRowsCount)
		{
			return oResult;
		}
		
		unsigned int i = 0;
		for(i = 1; i <= iRowsCount ; i++)
		{
			oResult.Set(i,1,oDiagonal.Get(1,i)*oVector.Get(i,1));
		}
		return oResult;
	}
	double Matrix::SumAllEntries()
	{
		unsigned int i = 0;
		unsigned int j = 0;
		double dSum = 0.0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = 0 ; j < m_iColumnsCount ; j++)
			{
				dSum = dSum + m_ppdData[i][j];
			}
		}
		return dSum;
	}
	void Matrix::ZeroValues()
	{
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = 0 ; j < m_iColumnsCount ; j++)
			{
				m_ppdData[i][j] = 0.0;
			}
		}
	}
	void Matrix::QRFactorize(Matrix& oQ,Matrix& oR) const
	{
		QRFactorize(oQ,oR,m_iRowsCount,m_iColumnsCount);
	}
	void Matrix::QRFactorize(Matrix& oQ,Matrix& oR,unsigned int iMaxRow,unsigned int iMaxColumn) const
	{
		// this function performs QR factorization for the submatrix between (1,1) and (max row,max xcol)
		if(iMaxRow > m_iRowsCount)
		{
			iMaxRow = m_iRowsCount;
		}
		if(iMaxColumn > m_iColumnsCount)
		{
			iMaxColumn = m_iColumnsCount;
		}
		oQ.SetSize(iMaxRow,iMaxColumn);
		oR.SetSize(iMaxColumn,iMaxColumn);
		Matrix oTempColumn(iMaxRow,1);
		unsigned int i = 0;
		unsigned int j = 0;
		unsigned int k = 0;
		double dProjection = 0.0;
		for(j = 0 ; j < iMaxColumn ; j++)
		{
			for(i = 0 ; i < iMaxRow ; i++)
			{
				oTempColumn.m_ppdData[i][0] = m_ppdData[i][j];
			}
			for(k = 0 ; k < j ; k++)
			{
				// calculate the projection length
				dProjection = 0.0;
				for(i = 0 ; i < iMaxRow ; i++)
				{
					dProjection = dProjection + oTempColumn.m_ppdData[i][0]*oQ.m_ppdData[i][k];
				}
				// set the R value
				oR.m_ppdData[k][j] = dProjection;
				// update the new vector
				for(i = 0 ; i < iMaxRow ; i++)
				{
					oTempColumn.AddToEntry(i + 1,1,-dProjection*oQ.m_ppdData[i][k]);
				}
			}
			// set the main entry in the R matrix
			oR.m_ppdData[j][j] = oTempColumn.GetNorm();
			oTempColumn.NormalizeColumns();
			for(i = 0 ; i < iMaxRow ; i++)
			{
				oQ.m_ppdData[i][j] = oTempColumn.m_ppdData[i][0];
			}
		}
	}
	Matrix Matrix::UpperTriangularSolve(const Matrix& oRHS) const
	{
		Matrix oX(oRHS.m_iRowsCount,oRHS.m_iColumnsCount);
		if(oRHS.m_iRowsCount != m_iRowsCount)
		{
			return oX;
		}
		int i = 0;
		int j = 0;
		int k = 0;
		for(k = 0 ; k < oRHS.m_iColumnsCount ; k++)
		{
			for(i = m_iRowsCount - 1 ; i >= 0 ; i--)
			{
				oX.m_ppdData[i][k] = oRHS.m_ppdData[i][k];
				for(j = i + 1 ; j < m_iColumnsCount ; j++)
				{
					oX.m_ppdData[i][k] = oX.m_ppdData[i][k] - m_ppdData[i][j]*oX.m_ppdData[j][k];
				}
				oX.m_ppdData[i][k] = oX.m_ppdData[i][k]/m_ppdData[i][i];
			}
		}
		
		return oX;
	}
	double Matrix::MultiplyColumnByRow(const unsigned int& iColumnIndex,const Matrix& oRow) const
	{
		double dResult = 0.0;
		unsigned int i = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			dResult = dResult + m_ppdData[i][iColumnIndex - 1]*oRow.m_ppdData[0][i];
		}
		return dResult;
	}
	double Matrix::MultiplyColumnByColumn(const unsigned int& iColumnIndex,const Matrix& oColumn) const
	{
		double dResult = 0.0;
		unsigned int i = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			dResult = dResult + m_ppdData[i][iColumnIndex - 1]*oColumn.m_ppdData[i][0];
		}
		return dResult;
	}
	double Matrix::MultiplyRowByRow(const unsigned int& iRowIndex,const Matrix& oRow) const
	{
		double dResult = 0.0;
		unsigned int i = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			dResult = dResult + m_ppdData[iRowIndex - 1][i]*oRow.m_ppdData[0][i];
		}
		return dResult;
	}
	double Matrix::MultiplyRowByColumn(const unsigned int& iRowIndex,const Matrix& oColumn) const
	{
		double dResult = 0.0;
		unsigned int i = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			dResult = dResult + m_ppdData[iRowIndex - 1][i]*oColumn.m_ppdData[i][0];
		}
		return dResult;
	}
	Matrix Matrix::LimitedColumnProduct(const unsigned int& iMaxColumn,const Matrix& oColumn)
	{
		// this function multiplies the matrix columns up to a certain column by the given column
		Matrix oResult;
		unsigned int i = 0;
		unsigned int j = 0;
		oResult.SetSize(m_iRowsCount,1);
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			oResult.m_ppdData[i][0] = 0.0;
			for(j = 0 ; j < iMaxColumn ; j++)
			{
				oResult.m_ppdData[i][0] = oResult.m_ppdData[i][0] + m_ppdData[i][j]*oColumn.m_ppdData[j][0];
			}
		}
		return oResult;
	}
	void Matrix::Print() const
	{
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = 0 ; j < m_iColumnsCount ; j++)
			{
				printf("%e\t",m_ppdData[i][j]);
			}
			printf("\n");
		}
	}
}



