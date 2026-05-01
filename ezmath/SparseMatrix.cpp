// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#include "SparseMatrix.h"
#include "float.h"
#include "map"
#include "math.h"

using namespace std;


namespace EZ
{
	SparseMatrix::SparseMatrix()
	{
		Initialize();
	}
	SparseMatrix::SparseMatrix(const unsigned int& iRowsCount)
	{
		Initialize();
		SetRowsCount(iRowsCount);
	}
	SparseMatrix::SparseMatrix(const SparseMatrix& oMatrix)
	{
		Initialize();
		*this = oMatrix;
	}
	SparseMatrix::~SparseMatrix()
	{
		Reset();
	}
	SparseMatrix& SparseMatrix::operator=(const SparseMatrix& oMatrix)
	{
		SetRowsCount(oMatrix.m_iRowsCount);
		unsigned int i = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			*m_ppoData[i] = *oMatrix.m_ppoData[i];
		}
		return *this;
	}
	void SparseMatrix::Reset()
	{
		unsigned int i = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			if(m_ppoData[i] != NULL)
			{
				m_ppoData[i]->clear();
				delete m_ppoData[i];
			}
		}
		if(m_ppoData != NULL)
		{
			delete[] m_ppoData;
		}
		m_ppoData = NULL;
		m_iRowsCount = 0;
		m_oPreconditioner.Reset();
		m_oGMRESQMatrix.Reset();
		m_oGMRESHessenbergMatrix.Reset();
		m_dGMRESBeta = 0.0;
	}
	void SparseMatrix::SetRowsCount(const unsigned int& iCount)
	{
		Reset();
		m_iRowsCount = iCount;
		m_ppoData = new map<unsigned int,double>*[m_iRowsCount];
		unsigned int i = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			m_ppoData[i] = new map<unsigned ,double>;
		}
	}
	unsigned int SparseMatrix::GetRowsCount() const
	{
		return m_iRowsCount;
	}
	bool SparseMatrix::Get(const unsigned int& iRowIndex,const unsigned int& iColumnIndex,double& dValue) const
	{
		map<unsigned int,double>::iterator poRowIterator = m_ppoData[iRowIndex - 1]->find(iColumnIndex);
		if(poRowIterator != m_ppoData[iRowIndex - 1]->end())
		{
			dValue = poRowIterator->second;
			return true;
		}
		return false;
	}
	void SparseMatrix::Set(const unsigned int& iRowIndex,const unsigned int& iColumnIndex,const double& dValue)
	{
		m_ppoData[iRowIndex - 1]->operator [](iColumnIndex) = dValue;
	}
	void SparseMatrix::DropEntry(const unsigned int& iRowIndex,const unsigned int& iColumnIndex)
	{
		map<unsigned int,double>::iterator poRowIterator = m_ppoData[iRowIndex - 1]->find(iColumnIndex);
		if(poRowIterator != m_ppoData[iRowIndex - 1]->end())
		{
			m_ppoData[iRowIndex - 1]->erase(poRowIterator);
		}
	}
	void SparseMatrix::AddToEntry(const unsigned int& iRowIndex,const unsigned int& iColumnIndex,const double& dValue)
	{
		double dEntry = dValue;
		if(Get(iRowIndex,iColumnIndex,dEntry))
		{
			dEntry = dEntry + dValue;
		}
		Set(iRowIndex,iColumnIndex,dEntry);
	}
	Matrix SparseMatrix::SolveConjugateGradient(const Matrix& oRHS) const
	{
		//////////////////////////////////////////////////
		// solution using the preconditioned conjugate gradient method
		//////////////////////////////////////////////////
		unsigned int iRHSCount = oRHS.GetColumnsCount();
 		Matrix oResult(oRHS.GetRowsCount(),iRHSCount);
 		if(oRHS.GetRowsCount() != m_iRowsCount)
 		{
 			printf("error: sparse matrix conjugate gradient solution size mismatch\n");
 			return oResult;
 		}
 		unsigned int i = 0;
 		Matrix oRk;
 		Matrix oRkp1;
 		Matrix oPk;
 		Matrix oPkp1;
 		double dAlpha = 0.0;
 		double dBeta = 0.0;
 		double dToleranceFactor = 1E-6;
 		double dTolerance = 0.0;
 		double dErrorNorm = 0.0;
 		Matrix oTemp;
 		double dNum = 0.0;
 		double dDen = 0.0;
 		Matrix oX;
 		Matrix oZ;
 		double dTemp = 0.0;
 		// solve the system
 		for(i = 1 ; i <= iRHSCount ; i++)
 		{
 			oRk = oRHS.GetColumn(i);
			if(i == 1)
			{
				oX.SetSize(m_iRowsCount,1);
			}
			else
			{
				oX.ZeroValues();
			}
 			dTolerance = dToleranceFactor*oRk.GetNorm();
 			oZ = Matrix::MultiplyDiagonalMatrixAndVector(m_oPreconditioner,oRk);
 			oPk = oZ;
 			dErrorNorm = 100*dTolerance;
 			while(dErrorNorm > dTolerance)
 			{
 				oTemp = (*this)*oPk;
 				dNum = (oRk.GetTranspose()*oZ).Get(1,1);
 				dDen = (oPk.GetTranspose()*oTemp).Get(1,1);
 				dAlpha = dNum/dDen;
 				oX = oX + oPk*dAlpha;
 				oRkp1 = oRk - oTemp*dAlpha;
 				dErrorNorm = oRkp1.GetNorm();
 				oZ = Matrix::MultiplyDiagonalMatrixAndVector(m_oPreconditioner,oRkp1);
 				dDen = dNum;
 				dNum = (oZ.GetTranspose()*oRkp1).Get(1,1);
 				dBeta = dNum/dDen;
 				oPkp1 = oZ + oPk*dBeta;
 				oPk = oPkp1;
 				oRk = oRkp1;
 			}
 			oResult.SetColumn(i,oX);
 		}
 		return oResult;
	}
// 	Matrix SparseMatrix::SolveGMRES(const Matrix& oRHSColumn,Matrix* poInitialSolution)
// 	{
// 		//////////////////////////////////////////////////
// 		// solution using the Generalized Minimum Residual method
// 		//////////////////////////////////////////////////
// 		// the GMRES matrices have already been initialized, zero out all of their entries
// 		m_oGMRESQMatrix.ZeroValues();
// 		m_oGMRESHessenbergMatrix.ZeroValues();
// 		unsigned int i = 0;
// 		Matrix oRHS(m_iRowsCount,1);
// 		for(i = 1 ; i <= m_iRowsCount ; i++)
// 		{
// 			oRHS.Set(i,1,oRHSColumn.Get(i,1)*m_oPreconditioner.Get(1,i));
// 		}
// 		// the right hand side can only be one column, because there is no performance gain 
// 		// when solving the system with multiple RHS columns over solving them one at a time
// 		InitializeArnodliIteration(oRHS);
//         unsigned int iCurrentStep = 1;
//         double dTolerance = 1.0E-8*m_dGMRESBeta;
// 		Matrix oX(m_iRowsCount,1);
// 		Matrix oQ;
// 		Matrix oR;
// 		Matrix oB;
// 		Matrix oY;
// 		
// 		double dError = 0.0;
// 		while(iCurrentStep <= m_iRowsCount)
// 		{
// 			dError = ((*this)*oX - oRHS).GetNorm();
// 			printf("%d : error : %e tol : %e norm X %e\n",iCurrentStep,dError,dTolerance,oX.GetNorm());
// 			if(dError <= dTolerance)
// 			{
// 				break;
// 			}
// 			PerformArnoldiIterationStep(iCurrentStep);
// 			m_oGMRESHessenbergMatrix.QRFactorize(oQ,oR,iCurrentStep + 1,iCurrentStep);
// 			oB.SetSize(iCurrentStep,1);
// 			for(i = 1 ; i <= iCurrentStep ; i++)
// 			{
// 				oB.Set(i,1,oQ.Get(1,i));
// 			}
// 			oB = oB*m_dGMRESBeta;
// 			oY = oR.UpperTriangularSolve(oB);
// 			oX = m_oGMRESQMatrix.LimitedColumnProduct(iCurrentStep,oY);
// 			iCurrentStep = iCurrentStep + 1;
// 		}
// 		return oX;
// 	}
	Matrix SparseMatrix::SolveGMRES(const Matrix& oRHS,Matrix* poInitialSolution)
	{
		//////////////////////////////////////////////////
		// solution using the Generalized Minimum Residual method
		//////////////////////////////////////////////////
		// the right hand side can only be one column, because there is no performance gain 
		// when solving the system with multiple RHS columns over solving them one at a time
		Matrix oX;
		Matrix oX0;
		if(poInitialSolution != NULL)
		{
			oX0 = *poInitialSolution;
		}
		else
		{
			oX0.SetSize(m_iRowsCount,1);
		}
		oX = oX0;
		unsigned int iRestartLimit = m_iRowsCount;
		unsigned int iCurrentStep = 0;
		Matrix oQ;
		Matrix oR;
		Matrix oB;
		Matrix oY;
		double dTolerance = 1.0E-8*oRHS.GetNorm();
		double dError = 100.0*dTolerance;
		unsigned int i = 0;
		while(dError > dTolerance)
		{
			// the GMRES matrices have already been initialized, zero out all of their entries
			m_oGMRESQMatrix.ZeroValues();
			m_oGMRESHessenbergMatrix.ZeroValues();
			oX0 = oX;
			oR = oRHS - (*this)*oX;
			InitializeArnodliIteration(oR);
			for(iCurrentStep = 1 ; iCurrentStep <= iRestartLimit ; iCurrentStep++)
			{
				dError = ((*this)*oX - oRHS).GetNorm();
				printf("%d : error : %e tol : %e norm X %e\n",iCurrentStep,dError,dTolerance,oX.GetNorm());
				if(dError <= dTolerance)
				{
					break;
				}
				PerformArnoldiIterationStep(iCurrentStep);
				m_oGMRESHessenbergMatrix.QRFactorize(oQ,oR,iCurrentStep + 1,iCurrentStep);
				oB.SetSize(iCurrentStep,1);
				for(i = 1 ; i <= iCurrentStep ; i++)
				{
					oB.Set(i,1,oQ.Get(1,i));
				}
				oB = oB*m_dGMRESBeta;
				oY = oR.UpperTriangularSolve(oB);
				oX = oX0 + Matrix::MultiplyDiagonalMatrixAndVector(m_oPreconditioner,m_oGMRESQMatrix.LimitedColumnProduct(iCurrentStep,oY));
			}
		}
		return oX;
	}
	Matrix SparseMatrix::SolveSteepestDescent(const Matrix& oRHS) const
	{
		//////////////////////////////////////////////////
		// solution using the steepest descent method
		//////////////////////////////////////////////////
		unsigned int iRHSCount = oRHS.GetColumnsCount();
 		Matrix oResult(oRHS.GetRowsCount(),iRHSCount);
 		if(oRHS.GetRowsCount() != m_iRowsCount)
 		{
 			printf("error: sparse matrix steepest descent solution size mismatch\n");
 			return oResult;
 		}
 		unsigned int i = 0;
 		Matrix oRk;
 		double dAlpha = 0.0;
 		double dToleranceFactor = 1E-12;
 		double dTolerance = 0.0;
 		double dErrorNorm = 0.0;
 		Matrix oTemp;
 		double dNum = 0.0;
 		double dDen = 0.0;
 		Matrix oX;
 		// solve the system
 		for(i = 1 ; i <= iRHSCount ; i++)
 		{
 			oRk = oRHS.GetColumn(i);
 			oX.SetSize(oRk.GetRowsCount(),1);
 			dTolerance = dToleranceFactor*oRk.GetNorm();
 			dErrorNorm = 100*dTolerance;
 			while(dErrorNorm > dTolerance)
 			{
 				oTemp = (*this)*oRk;
 				dNum = (oRk.GetTranspose()*oRk).Get(1,1);
 				dDen = (oRk.GetTranspose()*oTemp).Get(1,1);
 				dAlpha = dNum/dDen;
 				oX = oX + oRk*dAlpha;
 				oRk = oRk - oTemp*dAlpha;
 				dErrorNorm = oRk.GetNorm();
 			}
 			oResult.SetColumn(i,oX);
 		}
 		return oResult;
	}
	Matrix SparseMatrix::SolveGaussElimination(const Matrix& oRHS) const
	{		
		////////////////////////////////////////////
		// solution using Gauss elimination with partial pivoting
		////////////////////////////////////////////
		unsigned int iRHSCount = oRHS.GetColumnsCount();
		Matrix oResult(m_iRowsCount,iRHSCount);
		if(m_iRowsCount != oRHS.GetRowsCount())
		{
			printf("error: sparse matrix Gauss elimination solution size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		vector<double> vdScales;
		vdScales.resize(m_iRowsCount);
		double dTemp = 0.0;
		double dTemp2 = 0.0;
 		map<unsigned int,double>::iterator liRowIterator;
 		map<unsigned int,double>* poRow = NULL;
		// get the maximum for each row and place it in the scales vector
		for(i = 0; i < m_iRowsCount ; i++)
		{
			poRow = m_ppoData[i];
			if(poRow->empty())
			{
				continue;
			}
			dTemp = DBL_MIN;
			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
			{
				dTemp2 = fabs(liRowIterator->second);
				if(dTemp2 > dTemp)
				{
					dTemp = dTemp2;
				}
			}
			vdScales[i] = dTemp;
		}
		// set the temp matrix with the row wise scaled version
		// of the original system matrix
		SparseMatrix oTemp;
		oTemp.SetRowsCount(m_iRowsCount);
		Matrix oTempRHS;
		oTempRHS.SetSize(oRHS.GetRowsCount(),oRHS.GetColumnsCount());
		vector<unsigned int> viIndices;		// used for virtual row swapping
		viIndices.resize(m_iRowsCount);
		for(i = 0; i < m_iRowsCount ; i++)
		{
			poRow = m_ppoData[i];
			if(poRow->empty())
			{
				continue;
			}
			// scale matrix row
			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
			{
				oTemp.Set(i + 1,liRowIterator->first,liRowIterator->second/vdScales[i]);
			}
			// scale RHS
			for(j = 0; j < iRHSCount ; j++)
			{
				oTempRHS.Set(i + 1,j + 1,oRHS.Get(i + 1,j + 1)/vdScales[i]);
			}
			viIndices[i] = i;
		}
		
		unsigned int iPivotIndex = 0;
		unsigned int k = 0;
		double dDiagonalEntry = 0.0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			oTemp.Get(viIndices[i] + 1,i + 1,dTemp2);
			dTemp = fabs(dTemp2);
			iPivotIndex = i;
			for(j = i + 1; j < m_iRowsCount ; j++)
			{
				if(oTemp.Get(viIndices[j] + 1,i + 1,dTemp2))
				{
					dTemp2 = fabs(dTemp2);
					if(dTemp2 > dTemp)
					{
						dTemp = dTemp2;
						iPivotIndex = j;
					}
				}
			}
			if(dTemp < 1E-50)
			{
				return oResult;
			}
			
			j = viIndices[i];
			viIndices[i] = viIndices[iPivotIndex];
			viIndices[iPivotIndex] = j;			
 			oTemp.Get(viIndices[i] + 1,i + 1,dDiagonalEntry);
 			for(j = i + 1; j < m_iRowsCount; j++)
 			{
 				if(oTemp.Get(viIndices[j] + 1,i + 1,dTemp2))
 				{
 					poRow = oTemp.m_ppoData[viIndices[j]];
 					if(poRow->empty())
 					{
 						continue;
 					}
					dTemp = -dTemp2/dDiagonalEntry;
					// drop the entry which is guaranteed to be zero
					oTemp.DropEntry(viIndices[j] + 1,i + 1);
					for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
					{
						if(liRowIterator->first <= i)
						{
							continue;
						}
						if(oTemp.Get(viIndices[i] + 1,liRowIterator->first,dTemp2))
						{
							oTemp.AddToEntry(viIndices[j] + 1,liRowIterator->first,dTemp*dTemp2);
						}
					}
					
					poRow = oTemp.m_ppoData[viIndices[i]];
					for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
					{
						if(!(oTemp.Get(viIndices[j] + 1,liRowIterator->first,dTemp2)))
						{
							oTemp.Get(viIndices[i] + 1,liRowIterator->first,dTemp2);
							oTemp.AddToEntry(viIndices[j] + 1,liRowIterator->first,dTemp*dTemp2);
						}
					}
					
					for(k = 0; k < iRHSCount; k++)
					{
						oTempRHS.Set(viIndices[j] + 1,k + 1,oTempRHS.Get(viIndices[j] + 1,k + 1) + dTemp*oTempRHS.Get(viIndices[i] + 1,k + 1));
					}
 				}
 			}
		}

		// back substitution to get the final solution
		for(int i = m_iRowsCount - 1 ; i >= 0; i--)
		{
			oTemp.Get(viIndices[i] + 1,i + 1,dDiagonalEntry);
			for(k = 0; k < iRHSCount ; k++)
			{	
				dTemp = 0.0;
				poRow = oTemp.m_ppoData[viIndices[i]];
				for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
				{
					if(liRowIterator->first <= (i + 1))
					{
						continue;
					}
					if(oTemp.Get(viIndices[i] + 1,liRowIterator->first,dTemp2))
					{
						dTemp = dTemp + dTemp2*oTempRHS.Get(viIndices[liRowIterator->first - 1] + 1,k + 1);
					}
				}
				oTempRHS.Set(viIndices[i] + 1,k + 1,(oTempRHS.Get(viIndices[i] + 1,k + 1) - dTemp)/dDiagonalEntry);
			}
		}

		for(i = 0; i < m_iRowsCount; i++)
		{
			for(j = 0; j < iRHSCount ; j++)
			{
				oResult.Set(i + 1,j + 1,oTempRHS.Get(viIndices[i] + 1,j + 1));
			}
		}
		return oResult;
	}
 	Matrix SparseMatrix::SolveGaussSeidel(const Matrix& oRHS) const
 	{
 		//////////////////////////////////////////////////
 		// solution using the Gauss-Seidel method
 		//////////////////////////////////////////////////
 		unsigned int iRHSCount = oRHS.GetColumnsCount();
  		Matrix oResult(oRHS.GetRowsCount(),iRHSCount);
  		if(oRHS.GetRowsCount() != m_iRowsCount)
  		{
  			printf("error: sparse matrix Gauss-Seidel solution size mismatch\n");
  			return oResult;
  		}
  		unsigned int i = 0;
  		unsigned int j = 0;
  		double dToleranceFactor = 1E-12;
  		double dTolerance = 0.0;
  		double dErrorNorm = 0.0;
  		Matrix oX;
  		Matrix oB;
		double dTemp = 0.0;
		double dDiagonalElement = 0.0;
		map<unsigned int,double>::iterator liRowIteraor;
 		map<unsigned int,double>* poRow = NULL;
  		// solve the system
  		for(i = 1 ; i <= iRHSCount ; i++)
  		{
  			oB = oRHS.GetColumn(i);
  			oX = oB;
  			dTolerance = dToleranceFactor*oX.GetNorm();
   			dErrorNorm = 100*dTolerance;
   			while(dErrorNorm > dTolerance)
   			{
 				for(j = 1; j <= m_iRowsCount ; j++)
 				{
 					if(!Get(j,j,dDiagonalElement))
 					{
 						return oResult;
 					}
 					
 					poRow = m_ppoData[j - 1];
					if(poRow->empty())
					{
						continue;
					}
					dTemp = 0.0;
					for(liRowIteraor = poRow->begin() ; liRowIteraor != poRow->end() ; liRowIteraor++)
 					{
 						if(liRowIteraor->first == j)
 						{
 							continue;
 						}
 						dTemp = dTemp + liRowIteraor->second*oX.Get(liRowIteraor->first,1);
 					}
 					oX.Set(j,1,(oB.Get(j,1) - dTemp)/dDiagonalElement);
 				}
 				dErrorNorm = ((*this)*oX - oB).GetNorm();
 			}
 			oResult.SetColumn(i,oX);
  		}
  		return oResult;
 	}
 	 Matrix SparseMatrix::SolveJacobi(const Matrix& oRHS) const
 	{
 		//////////////////////////////////////////////////
 		// solution using the Jacobi method
 		//////////////////////////////////////////////////
 		unsigned int iRHSCount = oRHS.GetColumnsCount();
  		Matrix oResult(oRHS.GetRowsCount(),iRHSCount);
  		if(oRHS.GetRowsCount() != m_iRowsCount)
  		{
  			printf("error: sparse matrix Jacobi solution size mismatch\n");
  			return oResult;
  		}
  		unsigned int i = 0;
  		unsigned int j = 0;
  		map<unsigned int,double>::iterator liRowIteraor;
 		map<unsigned int,double>* poRow = NULL;
  		double dToleranceFactor = 1E-12;
  		double dTolerance = 0.0;
  		double dErrorNorm = 0.0;
  		Matrix oX;
  		Matrix oXTemp;
  		Matrix oB;
		double dTemp = 0.0;
		double dDiagonalElement = 0.0;

  		// solve the system
  		for(i = 1 ; i <= iRHSCount ; i++)
  		{
  			oB = oRHS.GetColumn(i);
  			oX = oB;
  			dTolerance = dToleranceFactor*oX.GetNorm();
   			dErrorNorm = 100*dTolerance;
   			while(dErrorNorm > dTolerance)
   			{
   				oXTemp.SetSize(m_iRowsCount,1);
 				for(j = 1; j <= m_iRowsCount ; j++)
 				{
 					if(!Get(j,j,dDiagonalElement))
 					{
 						return oResult;
 					}
 					
 					poRow = m_ppoData[j - 1];
					if(poRow->empty())
					{
						continue;
					}
					dTemp = 0.0;
					for(liRowIteraor = poRow->begin() ; liRowIteraor != poRow->end() ; liRowIteraor++)
 					{
 						if(liRowIteraor->first == j)
 						{
 							continue;
 						}
 						dTemp = dTemp + liRowIteraor->second*oX.Get(liRowIteraor->first,1);
 					}
 					oXTemp.Set(j,1,(oB.Get(j,1) - dTemp)/dDiagonalElement);
 				}
 				dErrorNorm = (oX - oXTemp).GetNorm();
 				oX = oXTemp;
 			}
 			oResult.SetColumn(i,oX);
  		}
  		return oResult;
 	}
	void SparseMatrix::StoreCRS(double*& pdValues,unsigned int*& piRowsStartingLocations,unsigned int*& piColumnsIndices) const
	{
		unsigned int iCRSSize = GetCRSSize();
		unsigned int i = 0;
		unsigned int iCount = 0;
		pdValues = new double[iCRSSize];
		piColumnsIndices = new unsigned int[iCRSSize];
		piRowsStartingLocations = new unsigned int[m_iRowsCount + 1];
		map<unsigned int,double>* poRow = NULL;
		map<unsigned int,double>::iterator liRowIterator;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			piRowsStartingLocations[i] = iCount;
			poRow = m_ppoData[i];
			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
			{
				piColumnsIndices[iCount] = liRowIterator->first;
				pdValues[iCount] = liRowIterator->second;
				iCount = iCount + 1;
			}
		}
		piRowsStartingLocations[m_iRowsCount] = iCRSSize;
	}
	void SparseMatrix::SetFromCRS(const unsigned int& iRowsCount,double* pdValues,unsigned int* piRowsStartingLocations,unsigned int* piColumnsIndices)
	{
		Reset();
		SetRowsCount(iRowsCount);
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = piRowsStartingLocations[i] ; j < piRowsStartingLocations[i + 1] ; j++)
			{
				Set(i + 1,piColumnsIndices[j],pdValues[j]);
			}
		}
	}
	void SparseMatrix::AddCRS(const unsigned int& iRowsCount,double* pdValues,unsigned int* piRowsStartingLocations,unsigned int* piColumnsIndices)
	{
		if(m_iRowsCount != iRowsCount)
		{
			return;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = piRowsStartingLocations[i] ; j < piRowsStartingLocations[i + 1] ; j++)
			{
				AddToEntry(i + 1,piColumnsIndices[j],pdValues[j]);
			}
		}
	}
	void SparseMatrix::SubtractCRS(const unsigned int& iRowsCount,double* pdValues,unsigned int* piRowsStartingLocations,unsigned int* piColumnsIndices)
	{
		if(m_iRowsCount != iRowsCount)
		{
			return;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(j = piRowsStartingLocations[i] ; j < piRowsStartingLocations[i + 1] ; j++)
			{
				AddToEntry(i + 1,piColumnsIndices[j],-pdValues[j]);
			}
		}
	}
	unsigned int SparseMatrix::GetMaximumColumnsCount() const
	{
		unsigned int i = 0;
		unsigned int iMaxCount = 0;
		unsigned int iTemp = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			iTemp = m_ppoData[i]->size();
			if(iMaxCount < iTemp)
			{
				iMaxCount = iTemp;
			}
		}
		return iMaxCount;
	}
	unsigned int SparseMatrix::GetMaximumColumnIndex() const
	{
		unsigned int i = 0;
		unsigned int iMaxIndex = 0;
		unsigned int iTemp = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			if(m_ppoData[i]->empty())
			{
				continue;
			}
			iTemp = m_ppoData[i]->rbegin()->first;
			if(iMaxIndex < iTemp)
			{
				iMaxIndex = iTemp;
			}
		}
		return iMaxIndex;
	}
	unsigned int SparseMatrix::GetEntriesCount() const
	{
		unsigned int i = 0;
		unsigned int iSum = 0;
		for(i = 0; i < m_iRowsCount ; i++)
		{
			iSum = iSum + m_ppoData[i]->size();
		}
		return iSum;
	}
	unsigned int SparseMatrix::GetCRSSize() const
	{
		unsigned int i = 0;
		unsigned int iSum = 0;
		vector<unsigned int> viSizes = GetRowsNonzeroCount();
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			iSum = iSum + viSizes[i];
		}
		return iSum;
	}
	SparseMatrix SparseMatrix::operator+(const SparseMatrix& oMatrix) const
	{
		SparseMatrix oResult = *this;
		if(m_iRowsCount != oMatrix.m_iRowsCount)
		{
			printf("error: sparse matrix-sparse matrix addition size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		map<unsigned int,double>* poRow = NULL;
		map<unsigned int,double>::iterator miRow;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			poRow = oMatrix.m_ppoData[i];
			for(miRow = poRow->begin() ; miRow != poRow->end() ; miRow++)
			{
				oResult.AddToEntry(i + 1,miRow->first,miRow->second);
			}
		}
		return oResult;
	}
	SparseMatrix SparseMatrix::operator-(const SparseMatrix& oMatrix) const
	{
		SparseMatrix oResult = *this;
		if(m_iRowsCount != oMatrix.m_iRowsCount)
		{
			printf("error: sparse matrix-sparse matrix subtraction size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		map<unsigned int,double>* poRow = NULL;
		map<unsigned int,double>::iterator miRow;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			poRow = oMatrix.m_ppoData[i];
			for(miRow = poRow->begin() ; miRow != poRow->end() ; miRow++)
			{
				oResult.AddToEntry(i + 1,miRow->first,-miRow->second);
			}
		}
		return oResult;
	}
	Matrix SparseMatrix::operator*(const Matrix& oMatrix) const
	{
		unsigned int iTemp = GetMaximumColumnIndex();
		unsigned int iColumnsCount = oMatrix.GetColumnsCount();
		Matrix oResult(m_iRowsCount,iColumnsCount);
		if(iTemp > oMatrix.GetRowsCount())
		{
			printf("error: sparse matrix-matrix multiplication size mismatch\n");
			return oResult;
		}
		unsigned int i = 0;
		unsigned int j = 0;
		double dTemp = 0.0;
		map<unsigned int,double>::iterator liRowIterator;
		map<unsigned int,double>* poRow = NULL;

		for(i = 0; i < m_iRowsCount ; i++)
		{
			poRow = m_ppoData[i];
			if(poRow->empty())
			{
				continue;
			}
			for(j = 0; j < iColumnsCount ; j++)
			{
				dTemp = 0.0;
				for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
				{
					iTemp = liRowIterator->first;
					dTemp = dTemp + liRowIterator->second*oMatrix.Get(iTemp,j + 1);
				}
				oResult.Set(i + 1,j + 1,dTemp);
			}
		}
		return oResult;
	}
	vector<unsigned int> SparseMatrix::GetRowsNonzeroCount() const
	{
		unsigned int i = 0;
		vector<unsigned int> viRowsNonzeroCount;
		viRowsNonzeroCount.resize(m_iRowsCount);
		for(i = 0; i < m_iRowsCount ; i++)
		{
			viRowsNonzeroCount[i] = (unsigned int)m_ppoData[i]->size();
		}
		return viRowsNonzeroCount;
	}
	double SparseMatrix::SumAllEntries() const
	{
		double dSum = 0.0;
		unsigned int i = 0;
		map<unsigned int,double>::iterator miRow;
		for(i = 0 ; i < m_iRowsCount ; i++)
		{
			for(miRow = m_ppoData[i]->begin() ; miRow != m_ppoData[i]->end() ; miRow++)
			{
				dSum = dSum + miRow->second;
			}
		}
		return dSum;
	}
	double SparseMatrix::GetEigenValueUpperBound() const
	{
		unsigned int i = 0;
		map<unsigned int,double>::iterator liRowIterator;
		map<unsigned int,double>* poRow = NULL;
		double dCenter = 0.0;
		double dRadius;
		double dUpperBound = -DBL_MAX;
		double dTemp = 0.0;
		for(i = 1 ; i <= m_iRowsCount ; i++)
		{
			poRow = m_ppoData[i - 1];
			if(poRow->empty())
			{
				continue;
			}
			Get(i,i,dCenter);
			dRadius = 0.0;
			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
			{
				if(liRowIterator->first != i)
				{
					dRadius = dRadius + fabs(liRowIterator->second);
				}
			}
			dTemp = dCenter + dRadius;
			if(dTemp > dUpperBound)
			{
				dUpperBound = dTemp;
			}
		}
		return dUpperBound;
	}
	void SparseMatrix::MultiplyRow(const unsigned int& iRowIndex,const double& dFactor)
	{
		if(iRowIndex > m_iRowsCount)
		{
			printf("error: sparse matrix row-factor multiplication index mismatch\n");
			return;
		}
		map<unsigned int,double>* poRow = poRow = m_ppoData[iRowIndex - 1];
		if(poRow->empty())
		{
			return;
		}
		map<unsigned int,double>::iterator liRowIterator;
		for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
		{
			Set(iRowIndex,liRowIterator->first,liRowIterator->second*dFactor);
		}
	}
	bool SparseMatrix::IsSymmetric(const double& dTolerance) const
	{
		unsigned int i = 0;
		map<unsigned int,double>::iterator liRowIterator;
		unsigned int j = 0;
		double dThisValue = 0.0;
		double dThatValue = 0.0;
		for(i = 1 ; i <= m_iRowsCount ; i++)
		{
			for(liRowIterator = m_ppoData[i - 1]->begin() ; liRowIterator != m_ppoData[i - 1]->end() ; liRowIterator++)
			{
				j = liRowIterator->first;
				dThisValue = liRowIterator->second;
				if(!Get(j,i,dThatValue))
				{
					return false;
				}
				if(fabs(dThisValue - dThatValue) > dTolerance)
				{
					return false;
				}
			}
		}
		return true;
	}
	/*SparseMatrix SparseMatrix::GetTranspose() const
	{
		unsigned int iTransposeRowsCount = GetMaximumColumnIndex();
		SparseMatrix oTranspose(iTransposeRowsCount);
		unsigned int i = 0;
		list<SparseMatrixEntry>* poRow = NULL;
		list<SparseMatrixEntry>::iterator liRowIterator;
		for(i = 1 ; i <= m_iRowsCount ; i++)
		{
			poRow = m_vpoData[i - 1];
			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
			{
				oTranspose.Set(liRowIterator->GetColumnIndex(),i,liRowIterator->GetValue());
			}
		}
		return oTranspose;
	}*/
	void SparseMatrix::BuildPreconditioner()
	{
		unsigned int i = 0;
		double dTemp = 0.0;
		m_oPreconditioner.SetSize(1,m_iRowsCount);
	 	for(i = 1; i <= m_iRowsCount ; i++)
 		{
 			if(Get(i,i,dTemp))
 			{
 				if(dTemp < 0.0)
 				{
 					dTemp = 1.0;
 				}
 			}
 			else
 			{
 				dTemp = 1.0;
 			}
 			m_oPreconditioner.Set(1,i,1.0/dTemp);
 		}
	}
	void SparseMatrix::InitializeGMRESMatrices()
	{
		m_oGMRESQMatrix.SetSize(m_iRowsCount,m_iRowsCount);
		m_oGMRESHessenbergMatrix.SetSize(m_iRowsCount,m_iRowsCount);
		BuildPreconditioner();
//		//multiply the preconditioner by the matrix
// 		unsigned int i = 0;
// 		map<unsigned int,double>::iterator liRowIterator;
// 		map<unsigned int,double>* poRow = NULL;
// 		double dTemp = 0.0;
// 		for(i = 0; i < m_iRowsCount ; i++)
// 		{
// 			poRow = m_ppoData[i];
// 			if(poRow->empty())
// 			{
// 				continue;
// 			}
// 			dTemp = m_oPreconditioner.Get(1,i + 1);
// 			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
// 			{
// 				liRowIterator->second = liRowIterator->second*dTemp;
// 			}
// 		}
	}
	void SparseMatrix::InitializeArnodliIteration(const Matrix& oB)
	{
		m_dGMRESBeta = oB.GetNorm();
		if(m_dGMRESBeta < 1.0E-18)
		{
			return;
		}
		unsigned int i = 0;
		for(i = 1 ; i <= m_iRowsCount ; i++)
		{
			m_oGMRESQMatrix.Set(i,1,oB.Get(i,1)/m_dGMRESBeta);
		}
	}
	void SparseMatrix::PerformArnoldiIterationStep(const unsigned int& iStepNumber)
	{
		if((iStepNumber < 1) || (iStepNumber > m_iRowsCount))
		{
			return;
		}
		
		unsigned int i = 0;
		unsigned int j = 0;
		double dTemp = 0.0;
		if(iStepNumber < m_iRowsCount)
		{
			Matrix oZ(m_iRowsCount,1);
			for(j = 1 ; j <= m_iRowsCount ; j++)
			{
				oZ.Set(j,1,m_oGMRESQMatrix.Get(j,iStepNumber)*m_oPreconditioner.Get(1,j));
			}
			Matrix oV = MultiplyByColumn(oZ);
			for(i = 1 ; i <= iStepNumber ; i++)
			{
				dTemp = m_oGMRESQMatrix.MultiplyColumnByColumn(i,oV);
				m_oGMRESHessenbergMatrix.Set(i,iStepNumber,dTemp);
				for(j = 1 ; j <= m_iRowsCount ; j++)
				{
					oV.AddToEntry(j,1,-dTemp*m_oGMRESQMatrix.Get(j,i));
				}
			}
			dTemp = oV.GetNorm();
			m_oGMRESHessenbergMatrix.Set(iStepNumber + 1,iStepNumber,dTemp);
			for(j = 1 ; j <= m_iRowsCount ; j++)
			{
				m_oGMRESQMatrix.Set(j,iStepNumber + 1,oV.Get(j,1)/dTemp);
			}
		}
		else
		{
			// generate the last column of H by from H = Q^t A Q
			unsigned int k = 0;
			double dValue = 0.0;
			for(k = 1 ; k <= m_iRowsCount ; k++)
			{
				dTemp = 0.0;
				for(i = 1 ; i <= m_iRowsCount ; i++)
				{
					for(j = 1 ; j <= m_iRowsCount ; j++)
					{
						if(Get(i,j,dValue))
						{
							dTemp = dTemp + m_oGMRESQMatrix.Get(i,k)*dValue*m_oGMRESQMatrix.Get(j,iStepNumber);
						}
					}
				}
				m_oGMRESHessenbergMatrix.Set(k,iStepNumber,dTemp);
			}
		}
	}
	Matrix SparseMatrix::MultiplyByColumn(const Matrix& oColumn) const
	{
		// this function is provided with no bound checks for efficiency, the passed arguments
		// have to meet the size requirements or else a segmentation fault will happen
		Matrix oResult(m_iRowsCount,1);

		unsigned int i = 0;
		double dTemp = 0.0;
		map<unsigned int,double>::iterator liRowIterator;
		map<unsigned int,double>* poRow = NULL;

		for(i = 0; i < m_iRowsCount ; i++)
		{
			poRow = m_ppoData[i];
			if(poRow->empty())
			{
				continue;
			}
			dTemp = 0.0;
			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
			{
				dTemp = dTemp + liRowIterator->second*oColumn.Get(liRowIterator->first,1);
			}
			oResult.Set(i + 1,1,dTemp);
		}
		return oResult;
	}
	Matrix SparseMatrix::MultiplyByMatrixColumn(const Matrix& oMatrix,const unsigned int& iColumnIndex) const
	{
		// this function is provided with no bound checks for efficiency, the passed arguments
		// have to meet the size requirements or else a segmentation fault will happen
		Matrix oResult(m_iRowsCount,1);
		if(iColumnIndex > oMatrix.GetColumnsCount())
		{
			printf("error: sparse matrix-matrix column multiplication index mismatch\n");
			return oResult;
		}

		unsigned int i = 0;
		double dTemp = 0.0;
		map<unsigned int,double>::iterator liRowIterator;
		map<unsigned int,double>* poRow = NULL;

		for(i = 0; i < m_iRowsCount ; i++)
		{
			poRow = m_ppoData[i];
			if(poRow->empty())
			{
				continue;
			}
			dTemp = 0.0;
			for(liRowIterator = poRow->begin() ; liRowIterator != poRow->end() ; liRowIterator++)
			{
				dTemp = dTemp + liRowIterator->second*oMatrix.Get(liRowIterator->first,iColumnIndex);
			}
			oResult.Set(i + 1,1,dTemp);
		}
		return oResult;
	}
	void SparseMatrix::Initialize()
	{
		m_iRowsCount = 0;
		m_ppoData = NULL;
		m_oPreconditioner.Reset();
		m_oGMRESQMatrix.Reset();
		m_oGMRESHessenbergMatrix.Reset();
		m_dGMRESBeta = 0.0;
	}
}


