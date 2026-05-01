// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef Matrix_H_
#define Matrix_H_

#include "vector"
#include "Vector.h"

using namespace std;

namespace EZ
{
	class Matrix
	{
	public:
		Matrix();
		Matrix(const Matrix& oMatrix);
		~Matrix();
		Matrix(const unsigned int& iRowsCount,const unsigned int& iColumnsCount);
		void SetSize(const unsigned int& iRowsCount,const unsigned int& iColumnsCount);
		void Set(const unsigned int& iRow,const unsigned int& iColumn,const double& dVal);
		double Get(const unsigned int& iRow,const unsigned int& iColumn) const;
		void Reset();
		Matrix& operator=(const Matrix& oMatrix);
		Matrix operator*(const Matrix& oMatrix) const;
		Vector operator*(const Vector& o3DVector) const;
		Matrix operator+(const Matrix& oMatrix) const;
		Matrix operator-(const Matrix& oMatrix) const;
		Matrix operator*(const double& dFactor) const;
		unsigned int GetRowsCount() const;
		unsigned int GetColumnsCount() const;
		Matrix Solve(const Matrix& oRHS) const;
		void FillMatrixFromVector(vector<double>* poVector);
		vector<double> Vectorize();
		void SetColumn(const unsigned int& iIndex,Matrix oColumn);
		void SetRow(const unsigned int& iIndex,Matrix oRow);
		void PlaceSubMatrix(const unsigned int& iEntryRowIndex,const unsigned int& iEntryColumnIndex,Matrix* poMatrix);
		Matrix GetInverse() const;
		double GetAbsoluteMaximum() const;
		Matrix GetRow(const unsigned int& iIndex) const;
		Matrix GetColumn(const unsigned int& iIndex) const;
		void Filter(const double& dTolerance = 1E-12);
		void SubtractFromIdentity();
		void AddToIdentity();
		void AddToEntry(const unsigned int& iRow,const unsigned int& iColumn,const double& dVal);
		Matrix GetTranspose() const;
		void ResetToZeros();
		static Matrix Invert3x3Matrix(const Matrix& oMatrix,double& dDeterminant);
		static Matrix Invert2x2Matrix(const Matrix& oMatrix,double& dDeterminant);
		static Matrix Solve3x3System(const Matrix& oSystem,const Matrix& oRHS);
		static Matrix Solve2x2System(const Matrix& oSystem,const Matrix& oRHS);
		void NormalizeColumns();
		void Randomize(const double& dMin,const double& dMax);
		void RandomizeNormal(const double& dMean,const double& dStandardDeviation);
		double GetNorm() const;
		bool IsSquare() const;
		Matrix GenerateGramSchmidtOrthogonalMatrix() const;
		static Matrix GenerateRandomOrthogonalMatrix(const unsigned int& iSize);
		static Matrix GenerateRandomSymmetricPositiveDefiniteMatrix(const unsigned int& iSize);
		static Matrix MultiplyDiagonalMatrixAndVector(const Matrix& oDiagonal,const Matrix& oVector);
		double SumAllEntries();
		void ZeroValues();
		void QRFactorize(Matrix& oQ,Matrix& oR) const;
		void QRFactorize(Matrix& oQ,Matrix& oR,unsigned int iMaxRow,unsigned int iMaxColumn) const;
		Matrix UpperTriangularSolve(const Matrix& oRHS) const;
		double MultiplyColumnByRow(const unsigned int& iColumnIndex,const Matrix& oRow) const;
		double MultiplyColumnByColumn(const unsigned int& iColumnIndex,const Matrix& oColumn) const;
		double MultiplyRowByRow(const unsigned int& iRowIndex,const Matrix& oRow) const;
		double MultiplyRowByColumn(const unsigned int& iRowIndex,const Matrix& oColumn) const;
		Matrix LimitedColumnProduct(const unsigned int& iMaxColumn,const Matrix& oColumn);
		void Print() const;
		
	private:
		unsigned int GetSSPDSize(const double& dTolerance = 1E-10) const;
		void StoreSSPD(double*& pdValues,unsigned int*& piRowsStartingLocations,unsigned int*& piColumnsIndices,const double& dTolerance = 1E-10) const;
		double** m_ppdData;
		unsigned int m_iRowsCount;
		unsigned int m_iColumnsCount;
	};
}

#endif

