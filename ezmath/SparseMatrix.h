// EZMath library by Ahmed M. Hussein
// July 2009
// Feel free to use and distribute, but please cite or acknowledge, thanks. 

#ifndef SPARSEMATRIX_H_
#define SPARSEMATRIX_H_

// this is a sparse matrix class, it supports optimized matrix storage and several
// iterative and direct solvers (direct solvers aren't recommended for sparse matrices)

// by Ahmed M. Hussein
// October 2011

#include "vector"
#include "map"
#include "Matrix.h"

using namespace std;

namespace EZ
{
	class SparseMatrix
	{
	public:
		SparseMatrix();
		SparseMatrix(const unsigned int& iRowsCount);
		SparseMatrix(const SparseMatrix& oMatrix);
		~SparseMatrix();
		SparseMatrix& operator=(const SparseMatrix& oMatrix);
		void Reset();
		void SetRowsCount(const unsigned int& iCount);
		unsigned int GetRowsCount() const;
		bool Get(const unsigned int& iRowIndex,const unsigned int& iColumnIndex,double& dValue) const;
		void Set(const unsigned int& iRowIndex,const unsigned int& iColumnIndex,const double& dValue);
		void DropEntry(const unsigned int& iRowIndex,const unsigned int& iColumnIndex);
		void AddToEntry(const unsigned int& iRowIndex,const unsigned int& iColumnIndex,const double& dValue);
		unsigned int GetMaximumColumnsCount() const;
		unsigned int GetMaximumColumnIndex() const;
		Matrix SolveConjugateGradient(const Matrix& oRHS) const;
		Matrix SolveGMRES(const Matrix& oRHSColumn,Matrix* poInitialSolution = NULL);
		Matrix SolveSteepestDescent(const Matrix& oRHS) const;
		Matrix SolveGaussElimination(const Matrix& oRHS) const;
		Matrix SolveGaussSeidel(const Matrix& oRHS) const;
		Matrix SolveJacobi(const Matrix& oRHS) const;
		void StoreCRS(double*& pdValues,unsigned int*& piRowsStartingLocations,unsigned int*& piColumnsIndices) const;
		void SetFromCRS(const unsigned int& iRowsCount,double* pdValues,unsigned int* piRowsStartingLocations,unsigned int* piColumnsIndices);
		void AddCRS(const unsigned int& iRowsCount,double* pdValues,unsigned int* piRowsStartingLocations,unsigned int* piColumnsIndices);
		void SubtractCRS(const unsigned int& iRowsCount,double* pdValues,unsigned int* piRowsStartingLocations,unsigned int* piColumnsIndices);
		unsigned int GetEntriesCount() const;
		unsigned int GetCRSSize() const;
		virtual SparseMatrix operator+(const SparseMatrix& oMatrix) const;
		virtual SparseMatrix operator-(const SparseMatrix& oMatrix) const;
		virtual Matrix operator*(const Matrix& oMatrix) const;
		vector<unsigned int> GetRowsNonzeroCount() const;
		double SumAllEntries() const;
		double GetEigenValueUpperBound() const;
		void MultiplyRow(const unsigned int& iRowIndex,const double& dFactor);
		//virtual SparseMatrix GetTranspose() const;
		bool IsSymmetric(const double& dTolerance = 1.0E-6) const;
		void BuildPreconditioner();
		void InitializeGMRESMatrices();
		void InitializeArnodliIteration(const Matrix& oB);
		void PerformArnoldiIterationStep(const unsigned int& iStepNumber);
		Matrix MultiplyByColumn(const Matrix& oColumn) const;
		Matrix MultiplyByMatrixColumn(const Matrix& oMatrix,const unsigned int& iColumnIndex) const;
		
	private:

	protected:
		void Initialize();
		map<unsigned int,double>** m_ppoData;
		unsigned int m_iRowsCount;
		Matrix m_oPreconditioner;
		Matrix m_oGMRESQMatrix;
		double m_dGMRESBeta;
		Matrix m_oGMRESHessenbergMatrix;
	};
}

#endif



