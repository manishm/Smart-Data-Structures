// Copyright (C) 2001-2009 Vivien Mallet
//
// This file is part of the linear-algebra library Seldon,
// http://seldon.sourceforge.net/.
//
// Seldon is free software; you can redistribute it and/or modify it under the
// terms of the GNU Lesser General Public License as published by the Free
// Software Foundation; either version 2.1 of the License, or (at your option)
// any later version.
//
// Seldon is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
// more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with Seldon. If not, see http://www.gnu.org/licenses/.


#define SELDON_DEBUG_LEVEL_4
#define SELDON_WITH_CBLAS

#if !defined(SELDON_WITH_UMFPACK) && !defined(SELDON_WITH_SUPERLU) \
  && !defined(SELDON_WITH_MUMPS)
#define SELDON_WITH_UMFPACK
//#define SELDON_WITH_SUPERLU
//#define SELDON_WITH_MUMPS
#endif

#ifdef SELDON_WITH_MPI
#include "mpi.h"
#endif

#include "Seldon.hxx"
#include "SeldonSolver.hxx"

using namespace Seldon;

typedef complex<double> cpx;
double epsilon(1e-14);


template<class T, class Prop, class Storage,
	 class Allocator1, class Allocator2, class Allocator3>
void Solve(Matrix<T, Prop, Storage, Allocator1>& A,
	   Vector<T, VectFull, Allocator2>& x,
	   const Vector<T, VectFull, Allocator3>& b)
{
#ifdef SELDON_WITH_UMFPACK
  MatrixUmfPack<T> mat_lu;
#endif
#ifdef SELDON_WITH_SUPERLU
  MatrixSuperLU<T> mat_lu;
#endif
#ifdef SELDON_WITH_MUMPS
  MatrixMumps<T> mat_lu;
#endif

  // The initial matrix is erased during the factorization process if you want
  // to keep it, you have to call 'GetLU(A, mat_lu, true)'.
  GetLU(A, mat_lu);
  x = b;
  SolveLU(mat_lu, x);
}


template<class Vector>
bool CheckSolution(Vector& x)
{
  Vector xt(x.GetM());
  xt.Fill();
  Add(-1, x, xt);
  return Norm2(xt) <= epsilon;
}


int main(int argc, char **argv)
{
#ifdef SELDON_WITH_MPI
  MPI_Init(&argc, &argv);
#endif
  
  // Number of rows.
  int n = 5;
  
  // Results.
  bool success_rowsym_real, success_rowsym_complex, success_colsym_real,
    success_colsym_complex, success_row_real, success_row_complex,
    success_col_real, success_col_complex;
  
  // We test RowSymSparse (real numbers).
  {
    // construction of a matrix
    Matrix<double, Symmetric, RowSymSparse> A;
    int colind_[] = {0, 1, 1, 2, 2, 3, 3, 4, 4};
    int rowptr_[] = {0, 2, 4, 6, 8, 9};
    double values_[] = {2.0, -1.5, 3.0, 1.0, 4.0, -2.0, 1.5, -0.5, 2.5};
    int nnz = rowptr_[n];
    IVect Col(nnz), Row(n+1);
    DVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];
    
    A.SetData(n, n, Values, Row, Col);
    DISP(A);
    // Computation of right hand side.
    DVect x_sol(n), b_vec(n); b_vec.Zero();
    x_sol.Fill();
    Mlt(1.0, A, x_sol, b_vec);
    DISP(b_vec);
    x_sol.Zero();
    
    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    cout << "Solution  "<< endl;
    DISP(x_sol);
    
    success_rowsym_real = CheckSolution(x_sol);
  }
  
  // We test RowSymSparse (complex numbers).
  {
    Matrix<complex<double>, Symmetric, RowSymSparse> A;
    int colind_[] = {0, 1, 1, 2, 2, 3, 3, 4, 4};
    int rowptr_[] = {0, 2, 4, 6, 8, 9};
    complex<double> values_[] = {cpx(2.0,-4.0), cpx(-1.5,0.5), cpx(3.0,-1.0),
				 cpx(1.0,0.0), cpx(4.0,2.0), cpx(-2.0,1.0),
				 cpx(1.5,-1.0), cpx(-0.5,0.0), cpx(2.5,5.0)};
    int nnz = rowptr_[n];
    IVect Col(nnz), Row(n+1);
    ZVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];
    
    A.SetData(n, n, Values, Row, Col);
    DISP(A);
    // Computation of right hand side
    ZVect x_sol(n), b_vec(n);
    x_sol.Fill(); b_vec.Zero();
    Mlt(1.0, A, x_sol, b_vec);
    DISP(b_vec);
    x_sol.Zero();
    
    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    cout << "Solution " << endl;
    DISP(x_sol);
    success_rowsym_complex = CheckSolution(x_sol);
  }
  
  // We test ColSymSparse (real numbers).
  {
    Matrix<double, Symmetric, ColSymSparse> A;
    int rowind_[] = {0, 0, 1, 1, 2, 2, 3, 3, 4};
    int colptr_[] = {0, 1, 3, 5, 7, 9};
    double values_[] = {2.0, -1.5, 3.0, 1.0, 4.0, -2.0, 1.5, -0.5, 2.5};
    
    int nnz = colptr_[n];
    IVect Col(n+1), Row(nnz);
    DVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Row(i) = rowind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Col(i) = colptr_[i];
    
    A.SetData(n, n, Values, Col, Row);
    DISP(A);
    // Computation of right hand side.
    DVect x_sol(n), b_vec(n);
    x_sol.Fill(); b_vec.Zero();
    Mlt(1.0, A, x_sol, b_vec);
    DISP(b_vec);
    x_sol.Zero();
    
    // we solve linear system
    Solve(A, x_sol, b_vec);
    cout << "Solution " << endl;
    DISP(x_sol);
    success_colsym_real = CheckSolution(x_sol);
  }
  
  // We test ColSymSparse (complex numbers).
  {
    Matrix<cpx, Symmetric, ColSymSparse> A;
    int rowind_[] = {0, 0, 1, 1, 2, 2, 3, 3, 4};
    int colptr_[] = {0, 1, 3, 5, 7, 9};
    cpx values_[] = {cpx(2.0,-4.0), cpx(-1.5,0.5), cpx(3.0,-1.0),
		     cpx(1.0,0.0), cpx(4.0,2.0), cpx(-2.0,1.0),
		     cpx(1.5,-1.0), cpx(-0.5,0.0), cpx(2.5,5.0)};
    
    int nnz = colptr_[n];
    IVect Col(n+1), Row(nnz);
    ZVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Row(i) = rowind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Col(i) = colptr_[i];
    
    A.SetData(n, n, Values, Col, Row);
    DISP(A);
    // Computation of right hand side.
    ZVect x_sol(n), b_vec(n); b_vec.Fill(0);
    x_sol.Fill(); b_vec.Zero();
    Mlt(1.0, A, x_sol, b_vec);
    x_sol.Zero();
    
    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    cout << "Solution " << endl;
    DISP(x_sol);
    success_colsym_complex = CheckSolution(x_sol);
  }
  
  // We test RowSparse (real numbers).
  {
    Matrix<double, General, RowSparse> A;
    int colind_[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3, 4, 3, 4};
    int rowptr_[] = {0, 3, 6, 10, 13, 15};
    double values_[] = {2.0, -1.5, 0.4, -1.0, 3.0, 1.0, -2.0,
			1.3, 4.0, -2.0, 0.5, 1.5, -0.5, -0.8, 2.5};
    
    int nnz = rowptr_[n];
    IVect Col(nnz), Row(n+1);
    DVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];
    
    A.SetData(n, n, Values, Row, Col);
    DISP(A);
    // computation of right hand side
    DVect x_sol(n), b_vec(n);
    x_sol.Fill();  b_vec.Zero();
    Mlt(1.0, A, x_sol, b_vec);
    DISP(b_vec);
    x_sol.Zero();
    
    // we solve linear system
    Solve(A, x_sol, b_vec);
    cout << "Solution " << endl;
    DISP(x_sol);
    success_row_real = CheckSolution(x_sol);
  }
  
  // We test RowSparse (complex numbers).
  {
    Matrix<cpx, General, RowSparse> A;
    int colind_[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3, 4, 3, 4};
    int rowptr_[] = {0, 3, 6, 10, 13, 15};
    cpx values_[] = {cpx(2.0,-4.0), cpx(-1.5,0.5), cpx(-0.4,0.2),
		     cpx(1.2,0.2), cpx(3.0,-1.0), cpx(1.0,0.0), cpx(-1.8,0.5),
		     cpx(0.6,0.8), cpx(4.0,2.0), cpx(-2.0,1.0), cpx(1.5,0.4),
		     cpx(1.5,-1.0), cpx(-0.5,0.0), cpx(0.7,-0.5),
		     cpx(2.5,5.0)};
    
    int nnz = rowptr_[n];
    IVect Col(nnz), Row(n+1);
    ZVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Col(i) = colind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Row(i) = rowptr_[i];
    
    A.SetData(n, n, Values, Row, Col);
    DISP(A);
    // Computation of right hand side.
    ZVect x_sol(n), b_vec(n);
    x_sol.Fill();  b_vec.Zero();
    Mlt(1.0, A, x_sol, b_vec);
    DISP(b_vec);
    x_sol.Zero();
    
    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    cout << "Solution " << endl;
    DISP(x_sol);
    success_row_complex = CheckSolution(x_sol);
  }
  
  // We test ColSparse (real numbers).
  {
    Matrix<double, General, ColSparse> A;
    int rowind_[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3, 4, 3, 4};
    int colptr_[] = {0, 3, 6, 10, 13, 15};
    double values_[] = {2.0, -1.0, -2.0, -1.5, 3.0, 1.3, 0.4, 1.0,
			4.0, 0.5, -2.0, 1.5, -0.8, -0.5, 2.5};
    int nnz = colptr_[n];
    IVect Col(n+1), Row(nnz);
    DVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Row(i) = rowind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Col(i) = colptr_[i];
    
    A.SetData(n, n, Values, Col, Row);
    DISP(A);
    // Computation of right hand side.
    DVect x_sol(n), b_vec(n);
    x_sol.Fill();  b_vec.Zero();
    Mlt(1.0, A, x_sol, b_vec);
    DISP(b_vec);
    x_sol.Zero();
    
    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    cout << "Solution " << endl;
    DISP(x_sol);
    success_col_real = CheckSolution(x_sol);
  }
  
  // We test ColSparse (complex numbers).
  {
    Matrix<cpx, General, ColSparse> A;
    int rowind_[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 2, 3, 4, 3, 4};
    int colptr_[] = {0, 3, 6, 10, 13, 15};
    cpx values_[] = {cpx(2.0,-4.0), cpx(1.2,0.2), cpx(-1.8,0.5),
		     cpx(-1.5,0.5), cpx(3.0,-1.0), cpx(0.6,0.8),
		     cpx(-0.4,0.2), cpx(1.0,0.0), cpx(4.0,2.0), cpx(1.5,0.4),
		     cpx(-2.0,1.0), cpx(1.5,-1.0), cpx(0.7,-0.5),
		     cpx(-0.5,0.0), cpx(2.5,5.0)};
    
    int nnz = colptr_[n];
    IVect Col(n+1), Row(nnz);
    ZVect Values(nnz);
    for (int i = 0; i < nnz; i++)
      {
	Values(i) = values_[i];
	Row(i) = rowind_[i];
      }
    
    for (int i = 0; i <= n; i++)
      Col(i) = colptr_[i];
    
    A.SetData(n, n, Values, Col, Row);
    DISP(A);
    // Computation of right hand side.
    ZVect x_sol(n), b_vec(n);
    x_sol.Fill();  b_vec.Zero();
    Mlt(1.0, A, x_sol, b_vec);
    DISP(b_vec);
    x_sol.Zero();
    
    // We solve the linear system.
    Solve(A, x_sol, b_vec);
    cout << "Solution " << endl;
    DISP(x_sol);
    success_col_complex = CheckSolution(x_sol);
  }
  
  bool overall_success = true;
  if (!success_rowsym_real)
    {
      cout << "Error during inversion of RowSymSparse real matrix" << endl;
      overall_success = false;
    }
  
  if (!success_rowsym_complex)
    {
      cout << "Error during inversion of RowSymSparse complex matrix" << endl;
      overall_success = false;
    }
  
  if (!success_colsym_real)
    {
      cout << "Error during inversion of ColSymSparse real matrix" << endl;
      overall_success = false;
    }
  
  if (!success_colsym_complex)
    {
      cout << "Error during inversion of ColSymSparse complex matrix" << endl;
      overall_success = false;
    }
  
  if (!success_row_real)
    {
      cout << "Error during inversion of RowSparse real matrix" << endl;
      overall_success = false;
    }
  
  if (!success_row_complex)
    {
      cout << "Error during inversion of RowSparse complex matrix" << endl;
      overall_success = false;
    }
  
  if (!success_col_real)
    {
      cout << "Error during inversion of ColSparse real matrix" << endl;
      overall_success = false;
    }
  
  if (!success_col_complex)
    {
      cout << "Error during inversion of ColSparse complex matrix" << endl;
      overall_success = false;
    }
  
  if (overall_success)
    cout << "All tests successfully completed" << endl;
  
  return 0;
}
