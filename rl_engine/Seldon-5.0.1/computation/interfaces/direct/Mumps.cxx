// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_MUMPS_CXX

#include "Mumps.hxx"

namespace Seldon
{
  
  //! Mumps is called in double precision
  template<>
  inline void MatrixMumps<double>::CallMumps()
  {
    dmumps_c(&struct_mumps);
  }
  
  
  //! Mumps is called in complex double precision
  template<>
  inline void MatrixMumps<complex<double> >::CallMumps()
  {
    zmumps_c(&struct_mumps);
  }
  
  
  //! initialization
  template<class T>
  inline MatrixMumps<T>::MatrixMumps()
  {
#ifdef SELDON_WITH_MPI
    // MPI initialization for parallel version
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    // struct_mumps.comm_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);
    struct_mumps.comm_fortran = -987654;
      
    // parameters for mumps
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    struct_mumps.sym = 0; // 0 -> unsymmetric matrix
    
    // mumps is called
    CallMumps();
    
    // other parameters
    struct_mumps.n = 0;
    type_ordering = 7; // default : we let Mumps choose the ordering
    print_level = 0;
    out_of_core = false;
  }
  
  
  //! initialization of the computation
  template<class T> template<class MatrixSparse>
  inline void MatrixMumps<T>::InitMatrix(const MatrixSparse& A)
  {
    // we clear previous factorization
    Clear();
    
    // the user has to infom the symmetry during the initialization stage
    struct_mumps.job = -1;
    struct_mumps.par = 1;
    if (IsSymmetricMatrix(A))
      struct_mumps.sym = 2; // general symmetric matrix
    else
      struct_mumps.sym = 0; // unsymmetric matrix
    
    // mumps is called
    CallMumps();
    
    struct_mumps.icntl[13] = 20;
    struct_mumps.icntl[6] = type_ordering;
    // setting out of core parameters
    if (out_of_core)
      struct_mumps.icntl[21] = 1;
    else
      struct_mumps.icntl[21] = 0;

    struct_mumps.icntl[17] = 0;
    
    // the print level is set in mumps
    if (print_level >= 0)
      {
	struct_mumps.icntl[0] = 6;
	struct_mumps.icntl[1] = 0;
	struct_mumps.icntl[2] = 6;
	struct_mumps.icntl[3] = 2;
      }
    else
      {
	struct_mumps.icntl[0] = -1;
	struct_mumps.icntl[1] = -1;
	struct_mumps.icntl[2] = -1;
	struct_mumps.icntl[3] = 0;
      }
  }
  
  
  //! selects another ordering scheme
  template<class T>
  inline void MatrixMumps<T>::SelectOrdering(int num_ordering)
  {
    type_ordering = num_ordering;
  }
  
  
  //! clears factorization
  template<class T>
  MatrixMumps<T>::~MatrixMumps()
  {
    Clear();
  }
  
  
  //! clears factorization
  template<class T>
  inline void MatrixMumps<T>::Clear()
  {
    if (struct_mumps.n > 0)
      {
	struct_mumps.job = -2;
	CallMumps(); /* Terminate instance */
	
	struct_mumps.n = 0;
      }
  }
  
  
  //! no display from Mumps
  template<class T>
  inline void MatrixMumps<T>::HideMessages()
  {
    print_level = -1;
    
    struct_mumps.icntl[0] = -1;
    struct_mumps.icntl[1] = -1;
    struct_mumps.icntl[2] = -1;
    struct_mumps.icntl[3] = 0;
    
  }
  
  
  //! standard display
  template<class T>
  inline void MatrixMumps<T>::ShowMessages()
  {
    print_level = 0;
    
    struct_mumps.icntl[0] = 6;
    struct_mumps.icntl[1] = 0;
    struct_mumps.icntl[2] = 6;
    struct_mumps.icntl[3] = 2;
    
  }
  
  
  template<class T>
  inline void MatrixMumps<T>::EnableOutOfCore()
  {
    out_of_core = true;
  }


  template<class T>
  inline void MatrixMumps<T>::DisableOutOfCore()
  {
    out_of_core = false;
  }
  
  
  //! computes row numbers
  /*!
    \param[in,out] mat matrix whose we want to find the ordering
    \param[out] numbers new row numbers
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class Prop,class Storage,class Allocator>
  void MatrixMumps<T>::FindOrdering(Matrix<T, Prop, Storage, Allocator> & mat,
				    IVect& numbers, bool keep_matrix)
  {
    InitMatrix(mat);
    
    int n = mat.GetM(), nnz = mat.GetNonZeros();
    // conversion in coordinate format
    IVect num_row, num_col; Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    // no values needed to renumber
    values.Clear();
    if (!keep_matrix)
      mat.Clear();
    
    /* Define the problem on the host */
    if (rank == 0)
      {
	struct_mumps.n = n; struct_mumps.nz = nnz;
	struct_mumps.irn = num_row.GetData();
	struct_mumps.jcn = num_col.GetData();
      }
    
    /* Call the MUMPS package. */
    struct_mumps.job = 1; // we analyse the system
    CallMumps();
    
    numbers.Reallocate(n);
    for (int i = 0; i < n; i++)
      numbers(i) = struct_mumps.sym_perm[i]-1;
  }
  
  
  //! factorization of a given matrix
  /*!
    \param[in,out] mat matrix to factorize
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class Prop, class Storage, class Allocator>
  void MatrixMumps<T>::FactorizeMatrix(Matrix<T,Prop,Storage,Allocator> & mat,
				       bool keep_matrix)
  {
    InitMatrix(mat);
    
    int n = mat.GetM(), nnz = mat.GetNonZeros();
    // conversion in coordinate format with fortran convention (1-index)
    IVect num_row, num_col; Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row, num_col, values, 1);
    if (!keep_matrix)
      mat.Clear();
    
    // DISP(num_row); DISP(num_col); DISP(values); DISP(nnz); DISP(values.GetM()); DISP(n);
    /* Define the problem on the host */
    if (rank == 0)
      {
	struct_mumps.n = n; struct_mumps.nz = nnz;
	struct_mumps.irn = num_row.GetData(); struct_mumps.jcn = num_col.GetData();
	struct_mumps.a = reinterpret_cast<pointer>(values.GetData());
      }
    
    /* Call the MUMPS package. */
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();
  }
  
  
  //! returns information about factorization performed
  template<class T>
  int MatrixMumps<T>::GetInfoFactorization() const
  {
    return struct_mumps.info[0];
  }
  
  
  //! Computation of Schur complement.
  /*!
    \param[in,out] mat initial matrix.
    \param[in] num numbers to keep in Schur complement.
    \param[out] mat_schur Schur matrix.
    \param[in] keep_matrix if false, \a mat is cleared.
  */
  template<class T> template<class Prop1, class Storage1, class Allocator,
			     class Prop2, class Storage2, class Allocator2>
  void MatrixMumps<T>::
  GetSchurMatrix(Matrix<T, Prop1, Storage1, Allocator>& mat, const IVect& num,
		 Matrix<T, Prop2, Storage2, Allocator2> & mat_schur,
		 bool keep_matrix)
  {
    InitMatrix(mat);
    
    int n_schur = num.GetM(), n = mat.GetM();
    // Subscripts are changed to respect fortran convention
    IVect index_schur(n_schur);
    for (int i = 0; i < n_schur; i++)
      index_schur(i) = num(i)+1;
    
    // array that will contain values of Schur matrix
    Vector<T, VectFull, Allocator2> vec_schur(n_schur*n_schur);
    
    struct_mumps.icntl[18] = n_schur;
    struct_mumps.size_schur = n_schur;
    struct_mumps.listvar_schur = index_schur.GetData();
    struct_mumps.schur = reinterpret_cast<pointer>(vec_schur.GetData());
    
    // factorization of the matrix
    FactorizeMatrix(mat, keep_matrix);
    
    // resetting parameters related to Schur complement
    struct_mumps.icntl[18] = 0;
    struct_mumps.size_schur = 0;
    struct_mumps.listvar_schur = NULL;
    struct_mumps.schur = NULL;
    
    // schur complement stored by rows
    int nb = 0;
    mat_schur.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n ;j++)
	mat_schur(i,j) = vec_schur(nb++);
    
    vec_schur.Clear(); index_schur.Clear();
  }
  
  
  //! resolution of a linear system using the computed factorization
  /*!
    \param[in,out] x right-hand-side on input, solution on output
    It is assumed that a call to FactorizeMatrix has been done before
  */
  template<class T> template<class Allocator2, class Transpose_status>
  void MatrixMumps<T>::Solve(const Transpose_status& TransA,
			     Vector<T, VectFull, Allocator2>& x)
  {
    if (TransA.Trans())
      struct_mumps.icntl[8] = 0;
    else
      struct_mumps.icntl[8] = 1;
    
    struct_mumps.rhs = reinterpret_cast<pointer>(x.GetData());
    struct_mumps.job = 3; // we solve system
    CallMumps();
  }



  template<class T> template<class Allocator2>
  void MatrixMumps<T>::Solve(Vector<T, VectFull, Allocator2>& x)
  {
    Solve(SeldonNoTrans, x);
  }
  
  
#ifdef SELDON_WITH_MPI
  //! factorization of a given matrix in distributed form (parallel execution)
  /*!
    \param[inout] mat columns of the matrix to factorize
    \param[in] sym Symmetric or General
    \param[in] glob_number row numbers (in the global matrix)
    \param[in] keep_matrix if false, the given matrix is cleared
  */
  template<class T> template<class Prop, class Allocator>
  void MatrixMumps<T>::
  FactorizeDistributedMatrix(Matrix<T, General, ColSparse, Allocator> & mat,
			     const Prop& sym, const IVect& glob_number,
			     bool keep_matrix)
  {
    // initialization depending on symmetric of the matrix
    Matrix<T, Prop, RowSparse, Allocator> Atest;
    InitMatrix(Atest);
    
    // distributed matrix
    struct_mumps.icntl[17] = 3;

    // global number of rows : mat.GetM()
    int N = mat.GetM();
    int nnz = mat.GetNonZeros();
    // conversion in coordinate format with C-convention (0-index)
    IVect num_row, num_col; Vector<T, VectFull, Allocator> values;
    ConvertMatrix_to_Coordinates(mat, num_row,
				 num_col, values, 0);
    
    // we replace num_col with global numbers
    for (int i = 0; i < num_row.GetM(); i++)
      {
	num_row(i)++;
	num_col(i) = glob_number(num_col(i)) + 1;
      }
    
    if (!keep_matrix)
      mat.Clear();
    
    /* Define the problem on the host */
    struct_mumps.n = N; struct_mumps.nz_loc = nnz;
    struct_mumps.irn_loc = num_row.GetData();
    struct_mumps.jcn_loc = num_col.GetData();
    struct_mumps.a_loc = reinterpret_cast<pointer>(values.GetData());
        
    /* Call the MUMPS package. */
    struct_mumps.job = 4; // we analyse and factorize the system
    CallMumps();
    cout<<"Factorization completed"<<endl;
  }
  
  
  //! solves linear system with parallel execution
  /*!
    \param[in] TransA we solve A x = b or A^T x = b
    \param[inout] x right-hand-side then solution
    \param[inout] glob_num global row numbers
  */
  template<class T> template<class Allocator2, class Transpose_status>
  void MatrixMumps<T>::SolveDistributed(const Transpose_status& TransA,
					Vector<T, VectFull, Allocator2>& x,
					const IVect& glob_num)
  {
    Vector<T, VectFull, Allocator2> rhs;
    int cplx = sizeof(T)/8;
    // allocating the global right hand side
    rhs.Reallocate(struct_mumps.n); rhs.Zero();
    
    if (rank == 0)
      {
	// on the host, we retrieve datas of all the other processors
	int nb_procs; MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &nb_procs);
	if (nb_procs > 1)
	  {
	    // assembling the right hand side
	    Vector<T, VectFull, Allocator2> xp;
	    IVect nump;
	    for (int i = 0; i < nb_procs; i++)
	      {
		
		if (i != 0)
		  {
		    int nb_dof;
		    MPI_Recv(&nb_dof, 1, MPI_INT, i, 34, MPI_COMM_WORLD, &status);
		    xp.Reallocate(nb_dof);
		    nump.Reallocate(nb_dof);
		    MPI_Recv(xp.GetDataVoid(), cplx*nb_dof, MPI_DOUBLE, i, 35, MPI_COMM_WORLD, &status);
		    MPI_Recv(nump.GetData(), nb_dof, MPI_INT, i, 36, MPI_COMM_WORLD, &status);
		  }
		else
		  {
		    xp = x; nump = glob_num;
		  }
		
		for (int j = 0; j < nump.GetM(); j++)
		  rhs(nump(j)) = xp(j);
	      }
	  }
	else
	  Copy(x, rhs);
	
	struct_mumps.rhs = reinterpret_cast<pointer>(rhs.GetData());
      }
    else
      {
	// on other processors, we send solution
	int nb = x.GetM();
	MPI_Send(&nb, 1, MPI_INT, 0, 34, MPI_COMM_WORLD);
	MPI_Send(x.GetDataVoid(), cplx*nb, MPI_DOUBLE, 0, 35, MPI_COMM_WORLD);
	MPI_Send(glob_num.GetData(), nb, MPI_INT, 0, 36, MPI_COMM_WORLD);
      }
    
    // we solve system
    if (TransA.Trans())
      struct_mumps.icntl[8] = 0;
    else
      struct_mumps.icntl[8] = 1;
    
    struct_mumps.job = 3;
    CallMumps();
    
    // we distribute solution on all the processors
    MPI_Bcast(rhs.GetDataVoid(), cplx*rhs.GetM(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // and we extract the solution on provided numbers
    for (int i = 0; i < x.GetM(); i++)
      x(i) = rhs(glob_num(i));
  }


  template<class T> template<class Allocator2>
  void MatrixMumps<T>::SolveDistributed(Vector<T, VectFull, Allocator2>& x,
					const IVect& glob_num)
  {
    SolveDistributed(SeldonNoTrans, x, glob_num);
  }
#endif
  
  
  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T,Symmetric,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  
  
  template<class T, class Storage, class Allocator>
  void GetLU(Matrix<T,General,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  
  
  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T,Symmetric,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
		      const IVect& num, MatrixFull& schur_matrix, bool keep_matrix = false)
  {
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }
  
  
  template<class T, class Storage, class Allocator, class MatrixFull>
  void GetSchurMatrix(Matrix<T,General,Storage,Allocator>& A, MatrixMumps<T>& mat_lu,
		      const IVect& num, MatrixFull& schur_matrix, bool keep_matrix = false)
  {
    mat_lu.GetSchurMatrix(A, num, schur_matrix, keep_matrix);
  }
  
  
  template<class T, class Allocator>
  void SolveLU(MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }
  
  
  template<class T, class Allocator, class Transpose_status>
  void SolveLU(const Transpose_status& TransA,
	       MatrixMumps<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(TransA, x);
  }
  
}

#define SELDON_FILE_MUMPS_CXX
#endif
