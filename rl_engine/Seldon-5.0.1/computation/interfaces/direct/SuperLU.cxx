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


#ifndef SELDON_FILE_SUPERLU_CXX

#include "SuperLU.hxx"

namespace Seldon
{
  //! default constructor
  template<class T>
  MatrixSuperLU_Base<T>::MatrixSuperLU_Base()
  {
    //   permc_spec = 0: use the natural ordering
    //   permc_spec = 1: use minimum degree ordering on structure of A'*A
    //   permc_spec = 2: use minimum degree ordering on structure of A'+A
    //   permc_spec = 3: use approximate mininum degree column ordering
    n = 0;
    permc_spec = 2;
    StatInit(&stat);
    set_default_options(&options);
    ShowMessages();
  }
  
  
  //! destructor
  template<class T>
  MatrixSuperLU_Base<T>::~MatrixSuperLU_Base()
  {
    Clear();
  }
  
  
  //! same effect as a call to the destructor
  template<class T>
  void MatrixSuperLU_Base<T>::Clear()
  {
    if (n > 0)
      {
	// SuperLU objects are cleared
	Destroy_CompCol_Matrix(&A);
	Destroy_SuperMatrix_Store(&B);
	Destroy_SuperNode_Matrix(&L);
	Destroy_CompCol_Matrix(&U);
	n = 0;
      }
  }
  
  
  //! no message from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::HideMessages()
  {
    display_info = false;
  }
  
  
  //! allows messages from SuperLU
  template<class T>
  void MatrixSuperLU_Base<T>::ShowMessages()
  {
    display_info = true;
  }
  
  
  //! factorization of matrix in double precision using SuperLU
  template<class Prop,class Storage,class Allocator>
  void MatrixSuperLU<double>::
  FactorizeMatrix(Matrix<double,Prop,Storage,Allocator> & mat,
		  bool keep_matrix)
  {
    // clearing previous factorization
    Clear();
    
    // conversion in CSR format
    n = mat.GetN();
    Matrix<double,General,ColSparse> Acsr;
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();
    
    // we get renumbering vectors perm_r and perm_c
    int nnz = Acsr.GetDataSize();
    dCreate_CompCol_Matrix(&A, n, n, nnz, Acsr.GetData(), Acsr.GetInd(),
			   Acsr.GetPtr(), SLU_NC, SLU_D, SLU_GE);
    
    perm_r.Reallocate(n);
    perm_c.Reallocate(n);
    
    // factorization -> no right hand side
    int nb_rhs = 0, info;
    dCreate_Dense_Matrix(&B, n, nb_rhs, NULL, n, SLU_DN, SLU_D, SLU_GE);
    
    dgssv(&options, &A, perm_c.GetData(), perm_r.GetData(), &L, &U, &B, &stat, &info);
    
    if ((info==0)&&(display_info))
      {
	mem_usage_t mem_usage;
	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
	cout<<"No of nonzeros in factor L = "<<Lstore->nnz<<endl;
	cout<<"No of nonzeros in factor U = "<<Ustore->nnz<<endl;
	cout<<"No of nonzeros in L+U     = "<<(Lstore->nnz+Ustore->nnz)<<endl;
	dQuerySpace(&L, &U, &mem_usage);
	cout<<"Memory used for factorisation in Mo "
	    <<mem_usage.total_needed/(1024*1024)<<endl;
      }
    
    Acsr.Nullify();
  }
  
  
  //! resolution of linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<double>::Solve(Vector<double,VectFull,Allocator2>& x)
  {
    trans_t trans = NOTRANS;
    int nb_rhs = 1, info;
    dCreate_Dense_Matrix(&B, x.GetM(), nb_rhs,
			 x.GetData(), x.GetM(), SLU_DN, SLU_D, SLU_GE);
    
    SuperLUStat_t stat;
    StatInit(&stat);
    dgstrs(trans, &L, &U, perm_r.GetData(),
	   perm_c.GetData(), &B, &stat, &info);
  }
  
  
  //! factorization of matrix in complex double precision using SuperLU
  template<class Prop, class Storage,class Allocator>
  void MatrixSuperLU<complex<double> >::
  FactorizeMatrix(Matrix<complex<double>,Prop,Storage,Allocator> & mat,
		  bool keep_matrix)
  {
    // clearing previous factorization
    Clear();
    
    // conversion in CSR format
    n = mat.GetN();
    Matrix<complex<double>,General,ColSparse> Acsr;
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();
    
    // we get renumbering vectors perm_r and perm_c
    int nnz = Acsr.GetDataSize();
    zCreate_CompCol_Matrix(&A, n, n, nnz,
			   reinterpret_cast<doublecomplex*>(Acsr.GetData()),
			   Acsr.GetInd(), Acsr.GetPtr(),
			   SLU_NC, SLU_Z, SLU_GE);
    
    perm_r.Reallocate(n);
    perm_c.Reallocate(n);
    
    int nb_rhs = 0, info;
    zCreate_Dense_Matrix(&B, n, nb_rhs, NULL, n, SLU_DN, SLU_Z, SLU_GE);

    zgssv(&options, &A, perm_c.GetData(), perm_r.GetData(),
	  &L, &U, &B, &stat, &info);
    
    if ((info==0)&&(display_info))
      {
	mem_usage_t mem_usage;
	Lstore = (SCformat *) L.Store;
	Ustore = (NCformat *) U.Store;
	cout<<"No of nonzeros in factor L = "<<Lstore->nnz<<endl;
	cout<<"No of nonzeros in factor U = "<<Ustore->nnz<<endl;
	cout<<"No of nonzeros in L+U     = "<<(Lstore->nnz+Ustore->nnz)<<endl;
	zQuerySpace(&L, &U, &mem_usage);
	cout<<"Memory used for factorisation in Mo "
	    <<mem_usage.total_needed/1e6<<endl;
      }
    
    Acsr.Nullify();
  }
  
  
  //! resolution of linear system A x = b
  template<class Allocator2>
  void MatrixSuperLU<complex<double> >::
  Solve(Vector<complex<double>,VectFull,Allocator2>& x)
  {
    trans_t trans = NOTRANS;
    int nb_rhs = 1, info;
    zCreate_Dense_Matrix(&B, x.GetM(), nb_rhs,
			 reinterpret_cast<doublecomplex*>(x.GetData()),
			 x.GetM(), SLU_DN, SLU_Z, SLU_GE);
    
    zgstrs (trans, &L, &U, perm_r.GetData(),
	    perm_c.GetData(), &B, &stat, &info);
  }
  
  
  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T,Prop,Storage,Allocator>& A, MatrixSuperLU<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  
  
  template<class T, class Allocator>
  void SolveLU(MatrixSuperLU<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }
  
}

#define SELDON_FILE_SUPERLU_CXX
#endif
