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


#ifndef SELDON_FILE_UMFPACK_CXX

#include "UmfPack.hxx"

namespace Seldon
{
  //! constructor
  template<class T>
  MatrixUmfPack_Base<T>::MatrixUmfPack_Base()
  {
    // first we clear arrays
    Clear();
    
    // allocation of arrays Control and Info
    Control.Reallocate(UMFPACK_CONTROL);
    Info.Reallocate(UMFPACK_INFO);
    
    display_info = true;
  }
  
  template<class T>
  MatrixUmfPack_Base<T>::~MatrixUmfPack_Base()
  {
    Clear();
  }
  
  template<class T>
  void MatrixUmfPack_Base<T>::Clear()
  {
    n = 0;
    Control.Clear(); Info.Clear();
  }
  
  //! no message will be displayed by UmfPack
  template<class T>
  void MatrixUmfPack_Base<T>::HideMessages()
  {
    display_info = false;
    Control(UMFPACK_PRL) = 0;
  }
  
  //! normal amount of message displayed by UmfPack
  template<class T>
  void MatrixUmfPack_Base<T>::ShowMessages()
  {
    display_info = true;
    Control(UMFPACK_PRL) = 2;
  }
  
  //! constructor
  MatrixUmfPack<double>::MatrixUmfPack() : MatrixUmfPack_Base<double>()
  {
    umfpack_di_defaults(this->Control.GetData());
    this->ShowMessages();
  }
  
  //! constructor
  MatrixUmfPack<complex<double> >::MatrixUmfPack()
    : MatrixUmfPack_Base<complex<double> >()
  {
    umfpack_zi_defaults(this->Control.GetData());
    this->ShowMessages();
  }
  
  //! destructor
  MatrixUmfPack<double>::~MatrixUmfPack()
  {
    if (this->n > 0)
      {
	// memory used by LU factorization is released
	umfpack_di_free_numeric(&this->Numeric) ;
	this->n = 0;
      }
  }
  
  //! destructor
  MatrixUmfPack<complex<double> >::~MatrixUmfPack()
  {
    if (this->n > 0)
      {
	// memory used by LU factorization is released
	umfpack_zi_free_numeric(&this->Numeric) ;
	this->n = 0;
      }
  }
  
  //! factorization of a matrix in double precision
  template<class Prop, class Storage, class Allocator>
  void MatrixUmfPack<double>::
  FactorizeMatrix(Matrix<double,Prop,Storage,Allocator> & mat,
		  bool keep_matrix)
  {
    // we clear previous factorization
    Clear();
    
    // conversion in CSR Format
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();
    
    // factorization with UmfPack
    this->n = mat.GetM();
    umfpack_di_symbolic(Acsr.GetM(), Acsr.GetN(), Acsr.GetPtr(),
			Acsr.GetInd(), Acsr.GetData(), &this->Symbolic,
			this->Control.GetData(), this->Info.GetData());
    
    int status =
      umfpack_di_numeric(Acsr.GetPtr(), Acsr.GetInd(), Acsr.GetData(),
			 this->Symbolic, &this->Numeric,
			 this->Control.GetData(), this->Info.GetData());
    
    // memory for numbering scheme is released
    umfpack_di_free_symbolic(&this->Symbolic) ;
    
    // we display informations about the performed operation
    if (this->display_info)
      {
	umfpack_di_report_status(this->Control.GetData(), status);
	umfpack_di_report_info(this->Control.GetData(),this->Info.GetData());
      }
  }
  
  template<class Allocator2>
  void MatrixUmfPack<double>::Solve(Vector<double,VectFull,Allocator2>& x)
  {
    // we call UmfPack
    Vector<double,VectFull,Allocator2> b(x);
    int status
      = umfpack_di_solve(UMFPACK_A, Acsr.GetPtr(), Acsr.GetInd(),
			 Acsr.GetData(), x.GetData(), b.GetData(),
			 this->Numeric, this->Control.GetData(),
			 this->Info.GetData());
    
    // we display informations about the performed operation
    if (this->display_info)
      umfpack_di_report_status(this->Control.GetData(), status);
  }
  
  //! LU factorization using UmfPack in double complex precision
  template<class Prop, class Storage,class Allocator>
  void MatrixUmfPack<complex<double> >::
  FactorizeMatrix(Matrix<complex<double>,Prop,Storage,Allocator> & mat,
		  bool keep_matrix)
  {
    Ptr.Clear(); Ind.Clear(); ValuesImag.Clear(); ValuesReal.Clear();
    
    this->n = mat.GetM();
    // conversion in CSR format
    Matrix<complex<double>, General, ColSparse, Allocator> Acsr;
    Copy(mat, Acsr);
    if (!keep_matrix)
      mat.Clear();
    
    int nnz = Acsr.GetDataSize();
    complex<double>* data = Acsr.GetData();
    int* ptr_ = Acsr.GetPtr(); int* ind_ = Acsr.GetInd();
    ValuesReal.Reallocate(nnz);
    ValuesImag.Reallocate(nnz);
    Ptr.Reallocate(this->n+1); Ind.Reallocate(nnz);
    
    for (int i = 0; i < nnz; i++)
      {
	ValuesReal(i) = real(data[i]);
	ValuesImag(i) = imag(data[i]);
	Ind(i) = ind_[i];
      }
    
    for (int i = 0; i <= this->n; i++)
      Ptr(i) = ptr_[i];
    
    // we clear intermediary matrix Acsr
    Acsr.Clear();
    
    // we call UmfPack
    umfpack_zi_symbolic(this->n, this->n, Ptr.GetData(), Ind.GetData(),
			ValuesReal.GetData(),ValuesImag.GetData(),
			&this->Symbolic, this->Control.GetData(),
			this->Info.GetData());
    
    int status
      = umfpack_zi_numeric(Ptr.GetData(), Ind.GetData(),
			   ValuesReal.GetData(), ValuesImag.GetData(),
			   this->Symbolic, &this->Numeric,
			   this->Control.GetData(), this->Info.GetData());
    
    umfpack_zi_free_symbolic(&this->Symbolic) ;
    if (this->display_info)
      {
	umfpack_zi_report_status(this->Control.GetData(),status);
	umfpack_zi_report_info(this->Control.GetData(),this->Info.GetData());
      }
  }
  
  //! solves linear system in complex double precision using UmfPack
  template<class Allocator2>
  void MatrixUmfPack<complex<double> >::
  Solve(Vector<complex<double>,VectFull,Allocator2>& x)
  {
    int m = x.GetM();
    // creation of vectors
    Vector<double> b_real(m),b_imag(m);
    for (int i = 0; i < m; i++)
      {
	b_real(i) = real(x(i));
	b_imag(i) = imag(x(i));
      }
    Vector<double> x_real(m),x_imag(m);
    x_real.Zero();
    x_imag.Zero();
    int status
      = umfpack_zi_solve(UMFPACK_A, Ptr.GetData(), Ind.GetData(),
			 ValuesReal.GetData(), ValuesImag.GetData(),
			 x_real.GetData(), x_imag.GetData(),
			 b_real.GetData(), b_imag.GetData(),
			 this->Numeric,
			 this->Control.GetData(), this->Info.GetData());
    b_real.Clear();
    b_imag.Clear();
    
    for (int i = 0; i < m; i++)
      x(i) = complex<double>(x_real(i), x_imag(i));
    
    if (this->display_info)
      umfpack_zi_report_status(this->Control.GetData(), status);
  }
  
  template<class T, class Prop, class Storage, class Allocator>
  void GetLU(Matrix<T,Prop,Storage,Allocator>& A, MatrixUmfPack<T>& mat_lu,
	     bool keep_matrix = false)
  {
    mat_lu.FactorizeMatrix(A, keep_matrix);
  }
  
  template<class T, class Allocator>
  void SolveLU(MatrixUmfPack<T>& mat_lu, Vector<T, VectFull, Allocator>& x)
  {
    mat_lu.Solve(x);
  }
  
}
#define SELDON_FILE_UMFPACK_CXX
#endif
