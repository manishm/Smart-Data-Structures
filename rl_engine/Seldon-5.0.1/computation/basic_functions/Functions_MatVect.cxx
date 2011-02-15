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


#ifndef SELDON_FILE_FUNCTIONS_MATVECT_CXX

/*
  Functions defined in this file:

  M*X -> Y
  Mlt(M, X, Y)

  alpha.M*X -> Y
  Mlt(alpha, M, X, Y)

  alpha.M*X + beta.Y -> Y
  MltAdd(alpha, M, X, beta, Y)
*/

namespace Seldon
{


  /////////
  // MLT //


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	   const Vector<T1, Storage1, Allocator1>& X,
	   Vector<T2, Storage2, Allocator2>& Y)
  {
    Y.Fill(T2(0));
    MltAdd(T2(1), M, X, T2(0), Y);
  }


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3, class Storage3, class Allocator3>
  void Mlt(const T0 alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	   const Vector<T2, Storage2, Allocator2>& X,
	   Vector<T3, Storage3, Allocator3>& Y)
  {
    Y.Fill(T2(0));
    MltAdd(alpha, M, X, T3(0), Y);
  }


  // MLT //
  /////////



  ////////////
  // MltAdd //


  /*** Sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    T4 zero(0);
    T4 temp;

    int* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSparse, Allocator1>::pointer
      data = M.GetData();

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = ptr[i]; j < ptr[i+1]; j++)
	  temp += data[j] * X(ind[j]);
	Y(i) += alpha * temp;
      }
  }


  /*** Complex sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    complex<T1> zero(0);
    complex<T1> temp;

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	Y(i) += alpha * temp;
      }

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T1>(T1(0), imag_data[j]) * X(imag_ind[j]);
	Y(i) += alpha * temp;
      }
  }


  /*** Symmetric sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i, j;
    T4 zero(0);
    T4 temp;

    int* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSymSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += data[j] * X(ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i < ma-1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += alpha * data[j] * X(i);
  }


  /*** Symmetric complex sparse matrices ***/


  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i, j;
    complex<T1> zero(0);
    complex<T1> temp;

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i<ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i<ma-1; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	if (real_ind[j] != i)
	  Y(real_ind[j]) += alpha * real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T1>(T1(0), imag_data[j]) * X(imag_ind[j]);
	Y(i) += alpha * temp;
      }
    for (i = 0; i<ma-1; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	if (imag_ind[j] != i)
	  Y(imag_ind[j]) += alpha * complex<T1>(T1(0), imag_data[j]) * X(i);
  }


  /*** Sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<T1, Prop1, RowSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	Y(ind[j]) += alpha * data[j] * X(i);
  }


  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, RowSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<complex<T1>, Prop1, RowSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	Y(ind[j]) += alpha * conj(data[j]) * X(i);
  }


  /*** Complex sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	Y(real_ind[j]) += alpha * real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	Y(imag_ind[j]) += alpha * complex<T1>(T1(0), imag_data[j]) * X(i);
  }


  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, Prop1, RowComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int i, j;

    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	Y(real_ind[j]) += alpha * real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	Y(imag_ind[j]) += alpha * complex<T1>(T1(0), - imag_data[j]) * X(i);
  }


  /*** Symmetric sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<complex<T1>, Prop1, RowSymSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i, j;
    complex<T1> zero(0);
    complex<T1> temp;

    int* ptr = M.GetPtr();
    int* ind = M.GetInd();
    typename Matrix<complex<T1>, Prop1, RowSymSparse, Allocator1>::pointer
      data = M.GetData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = ptr[i]; j < ptr[i + 1]; j++)
	  temp += conj(data[j]) * X(ind[j]);
	Y(i) += temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = ptr[i]; j < ptr[i + 1]; j++)
	if (ind[j] != i)
	  Y(ind[j]) += conj(data[j]) * X(i);
  }


  /*** Symmetric complex sparse matrices, *Trans ***/


  // NoTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonNoTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // Trans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    MltAdd(alpha, M, X, beta, Y);
  }


  // ConjTrans.
  template <class T0,
	    class T1, class Prop1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const class_SeldonConjTrans& Trans,
	      const Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta, Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, SeldonConjTrans, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    int i, j;
    complex<T1> zero(0);
    complex<T1> temp;

    int* real_ptr = M.GetRealPtr();
    int* imag_ptr = M.GetImagPtr();
    int* real_ind = M.GetRealInd();
    int* imag_ind = M.GetImagInd();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      real_data = M.GetRealData();
    typename Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>::pointer
      imag_data = M.GetImagData();

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	  temp += real_data[j] * X(real_ind[j]);
	Y(i) += temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = real_ptr[i]; j < real_ptr[i + 1]; j++)
	if (real_ind[j] != i)
	  Y(real_ind[j]) += real_data[j] * X(i);

    for (i = 0; i < ma; i++)
      {
	temp = zero;
	for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	  temp += complex<T1>(T1(0), - imag_data[j]) * X(imag_ind[j]);
	Y(i) += temp;
      }
    for (i = 0; i < ma - 1; i++)
      for (j = imag_ptr[i]; j < imag_ptr[i + 1]; j++)
	if (imag_ind[j] != i)
	  Y(imag_ind[j]) += complex<T1>(T1(0), - imag_data[j]) * X(i);
  }


  // MltAdd //
  ////////////



  ////////////
  // MltAdd //


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y)
  {
    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(M, X, Y, "MltAdd(alpha, M, X, beta, Y)");
#endif

    Mlt(beta, Y);

    T4 zero(0);
    T4 temp;
    T4 alpha_(alpha);

    for (int i = 0; i < ma; i++)
      {
	temp = zero;
	for (int j = 0; j < na; j++)
	  temp += M(i, j) * X(j);
	Y(i) += alpha_ * temp;
      }
  }
  
  
  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const SeldonTranspose& Trans,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& M,
	      const Vector<T2, Storage2, Allocator2>& X,
	      const T3 beta,
	      Vector<T4, Storage4, Allocator4>& Y)
  {
    if (Trans.NoTrans())
      {
        MltAdd(alpha, M, X, beta, Y);
        return;
      }
    else if (Trans.ConjTrans())
      throw WrongArgument("MltAdd(alpha, trans, M, X, beta, Y)",
                          "Complex conjugation not supported.");

    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(Trans, M, X, Y, "MltAdd(alpha, trans, M, X, beta, Y)");
#endif

    if (beta == T3(0))
      Y.Fill(T4(0));
    else
      Mlt(beta, Y);

    T4 zero(0);
    T4 temp;
    T4 alpha_(alpha);

    for (int i = 0; i < na; i++)
      {
	temp = zero;
	for (int j = 0; j < ma; j++)
	  temp += M(j, i) * X(j);
	Y(i) += alpha_ * temp;
      }
  }
  

  // MltAdd //
  ////////////



  ///////////
  // Gauss //


  // Solve X = M*Y with Gauss method.
  // Warning: M is modified. The results are stored in X.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  inline void Gauss(Matrix<T0, Prop0, Storage0, Allocator0>& M,
		    Vector<T1, Storage1, Allocator1>& X)
  {
    int i, j, k;
    T1 r, S;

    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_BOUNDS
    if (na != ma)
      throw WrongDim("Gauss(M, X)",
		     "The matrix must be squared.");

    CheckDim(M, X, "Gauss(M, X)");
#endif

    for (k = 0; k < ma - 1; k++)
      for (i = k + 1; i < ma; i++)
	{
	  r = M(i, k) / M(k, k);
	  for (j = k + 1; j < ma; j++)
	    M(i, j) -= r * M(k, j);
	  X(i) -= r *= X(k);
	}

    X(ma - 1) = X(ma - 1) / M(ma - 1, ma - 1);
    for (k = ma - 2; k > -1; k--)
      {
	S = X(k);
	for (j = k + 1; j < ma; j++)
	  S -= M(k, j) * X(j);
	X(k) = S / M(k, k);
      }
  }


  // Gauss //
  ///////////



  ////////////////////
  // Gauss - Seidel //


  // Solve X = M*Y with Gauss-Seidel method.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  inline void GaussSeidel(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
			  const Vector<T1, Storage1, Allocator1>& X,
			  Vector<T2, Storage2, Allocator2>& Y,
			  int iter)
  {
    int i, j, k;
    T1 temp;

    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_BOUNDS
    if (na != ma)
      throw WrongDim("GaussSeidel(M, X, Y, iter)",
		     "The matrix must be squared.");

    CheckDim(M, X, Y, "GaussSeidel(M, X, Y, iter)");
#endif
    
    for (i = 0; i < iter; i++)
      for (j = 0; j < na; j++)
	{
	  temp = 0;
	  for (k = 0; k < j; k++)
	    temp -= M(j, k) * Y(k);
	  for (k = j + 1; k < na; k++)
	    temp -= M(j, k) * Y(k);
	  Y(j) = (X(j) + temp) / M(j, j);
	}
  }


  // Gauss-Seidel //
  //////////////////



  ///////////////////
  // S.O.R. method //


  // Solve X = M*Y with S.O.R. method.
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2,
	    class T3>
  inline void SOR(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		  const Vector<T1, Storage1, Allocator1>& X,
		  Vector<T2, Storage2, Allocator2>& Y,
		  T3 omega,
		  int iter)
  {
    int i, j, k;
    T1 temp;

    int ma = M.GetM();
    int na = M.GetN();

#ifdef SELDON_CHECK_BOUNDS
    if (na != ma)
      throw WrongDim("SOR(M, X, Y, omega, iter)",
		     "The matrix must be squared.");

    CheckDim(M, X, Y, "SOR(M, X, Y, omega, iter)");
#endif
    
    for (i = 0; i < iter; i++)
      for (j = 0; j < na; j++)
	{
	  temp = 0;
	  for (k = 0; k < j; k++)
	    temp -= M(j, k) * Y(k);
	  for (k = j + 1; k < na; k++)
	    temp -= M(j, k) * Y(k);
	  Y(j) = (T3(1) - omega) * Y(j) + omega * (X(j) + temp) / M(j, j);
	}
  }


  // Gauss-Seidel //
  //////////////////



  /////////////
  // SolveLU //


  // Solves M.X = Y where A has been decomposed in a LU form.
  // Y is overwritten (Y <- X).
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SolveLU(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	       Vector<T1, Storage1, Allocator1>& Y)
  {
    int i, k;
    T1 temp;

    int ma = M.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = M.GetN();
    if (na != ma)
      throw WrongDim("SolveLU(M, Y)",
		     "The matrix must be squared.");

    CheckDim(M, Y, "SolveLU(M, Y)");
#endif

    // Forward substitution.
    for (i = 0; i < ma; i++)
      {
	temp = 0;
	for (k = 0; k < i; k++)
	  temp += M(i, k) * Y(k);
	Y(i) = (Y(i) - temp) / M(i, i);
      }
    // Back substitution.
    for (i = ma - 2; i > -1; i--)
      {
	temp = 0;
	for (k = i + 1; k < ma; k++)
	  temp += M(i, k) * Y(k);
	Y(i) -= temp;
      }
  }


  // SolveLU //
  /////////////



  //////////////
  // CHECKDIM //


  //! Checks the compatibility of the dimensions.
  /*! Checks that M.X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "")
  {
    if (X.GetLength() != M.GetN() || Y.GetLength() != M.GetM())
      throw WrongDim(function, string("Operation M.X + Y -> Y not permitted:")
		     + string("\n     M (") + to_str(&M) + string(") is a ")
		     + to_str(M.GetM()) + string(" x ") + to_str(M.GetN())
		     + string(" matrix;\n     X (") + to_str(&X)
		     + string(") is vector of length ")
		     + to_str(X.GetLength()) + string(";\n     Y (")
		     + to_str(&Y) + string(") is vector of length ")
		     + to_str(Y.GetLength()) + string("."));
  }


  //! Checks the compatibility of the dimensions.
  /*! Checks that M.X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param status status of the matrix M, e.g. transposed.
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "M.X + Y -> Y".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const SeldonTranspose& status,
		const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "", string op = "M.X + Y -> Y")
  {
    if (op == "M.X + Y -> Y")
      if (status.Trans())
	op = string("Operation M'.X + Y -> Y not permitted:");
      else if (status.ConjTrans())
	op = string("Operation M*.X + Y -> Y not permitted:");
      else
	op = string("Operation M.X + Y -> Y not permitted:");
    else
      op = string("Operation ") + op + string(" not permitted:");
    if (X.GetLength() != M.GetN(status) || Y.GetLength() != M.GetM(status))
      throw WrongDim(function, op + string("\n     M (") + to_str(&M)
		     + string(") is a ") + to_str(M.GetM()) + string(" x ")
		     + to_str(M.GetN()) + string(" matrix;\n     X (")
		     + to_str(&X) + string(") is vector of length ")
		     + to_str(X.GetLength()) + string(";\n     Y (")
		     + to_str(&Y) + string(") is vector of length ")
		     + to_str(Y.GetLength()) + string("."));
  }


#ifdef SELDON_WITH_CBLAS
  //! Checks the compatibility of the dimensions.
  /*! Checks that M.X + Y -> Y is possible according to the dimensions of
    the matrix M and the vectors X and Y. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param trans status of the matrix M, e.g. transposed.
    \param M matrix.
    \param X vector.
    \param Y vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed on the vectors.
    Default: "M.X + Y -> Y".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1,
	    class T2, class Storage2, class Allocator2>
  void CheckDim(const enum CBLAS_TRANSPOSE trans,
		const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		const Vector<T2, Storage2, Allocator2>& Y,
		string function = "", string op = "M.X + Y -> Y")
  {
    SeldonTranspose status(trans);
    if (op == "M.X + Y -> Y")
      if (status.Trans())
	op = string("Operation M'.X + Y -> Y not permitted:");
      else if (status.ConjTrans())
	op = string("Operation M*.X + Y -> Y not permitted:");
      else
	op = string("Operation M.X + Y -> Y not permitted:");
    else
      op = string("Operation ") + op + string(" not permitted:");
    if (X.GetLength() != M.GetN(status) || Y.GetLength() != M.GetM(status))
      throw WrongDim(function, op + string("\n     M (") + to_str(&M)
		     + string(") is a ") + to_str(M.GetM()) + string(" x ")
		     + to_str(M.GetN()) + string(" matrix;\n     X (")
		     + to_str(&X) + string(") is vector of length ")
		     + to_str(X.GetLength()) + string(";\n     Y (")
		     + to_str(&Y) + string(") is vector of length ")
		     + to_str(Y.GetLength()) + string("."));
  }
#endif


  //! Checks the compatibility of the dimensions.
  /*! Checks that M.X is possible according to the dimensions of
    the matrix M and the vector X. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param M matrix.
    \param X vector.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
    \param op (optional) operation to be performed. Default: "M.X".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
		const Vector<T1, Storage1, Allocator1>& X,
		string function = "", string op = "M.X")
  {
    if (X.GetLength() != M.GetN())
      throw WrongDim(function, string("Operation ") + op + " not permitted:"
		     + string("\n     M (") + to_str(&M) + string(") is a ")
		     + to_str(M.GetM()) + string(" x ") + to_str(M.GetN())
		     + string(" matrix;\n     X (") + to_str(&X)
		     + string(") is vector of length ")
		     + to_str(X.GetLength()) + string("."));
  }


  // CHECKDIM //
  //////////////


}  // namespace Seldon.

#define SELDON_FUNCTIONS_MATVECT_CXX
#endif
