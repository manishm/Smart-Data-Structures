// Copyright (C) 2001-2009 Vivien Mallet
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


#ifndef SELDON_FILE_FUNCTIONS_MATRIX_CXX

/*
  Function defined in this file:

  alpha.A -> A
  Mlt(alpha, A)

  A*B -> C
  Mlt(A, B, C)

  alpha.A*B -> C
  Mlt(alpha, A, B, C)

  alpha.A*B + beta.C -> C
  MltAdd(alpha, A, B, beta, C)
  
  alpha*A + B -> B
  Add(alpha, A, B)
  
  LU factorization of matrix A without pivoting.
  GetLU(A)
  
  Highest absolute value of A.
  MaxAbs(A)
  
  1-norm of matrix A.
  Norm1(A)
  
  infinity norm of matrix A.
  NormInf(A)
  
  Transpose(A)
*/

namespace Seldon
{


  /////////
  // MLT //


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void Mlt(const T0 alpha,
	   Matrix<T1, Prop1, Storage1, Allocator1>& A)  throw()
  {
    T1 alpha_ = alpha;

    typename Matrix<T1, Prop1, Storage1, Allocator1>::pointer
      data = A.GetData();

    for (int i = 0; i < A.GetDataSize(); i++)
      data[i] = alpha_ * data[i];
  }


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2,
	    class T3, class Prop3, class Storage3, class Allocator3>
  void Mlt(const T0 alpha,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	   const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	   Matrix<T3, Prop3, Storage3, Allocator3>& C)
  {
    C.Fill(T3(0));
    MltAdd(alpha, A, B, T3(0), C);
  }


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void Mlt(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
	   const Matrix<T1, Prop1, Storage1, Allocator1>& B,
	   Matrix<T2, Prop2, Storage2, Allocator2>& C)
  {
    C.Fill(T2(0));
    MltAdd(T0(1), A, B, T2(0), C);
  }


  // MLT //
  /////////



  ////////////
  // MLTADD //


  template <class T0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2,
	    class T3,
	    class T4, class Prop4, class Storage4, class Allocator4>
  void MltAdd(const T0 alpha,
	      const Matrix<T1, Prop1, Storage1, Allocator1>& A,
	      const Matrix<T2, Prop2, Storage2, Allocator2>& B,
	      const T3 beta,
	      Matrix<T4, Prop4, Storage4, Allocator4>& C)
  {
    int na = A.GetN();
    int mc = C.GetM();
    int nc = C.GetN();

#ifdef SELDON_CHECK_BOUNDS
    CheckDim(A, B, C, "MltAdd(alpha, A, B, beta, C)");
#endif

    T4 temp;
    T4 alpha_(alpha);
    T4 beta_(beta);
    
    if (beta_ != T4(0))
      for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
	for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
	  C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j))
	    *= beta_;
    else
      for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
	for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
	  C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j)) = T4(0);

    for (int i = 0; i < Storage4::GetFirst(mc, nc); i++)
      for (int j = 0; j < Storage4::GetSecond(mc, nc); j++)
	{
	  temp = T4(0);
	  for (int k = 0; k < na; k++)
	    temp += A(Storage4::GetFirst(i, j), k)
	      * B(k, Storage4::GetSecond(i, j));
	  C(Storage4::GetFirst(i, j), Storage4::GetSecond(i, j))
	    += alpha_ * temp;
	}
  }


  // MLTADD //
  ////////////
  

  
  /////////
  // ADD //
  

  template<class T0, class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, General, Storage1, Allocator1>& A,
	   Matrix<T2, General, Storage2, Allocator2>& B)
  {
    int i, j;
    for (i = 0; i < A.GetM(); i++)
      for (j = 0; j < A.GetN(); j++)
	B(i, j) += alpha * A(i, j);
  }
  

  template<class T0, class T1, class Storage1, class Allocator1,
	   class T2, class Storage2, class Allocator2>
  void Add(const T0& alpha,
	   const Matrix<T1, Symmetric, Storage1, Allocator1>& A,
	   Matrix<T2, Symmetric, Storage2, Allocator2>& B)
  {
    int i, j;
    for (i = 0; i < A.GetM(); i++)
      for (j = i; j < A.GetN(); j++)
	B(i, j) += alpha * A(i, j);
  }

  
  // ADD //
  /////////

  
  
  ///////////
  // GETLU //


  // Returns the LU decomposition of A = LU (in A)
  // where L diagonal elements are set to unit value.
  template <class T0, class Prop0, class Storage0, class Allocator0>
  void GetLU(Matrix<T0, Prop0, Storage0, Allocator0>& A)
  {
    int i, p, q, k;
    T0 temp;

    int ma = A.GetM();

#ifdef SELDON_CHECK_BOUNDS
    int na = A.GetN();
    if (na != ma)
      throw WrongDim("GetLU(A)", "The matrix must be squared.");
#endif

    for (i = 0; i < ma; i++)
      {
	for (p = i; p < ma; p++)
	  {
	    temp = 0;
	    for (k = 0; k < i; k++)
	      temp += A(p, k) * A(k, i);
	    A(p, i) -= temp;
	  }
	for (q = i+1; q < ma; q++)
	  {
	    temp = 0;
	    for (k = 0; k < i; k++)
	      temp += A(i, k) * A(k, q);
	    A(i, q) = (A(i,q ) - temp) / A(i, i);
	  }
      }
  }
  
  
  // GETLU //
  ///////////



  //////////////
  // CHECKDIM //


  //! Checks the compatibility of the dimensions.
  /*! Checks that A.B + C -> C is possible according to the dimensions of
    the matrices A, B and C. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param A matrix.
    \param B matrix.
    \param C matrix.
    \param function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function = "")
  {
    if (B.GetM() != A.GetN() || C.GetM() != A.GetM() || B.GetN() != C.GetN())
      throw WrongDim(function, string("Operation A.B + C -> C not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
  }


#ifdef SELDON_WITH_CBLAS
  //! Checks the compatibility of the dimensions.
  /*! Checks that A.B + C -> C or B.A + C -> C is possible according to the
    dimensions of the matrices A, B and C. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param side side by which A is multiplied by B.
    \param A matrix.
    \param B matrix.
    \param C matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const enum CBLAS_SIDE side,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function = "")
  {
    if ( SeldonSide(side).Left() &&
	 (B.GetM() != A.GetN() || C.GetM() != A.GetM()
	  || B.GetN() != C.GetN()) )
      throw WrongDim(function, string("Operation A.B + C -> C not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
    else if ( SeldonSide(side).Right() &&
	      (B.GetN() != A.GetM() || C.GetM() != B.GetM()
	       || A.GetN() != C.GetN()) )
      throw WrongDim(function, string("Operation B.A + C -> C not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
  }
#endif


#ifdef SELDON_WITH_CBLAS
  //! Checks the compatibility of the dimensions.
  /*! Checks that A.B + C -> C is possible according to the dimensions of
    the matrices A, B and C. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param TransA status of A, e.g. transposed.
    \param A matrix.
    \param TransB status of B, e.g. transposed.
    \param B matrix.
    \param C matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1,
	    class T2, class Prop2, class Storage2, class Allocator2>
  void CheckDim(const enum CBLAS_TRANSPOSE TransA,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const enum CBLAS_TRANSPOSE TransB,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		const Matrix<T2, Prop2, Storage2, Allocator2>& C,
		string function = "")
  {
    SeldonTranspose status_A(TransA);
    SeldonTranspose status_B(TransB);
    string op;
    if (status_A.Trans())
      op = string("A'");
    else if (status_A.ConjTrans())
      op = string("A*");
    else
      op = string("A");
    if (status_B.Trans())
      op += string(".B' + C");
    else if (status_B.ConjTrans())
      op += string(".B* + C");
    else
      op += string(".B + C");
    op = string("Operation ") + op + string(" not permitted:");
    if (B.GetM(status_B) != A.GetN(status_A) || C.GetM() != A.GetM(status_A)
	|| B.GetN(status_B) != C.GetN())
      throw WrongDim(function, op
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix;\n     C (")
		     + to_str(&C) + string(") is a ") + to_str(C.GetM())
		     + string(" x ") + to_str(C.GetN()) + string(" matrix."));
  }
#endif


#ifdef SELDON_WITH_CBLAS
  //! Checks the compatibility of the dimensions.
  /*! Checks that A.B or B.A is possible according to the dimensions of
    the matrices A and B. If the dimensions are incompatible,
    an exception is raised (a WrongDim object is thrown).
    \param side side by which A is multiplied by B.
    \param A matrix.
    \param B matrix.
    \function (optional) function in which the compatibility is checked.
    Default: "".
  */
  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Prop1, class Storage1, class Allocator1>
  void CheckDim(const enum CBLAS_SIDE side,
		const Matrix<T0, Prop0, Storage0, Allocator0>& A,
		const Matrix<T1, Prop1, Storage1, Allocator1>& B,
		string function = "")
  {
    if (SeldonSide(side).Left() && B.GetM() != A.GetN())
      throw WrongDim(function, string("Operation A.B not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix."));
    else if (SeldonSide(side).Right() && B.GetN() != A.GetM())
      throw WrongDim(function, string("Operation B.A not permitted:")
		     + string("\n     A (") + to_str(&A) + string(") is a ")
		     + to_str(A.GetM()) + string(" x ") + to_str(A.GetN())
		     + string(" matrix;\n     B (") + to_str(&B)
		     + string(") is a ") + to_str(B.GetM())  + string(" x ")
		     + to_str(B.GetN()) + string(" matrix."));
  }
#endif


  // CHECKDIM //
  //////////////


  ///////////
  // NORMS //
  
  
  //! Returns the maximum (in absolute value) of a matrix.
  /*!
    \param A matrix.
    \return The maximum (in absolute value) of matrix A.
  */
  template <class T, class Prop, class Storage, class Allocator>
  T MaxAbs(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    T res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
	res = max(res, abs(A(i, j)) );

    return res;
  }
  
  
  //! Returns the 1-norm of a matrix.
  /*!
    \param A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T Norm1(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int j = 0; j < A.GetN(); j++)
      {
	sum = T(0);
	for (int i = 0; i < A.GetM(); i++)
	  sum += abs( A(i, j) );
	
	res = max(res, sum);
      }
    
    return res;
  }
  
  
  //! Returns the infinity-norm of a matrix.
  /*!
    \param A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T NormInf(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
	sum = T(0);
	for (int j = 0; j < A.GetN(); j++)
	  sum += abs( A(i, j) );
	
	res = max(res, sum);
      }

    return res;
  }
  
  
  //! Returns the maximum (in modulus) of a matrix.
  /*!
    \param A matrix.
    \return The maximum (in modulus) of matrix A.
  */
  template <class T, class Prop, class Storage, class Allocator>
  T MaxAbs(const Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    T res(0);
    for (int i = 0; i < A.GetM(); i++)
      for (int j = 0; j < A.GetN(); j++)
	{
	  res = max(res, abs(A(i, j)) );
	}
    
    return res;
  }
  
  
  //! Returns the 1-norm of a matrix.
  /*!
    \param A matrix.
    \return \f$ max_j \sum_i |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T Norm1(const Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int j = 0; j < A.GetN(); j++)
      {
	sum = T(0);
	for (int i = 0; i < A.GetM(); i++)
	  sum += abs( A(i, j) );
	
	res = max(res, sum);
      }
    
    return res;
  }
  
  
  //! Returns the infinity-norm of a matrix.
  /*!
    \param A matrix.
    \return \f$ max_i \sum_j |A_{ij}| \f$
  */
  template <class T, class Prop, class Storage, class Allocator>
  T NormInf(const Matrix<complex<T>, Prop, Storage, Allocator>& A)
  {
    T res(0), sum;
    for (int i = 0; i < A.GetM(); i++)
      {
	sum = T(0);
	for (int j = 0; j < A.GetN(); j++)
	  sum += abs( A(i, j) );
	
	res = max(res, sum);
      }

    return res;
  }
  
  
  // NORMS //
  ///////////
  
  
  
  ///////////////
  // TRANSPOSE //
  
  
  //! Matrix transposition.
  template<class T, class Prop, class Storage, class Allocator>
  void Transpose(Matrix<T, Prop, Storage, Allocator>& A)
  {
    int m = A.GetM();
    int n = A.GetN();
    
    if (m == n)
      {
	T tmp;
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < i; j++)
	    {
	      tmp = A(i,j);
	      A(i,j) = A(j,i);
	      A(j,i) = tmp;
	    }
      }
    else
      {
	Matrix<T, Prop, Storage, Allocator> B;
	B = A;
	A.Reallocate(n,m);
	for (int i = 0; i < m; i++)
	  for (int j = 0; j < n; j++)
	    A(j,i) = B(i,j);
      }
  }
  
  
  //! Matrix transposition and conjugation.
  template<class T, class Prop, class Storage, class Allocator>
  void TransposeConj(Matrix<T, Prop, Storage, Allocator>& A)
  {
    int i, j;

    int m = A.GetM();
    int n = A.GetN();
    
    if (m == n)
      {
	T tmp;
	for (i = 0; i < m; i++)
	  for (j = 0; j < i; j++)
	    {
	      tmp = A(i, j);
	      A(i, j) = conj(A(j, i));
	      A(j, i) = conj(tmp);
	    }
      }
    else
      {
	Matrix<T, Prop, Storage, Allocator> B;
	B = A;
	A.Reallocate(n, m);
	for (i = 0; i < m; i++)
	  for (j = 0; j < n; j++)
	    A(j, i) = conj(B(i, j));
      }
  }
  
  
  // TRANSPOSE //
  ///////////////

  
  ///////////////////////
  // ISSYMMETRICMATRIX //
  
  
  //! returns true if the matrix is symmetric
  template<class T, class Prop, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Prop, Storage, Allocator>& A)
  {
    return false;
  }

  
  //! returns true if the matrix is symmetric
  template<class T, class Storage, class Allocator>
  bool IsSymmetricMatrix(const Matrix<T, Symmetric, Storage, Allocator>& A)
  {
    return true;
  }


  // ISSYMMETRICMATRIX //
  ///////////////////////
  
} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_MATRIX_CXX
#endif
