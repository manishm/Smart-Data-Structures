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


#ifndef SELDON_FILE_SELDON_HXX


#include "SeldonHeader.hxx"

namespace Seldon
{
  

  // Vector class - specialized for each used type.
  template <class T, class Storage, class Allocator>
  class Vector
  {
    // Nothing in it: no default vector is supplied so as to avoid suprises!
  };


  // Matrix class - specialized for each used type.
  template <class T, class Prop, class Storage, class Allocator>
  class Matrix
  {
    // Nothing in it: no default matrix is supplied so as to avoid suprises!
  };


  class SeldonTranspose
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_TRANSPOSE cblas_status_;
#endif
  protected:
    // 0: Trans, 1: NoTrans, 2: ConjTrans.
    int status_;
  public:
    SeldonTranspose(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasTrans;
      else if (status_ == 1)
	cblas_status_ = CblasNoTrans;
      else
	cblas_status_ = CblasConjTrans;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    SeldonTranspose(const enum CBLAS_TRANSPOSE status):
      cblas_status_(status)
    {
      if (cblas_status_ == CblasTrans)
	status_ = 0;
      else if (cblas_status_ == CblasNoTrans)
	status_ = 1;
      else
	status_ = 2;
    }
#endif
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_TRANSPOSE() const
    {
      return cblas_status_;
    }
#endif
    
    char Char() const
    {
      if (status_ == 0)
	return 'T';
      else if (status_ == 1)
	return 'N';
      else
	return 'C';
    }
    
    char RevChar() const
    {
      if (status_ == 0)
	return 'N';
      else if (status_ == 1)
	return 'T';
      else
	return 'N';
    }
    
    bool Trans() const {return (status_ == 0);}
    bool NoTrans() const {return (status_ == 1);}
    bool ConjTrans() const {return (status_ == 2);}
  };

  class class_SeldonTrans: public SeldonTranspose
  {
  public:
    class_SeldonTrans(): SeldonTranspose(0) {};
  };

  class class_SeldonNoTrans: public SeldonTranspose
  {
  public:
    class_SeldonNoTrans(): SeldonTranspose(1) {};
  };

  class class_SeldonConjTrans: public SeldonTranspose
  {
  public:
    class_SeldonConjTrans(): SeldonTranspose(2) {};
  };

  class_SeldonTrans SeldonTrans;
  class_SeldonNoTrans SeldonNoTrans;
  class_SeldonConjTrans SeldonConjTrans;

  //

  class SeldonDiag
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_DIAG cblas_status_;
#endif
  protected:
    // 0: NonUnit, 1: Unit.
    int status_;
  public:
    SeldonDiag(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasNonUnit;
      else
	cblas_status_ = CblasUnit;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_DIAG() const
    {
      return cblas_status_;
    }
#endif
    
    char Char() const { return (status_ == 0) ? 'N' : 'U'; }
    bool NonUnit() const {return (status_ == 0);}
    bool Unit() const {return (status_ == 1);}
  };

  class class_SeldonNonUnit: public SeldonDiag
  {
  public:
    class_SeldonNonUnit(): SeldonDiag(0) {};
  };

  class class_SeldonUnit: public SeldonDiag
  {
  public:
    class_SeldonUnit(): SeldonDiag(1) {};
  };

  class_SeldonNonUnit SeldonNonUnit;
  class_SeldonUnit SeldonUnit;

  //

  class SeldonUplo
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_UPLO cblas_status_;
#endif
  protected:
    // 0: Upper, 1: Lower.
    int status_;
  public:
    SeldonUplo(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasUpper;
      else
	cblas_status_ = CblasLower;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_UPLO() const
    {
      return cblas_status_;
    }
#endif
    
    operator char() const
    {
      return (status_ == 0) ? 'U' : 'L';
    }

    bool Upper() const {return (status_ == 0);}
    bool Lower() const {return (status_ == 1);}

    char Char() const
    {
      return (status_ == 0) ? 'U' : 'L';
    }
    char RevChar() const
    {
      return (status_ == 0) ? 'L' : 'U';
    }

  };

  SeldonUplo SeldonUpper(0);
  SeldonUplo SeldonLower(1);

  //

  class SeldonNorm
  {
  protected:
    // 0: Infinity-norm, 1: 1-norm.
    int status_;
  public:
    SeldonNorm(int status)
    {
      status_ = status;
    }
    
    operator char() const
    {
      return (status_ == 0) ? 'I' : '1';
    }
    
    char Char() const
    {
      return (status_ == 0) ? 'I' : '1';
    }
    char RevChar() const
    {
      return (status_ == 0) ? '1' : 'I';
    }
    
  };
  
  SeldonNorm SeldonNormInf(0);
  SeldonNorm SeldonNorm1(1);

  //

  class SeldonConjugate
  {
  protected:
    // false: Unconj, true: Conj.
    bool status_;
  public:
    SeldonConjugate(bool status)
    {
      status_ = status;
    }
    inline bool Conj() const
    {
      return status_;
    }
  };

  SeldonConjugate SeldonUnconj(false);
  SeldonConjugate SeldonConj(true);
  

  //

  class SeldonSide
  {
#ifdef SELDON_WITH_CBLAS
  protected:
    CBLAS_SIDE cblas_status_;
#endif
  protected:
    // 0: Left, 1: Right.
    int status_;
  public:
    SeldonSide(int status)
    {
      status_ = status;
#ifdef SELDON_WITH_CBLAS
      if (status_ == 0)
	cblas_status_ = CblasLeft;
      else
	cblas_status_ = CblasRight;
#endif
    }
#ifdef SELDON_WITH_CBLAS
    SeldonSide(const enum CBLAS_SIDE status):
      cblas_status_(status)
    {
      if (cblas_status_ == CblasLeft)
	status_ = 0;
      else
	status_ = 1;
    }
#endif
#ifdef SELDON_WITH_CBLAS
    operator CBLAS_SIDE() const
    {
      return cblas_status_;
    }
#endif
    
    bool Left() const {return (status_ == 0);}
    bool Right() const {return (status_ == 1);}
  };

  class class_SeldonLeft: public SeldonSide
  {
  public:
    class_SeldonLeft(): SeldonSide(0) {};
  };

  class class_SeldonRight: public SeldonSide
  {
  public:
    class_SeldonRight(): SeldonSide(1) {};
  };

  class_SeldonLeft SeldonLeft;
  class_SeldonRight SeldonRight;


} // namespace Seldon.


#include "share/Common.cxx"

// Memory management.
#include "share/Allocator.cxx"

// Storage type.
#include "share/Storage.cxx"

#include "share/Errors.cxx"

#include "array3d/Array3D.cxx"
#include "matrix/Matrix_Base.cxx"
#include "matrix/Matrix_Pointers.cxx"
#include "matrix/Matrix_Triangular.cxx"
#include "matrix/Matrix_Symmetric.cxx"
#include "matrix/Matrix_Hermitian.cxx"
#include "matrix_sparse/Matrix_Sparse.cxx"
#include "matrix_sparse/Matrix_ComplexSparse.cxx"
#include "matrix_sparse/Matrix_SymSparse.cxx"
#include "matrix_sparse/Matrix_SymComplexSparse.cxx"
#include "matrix/Matrix_SymPacked.cxx"
#include "matrix/Matrix_HermPacked.cxx"
#include "matrix/Matrix_TriangPacked.cxx"
#include "vector/Vector.cxx"
#include "vector/Functions_Arrays.cxx"
#include "vector/SparseVector.cxx"
#include "matrix/Functions.cxx"
#include "computation/basic_functions/Functions_Matrix.cxx"
#include "computation/basic_functions/Functions_Vector.cxx"
#include "computation/basic_functions/Functions_MatVect.cxx"

// Blas interface.
#ifdef SELDON_WITH_CBLAS
#include "computation/interfaces/Blas_1.cxx"
#include "computation/interfaces/Blas_2.cxx"
#include "computation/interfaces/Blas_3.cxx"
#endif

// Lapack interface.
#ifdef SELDON_WITH_LAPACK
#include "computation/interfaces/Lapack_LinearEquations.cxx"
#include "computation/interfaces/Lapack_LeastSquares.cxx"
#include "computation/interfaces/Lapack_Eigenvalues.cxx"
#endif // SELDON_WITH_LAPACK.


#define SELDON_FILE_SELDON_HXX
#endif
