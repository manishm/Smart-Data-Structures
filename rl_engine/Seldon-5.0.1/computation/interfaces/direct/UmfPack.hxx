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


#ifndef SELDON_FILE_UMFPACK_HXX

extern "C"
{
#include "umfpack.h"
}

namespace Seldon
{
  //!< base class to solve linear system by using UmfPack
  template<class T>
  class MatrixUmfPack_Base
  {
  public :
    Vector<double> Control, Info; //!< parameters for UmfPack
    void *Symbolic, *Numeric ; //!< pointers of UmfPack objects
    int n; //!< number of rows in the matrix
    bool display_info; //!< true if display is allowed
    
  public :
    MatrixUmfPack_Base();
    ~MatrixUmfPack_Base();
    
    void HideMessages();
    void ShowMessages();
    void Clear();
    
  };
  
  //! empty class
  template<class T>
  class MatrixUmfPack : public MatrixUmfPack_Base<T>
  {
  };
  
  //! class to solve linear system in double precision with UmfPack
  template<>
  class MatrixUmfPack<double> : public MatrixUmfPack_Base<double>
  {
  public :
    //! unsymmetric matrix in Column Sparse Row Format
    Matrix<double, General, ColSparse> Acsr;
    
    MatrixUmfPack();
    ~MatrixUmfPack();
    
    void Clear(){this->~MatrixUmfPack();}
    
    template<class Prop, class Storage,class Allocator>
    void FactorizeMatrix(Matrix<double,Prop,Storage,Allocator> & mat,
			 bool keep_matrix = false);
    
    template<class Allocator2>
    void Solve(Vector<double, VectFull, Allocator2>& x);
    
  };
  
  //! class to solve linear system in complex double precision with UmfPack
  template<>
  class MatrixUmfPack<complex<double> >
    : public MatrixUmfPack_Base<complex<double> >
  {
  public:
    //! Index of unsymmetric matrix in Column Sparse Row Format
    IVect Ptr, Ind;
    
    //! imaginary part of unsymmetric matrix in Column Sparse Row Format
    Vector<double> ValuesImag;
    
    //! real part of unsymmetric matrix in Column Sparse Row Format
    Vector<double> ValuesReal;
    
    MatrixUmfPack();
    ~MatrixUmfPack();
    
    void Clear(){this->~MatrixUmfPack();}
    
    template<class Prop, class Storage,class Allocator>
    void FactorizeMatrix(Matrix<complex<double>,Prop,Storage,Allocator> & mat,
			 bool keep_matrix = false);
    
    template<class Allocator2>
    void Solve(Vector<complex<double>,VectFull,Allocator2>& x);
    
  };

}

#define SELDON_FILE_UMFPACK_HXX
#endif
