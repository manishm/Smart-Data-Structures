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


#ifndef SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_CXX

#include "Matrix_SymComplexSparse.hxx"

namespace Seldon
{


  /****************
   * CONSTRUCTORS *
   ****************/


  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(): Matrix_Base<T, Allocator>()
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
  }


  //! Constructor.
  /*!
    Builds an empty i by j sparse matrix.
    \param i number of rows.
    \param j number of columns.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(int i, int j): Matrix_Base<T, Allocator>(i, i)
  {
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
  }


  //! Constructor.
  /*! Builds a sparse matrix of size i by j , with real_nz
    non-zero (stored) elements in the real part of the matrix and imag_nz
    non-zero elements in the imaginary part of the matrix.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements that are stored
    for the real part.
    \param imag_nz number of non-zero elements that are stored
    for the imaginary part.
    \note Matrix values are not initialized. Indices of non-zero entries
    are not initialized either.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(int i, int j, int real_nz, int imag_nz):
    Matrix_Base<T, Allocator>(i, i)
  {
    this->real_nz_ = real_nz;
    this->imag_nz_ = imag_nz;

#ifdef SELDON_CHECK_DIMENSIONS
    if ( (static_cast<long int>(2 * real_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + "Matrix_SymComplexSparse(int, int, int, int)",
		       string("There are more values to be stored (")
		       + to_str(real_nz) + " values for the real part and "
		       + to_str(imag_nz) + string(" values for the imaginary")
		       + " part) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ptr_ = reinterpret_cast<int*>( calloc(i + 1, sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	imag_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1)) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	imag_ptr_ = reinterpret_cast<int*>( calloc(i + 1, sizeof(int)) );

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	real_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1) ) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	real_ind_ = reinterpret_cast<int*>( calloc(real_nz_, sizeof(int)) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz) + " bytes to store "
		     + to_str(real_nz)
		     + " row or column indices (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	imag_ind_ = reinterpret_cast<int*>( calloc(imag_nz_, sizeof(int)) );
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(imag_ind_);
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz) + " bytes to store "
		     + to_str(imag_nz)
		     + " row or column indices (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->real_data_ = this->allocator_.allocate(real_nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	imag_data_ = NULL;
      }
    if (real_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz) + " bytes to store "
		     + to_str(real_nz) + " values (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->imag_data_ = this->allocator_.allocate(imag_nz_, this);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->allocator_.deallocate(this->real_data_, real_nz_);
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->allocator_.deallocate(this->real_data_, real_nz_);
	real_data_ = NULL;
      }
    if (imag_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz) + " bytes to store "
		     + to_str(imag_nz) + " values (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  Matrix_SymComplexSparse(int i, int j,
			  Vector<T, Storage0, Allocator0>& real_values,
			  Vector<int, Storage1, Allocator1>& real_ptr,
			  Vector<int, Storage2, Allocator2>& real_ind,
			  Vector<T, Storage0, Allocator0>& imag_values,
			  Vector<int, Storage1, Allocator1>& imag_ptr,
			  Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_Base<T, Allocator>(i, j)
  {
    real_nz_ = real_values.GetLength();
    imag_nz_ = imag_values.GetLength();
    
#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.
    
    if (real_ind.GetLength() != real_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(real_nz_)
		       + " values (real part) but "
		       + to_str(real_ind.GetLength())
		       + " row or column indices.");
      }

    if (imag_ind.GetLength() != imag_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(imag_nz_)
		       + " values (imaginary part) but "
		       + to_str(imag_ind.GetLength())
		       + " row or column indices.");
      }

    if (real_ptr.GetLength() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (real part)")
		       + " contains " + to_str(real_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (imag_ptr.GetLength() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (imaginary part)")
		       + " contains " + to_str(imag_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if ( (static_cast<long int>(2 * real_nz_ - 2)
	  / static_cast<long int>(i + 1)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2)
	  / static_cast<long int>(i + 1)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + string("Matrix_SymComplexSparse(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(real_values.GetLength())
		       + " values for the real part and "
		       + to_str(real_values.GetLength()) + " values for"
		       + string(" the imaginary part) than elements in the")
		       + " matrix (" + to_str(i) + " by " + to_str(i) + ").");
      }
#endif

    this->real_ptr_ = real_ptr.GetData();
    this->imag_ptr_ = imag_ptr.GetData();
    this->real_ind_ = real_ind.GetData();
    this->imag_ind_ = imag_ind.GetData();
    this->real_data_ = real_values.GetData();
    this->imag_data_ = imag_values.GetData();

    real_ptr.Nullify();
    imag_ptr.Nullify();
    real_ind.Nullify();
    imag_ind.Nullify();
    real_values.Nullify();
    imag_values.Nullify();
  }


  //! Copy constructor
  template <class T, class Prop, class Storage, class Allocator>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::Matrix_SymComplexSparse(const Matrix_SymComplexSparse<T, Prop,
			    Storage, Allocator>& A)
  {
    this->m_ = 0;
    this->n_ = 0;
    real_nz_ = 0;
    imag_nz_ = 0;
    real_ptr_ = NULL;
    imag_ptr_ = NULL;
    real_ind_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
    
    this->Copy(A);
  }
  
  
  /**************
   * DESTRUCTOR *
   **************/


  //! Destructor.
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::~Matrix_SymComplexSparse()
  {
    this->m_ = 0;
    this->n_ = 0;

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (real_ptr_ != NULL)
	  {
	    free(real_ptr_);
	    real_ptr_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (imag_ptr_ != NULL)
	  {
	    free(imag_ptr_);
	    imag_ptr_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ptr_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (real_ind_ != NULL)
	  {
	    free(real_ind_);
	    real_ind_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	real_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (imag_ind_ != NULL)
	  {
	    free(imag_ind_);
	    imag_ind_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	imag_ind_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (this->real_data_ != NULL)
	  {
	    this->allocator_.deallocate(this->real_data_, real_nz_);
	    this->real_data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->real_nz_ = 0;
	this->real_data_ = NULL;
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	if (this->imag_data_ != NULL)
	  {
	    this->allocator_.deallocate(this->imag_data_, imag_nz_);
	    this->imag_data_ = NULL;
	  }
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->imag_nz_ = 0;
	this->imag_data_ = NULL;
      }
#endif

    this->real_nz_ = 0;
    this->imag_nz_ = 0;
  }
  

  //! Clears the matrix.
  /*! This methods is equivalent to the destructor. On exit, the matrix
    is empty (0x0).
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Clear()
  {
    this->~Matrix_SymComplexSparse();
  }


  /*********************
   * MEMORY MANAGEMENT *
   *********************/


  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by
    'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j,
	  Vector<T, Storage0, Allocator0>& real_values,
	  Vector<int, Storage1, Allocator1>& real_ptr,
	  Vector<int, Storage2, Allocator2>& real_ind,
	  Vector<T, Storage0, Allocator0>& imag_values,
	  Vector<int, Storage1, Allocator1>& imag_ptr,
	  Vector<int, Storage2, Allocator2>& imag_ind)
  {
    this->Clear();
    this->m_ = i;
    this->n_ = i;
    real_nz_ = real_values.GetLength();
    imag_nz_ = imag_values.GetLength();
    
#ifdef SELDON_CHECK_DIMENSIONS
    // Checks whether vector sizes are acceptable.
    
    if (real_ind.GetLength() != real_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(real_nz_)
		       + " values (real part) but "
		       + to_str(real_ind.GetLength())
		       + " row or column indices.");
      }

    if (imag_ind.GetLength() != imag_nz_)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are ") + to_str(imag_nz_)
		       + " values (imaginary part) but "
		       + to_str(imag_ind.GetLength())
		       + " row or column indices.");
      }

    if (real_ptr.GetLength() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (real part)")
		       + " contains " + to_str(real_ptr.GetLength() - 1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if (imag_ptr.GetLength() - 1 != i)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("The vector of start indices (imaginary part)")
		       + " contains " + to_str(imag_ptr.GetLength()-1)
		       + string(" row or column start indices (plus the")
		       + " number of non-zero entries) but there are "
		       + to_str(i) + " rows or columns ("
		       + to_str(i) + " by " + to_str(i) + " matrix).");
      }

    if ( (static_cast<long int>(2 * real_nz_ - 2) / static_cast<long int>(i)
	  >= static_cast<long int>(i + 1)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2) / static_cast<long int>(i)
	  >= static_cast<long int>(i + 1)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::SetData(int, int, ")
		       + string("const Vector&, const Vector&, const Vector&")
		       + ", const Vector&, const Vector&, const Vector&)",
		       string("There are more values (")
		       + to_str(real_values.GetLength())
		       + " values for the real part and "
		       + to_str(real_values.GetLength()) + string(" values")
		       + string(" for the imaginary part) than elements in")
		       + " the matrix (" + to_str(i)
		       + " by " + to_str(i) + ").");
      }
#endif

    this->real_ptr_ = real_ptr.GetData();
    this->imag_ptr_ = imag_ptr.GetData();
    this->real_ind_ = real_ind.GetData();
    this->imag_ind_ = imag_ind.GetData();
    this->real_data_ = real_values.GetData();
    this->imag_data_ = imag_values.GetData();

    real_ptr.Nullify();
    imag_ptr.Nullify();
    real_ind.Nullify();
    imag_ind.Nullify();
    real_values.Nullify();
    imag_values.Nullify();
  }

  
  //! Redefines the matrix.
  /*! It clears the matrix and sets it to a new matrix defined by arrays
    'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part).
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero entries (real part).
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_nz number of non-zero entries (imaginary part).
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning On exit, arrays 'real_values', 'real_ptr', 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are managed by the matrix.
    For example, it means that the destructor will release those arrays;
    therefore, the user mustn't release those arrays.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  SetData(int i, int j, int real_nz,
	  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
	  ::pointer real_values,
	  int* real_ptr, int* real_ind, int imag_nz,
	  typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
	  ::pointer imag_values,
	  int* imag_ptr, int* imag_ind)
  {
    this->Clear();

    this->m_ = i;
    this->n_ = i;

    this->real_nz_ = real_nz;
    this->imag_nz_ = imag_nz;

    real_data_ = real_values;
    imag_data_ = imag_values;
    real_ind_ = real_ind;
    imag_ind_ = imag_ind;
    real_ptr_ = real_ptr;
    imag_ptr_ = imag_ptr;
  }


  //! Clears the matrix without releasing memory.
  /*!
    On exit, the matrix is empty and the memory has not been released.
    It is useful for low level manipulations on a Matrix instance.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Nullify()
  {
    this->data_ = NULL;
    this->m_ = 0;
    this->n_ = 0;
    real_nz_ = 0;
    real_ptr_ = NULL;
    real_ind_ = NULL;
    imag_nz_ = 0;
    imag_ptr_ = NULL;
    imag_ind_ = NULL;
    real_data_ = NULL;
    imag_data_ = NULL;
  }


  //! Copies a matrix
  template <class T, class Prop, class Storage, class Allocator>
  inline void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::
  Copy(const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Clear();
    int i = A.m_;
    int j = A.n_;
    real_nz_ = A.real_nz_;
    imag_nz_ = A.imag_nz_;
    this->m_ = i;
    this->n_ = j;
    if ((i == 0)||(j == 0))
      {
	this->m_ = 0;
	this->n_ = 0;
	this->real_nz_ = 0;
	this->imag_nz_ = 0;
	return;
      }
    
#ifdef SELDON_CHECK_DIMENSIONS
    if ( (static_cast<long int>(2 * real_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) ||
	 (static_cast<long int>(2 * imag_nz_ - 2) / static_cast<long int>(i+1)
	  >= static_cast<long int>(i)) )
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
	throw WrongDim(string("Matrix_SymComplexSparse::")
		       + "Matrix_SymComplexSparse(int, int, int, int)",
		       string("There are more values to be stored (")
		       + to_str(real_nz_) + " values for the real part and "
		       + to_str(imag_nz_) + string(" values for the imaginary")
		       + " part) than elements in the matrix ("
		       + to_str(i) + " by " + to_str(j) + ").");
      }
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	real_ptr_ = reinterpret_cast<int*>( calloc(i + 1, sizeof(int)) );
	memcpy(this->real_ptr_, A.real_ptr_, i+1);
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	imag_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1)) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	imag_ptr_ = reinterpret_cast<int*>( calloc(i + 1, sizeof(int)) );
	memcpy(this->imag_ptr_, A.imag_ptr_, i+1);

#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	real_ptr_ = 0;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ptr_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * (i + 1) ) + " bytes to store "
		     + to_str(i + 1) + string(" row or column")
		     + " start indices (for the imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	real_ind_ = reinterpret_cast<int*>( calloc(real_nz_, sizeof(int)) );
	memcpy(this->real_ind_, A.real_ind_, real_nz_);
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz_) + " bytes to store "
		     + to_str(real_nz_)
		     + " row or column indices (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif
	
	imag_ind_ = reinterpret_cast<int*>( calloc(imag_nz_, sizeof(int)) );
	memcpy(this->imag_ind_, A.imag_ind_, imag_nz_);
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_ind_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	real_nz_ = 0;
	imag_nz_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(imag_ind_);
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (imag_ind_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz_) + " bytes to store "
		     + to_str(imag_nz_)
		     + " row or column indices (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->real_data_ = this->allocator_.allocate(real_nz_, this);
	this->allocator_.memorycpy(this->real_data_, A.real_data_, real_nz_);
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	imag_data_ = NULL;
      }
    if (real_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * real_nz_) + " bytes to store "
		     + to_str(real_nz_) + " values (real part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif

#ifdef SELDON_CHECK_MEMORY
    try
      {
#endif

	this->imag_data_ = this->allocator_.allocate(imag_nz_, this);
	this->allocator_.memorycpy(this->imag_data_, A.imag_data_, imag_nz_);
	
#ifdef SELDON_CHECK_MEMORY
      }
    catch (...)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->allocator_.deallocate(this->real_data_, real_nz_);
	this->real_data_ = NULL;
	this->imag_data_ = NULL;
      }
    if (real_data_ == NULL)
      {
	this->m_ = 0;
	this->n_ = 0;
	free(real_ptr_);
	free(imag_ptr_);
	real_ptr_ = NULL;
	imag_ptr_ = NULL;
	free(real_ind_);
	free(imag_ind_);
	real_ind_ = NULL;
	imag_ind_ = NULL;
	this->allocator_.deallocate(this->real_data_, real_nz_);
	real_data_ = NULL;
      }
    if (imag_data_ == NULL && i != 0)
      throw NoMemory(string("Matrix_SymComplexSparse::")
		     + "Matrix_SymComplexSparse(int, int, int, int)",
		     string("Unable to allocate ")
		     + to_str(sizeof(int) * imag_nz_) + " bytes to store "
		     + to_str(imag_nz_) + " values (imaginary part), for a "
		     + to_str(i) + " by " + to_str(i) + " matrix.");
#endif
    
  }
  
  
  /*******************
   * BASIC FUNCTIONS *
   *******************/
  

  //! Returns the number of elements stored in memory.
  /*!
    Returns the number of elements stored in memory, i.e.
    the cumulated number of non-zero entries of both the real and
    the imaginary part.
    \return The number of elements stored in memory.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetDataSize() const
  {
    return real_nz_ + imag_nz_;
  }


  //! Returns (row or column) start indices for the real part.
  /*!
    Returns the array ('ptr_') of start indices for the real part.
    \return The array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealPtr() const
  {
    return real_ptr_;
  }


  //! Returns (row or column) start indices for the imaginary part.
  /*!
    Returns the array ('ptr_') of start indices for the imaginary part.
    \return The array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagPtr() const
  {
    return imag_ptr_;
  }


  //! Returns (row or column) indices of non-zero entries for the real part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the real part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealInd() const
  {
    return real_ind_;
  }


  //! Returns (row or column) indices of non-zero entries for
  //! the imaginary part.
  /*!
    Returns the array ('ind_') of (row or column) indices
    of non-zero entries for the imaginary part. This array defines non-zero
    entries indices if coupled with (column or row) start indices.
    \return The array of (row or column) indices of
    non-zero entries for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagInd() const
  {
    return imag_ind_;
  }


  //! Returns the length of the array of start indices for the real part.
  /*!
    \return The length of the array of start indices for the real part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealPtrSize() const
  {
    return (this->m_ + 1);
  }

  
  //! Returns the length of the array of start indices for the imaginary part.
  /*!
    \return The length of the array of start indices for the imaginary part.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagPtrSize() const
  {
    return (this->m_ + 1);
  }

  
  //! Returns the length of the array of (column or row) indices for
  //! the real part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored) for the real part. This array
    defines non-zero entries indices if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices for
    the real part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetRealIndSize() const
  {
    return real_nz_;
  }


  //! Returns the length of the array of (column or row) indices
  //! for the imaginary part.
  /*!
    Returns the length of the array ('ind_') of (row or column) indices
    of non-zero entries (that are stored) for the imaginary part. This array
    defines non-zero entries indices if coupled with (column or row)
    start indices.
    \return The length of the array of (column or row) indices
    for the imaginary part.
    \note The length of the array of (column or row) indices is the
    number of non-zero entries that are stored.
  */
  template <class T, class Prop, class Storage, class Allocator>
  int Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::GetImagIndSize() const
  {
    return imag_nz_;
  }


  //! Returns the array of values of the real part.
  /*!
    \return The array 'real_data_' of values of the real part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  T* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetRealData() const
  {
    return real_data_;
  }


  //! Returns the array of values of the imaginary part.
  /*!
    \return The array 'imag_data_' of values of the imaginary part..
  */
  template <class T, class Prop, class Storage, class Allocator>
  T* Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::GetImagData() const
  {
    return imag_data_;
  }


  /**********************************
   * ELEMENT ACCESS AND AFFECTATION *
   **********************************/


  //! Access operator.
  /*!
    Returns the value of element (i, j).
    \param i row index.
    \param j column index.
    \return Element (i, j) of the matrix.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline complex<typename Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::value_type>
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::operator() (int i, int j) const
  {

#ifdef SELDON_CHECK_BOUNDS
    if (i < 0 || i >= this->m_)
      throw WrongRow("Matrix_SymComplexSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->m_-1)
		     + "], but is equal to " + to_str(i) + ".");
    if (j < 0 || j >= this->n_)
      throw WrongCol("Matrix_SymComplexSparse::operator()",
		     string("Index should be in [0, ") + to_str(this->n_-1)
		     + "], but is equal to " + to_str(j) + ".");
#endif

    int real_k, imag_k, l;
    int real_a, real_b;
    int imag_a, imag_b;

    // Only the upper part is stored.
    if (i > j)
      {
	l = i;
	i = j;
	j = l;
      }

    real_a = real_ptr_[Storage::GetFirst(i, j)];
    real_b = real_ptr_[Storage::GetFirst(i, j) + 1];

    imag_a = imag_ptr_[Storage::GetFirst(i, j)];
    imag_b = imag_ptr_[Storage::GetFirst(i, j) + 1];
    
    if (real_a != real_b)
      {
	l = Storage::GetSecond(i, j);
	for (real_k = real_a;
	     (real_k < real_b - 1) && (real_ind_[real_k] < l);
	     real_k++);
	if (imag_a != imag_b)
	  {
	    for (imag_k = imag_a;
		 (imag_k < imag_b - 1) && (imag_ind_[imag_k] < l);
		 imag_k++);
	    if (real_ind_[real_k] == l)
	      {
		if (imag_ind_[imag_k] == l)
		  return complex<T>(real_data_[real_k], imag_data_[imag_k]);
		else
		  return complex<T>(real_data_[real_k], T(0));
	      }
	    else
	      if (imag_ind_[imag_k] == l)
		return complex<T>(T(0), imag_data_[imag_k]);
	      else
		return complex<T>(T(0), T(0));
	  }
	else
	  {
	    if (real_ind_[real_k] == l)
	      return complex<T>(real_data_[real_k], T(0));
	    else
	      return complex<T>(T(0), T(0));
	  }
      }
    else
      {
	if (imag_a != imag_b)
	  {
	    l = Storage::GetSecond(i, j);
	    for (imag_k = imag_a;
		 (imag_k < imag_b - 1) && (imag_ind_[imag_k] < l);
		 imag_k++);
	    if (imag_ind_[imag_k] == l)
	      return complex<T>(T(0), imag_data_[imag_k]);
	    else
	      return complex<T>(T(0), T(0));
	  }
	else
	  return complex<T>(T(0), T(0));
      }
  }


  //! Duplicates a matrix (assignment operator).
  /*!
    \param A matrix to be copied.
    \note Memory is duplicated: 'A' is therefore independent from the current
    instance after the copy.
  */
  template <class T, class Prop, class Storage, class Allocator>
  inline Matrix_SymComplexSparse<T, Prop, Storage, Allocator>&
  Matrix_SymComplexSparse<T, Prop, Storage, Allocator>
  ::operator= (const Matrix_SymComplexSparse<T, Prop, Storage, Allocator>& A)
  {
    this->Copy(A);

    return *this;
  }
  
  
  /************************
   * CONVENIENT FUNCTIONS *
   ************************/


  //! Displays the matrix on the standard output.
  /*!
    Displays elements on the standard output, in text format.
    Each row is displayed on a single line and elements of
    a row are delimited by tabulations.
  */
  template <class T, class Prop, class Storage, class Allocator>
  void Matrix_SymComplexSparse<T, Prop, Storage, Allocator>::Print() const
  {
    for (int i = 0; i < this->m_; i++)
      {
	for (int j = 0; j < this->n_; j++)
	  cout << (*this)(i, j) << "\t";
	cout << endl;
      }
  }


 
  /////////////////////////////////
  // MATRIX<COLSYMCOMPLEXSPARSE> //
  /////////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSymComplexSparse, Allocator>::Matrix()  throw():
    Matrix_SymComplexSparse<T, Prop, ColSymComplexSparse, Allocator>()
  {
  }


  //! Builds a i by j matrix.
  /*!
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSymComplexSparse, Allocator>
  ::Matrix(int i, int j):
    Matrix_SymComplexSparse<T, Prop, ColSymComplexSparse, Allocator>(i, j,
                                                                     0, 0)
  {
  }


  //! Constructor.
  /*! Builds a i by j matrix with real_nz and imag_nz non-zero (stored)
    elements for the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements that are stored
    for the real part.
    \param imag_nz number of non-zero elements that are stored
    for the imaginary part.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, ColSymComplexSparse, Allocator>::Matrix(int i, int j,
							  int real_nz,
							  int imag_nz):
    Matrix_SymComplexSparse<T, Prop,
			    ColSymComplexSparse, Allocator>(i, j,
							    real_nz, imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, ColSymComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& real_values,
	 Vector<int, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<T, Storage0, Allocator0>& imag_values,
	 Vector<int, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_SymComplexSparse<T, Prop,
			    ColSymComplexSparse, Allocator>(i, j,
							    real_values,
							    real_ptr,
							    real_ind,
							    imag_values,
							    imag_ptr,
							    imag_ind)
  {
  }



  /////////////////////////////////
  // MATRIX<ROWSYMCOMPLEXSPARSE> //
  /////////////////////////////////


  /****************
   * CONSTRUCTORS *
   ****************/

  //! Default constructor.
  /*!
    Builds an empty 0x0 matrix.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSymComplexSparse, Allocator>::Matrix()  throw():
    Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>()
  {
  }


  //! Builds a i by j matrix.
  /*!
    \param i number of rows.
    \param j number of columns.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSymComplexSparse, Allocator>
  ::Matrix(int i, int j):
    Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>(i, j,
                                                                     0, 0)
  {
  }


  /*! Builds a i by j matrix with real_nz and imag_nz non-zero (stored)
    elements for the real part and the imaginary part respectively.
    \param i number of rows.
    \param j number of columns.
    \param real_nz number of non-zero elements that are stored
    for the real part.
    \param imag_nz number of non-zero elements that are store
    for the imaginary part.
    \note Matrix values are not initialized.
    \warning 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  Matrix<T, Prop, RowSymComplexSparse, Allocator>
  ::Matrix(int i, int j, int real_nz, int imag_nz):
    Matrix_SymComplexSparse<T, Prop, RowSymComplexSparse, Allocator>(i, j,
								     real_nz,
								     imag_nz)
  {
  }


  //! Constructor.
  /*!
    Builds a i by j sparse matrix with non-zero values and indices
    provided by 'real_values' (values of the real part), 'real_ptr'
    (pointers for the real part), 'real_ind' (indices for the real part),
    'imag_values' (values of the imaginary part), 'imag_ptr'
    (pointers for the imaginary part) and 'imag_ind' (indices for the
    imaginary part). Input vectors are released and are empty on exit.
    \param i number of rows.
    \param j number of columns.
    \param real_values values of non-zero entries for the real part.
    \param real_ptr row or column start indices for the real part.
    \param real_ind row or column indices for the real part.
    \param imag_values values of non-zero entries for the imaginary part.
    \param imag_ptr row or column start indices for the imaginary part.
    \param imag_ind row or column indices for the imaginary part.
    \warning Input vectors 'real_values', 'real_ptr' and 'real_ind',
    'imag_values', 'imag_ptr' and 'imag_ind' are empty on exit.
    Moreover 'j' is assumed to be equal to 'i' so that 'j' is discarded.
  */
  template <class T, class Prop, class Allocator>
  template <class Storage0, class Allocator0,
	    class Storage1, class Allocator1,
	    class Storage2, class Allocator2>
  Matrix<T, Prop, RowSymComplexSparse, Allocator>::
  Matrix(int i, int j,
	 Vector<T, Storage0, Allocator0>& real_values,
	 Vector<int, Storage1, Allocator1>& real_ptr,
	 Vector<int, Storage2, Allocator2>& real_ind,
	 Vector<T, Storage0, Allocator0>& imag_values,
	 Vector<int, Storage1, Allocator1>& imag_ptr,
	 Vector<int, Storage2, Allocator2>& imag_ind):
    Matrix_SymComplexSparse<T, Prop,
			    RowSymComplexSparse, Allocator>(i, j,
							    real_values,
							    real_ptr,
							    real_ind,
							    imag_values,
							    imag_ptr,
							    imag_ind)
  {
  }


} // namespace Seldon.

#define SELDON_FILE_MATRIX_SYMCOMPLEXSPARSE_CXX
#endif
