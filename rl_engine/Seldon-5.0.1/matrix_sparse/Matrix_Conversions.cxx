// Copyright (C) 2003-2009 Marc Duruflé
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


#ifndef SELDON_FILE_MATRIX_CONVERSIONS_CXX

namespace Seldon
{
  
  /*
    From CSR formats to "Matlab" coordinate format.
  */
  
  
  //! Conversion from RowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, RowSparse,
			       Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i < m; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndRow(j) = i + index;
	  IndCol(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }
  
  
  //! Conversion from ColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop, ColSparse,
			       Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    IndCol.Reallocate(nnz);
    IndRow.Reallocate(nnz);
    Val.Reallocate(nnz);
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    for (i = 0; i< n; i++)
      for (j = ptr[i]; j< ptr[i+1]; j++)
	{
	  IndCol(j) = i + index;
	  IndRow(j) = ind[j] + index;
	  Val(j) = val[j];
	}
  }
  
  
  //! Conversion from RowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       RowSymSparse, Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  if (ind[ptr[i]] == i)
	    nnz--;
	
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	IVect Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;
	      
	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }
	
	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);
	
	// ... and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
	
      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j< ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }
  
  
  //! Conversion from ColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ColSymSparse, Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    int* ptr = A.GetPtr();
    int* ind = A.GetInd();
    T* val = A.GetData();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; i++)
	    if (ind[j] == i)
	      nnz--;
	
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	IVect Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j < ptr[i + 1]; j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = ind[j] + index;
	      Val(nb) = val[j];
	      Ptr(ind[j])++;
	      nb++;
	      
	      if (ind[j] != i)
		{
		  IndRow(nb) = ind[j] + index;
		  IndCol(nb) = i + index;
		  Val(nb) = val[j];
		  Ptr(i)++;
		  nb++;
		}
	    }
	
	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);
	
	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
	
      }
    else
      {
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	for (i = 0; i < m; i++)
	  for (j = ptr[i]; j< ptr[i + 1]; j++)
	    {
	      IndRow(j) = i + index;
	      IndCol(j) = ind[j] + index;
	      Val(j) = val[j];
	    }
      }
  }
  
  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */
  
  
  //! Conversion from ArrayRowSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayRowSparse, Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < m; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	{
	  IndRow(nb) = i + index;
	  IndCol(nb) = A.Index(i, j) + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }
  
  
  //! Conversion from ArrayColSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayColSparse, Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int n = A.GetN();
    int nnz = A.GetDataSize();
    // Allocating arrays.
    IndRow.Reallocate(nnz);
    IndCol.Reallocate(nnz);
    Val.Reallocate(nnz);
    int nb = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetColumnSize(i); j++)
	{
	  IndRow(nb) = A.Index(i, j) + index;
	  IndCol(nb) = i + index;
	  Val(nb) = A.Value(i, j);
	  nb++;
	}
  }
  
  
  //! Conversion from ArrayRowSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayRowSymSparse, Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;
	
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	IVect Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;
	      
	      if (A.Index(i, j) != i)
		{
		  IndRow(nb) = A.Index(i, j) + index;
		  IndCol(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }
	
	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);
	
	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }
  
  
  //! Conversion from ArrayColSymSparse to coordinate format.
  template<class T, class Prop, class Allocator1, class Allocator2>
  void
  ConvertMatrix_to_Coordinates(const Matrix<T, Prop,
			       ArrayColSymSparse, Allocator1>& A,
			       IVect& IndRow, IVect& IndCol,
			       Vector<T, VectFull, Allocator2>& Val,
			       int index = 0, bool sym = false)
  {
    int i, j;
    int m = A.GetM();
    int nnz = A.GetDataSize();
    if (sym)
      {
	nnz *= 2;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    if (A.Index(i, j) == i)
	      nnz--;
	
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	IVect Ptr(m);
	Ptr.Zero();
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetRowSize(i); j++)
	    {
	      IndRow(nb) = i + index;
	      IndCol(nb) = A.Index(i, j) + index;
	      Val(nb) = A.Value(i, j);
	      Ptr(A.Index(i, j))++;
	      nb++;
	      
	      if (A.Index(i, j) != i)
		{
		  IndRow(nb) = A.Index(i, j) + index;
		  IndCol(nb) = i + index;
		  Val(nb) = A.Value(i, j);
		  Ptr(i)++;
		  nb++;
		}
	    }
	
	// Sorting the row numbers...
	Sort(IndRow, IndCol, Val);
	
	// ...and the column numbers.
	int offset = 0;
	for (i = 0; i < m; i++)
	  {
	    Sort(offset, offset + Ptr(i) - 1, IndCol, Val);
	    offset += Ptr(i);
	  }
      }
    else
      {
	// Allocating arrays.
	IndRow.Reallocate(nnz);
	IndCol.Reallocate(nnz);
	Val.Reallocate(nnz);
	int nb = 0;
	for (i = 0; i < m; i++)
	  for (j = 0; j < A.GetColumnSize(i); j++)
	    {
	      IndRow(nb) = A.Index(i, j) + index;
	      IndCol(nb) = i + index;
	      Val(nb) = A.Value(i, j);
	      nb++;
	    }
      }
  }
  
  
  /*
    From "Matlab" coordinate format to CSR formats.
  */
  
  
  //! Conversion from coordinate format to RowSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, RowSparse, Allocator>& A,
				 int index = 0)
  {
    // Assuming that there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;
    
    int i;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    IVect Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i) + 1)++;
      }
    
    for (i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
    // Sorts 'IndCol'.
    for (i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);
    
    A.SetData(m, n, nnz, Val, Ptr, IndCol);
  }
  
  
  //! Conversion from coordinate format to ColSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ColSparse, Allocator>& A,
				 int index = 0)
  {
    // Assuming that there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;
    
    int i;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // Sorts the array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    IVect Ptr(n + 1);
    Ptr.Zero();
    for (i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i) + 1)++;
      }
    
    for (i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
    // Sorts 'IndRow'
    for (i = 0; i < n; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndRow, Val);
    
    A.SetData(m, n, nnz, Val, Ptr, IndRow);
  }
  
  
  //! Conversion from coordinate format to RowSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, RowSymSparse, Allocator>& A,
				 int index = 0)
  {
    // Assuming there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;
    
    int i;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;
    
    if (nb_low > 0)
      {
	int nb = 0;
	for (i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }
	
	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    IVect Ptr(m + 1);
    Ptr.Zero();
    for (int i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i) + 1)++;
      }
    
    for (int i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
    // Sorts 'IndCol'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);
    
    A.SetData(m, n, nnz, Val, Ptr, IndCol);
  }
  
  
  //! Conversion from coordinate format to ColSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ColSymSparse, Allocator>& A,
				 int index = 0)
  {
    // Assuming there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;
    
    int i;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;
    
    if (nb_low > 0)
      {
	int nb = 0;
	for (i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }
	
	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    IVect Ptr(m + 1);
    Ptr.Zero();
    for (i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i) + 1)++;
      }
    
    for (i = 0; i < m; i++)
      Ptr(i + 1) += Ptr(i);
    
    // Sorts 'IndCol'.
    for (int i = 0; i < m; i++)
      Sort(Ptr(i), Ptr(i + 1) - 1, IndCol, Val);
    
    A.SetData(m, n, nnz, Val, Ptr, IndCol);
  }
  
  
  /*
    From Sparse Array formats to "Matlab" coordinate format.
  */
  
  
  //! Conversion from coordinate format to ArrayRowSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ArrayRowSparse,
				 Allocator>& A,
				 int index = 0)
  {
    if (IndRow.GetM() <= 0)
      return;
    
    int i, j;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Number of elements per row.
    IVect Ptr(m);
    Ptr.Zero();
    for (i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }
    
    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateRow(i, Ptr(i));
	  for (j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}
    
    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }
  
  
  //! Conversion from coordinate format to ArrayColSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop, ArrayColSparse,
				 Allocator>& A,
				 int index = 0)
  {
    // Assuming there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;
    
    int i, j;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // Sorts array 'IndCol'.
    Sort(IndCol, IndRow, Val);

    // Construction of array 'Ptr'.
    IVect Ptr(n);
    Ptr.Zero();
    for (i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndCol(i))++;
      }
    
    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (i = 0; i < n; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndRow(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}
    
    // Assembles 'A' to sort row numbers.
    A.Assemble();
  }
  
  
  //! Conversion from coordinate format to ArrayRowSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop,
				 ArrayRowSymSparse, Allocator>& A,
				 int index = 0)
  {
    // Assuming that there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;
    
    int i, j;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;
    
    if (nb_low > 0)
      {
	int nb = 0;
	for (i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }
	
	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    IVect Ptr(m);
    Ptr.Zero();
    for (i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }
    
    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateRow(i, Ptr(i));
	  for (j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}
    
    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }
  
  
  //! Conversion from coordinate format to ArrayColSymSparse.
  template<class T, class Prop, class Allocator>
  void
  ConvertMatrix_from_Coordinates(IVect& IndRow, IVect& IndCol,
				 Vector<T, VectFull, Allocator>& Val,
				 Matrix<T, Prop,
				 ArrayColSymSparse, Allocator>& A,
				 int index = 0)
  {
    // Assuming that there is no duplicate value.
    if (IndRow.GetM() <= 0)
      return;
    
    int i, j;

    int row_max = IndRow.GetNormInf();
    int col_max = IndCol.GetNormInf();
    int m = row_max - index + 1;
    int n = col_max - index + 1;
    int nnz = IndRow.GetM();
    
    // First, removing the lower part of the matrix (if present).
    int nb_low = 0;
    for (i = 0; i < nnz; i++)
      if (IndRow(i) > IndCol(i))
	nb_low++;
    
    if (nb_low > 0)
      {
	int nb = 0;
	for (i = 0; i < nnz; i++)
	  if (IndRow(i) <= IndCol(i))
	    {
	      IndRow(nb) = IndRow(i);
	      IndCol(nb) = IndCol(i);
	      Val(nb) = Val(i);
	      nb++;
	    }
	
	IndRow.Resize(nb);
	IndCol.Resize(nb);
	Val.Resize(nb);
	nnz = nb;
      }

    // Sorts the array 'IndRow'.
    Sort(IndRow, IndCol, Val);

    // Construction of array 'Ptr'.
    IVect Ptr(m);
    Ptr.Zero();
    for (i = 0; i < nnz; i++)
      {
	IndRow(i) -= index;
	IndCol(i) -= index;
	Ptr(IndRow(i))++;
      }
    
    // Fills matrix 'A'.
    A.Reallocate(m, n);
    int offset = 0;
    for (i = 0; i < m; i++)
      if (Ptr(i) > 0)
	{
	  A.ReallocateColumn(i, Ptr(i));
	  for (j = 0; j < Ptr(i); j++)
	    {
	      A.Index(i, j) = IndCol(offset + j);
	      A.Value(i, j) = Val(offset + j);
	    }
	  offset += Ptr(i);
	}
    
    // Assembles 'A' to sort column numbers.
    A.Assemble();
  }
  
  
  /*
    From CSR to other CSR formats.
  */
  
  
  //! B = A.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, ColSparse, Alloc1>& A,
	    Matrix<T, Prop, ColSparse, Alloc2>& B)
  {
    int i;

    int m = A.GetM(), n = A.GetN(), nnz = A.GetDataSize();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    IVect Ptr(n+1), Ind(nnz);
    Vector<T, VectFull, Alloc2> Val(nnz);
    for (i = 0; i <= n; i++)
      Ptr(i) = ptr_[i];
    
    for (i = 0; i < nnz; i++)
      {
	Ind(i) = ind_[i];
	Val(i) = data_[i];
      }
    
    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop, RowSparse, Alloc1>& A,
	    Matrix<T, Prop, ColSparse, Alloc2>& B)
  {
    int i, j;

    int m = A.GetM(), n = A.GetN(), nnz = A.GetDataSize();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    // Computation of the indices of the beginning of columns.
    IVect Ptr(n + 1);
    Ptr.Fill(0);
    // Counting the number of entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(ind_[i])++;
    
    // Incrementing in order to get indices.
    int increment = 0, size, num_col;
    for (i = 0; i < n; i++)
      {
	size = Ptr(i);
	Ptr(i) = increment;
	increment += size;
      }
    // Last index.
    Ptr(n) = increment;
    
    // 'Offset' will be used to get current positions of new entries.
    IVect Offset = Ptr;
    IVect Ind(nnz);
    Vector<T, VectFull, Alloc2> Val(nnz);
    
    // Loop over rows.
    for (i = 0; i < m; i++)
      for (j = ptr_[i]; j < ptr_[i + 1]; j++)
	{
	  num_col = ind_[j];
	  Ind(Offset(num_col)) = i;
	  Val(Offset(num_col)) = data_[j];
	  Offset(num_col)++;
	}
    
    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  //! Conversion from row-sparse to column-sparse.
  template<class T, class Prop1, class Prop2,
	   class Storage, class Alloc1, class Alloc2>
  void
  ConvertMatrixSymSparse_to_ColSparse(const Matrix<T, Prop1,
				      Storage, Alloc1>& A,
				      Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    int i, j;

    int m = A.GetM(), n = A.GetN();
    int* ptr_ = A.GetPtr();
    int* ind_ = A.GetInd();
    T* data_ = A.GetData();
    
    // Computation of the indices of the beginning of columns.
    IVect Ptr(n + 1);
    Ptr.Fill(0);
    // Counting the number of entries per column.
    for (i = 0; i < m; i++)
      for (j = ptr_[i]; j < ptr_[i + 1]; j++)
	{
	  Ptr(ind_[j])++;
	  if (i != ind_[j])
	    Ptr(i)++;
	}
    
    // Incrementing in order to get indices.
    int nnz = 0, size, num_col;
    for (i = 0; i < n; i++)
      {
	size = Ptr(i);
	Ptr(i) = nnz;
	nnz += size;
      }
    // Last index.
    Ptr(n) = nnz;
    
    // 'Offset' will be used to get current positions of new entries.
    IVect Offset = Ptr;
    IVect Ind(nnz);
    Vector<T, VectFull, Alloc2> Val(nnz);
    
    // Loop over rows.
    for (i = 0; i < m; i++)
      for (j = ptr_[i]; j < ptr_[i + 1]; j++)
	{
	  num_col = ind_[j];
	  Ind(Offset(num_col)) = i;
	  Val(Offset(num_col)) = data_[j];
	  Offset(num_col)++;
	  if (i != ind_[j])
	    {
	      Ind(Offset(i)) = num_col;
	      Val(Offset(i)) = data_[j];
	      Offset(i)++;
	    }
	}
    
    B.SetData(m, n, Val, Ptr, Ind);
  }
  
  
  //! Conversion from RowSymSparse to column-sparse.
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, RowSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    ConvertMatrixSymSparse_to_ColSparse(A, B);
  }
  
  
  //! Conversion from ColSymSparse to column-sparse.
  template<class T, class Prop1, class Prop2, class Alloc1, class Alloc2>
  void Copy(const Matrix<T, Prop1, ColSymSparse, Alloc1>& A,
	    Matrix<T, Prop2, ColSparse, Alloc2>& B)
  {
    ConvertMatrixSymSparse_to_ColSparse(A, B);
  }
  
  
  /*
    From ArraySparse matrices to CSR matrices.
  */
  
  
  //! Conversion from ArrayRowSparse to RowSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Matrix (m,n) with 'nnz' entries.
    int nnz = mat_array.GetDataSize();
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    
    // Allocating arrays needed for CSR format.
    Vector<T1> Val(nnz);
    IVect IndRow(m + 1);
    IVect IndCol(nnz);
    
    // Filling the arrays.
    int ind = 0;
    IndRow(0) = 0;
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRowSize(i); k++)
	  {
	    IndCol(ind) = mat_array.Index(i, k);
	    Val(ind) = mat_array.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
    
    mat_csr.SetData(m, n, Val, IndRow, IndCol);
  }
  
  
  //! Conversion from ArrayRowSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csr)
  {
    int i;

    // Matrix (m,n) with nnz entries.
    int nnz = mat_array.GetDataSize();
    int m = mat_array.GetM();
    int n = mat_array.GetN();
   
    // Conversion in coordinate format.
    Vector<T1> Val;
    IVect IndRow, IndCol;
    ConvertMatrix_to_Coordinates(mat_array, IndRow, IndCol, Val);
    
    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);
    
    // Constructing pointer array 'Ptr'.
    IVect Ptr(n + 1);
    
    // Counting non-zero entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;
    
    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
    mat_csr.SetData(m, n, Val, Ptr, IndRow);
  }
  
  
  //! Conversion from ArrayRowSymSparse to RowSymSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse,Allocator0>& mat_array,
	    Matrix<T1, Prop1, RowSymSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Number of rows and non-zero entries.
    int nnz = mat_array.GetDataSize();
    int m = mat_array.GetM();
    
    // Allocation of arrays for CSR format.
    Vector<T1> Val(nnz);
    Vector<int, VectFull, CallocAlloc<int> > IndRow(m + 1);
    Vector<int, VectFull, CallocAlloc<int> > IndCol(nnz);
    
    int ind = 0;
    IndRow(0) = 0;
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRowSize(i); k++)
	  {
	    IndCol(ind) = mat_array.Index(i, k);
	    Val(ind) = mat_array.Value(i, k);
	    ind++;
	  }
	IndRow(i + 1) = ind;
      }
    
    mat_csr.SetData(m, m, Val, IndRow, IndCol);
  }
  
  
  //! Conversion from ArrayRowComplexSparse to RowComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0, ArrayRowComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();
    
    // Allocation of arrays for CSR.
    Vector<T0> Val_real(nnz_real), Val_imag(nnz_imag);
    IVect IndRow_real(m + 1), IndRow_imag(m + 1);
    IVect IndCol_real(nnz_real), IndCol_imag(nnz_imag);
    
    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }
	
	IndRow_real(i + 1) = ind_real;
	for (k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }
	
	IndRow_imag(i + 1) = ind_imag;
      }
    
    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }
  
  
  //! Conversion from ArrayRowSymComplexSparse to RowSymComplexSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void
  Copy(const Matrix<T0, Prop0,
       ArrayRowSymComplexSparse, Allocator0>& mat_array,
       Matrix<T1, Prop1, RowSymComplexSparse, Allocator1>& mat_csr)
  {
    int i, k;

    // Non-zero entries (real and imaginary part).
    int nnz_real = mat_array.GetRealDataSize();
    int nnz_imag = mat_array.GetImagDataSize();
    int m = mat_array.GetM();
    
    // Allocation of arrays for CSR.
    Vector<T0> Val_real(nnz_real), Val_imag(nnz_imag);
    IVect IndRow_real(m + 1), IndRow_imag(m + 1);
    IVect IndCol_real(nnz_real), IndCol_imag(nnz_imag);
    
    int ind_real = 0, ind_imag = 0;
    IndRow_real(0) = 0;
    IndRow_imag(0) = 0;
    // Loop over rows.
    for (i = 0; i < m; i++)
      {
	for (k = 0; k < mat_array.GetRealRowSize(i); k++)
	  {
	    IndCol_real(ind_real) = mat_array.IndexReal(i, k);
	    Val_real(ind_real) = mat_array.ValueReal(i, k);
	    ind_real++;
	  }
	
	IndRow_real(i + 1) = ind_real;
	for (int k = 0; k < mat_array.GetImagRowSize(i); k++)
	  {
	    IndCol_imag(ind_imag) = mat_array.IndexImag(i, k);
	    Val_imag(ind_imag) = mat_array.ValueImag(i, k);
	    ind_imag++;
	  }
	
	IndRow_imag(i + 1) = ind_imag;
      }
    
    mat_csr.SetData(m, m, Val_real, IndRow_real, IndCol_real,
		    Val_imag, IndRow_imag, IndCol_imag);
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayRowSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    IVect IndRow(nnz),IndCol(nnz);
    Vector<T1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;
    
    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeRow(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleRow(i);
      }
  }
  
  
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSymSparse, Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i, j;

    int nnz = A.GetDataSize();
    int n = A.GetM();
    IVect IndRow(nnz), IndCol(nnz);
    Vector<T1> Val(nnz);
    int ind = 0;
    for (i = 0; i < n; i++)
      for (j = 0; j < A.GetRowSize(i); j++)
	if (A.Index(i, j) != i)
	  {
	    IndRow(ind) = i;
	    IndCol(ind) = A.Index(i, j);
	    Val(ind) = A.Value(i, j);
	    ind++;
	  }
    Sort(ind, IndCol, IndRow, Val);
    nnz = ind;
    ind = 0;
    
    B.Reallocate(n, n);
    for (i = 0; i < n; i++)
      {
	int first_index = ind;
	while (ind < nnz && IndCol(ind) <= i)
	  ind++;
	int size_lower = ind - first_index;
	int size_upper = A.GetRowSize(i);
	int size_row = size_lower + size_upper;
	B.ResizeColumn(i, size_row);
	ind = first_index;
	for (j = 0; j < size_lower; j++)
	  {
	    B.Index(i, j) = IndRow(ind);
	    B.Value(i, j) = Val(ind);
	    ind++;
	  }
	for (j = 0; j < size_upper; j++)
	  {
	    B.Index(i, size_lower + j) = A.Index(i, j);
	    B.Value(i, size_lower + j) = A.Value(i, j);
	  }
	B.AssembleColumn(i);
      }
  }
  
  
  //! Conversion from ArrayColSparse to ColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayColSparse, Allocator0>& mat_array,
	    Matrix<T1, Prop1, ColSparse, Allocator1>& mat_csc)
  {
    int i, k;

    // Matrix (m,n) with 'nnz' entries.
    int nnz = mat_array.GetDataSize();
    int m = mat_array.GetM();
    int n = mat_array.GetN();
    
    // Allocating arrays needed for CSC format.
    Vector<T1> Val(nnz);
    IVect IndRow(nnz);
    IVect IndCol(n+1);
    
    // Filling the arrays.
    int ind = 0;
    IndCol(0) = 0;
    for (i = 0; i < n; i++)
      {
	for (k = 0; k < mat_array.GetColumnSize(i); k++)
	  {
	    IndRow(ind) = mat_array.Index(i, k);
	    Val(ind) = mat_array.Value(i, k);
	    ind++;
	  }
	IndCol(i + 1) = ind;
      }
    
    mat_csc.SetData(m, n, Val, IndCol, IndRow);
  }
  
  
  //! Conversion from ArrayRowSparse to ArrayColSparse.
  template<class T0, class Prop0, class Allocator0,
	   class T1, class Prop1, class Allocator1>
  void Copy(const Matrix<T0, Prop0, ArrayRowSparse,Allocator0>& A,
	    Matrix<T1, Prop1, ArrayColSparse, Allocator1>& B)
  {
    int i;
    
    // Matrix (m,n) with nnz entries.
    int nnz = A.GetDataSize();
    int n = A.GetN();
   
    // Conversion in coordinate format.
    Vector<T1> Val;
    IVect IndRow, IndCol;
    ConvertMatrix_to_Coordinates(A, IndRow, IndCol, Val);
    
    // Sorting with respect to column numbers.
    Sort(IndCol, IndRow, Val);
    
    // Constructing pointer array 'Ptr'.
    IVect Ptr(n + 1);
    
    // Counting non-zero entries per column.
    for (i = 0; i < nnz; i++)
      Ptr(IndCol(i) + 1)++;
    
    // Accumulation to get pointer array.
    Ptr(0) = 0;
    for (i = 0; i < n; i++)
      Ptr(i + 1) += Ptr(i);
    
    // we fill matrix B
    for (int i = 0; i < n; i++)
      {
	int size_col = Ptr(i+1) - Ptr(i);
	if (size_col > 0)
	  {
	    B.ReallocateColumn(i, size_col);
	    for (int j = Ptr(i); j < Ptr(i+1); j++)
	      {
		B.Index(i, j-Ptr(i)) = IndRow(j);
		B.Value(i, j-Ptr(i)) = Val(j);
	      }
	  }
      }
  }
  
} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_MATRIX_CXX
#endif
