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


#ifndef SELDON_FILE_FUNCTIONS_CXX

#include "Functions.hxx"

namespace Seldon
{


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetRow(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int i, Vector<T1, Storage1, Allocator1>& X)
  {
    X.Reallocate(M.GetN());
    for (int j = 0; j < M.GetN(); j++)
      X(j) = M(i, j);
  }


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void GetCol(const Matrix<T0, Prop0, Storage0, Allocator0>& M,
	      int j, Vector<T1, Storage1, Allocator1>& X)
  {
    X.Reallocate(M.GetM());
    for (int i = 0; i < M.GetM(); i++)
      X(i) = M(i, j);
  }


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetRow(const Vector<T1, Storage1, Allocator1>& X,
	      int i, Matrix<T0, Prop0, Storage0, Allocator0>& M)
  {
    for (int j = 0; j < M.GetN(); j++)
      M(i, j) = X(j);
  }


  template <class T0, class Prop0, class Storage0, class Allocator0,
	    class T1, class Storage1, class Allocator1>
  void SetCol(const Vector<T1, Storage1, Allocator1>& X,
	      int j, Matrix<T0, Prop0, Storage0, Allocator0>& M)
  {
    for (int i = 0; i < M.GetM(); i++)
      M(i, j) = X(i);
  }


} // namespace Seldon.

#define SELDON_FILE_FUNCTIONS_CXX
#endif
