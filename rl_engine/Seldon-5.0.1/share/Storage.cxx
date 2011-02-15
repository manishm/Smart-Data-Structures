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


#ifndef SELDON_FILE_STORAGE_CXX

#include "Storage.hxx"

namespace Seldon
{


  //////////////////////
  // GENERAL MATRICES //
  //////////////////////


  int ColMajor::GetFirst(int i, int j)
  {
    return j;
  }
  int ColMajor::GetSecond(int i, int j)
  {
    return i;
  }


  int RowMajor::GetFirst(int i, int j)
  {
    return i;
  }
  int RowMajor::GetSecond(int i, int j)
  {
    return j;
  }



  /////////////
  // VECTORS //
  /////////////


  class VectFull
  {
  };
  
  
  class VectSparse
  {
  };


  ////////////
  // SPARSE //
  ////////////


  int ColSparse::GetFirst(int i, int j)
  {
    return j;
  }
  int ColSparse::GetSecond(int i, int j)
  {
    return i;
  }


  int RowSparse::GetFirst(int i, int j)
  {
    return i;
  }
  int RowSparse::GetSecond(int i, int j)
  {
    return j;
  }


  int ColComplexSparse::GetFirst(int i, int j)
  {
    return j;
  }
  int ColComplexSparse::GetSecond(int i, int j)
  {
    return i;
  }


  int RowComplexSparse::GetFirst(int i, int j)
  {
    return i;
  }
  int RowComplexSparse::GetSecond(int i, int j)
  {
    return j;
  }


  int ColSymSparse::GetFirst(int i, int j)
  {
    return j;
  }
  int ColSymSparse::GetSecond(int i, int j)
  {
    return i;
  }


  int RowSymSparse::GetFirst(int i, int j)
  {
    return i;
  }
  int RowSymSparse::GetSecond(int i, int j)
  {
    return j;
  }


  int ColSymComplexSparse::GetFirst(int i, int j)
  {
    return j;
  }
  int ColSymComplexSparse::GetSecond(int i, int j)
  {
    return i;
  }


  int RowSymComplexSparse::GetFirst(int i, int j)
  {
    return i;
  }
  int RowSymComplexSparse::GetSecond(int i, int j)
  {
    return j;
  }



  ///////////////
  // SYMMETRIC //
  ///////////////


  int ColSymPacked::GetFirst(int i, int j)
  {
    return j;
  }
  int ColSymPacked::GetSecond(int i, int j)
  {
    return i;
  }


  int RowSymPacked::GetFirst(int i, int j)
  {
    return i;
  }
  int RowSymPacked::GetSecond(int i, int j)
  {
    return j;
  }


  int ColSym::GetFirst(int i, int j)
  {
    return j;
  }
  int ColSym::GetSecond(int i, int j)
  {
    return i;
  }


  int RowSym::GetFirst(int i, int j)
  {
    return i;
  }
  int RowSym::GetSecond(int i, int j)
  {
    return j;
  }



  ///////////////
  // HERMITIAN //
  ///////////////


  int ColHerm::GetFirst(int i, int j)
  {
    return j;
  }
  int ColHerm::GetSecond(int i, int j)
  {
    return i;
  }


  int RowHerm::GetFirst(int i, int j)
  {
    return i;
  }
  int RowHerm::GetSecond(int i, int j)
  {
    return j;
  }


  int ColHermPacked::GetFirst(int i, int j)
  {
    return j;
  }
  int ColHermPacked::GetSecond(int i, int j)
  {
    return i;
  }


  int RowHermPacked::GetFirst(int i, int j)
  {
    return i;
  }
  int RowHermPacked::GetSecond(int i, int j)
  {
    return j;
  }



  ////////////////
  // TRIANGULAR //
  ////////////////


  int ColUpTriang::GetFirst(int i, int j)
  {
    return j;
  }
  int ColUpTriang::GetSecond(int i, int j)
  {
    return i;
  }
  bool ColUpTriang::UpLo()
  {
    return true;
  }


  int ColLoTriang::GetFirst(int i, int j)
  {
    return j;
  }
  int ColLoTriang::GetSecond(int i, int j)
  {
    return i;
  }
  bool ColLoTriang::UpLo()
  {
    return false;
  }


  int RowUpTriang::GetFirst(int i, int j)
  {
    return i;
  }
  int RowUpTriang::GetSecond(int i, int j)
  {
    return j;
  }
  bool RowUpTriang::UpLo()
  {
    return true;
  }


  int RowLoTriang::GetFirst(int i, int j)
  {
    return i;
  }
  int RowLoTriang::GetSecond(int i, int j)
  {
    return j;
  }
  bool RowLoTriang::UpLo()
  {
    return false;
  }


  int ColUpTriangPacked::GetFirst(int i, int j)
  {
    return j;
  }
  int ColUpTriangPacked::GetSecond(int i, int j)
  {
    return i;
  }
  bool ColUpTriangPacked::UpLo()
  {
    return true;
  }


  int ColLoTriangPacked::GetFirst(int i, int j)
  {
    return j;
  }
  int ColLoTriangPacked::GetSecond(int i, int j)
  {
    return i;
  }
  bool ColLoTriangPacked::UpLo()
  {
    return false;
  }


  int RowUpTriangPacked::GetFirst(int i, int j)
  {
    return i;
  }
  int RowUpTriangPacked::GetSecond(int i, int j)
  {
    return j;
  }
  bool RowUpTriangPacked::UpLo()
  {
    return true;
  }


  int RowLoTriangPacked::GetFirst(int i, int j)
  {
    return i;
  }
  int RowLoTriangPacked::GetSecond(int i, int j)
  {
    return j;
  }
  bool RowLoTriangPacked::UpLo()
  {
    return false;
  }


} // namespace Seldon.

#define SELDON_FILE_STORAGE_CXX
#endif
