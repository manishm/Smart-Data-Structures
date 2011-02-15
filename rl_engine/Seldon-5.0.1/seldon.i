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


%module seldon
%{
#include "SeldonHeader.hxx"
  %}

%include "std_string.i"

using namespace std;

%include "share/Errors.hxx"
%exception
{
  try
    {
      $action
	}
  catch(Seldon::Error& e)
    {
      PyErr_SetString(PyExc_Exception, e.What().c_str());
      return NULL;
    }
  catch(std::exception& e)
    {
      PyErr_SetString(PyExc_Exception, e.what());
      return NULL;
    }
  catch(std::string& s)
    {
      PyErr_SetString(PyExc_Exception, s.c_str());
      return NULL;
    }
  catch(const char* s)
    {
      PyErr_SetString(PyExc_Exception, s);
      return NULL;
    }
  catch(...)
    {
      PyErr_SetString(PyExc_Exception, "Unknown exception...");
      return NULL;
    }
}

%include "SeldonHeader.hxx"
%include "share/Common.hxx"
%include "share/Storage.hxx"
%include "share/Properties.hxx"
%include "vector/Vector.hxx"
%include "matrix/Matrix_Base.hxx"
%include "matrix/Matrix_Pointers.hxx"
%include "share/Allocator.hxx"

namespace Seldon
{
  %extend Vector<int, VectFull, MallocAlloc<int> >
  {
    int __getitem__(int index) {
      return (*self)(index);
    }
    void __setitem__(int index, int value) {
      (*self)(index) = value;
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }
  %extend Vector<double, VectFull, MallocAlloc<double> >
  {
    double __getitem__(int index) {
      return (*self)(index);
    }
    void __setitem__(int index, double value) {
      (*self)(index) = value;
    }
    unsigned long __len__() {
      return self->GetM();
    }
  }

  %extend Matrix<int, General, RowMajor, MallocAlloc<int> >
  {
    int __getitem__(PyObject* args)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j);
    }
    Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > __getitem__(int i)
    {
      Seldon::Vector<int, Seldon::VectFull, Seldon::MallocAlloc<int> > v(self->GetN());
      for (int j = 0; j < self->GetN(); j++)
	v(j) = (*self)(i, j);
      return v;
    }
    void __setitem__(PyObject* args, int value)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j) = value;
    }
    unsigned long __len__()
    {
      return self->GetM();
    }
  }
  %extend Matrix<double, General, RowMajor, MallocAlloc<double> >
  {
    double __getitem__(PyObject* args)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      return (*self)(i, j);
    }
    Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > __getitem__(int i)
    {
      Seldon::Vector<double, Seldon::VectFull, Seldon::MallocAlloc<double> > v(self->GetN());
      for (int j = 0; j < self->GetN(); j++)
	v(j) = (*self)(i, j);
      return v;
    }
    void __setitem__(PyObject* args, double value)
    {
      int i, j;
      int success = PyArg_ParseTuple(args, "ii", &i, &j);
      if (!success)
	throw std::out_of_range("Failed!");
      (*self)(i, j) = value;
    }
    unsigned long __len__()
    {
      return self->GetM();
    }
  }

  %template(IntMalloc) MallocAlloc<int>;
  %template(BaseSeldonVectorInt) Vector_Base<int, MallocAlloc<int> >;
  %template(VectorInt) Vector<int, VectFull, MallocAlloc<int> >;
  %template(DoubleMalloc) MallocAlloc<double>;
  %template(BaseSeldonVectorDouble) Vector_Base<double, MallocAlloc<double> >;
  %template(VectorDouble) Vector<double, VectFull, MallocAlloc<double> >;

  %template(MatrixBaseInt) Matrix_Base<int, MallocAlloc<int> >;
  %template(MatrixPointersInt) Matrix_Pointers<int, General, RowMajor, MallocAlloc<int> >;
  %template(MatrixInt) Matrix<int, General, RowMajor, MallocAlloc<int> >;
  %template(MatrixBaseDouble) Matrix_Base<double, MallocAlloc<double> >;
  %template(MatrixPointersDouble) Matrix_Pointers<double, General, RowMajor, MallocAlloc<double> >;
  %template(MatrixDouble) Matrix<double, General, RowMajor, MallocAlloc<double> >;
}
