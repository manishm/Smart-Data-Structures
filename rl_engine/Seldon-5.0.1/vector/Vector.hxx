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


#ifndef SELDON_FILE_VECTOR_HXX

#include "../share/Common.hxx"
#include "../share/Properties.hxx"
#include "../share/Storage.hxx"
#include "../share/Errors.hxx"
#include "../share/Allocator.hxx"

namespace Seldon
{


  //! Base structure for all vectors.
  /*!
    It stores data and the vector size. It defines basic
    methods as well.
  */
  template <class T, class Allocator = SELDON_DEFAULT_ALLOCATOR<T> >
  class Vector_Base
  {
    // typdef declarations.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    // Static attributes.
  protected:
    static Allocator vect_allocator_;

    // Attributes.
  protected:
    // Number of elements.
    int m_;
    // Pointer to stored elements.
    pointer data_;

    // Methods.
  public:
    // Constructors.
    Vector_Base();
    explicit Vector_Base(int i);
    Vector_Base(const Vector_Base<T, Allocator>& A);
  
    // Destructor.
    ~Vector_Base();

    // Basic methods.
    int GetM() const;
    int GetLength() const;
    int GetSize() const;
    pointer GetData() const;
    const_pointer GetDataConst() const;
    void* GetDataVoid() const;
    const void* GetDataConstVoid() const;

  };


  // Vector allocator.
  template <class T, class Allocator>
  Allocator Vector_Base<T, Allocator>::vect_allocator_;


  //! Full vector class.
  /*!
    Basic vector class (i.e. not sparse).
  */
  template <class T, class Allocator>
  class Vector<T, VectFull, Allocator>: public Vector_Base<T, Allocator>
  {
    // typedef declaration.
  public:
    typedef typename Allocator::value_type value_type;
    typedef typename Allocator::pointer pointer;
    typedef typename Allocator::const_pointer const_pointer;
    typedef typename Allocator::reference reference;
    typedef typename Allocator::const_reference const_reference;

    // Attributes.
  private:

    // Methods.
  public:
    // Constructor.
    explicit Vector()  throw();
    explicit Vector(int i);
    Vector(const Vector<T, VectFull, Allocator>& A);
    
    // Destructor.
    ~Vector();
    void Clear();

    // Memory management.
    void Reallocate(int i);
    void Resize(int i);
    void SetData(int i, pointer data);
    void Nullify();

    // Element access and affectation.
    reference operator() (int i);
#ifndef SWIG
    const_reference operator() (int i) const;
    Vector<T, VectFull, Allocator>& operator= (const Vector<T, VectFull,
					       Allocator>& X);
#endif
    void Copy(const Vector<T, VectFull, Allocator>& X);
    void Append(const T& x);
    template<class T0>
    void PushBack(const T0& x);
    template<class Allocator0>
    void PushBack(const Vector<T, VectFull, Allocator0>& X);
    
    // Basic functions.
    int GetDataSize();

    // Convenient functions.
    void Zero();
    void Fill();
    template <class T0>
    void Fill(const T0& x);
    template <class T0>
#ifndef SWIG
    Vector<T, VectFull, Allocator>& operator= (const T0& X);
#endif
    template <class T0>
    Vector<T, VectFull, Allocator>& operator*= (const T0& X);
    void FillRand();
    void Print() const;

#ifndef SELDON_EXCLUDE_FROM_EXPLICIT_INSTANTIATION
#ifndef SWIG
    // Norms.
    value_type GetNormInf() const;
    int GetNormInfIndex() const;
#endif
#endif

    // Input/output functions.
    void Write(string FileName, bool with_size) const;
#ifndef SWIG
    void Write(ostream& FileStream, bool with_size) const;
#endif
    void WriteText(string FileName) const;
#ifndef SWIG
    void WriteText(ostream& FileStream) const;
#endif
    void Read(string FileName);
#ifndef SWIG
    void Read(istream& FileStream);
#endif
#ifndef SELDON_EXCLUDE_FROM_EXPLICIT_INSTANTIATION
#ifndef SWIG
    void ReadText(string FileName);
    void ReadText(istream& FileStream);
#endif
#endif
    
  };

#ifndef SWIG
  template <class T, class Storage, class Allocator>
  ostream& operator << (ostream& out,
			const Vector<T, Storage, Allocator>& V);
#endif


} // namespace Seldon.

#define SELDON_FILE_VECTOR_HXX
#endif
