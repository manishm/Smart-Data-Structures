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


#ifndef SELDON_FILE_ERRORS_HXX

namespace Seldon
{

  
  ///////////
  // ERROR //
  ///////////

  class Error
  {
  protected:
    //! Message describing the exception type.
    string description_;
    //! Name of the function in which the error occurred.
    string function_;
    //! A comment about the error.
    string comment_;

  public:
    Error(string function, string comment)  throw();
    Error(string description, string function, string comment)  throw();
    virtual ~Error()  throw();

    virtual string What();
    void CoutWhat();
  };

  
  ///////////////////
  // WRONGARGUMENT //
  ///////////////////

  class WrongArgument: public Error
  {
  public:
    WrongArgument(string function, string comment)  throw();

    virtual string What();
  };
  

  //////////////
  // NOMEMORY //
  //////////////

  class NoMemory: public Error
  {
  public:
    NoMemory(string function, string comment)  throw();
  };
  

  //////////////
  // WRONGDIM //
  //////////////

  class WrongDim: public Error
  {
  public:
    WrongDim(string function, string comment)  throw();
  };
  

  ////////////////
  // WRONGINDEX //
  ////////////////

  class WrongIndex: public Error
  {
  public:
    WrongIndex(string function, string comment)  throw();
  };
  

  //////////////
  // WRONGROW //
  //////////////

  class WrongRow: public Error
  {
  public:
    WrongRow(string function, string comment)  throw();
  };
  

  //////////////
  // WRONGCOL //
  //////////////

  class WrongCol: public Error
  {
  public:
    WrongCol(string function, string comment)  throw();
  };
  

  /////////////
  // IOERROR //
  /////////////

  class IOError: public Error
  {
  public:
    IOError(string function, string comment)  throw();
  };
  

  /////////////////
  // LAPACKERROR //
  /////////////////

  class LapackError: public Error
  {
  protected:
    int info_;

  public:
    LapackError(int info, string function, string comment)  throw();

    virtual string What();
  };
  

} // namespace Seldon.

#define SELDON_FILE_ERRORS_HXX
#endif
