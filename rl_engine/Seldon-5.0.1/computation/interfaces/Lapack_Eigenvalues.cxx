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


#ifndef SELDON_FILE_LAPACK_EIGENVALUES_CXX

/*
  Functions included in this file:

  xGEEV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xSYEV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xHEEV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xSPEV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xHPEV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xSYGV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xGGEV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xHEGV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xSPGV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xHPGV   (GetEigenvalues, GetEigenvaluesEigenvectors)
  xGESVD  (GetSVD)
  xGEQRF  (GetHessian)
  ZGEQRF + ZUNGQR + ZUNMQR + ZGGHRD   (GetHessian)
  ZGEQRF + ZUNGQR + ZUNMQR + ZGGHRD + ZHGEQZ   (GetQZ)
  (SolveSylvester)
*/

namespace Seldon
{
  
    
  /////////////////////////////////
  // STANDARD EIGENVALUE PROBLEM //
  
  
  /* RowMajor */
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<float, Prop, RowMajor, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& wr,
		      Vector<float, VectFull, Allocator3>& wi,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('N');
    char jobvr('N');
    Vector<float> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    sgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   A.GetData(), &n, A.GetData(), &n, work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop, RowMajor, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& wr,
				  Vector<float, VectFull, Allocator3>& wi,
				  Matrix<float, General, RowMajor, Allocator4>& zr,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('V');
    char jobvr('N');
    Vector<float> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    zr.Reallocate(n, n);
    sgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   zr.GetData(), &n, zr.GetData(), &n, work.GetData(), &lwork,
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(zr);
    // conjugate if necessary
    int i = 0;
    while (i < n)
      {
	if (i < (n-1))
	  if (wr(i) == wr(i+1))
	    {
	      for (int j = 0; j < n; j++)
		zr(j,i+1) = -zr(j,i+1);
	      
	      i++;
	    }
	
	i++;
      }
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>, Prop, RowMajor, Allocator1>& A,
		      Vector<complex<float>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('N'), jobr('N'); int lwork = 3*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(2*n);
    w.Reallocate(n);
    cgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   A.GetDataVoid(), &n, A.GetData(), &n, work.GetDataVoid(), &lwork,
	   rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, RowMajor, Allocator1>& A,
				  Vector<complex<float>, VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('V'), jobr('N'); int lwork = 3*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(2*n);
    w.Reallocate(n);
    z.Reallocate(n, n);
    cgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   z.GetDataVoid(), &n, z.GetData(), &n, work.GetDataVoid(), &lwork,
	   rwork.GetData(), &info.GetInfoRef());
    
#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

    TransposeConj(z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<double, Prop, RowMajor, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& wr,
		      Vector<double, VectFull, Allocator3>& wi,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('N'), jobvr('N'); Vector<double> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    dgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   A.GetData(), &n, A.GetData(), &n, work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop, RowMajor, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& wr,
				  Vector<double, VectFull, Allocator3>& wi,
				  Matrix<double, General, RowMajor, Allocator4>& zr,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('V'), jobvr('N'); Vector<double> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    zr.Reallocate(n, n);
    dgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   zr.GetData(), &n, zr.GetData(), &n, work.GetData(), &lwork,
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(zr);
    // conjugate if necessary
    int i = 0;
    while (i < n)
      {
	if (i < (n-1))
	  if (wr(i) == wr(i+1))
	    {
	      for (int j = 0; j < n; j++)
		zr(j,i+1) = -zr(j,i+1);
	      
	      i++;
	    }
	
	i++;
      }
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>, Prop, RowMajor, Allocator1>& A,
		      Vector<complex<double>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('N'), jobr('N'); int lwork = 3*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(2*n);
    w.Reallocate(n);
    zgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   A.GetDataVoid(), &n, A.GetData(), &n, work.GetDataVoid(), &lwork,
	   rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, RowMajor, Allocator1>& A,
				  Vector<complex<double>,
				  VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('V'), jobr('N'); int lwork = 3*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(2*n);
    w.Reallocate(n);
    z.Reallocate(n, n);
    zgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   z.GetDataVoid(), &n, z.GetData(), &n, work.GetDataVoid(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    TransposeConj(z);
  }
  
  
  /* ColMajor */
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<float, Prop, ColMajor, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& wr,
		      Vector<float, VectFull, Allocator3>& wi,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('N'), jobvr('N'); Vector<float> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    sgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   A.GetData(), &n, A.GetData(), &n, work.GetData(), &lwork,
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  template<class Prop, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop, ColMajor, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& wr,
				  Vector<float, VectFull, Allocator3>& wi,
				  Matrix<float, General, ColMajor, Allocator4>&zr,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('N'), jobvr('V');
    Vector<float> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    zr.Reallocate(n, n);
    sgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   zr.GetData(), &n, zr.GetData(), &n, work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>, Prop, ColMajor, Allocator1>& A,
		      Vector<complex<float>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('N'), jobr('N'); int lwork = 3*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(2*n);
    w.Reallocate(n);
    cgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   A.GetDataVoid(), &n, A.GetData(), &n, work.GetDataVoid(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, ColMajor, Allocator1>& A,
				  Vector<complex<float>,
				  VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('N'), jobr('V'); int lwork = 3*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(2*n);
    w.Reallocate(n);
    z.Reallocate(n, n);
    cgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   z.GetDataVoid(), &n, z.GetData(), &n, work.GetDataVoid(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<double, Prop, ColMajor, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& wr,
		      Vector<double, VectFull, Allocator3>& wi,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('N'), jobvr('N'); Vector<double> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    dgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   A.GetData(), &n, A.GetData(), &n, work.GetData(), &lwork,
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop, ColMajor, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& wr,
				  Vector<double, VectFull, Allocator3>& wi,
				  Matrix<double, General, ColMajor, Allocator4>&zr,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM(), lwork = 6*n;
    char jobvl('N'), jobvr('V');
    Vector<double> work(lwork);
    wr.Reallocate(n);
    wi.Reallocate(n);
    zr.Reallocate(n, n);
    dgeev_(&jobvl, &jobvr, &n, A.GetData(), &n, wr.GetData(), wi.GetData(),
	   zr.GetData(), &n, zr.GetData(), &n, work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>, Prop, ColMajor, Allocator1>& A,
		      Vector<complex<double>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('N'), jobr('N'); int lwork = 3*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(2*n);
    w.Reallocate(n);
    zgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   A.GetDataVoid(), &n, A.GetData(), &n, work.GetDataVoid(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, ColMajor, Allocator1>& A,
				  Vector<complex<double>,
				  VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobl('N'), jobr('V'); int lwork = 3*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(2*n);
    w.Reallocate(n);
    z.Reallocate(n, n);
    zgeev_(&jobl, &jobr, &n, A.GetDataVoid(), &n, w.GetDataVoid(),
	   z.GetDataVoid(), &n, z.GetData(), &n, work.GetDataVoid(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  /* RowSym */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<float, Prop, RowSym, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    ssyev_(&job, &uplo, &n, A.GetData(), &n, w.GetData(), work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop, RowSym, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<float, General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
 
    ssyev_(&job, &uplo, &n, z.GetData(), &n, w.GetData(), work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(z);
  }
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>, Prop, RowSym, Allocator1>& A,
		      Vector<complex<float>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    w.Reallocate(n);
    Matrix<complex<float>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    GetEigenvalues(B, w);
  }
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, RowSym, Allocator1>& A,
				  Vector<complex<float>,
				  VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    w.Reallocate(n);
    z.Reallocate(n, n);
    Matrix<complex<float>, General, RowMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<double, Prop, RowSym, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    dsyev_(&job, &uplo, &n, A.GetData(), &n, w.GetData(), work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop, RowSym, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<double, General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
 
    dsyev_(&job, &uplo, &n, z.GetData(), &n, w.GetData(), work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(z);
  }
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>, Prop, RowSym, Allocator1>& A,
		      Vector<complex<double>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    w.Reallocate(n);
    Matrix<complex<double>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    GetEigenvalues(B, w);
  }
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, RowSym, Allocator1>& A,
				  Vector<complex<double>,
				  VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, RowMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    z.Reallocate(n, n);
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  /* ColSym */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<float, Prop, ColSym, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('N'); int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    ssyev_(&job, &uplo, &n, A.GetData(), &n, w.GetData(), work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop, ColSym, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<float, General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    ssyev_(&job, &uplo, &n, z.GetData(), &n, w.GetData(),
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>, Prop, ColSym, Allocator1>& A,
		      Vector<complex<float>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    w.Reallocate(n);
    Matrix<complex<float>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    GetEigenvalues(B, w);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, ColSym, Allocator1>& A,
				  Vector<complex<float>,
				  VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    z.Reallocate(n, n);
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<double, Prop, ColSym, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('N'); int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    dsyev_(&job, &uplo, &n, A.GetData(), &n, w.GetData(), work.GetData(),
	   &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop, ColSym, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<double, General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    dsyev_(&job, &uplo, &n, z.GetData(), &n, w.GetData(),
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>, Prop, ColSym, Allocator1>& A,
		      Vector<complex<double>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> B(n,n);
    w.Reallocate(n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    GetEigenvalues(B, w);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, ColSym, Allocator1>& A,
				  Vector<complex<double>,
				  VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    z.Reallocate(n, n);
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  /* RowHerm */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>, Prop, RowHerm, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 2*n; Vector<complex<float> > work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    cheev_(&job, &uplo, &n, A.GetDataVoid(), &n, w.GetData(),
	   work.GetDataVoid(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, RowHerm, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 2*n; Vector<complex<float> > work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    cheev_(&job, &uplo,&n, z.GetDataVoid(),&n, w.GetData(), work.GetDataVoid(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());
    
#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

    Transpose(z);
  }
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>, Prop, RowHerm, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 2*n; Vector<complex<double> > work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zheev_(&job, &uplo, &n, A.GetDataVoid(), &n, w.GetData(),
	   work.GetDataVoid(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, RowHerm, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 2*n; Vector<complex<double> > work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    zheev_(&job, &uplo,&n, z.GetDataVoid(),&n, w.GetData(), work.GetDataVoid(),
	   &lwork, rwork.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(z);
  }
  
  
  /* ColHerm */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>, Prop, ColHerm, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('N'); int lwork = 2*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    cheev_(&job, &uplo, &n, A.GetDataVoid(), &n, w.GetData(),
	   work.GetDataVoid(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, ColHerm, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V'); int lwork = 2*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    cheev_(&job, &uplo, &n, z.GetDataVoid(), &n,
	   w.GetData(), work.GetDataVoid(),
	   &lwork, rwork.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>, Prop, ColHerm, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('N'); int lwork = 2*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zheev_(&job, &uplo, &n, A.GetDataVoid(), &n, w.GetData(),
	   work.GetDataVoid(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, ColHerm, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V'); int lwork = 2*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    zheev_(&job, &uplo, &n, z.GetDataVoid(), &n,
	   w.GetData(), work.GetDataVoid(),
	   &lwork, rwork.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  /* RowSymPacked */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<float, Prop, RowSymPacked, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N');
    Vector<float> work(3*n);
    w.Reallocate(n);
    sspev_(&job, &uplo, &n, A.GetData(), w.GetData(), A.GetData(), &n,
	   work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop,RowSymPacked, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<float, General, RowMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V');
    Vector<float> work(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    sspev_(&job, &uplo, &n, A.GetData(), w.GetData(), z.GetData(), &n,
	   work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>, Prop, RowSymPacked, Allocator1>& A,
		      Vector<complex<float>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    GetEigenvalues(B, w);
  }
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>, Prop,
				  RowSymPacked, Allocator1>& A,
				  Vector<complex<float>,VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, RowMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<double, Prop, RowSymPacked, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N'); Vector<double> work(3*n);
    w.Reallocate(n);
    dspev_(&job, &uplo, &n, A.GetData(), w.GetData(), A.GetData(), &n,
	   work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<double,Prop,RowSymPacked, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<double, General, RowMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V'); Vector<double> work(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    dspev_(&job, &uplo, &n, A.GetData(), w.GetData(), z.GetData(), &n,
	   work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>, Prop, RowSymPacked, Allocator1>& A,
		      Vector<complex<double>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    GetEigenvalues(B, w);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>, Prop,
				  RowSymPacked, Allocator1>& A,
				  Vector<complex<double>,VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, RowMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  /* ColSymPacked */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<float, Prop, ColSymPacked, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('N');
    Vector<float> work(3*n);
    w.Reallocate(n);
    sspev_(&job, &uplo, &n, A.GetData(), w.GetData(), A.GetData(),
	   &n, work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop,ColSymPacked, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<float, General, ColMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V');
    Vector<float> work(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    sspev_(&job, &uplo, &n, A.GetData(), w.GetData(), z.GetData(),
	   &n, work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>,
		      Prop, ColSymPacked, Allocator1>& A,
		      Vector<complex<float>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    GetEigenvalues(B, w);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, ColSymPacked, Allocator1>& A,
				  Vector<complex<float>, VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<double, Prop, ColSymPacked, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('N'); Vector<double> work(3*n);
    w.Reallocate(n);
    dspev_(&job, &uplo, &n, A.GetData(), w.GetData(), A.GetData(),
	   &n, work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<double,Prop,ColSymPacked, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<double, General, ColMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V');
    Vector<double> work(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    dspev_(&job, &uplo, &n, A.GetData(), w.GetData(), z.GetData(),
	   &n, work.GetData() , &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>,
		      Prop, ColSymPacked, Allocator1>& A,
		      Vector<complex<double>, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    GetEigenvalues(B, w);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, ColSymPacked, Allocator1>& A,
				  Vector<complex<double>, VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator3>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> B(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	B(i,j) = A(i,j);
    
    w.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(B, w, z);
  }
  
  
  /* RowHermPacked */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>,
		      Prop, RowHermPacked, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N');
    Vector<complex<float> > work(2*n);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    chpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), A.GetDataVoid(), &n,
	   work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, RowHermPacked, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V');
    Vector<complex<float> > work(2*n);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    chpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), z.GetDataVoid(),
	   &n, work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(z);
  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>,
		      Prop, RowHermPacked, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('N');
    Vector<complex<double> > work(2*n);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zhpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), A.GetDataVoid(), &n,
	   work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, RowHermPacked, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('L'); char job('V');
    Vector<complex<double> > work(2*n);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    zhpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), z.GetDataVoid(),
	   &n, work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(z);
  }
  
  
  /* ColHermPacked */
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<float>,
		      Prop, ColHermPacked, Allocator1>& A,
		      Vector<float, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U');
    char job('N');
    Vector<complex<float> > work(2*n);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    chpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), A.GetDataVoid(), &n,
	   work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop, ColHermPacked, Allocator1>& A,
				  Vector<float, VectFull, Allocator2>& w,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V');
    Vector<complex<float> > work(2*n);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    chpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), z.GetDataVoid(), &n,
	   work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2>
  void GetEigenvalues(Matrix<complex<double>,
		      Prop, ColHermPacked, Allocator1>& A,
		      Vector<double, VectFull, Allocator2>& w,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U');
    char job('N');
    Vector<complex<double> > work(2*n);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zhpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), A.GetDataVoid(), &n,
	   work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop, class Allocator1, class Allocator2, class Allocator3>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop, ColHermPacked, Allocator1>& A,
				  Vector<double, VectFull, Allocator2>& w,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator3>&z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char uplo('U'); char job('V');
    Vector<complex<double> > work(2*n);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    zhpev_(&job, &uplo, &n, A.GetDataVoid(), w.GetData(), z.GetDataVoid(), &n,
	   work.GetDataVoid(), rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  // STANDARD EIGENVALUE PROBLEM //
  /////////////////////////////////
  
  
  ////////////////////////////////////
  // GENERALIZED EIGENVALUE PROBLEM //
  
  
  /* RowSym */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<float, Prop1, RowSym, Allocator1>& A,
		      Matrix<float, Prop2, RowSym, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    ssygv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop1, RowSym, Allocator1>& A,
				  Matrix<float, Prop2, RowSym, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<float, General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    ssygv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<float>, Prop1, RowSym, Allocator1>& A,
		      Matrix<complex<float>, Prop2, RowSym, Allocator2>& B,
		      Vector<complex<float>, VectFull, Allocator4>& alpha,
		      Vector<complex<float>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N');
    int lwork = 2*n; Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, RowSym, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, RowSym, Allocator2>& B,
				  Vector<complex<float>,
				  VectFull, Allocator4>& alpha,
				  Vector<complex<float>,
				  VectFull, Allocator5>& beta,
				  Matrix<complex<float>,
				  Prop3, RowMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 2*n; Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n);
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    TransposeConj(V);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<double, Prop1, RowSym, Allocator1>& A,
		      Matrix<double, Prop2, RowSym, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    dsygv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop1, RowSym, Allocator1>& A,
				  Matrix<double, Prop2, RowSym, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<double, General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    dsygv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<double>, Prop1, RowSym, Allocator1>& A,
		      Matrix<complex<double>, Prop2, RowSym, Allocator2>& B,
		      Vector<complex<double>, VectFull, Allocator4>& alpha,
		      Vector<complex<double>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N');
    int lwork = 2*n; Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, RowSym, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, RowSym, Allocator2>& B,
				  Vector<complex<double>,
				  VectFull, Allocator4>& alpha,
				  Vector<complex<double>,
				  VectFull, Allocator5>& beta,
				  Matrix<complex<double>,
				  Prop3, RowMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 2*n; Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    TransposeConj(V);
  }
  
  
  /* ColSym */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<float, Prop1, ColSym, Allocator1>& A,
		      Matrix<float, Prop2, ColSym, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    ssygv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop1, ColSym, Allocator1>& A,
				  Matrix<float, Prop2, ColSym, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<float, General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    ssygv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<float>, Prop1, ColSym, Allocator1>& A,
		      Matrix<complex<float>, Prop2, ColSym, Allocator2>& B,
		      Vector<complex<float>, VectFull, Allocator4>& alpha,
		      Vector<complex<float>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N'); int lwork = 2*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Alloc1,
	   class Alloc2, class Alloc4, class Alloc5, class Alloc6>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, ColSym, Alloc1>& A,
				  Matrix<complex<float>,
				  Prop2, ColSym, Alloc2>& B,
				  Vector<complex<float>, VectFull, Alloc4>& alpha,
				  Vector<complex<float>, VectFull, Alloc5>& beta,
				  Matrix<complex<float>,
				  Prop3, ColMajor, Alloc6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('V');
    int lwork = 2*n; Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<double, Prop1, ColSym, Allocator1>& A,
		      Matrix<double, Prop2, ColSym, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    dsygv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop1, ColSym, Allocator1>& A,
				  Matrix<double, Prop2, ColSym, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<double, General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    dsygv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<double>, Prop1, ColSym, Allocator1>& A,
		      Matrix<complex<double>, Prop2, ColSym, Allocator2>& B,
		      Vector<complex<double>, VectFull, Allocator4>& alpha,
		      Vector<complex<double>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N'); int lwork = 2*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Alloc1,
	   class Alloc2, class Alloc4, class Alloc5, class Alloc6>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, ColSym, Alloc1>& A,
				  Matrix<complex<double>,
				  Prop2, ColSym, Alloc2>& B,
				  Vector<complex<double>, VectFull, Alloc4>& alpha,
				  Vector<complex<double>, VectFull, Alloc5>& beta,
				  Matrix<complex<double>,
				  Prop3, ColMajor, Alloc6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('V');
    int lwork = 2*n; Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  /* RowHerm */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<float>, Prop1, RowHerm, Allocator1>& A,
		      Matrix<complex<float>, Prop2, RowHerm, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 2*n; Vector<float> work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    chegv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, RowHerm, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, RowHerm, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 2*n; Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    Vector<float> rwork(3*n);
    chegv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<double>, Prop1, RowHerm, Allocator1>& A,
		      Matrix<complex<double>, Prop2, RowHerm, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 3*n; Vector<double> work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zhegv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, RowHerm, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, RowHerm, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    Vector<double> rwork(3*n);
    zhegv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  /* ColHerm */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<float>, Prop1, ColHerm, Allocator1>& A,
		      Matrix<complex<float>, Prop2, ColHerm, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 2*n; Vector<float> work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    chegv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, ColHerm, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, ColHerm, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    Vector<float> rwork(3*n);
    chegv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork,
	   rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<double>, Prop1, ColHerm, Allocator1>& A,
		      Matrix<complex<double>, Prop2, ColHerm, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 3*n; Vector<double> work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zhegv_(&itype, &job, &uplo, &n, A.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork, rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, ColHerm, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, ColHerm, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	z(i,j) = A(i,j);
    
    Vector<double> rwork(3*n);
    zhegv_(&itype, &job, &uplo, &n, z.GetData(), &n, B.GetData(), &n,
	   w.GetData(), work.GetData(), &lwork,
	   rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  /* RowSymPacked */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<float, Prop1, RowSymPacked, Allocator1>& A,
		      Matrix<float, Prop2, RowSymPacked, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    sspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   A.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<float,
				  Prop1, RowSymPacked, Allocator1>& A,
				  Matrix<float,
				  Prop2, RowSymPacked, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<float,
				  General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    sspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvalues(Matrix<complex<float>,
		      Prop1, RowSymPacked, Allocator1>& A,
		      Matrix<complex<float>,
		      Prop2, RowSymPacked, Allocator2>& B,
		      Vector<complex<float>, VectFull, Allocator3>& alpha,
		      Vector<complex<float>, VectFull, Allocator4>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, ColMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    GetEigenvalues(C, D, alpha, beta);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Allocator5>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, RowSymPacked, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, RowSymPacked, Allocator2>& B,
				  Vector<complex<float>,
				  VectFull, Allocator3>& alpha,
				  Vector<complex<float>,
				  VectFull, Allocator4>& beta,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator5>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, RowMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(C, D, alpha, beta, z);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<double, Prop1, RowSymPacked, Allocator1>& A,
		      Matrix<double, Prop2, RowSymPacked, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    dspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   A.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<double,
				  Prop1, RowSymPacked, Allocator1>& A,
				  Matrix<double,
				  Prop2, RowSymPacked, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<double,
				  General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    dspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvalues(Matrix<complex<double>,
		      Prop1, RowSymPacked, Allocator1>& A,
		      Matrix<complex<double>,
		      Prop2, RowSymPacked, Allocator2>& B,
		      Vector<complex<double>, VectFull, Allocator3>& alpha,
		      Vector<complex<double>, VectFull, Allocator4>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    GetEigenvalues(C, D, alpha, beta);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Allocator5>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, RowSymPacked, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, RowSymPacked, Allocator2>& B,
				  Vector<complex<double>,
				  VectFull, Allocator3>& alpha,
				  Vector<complex<double>,
				  VectFull, Allocator4>& beta,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator5>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, RowMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(C, D, alpha, beta, z);
  }
  
  
  /* ColSymPacked */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<float, Prop1, ColSymPacked, Allocator1>& A,
		      Matrix<float, Prop2, ColSymPacked, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 3*n; Vector<float> work(lwork);
    w.Reallocate(n);
    sspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   A.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<float,
				  Prop1, ColSymPacked, Allocator1>& A,
				  Matrix<float,
				  Prop2, ColSymPacked, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<float,
				  General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V'); int lwork = 3*n;
    Vector<float> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    sspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvalues(Matrix<complex<float>,
		      Prop1, ColSymPacked, Allocator1>& A,
		      Matrix<complex<float>,
		      Prop2, ColSymPacked, Allocator2>& B,
		      Vector<complex<float>, VectFull, Allocator3>& alpha,
		      Vector<complex<float>, VectFull, Allocator4>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, ColMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    GetEigenvalues(C, D, alpha, beta);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Allocator5>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, ColSymPacked, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, ColSymPacked, Allocator2>& B,
				  Vector<complex<float>,
				  VectFull, Allocator3>& alpha,
				  Vector<complex<float>,
				  VectFull, Allocator4>& beta,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator5>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<float>, General, ColMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(C, D, alpha, beta, z);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<double, Prop1, ColSymPacked, Allocator1>& A,
		      Matrix<double, Prop2, ColSymPacked, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 3*n; Vector<double> work(lwork);
    w.Reallocate(n);
    dspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   A.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<double,
				  Prop1, ColSymPacked, Allocator1>& A,
				  Matrix<double,
				  Prop2, ColSymPacked, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<double,
				  General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V'); int lwork = 3*n;
    Vector<double> work(lwork);
    w.Reallocate(n);
    z.Reallocate(n,n);
    dspgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvalues(Matrix<complex<double>,
		      Prop1, ColSymPacked, Allocator1>& A,
		      Matrix<complex<double>,
		      Prop2, ColSymPacked, Allocator2>& B,
		      Vector<complex<double>, VectFull, Allocator3>& alpha,
		      Vector<complex<double>, VectFull, Allocator4>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    GetEigenvalues(C, D, alpha, beta);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4, class Allocator5>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, ColSymPacked, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, ColSymPacked, Allocator2>& B,
				  Vector<complex<double>,
				  VectFull, Allocator3>& alpha,
				  Vector<complex<double>,
				  VectFull, Allocator4>& beta,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator5>& z,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> C(n,n), D(n,n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
	{
	  C(i,j) = A(i,j);
	  D(i,j) = B(i,j);
	}
    
    alpha.Reallocate(n);
    beta.Reallocate(n);
    z.Reallocate(n,n);
    GetEigenvaluesEigenvectors(C, D, alpha, beta, z);
  }
  
  
  /* RowHermPacked */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<float>,
		      Prop1, RowHermPacked, Allocator1>& A,
		      Matrix<complex<float>,
		      Prop2, RowHermPacked, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 2*n;
    Vector<float> work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    chpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   A.GetData(), &n, work.GetData(), rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, RowHermPacked, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, RowHermPacked, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<complex<float>,
				  General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V'); int lwork = 2*n;
    Vector<float> work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    chpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<double>,
		      Prop1, RowHermPacked, Allocator1>& A,
		      Matrix<complex<double>,
		      Prop2, RowHermPacked, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('N');
    int lwork = 2*n;
    Vector<double> work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zhpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   A.GetData(), &n, work.GetData(), rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, RowHermPacked, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, RowHermPacked, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<complex<double>,
				  General, RowMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('L'); char job('V'); int lwork = 2*n;
    Vector<double> work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    zhpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  /* ColHermPacked */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<float>,
		      Prop1, ColHermPacked, Allocator1>& A,
		      Matrix<complex<float>,
		      Prop2, ColHermPacked, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 2*n;
    Vector<float> work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    chpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(),
	   w.GetData(), A.GetData(), &n, work.GetData(),
	   rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, ColHermPacked, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, ColHermPacked, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& w,
				  Matrix<complex<float>,
				  General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n;
    Vector<float> work(lwork);
    Vector<float> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    chpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3>
  void GetEigenvalues(Matrix<complex<double>,
		      Prop1, ColHermPacked, Allocator1>& A,
		      Matrix<complex<double>,
		      Prop2, ColHermPacked, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& w,
		      LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('N');
    int lwork = 2*n;
    Vector<double> work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    zhpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(),
	   w.GetData(), A.GetData(), &n, work.GetData(),
	   rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, ColHermPacked, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, ColHermPacked, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& w,
				  Matrix<complex<double>,
				  General, ColMajor, Allocator4>& z,
				  LapackInfo& info = lapack_info)
  {
    int itype = 1;
    int n = A.GetM();
    char uplo('U'); char job('V');
    int lwork = 3*n;
    Vector<double> work(lwork);
    Vector<double> rwork(3*n);
    w.Reallocate(n);
    z.Reallocate(n,n);
    zhpgv_(&itype, &job, &uplo, &n, A.GetData(), B.GetData(), w.GetData(),
	   z.GetData(), &n, work.GetData(), rwork.GetData(),
	   &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  /* RowMajor */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3,
	   class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<float, Prop1, RowMajor, Allocator1>& A,
		      Matrix<float, Prop2, RowMajor, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& alpha_real,
		      Vector<float, VectFull, Allocator4>& alpha_imag,
		      Vector<float, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N');
    int lwork = 8*n+16; Vector<float> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);
    sggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha_real.GetData(), alpha_imag.GetData(), beta.GetData(),
	   A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop1, RowMajor, Allocator1>& A,
				  Matrix<float, Prop2, RowMajor, Allocator2>& B,
				  Vector<float, VectFull, Allocator3>& alpha_real,
				  Vector<float, VectFull, Allocator4>& alpha_imag,
				  Vector<float, VectFull, Allocator5>& beta,
				  Matrix<float, Prop3, RowMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 8*n+16; Vector<float> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    sggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha_real.GetData(), alpha_imag.GetData(), beta.GetData(),
	   V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());


#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(V);
    // conjugate if necessary
    int i = 0;
    while (i < n)
      {
	if (i < (n-1))
	  if (alpha_real(i) == alpha_real(i+1))
	    {
	      for (int j = 0; j < n; j++)
		V(j,i+1) = -V(j,i+1);
	      
	      i++;
	    }
	
	i++;
      }
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<float>, Prop1, RowMajor, Allocator1>& A,
		      Matrix<complex<float>, Prop2, RowMajor, Allocator2>& B,
		      Vector<complex<float>, VectFull, Allocator4>& alpha,
		      Vector<complex<float>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N'); int lwork = 2*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, RowMajor, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, RowMajor, Allocator2>& B,
				  Vector<complex<float>,
				  VectFull, Allocator4>& alpha,
				  Vector<complex<float>,
				  VectFull, Allocator5>& beta,
				  Matrix<complex<float>,
				  Prop3, RowMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 2*n; Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n, alpha.GetData(),
	   beta.GetData(), V.GetData(), &n, V.GetData(), &n, work.GetData(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());


#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    TransposeConj(V);
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3,
	   class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<double, Prop1, RowMajor, Allocator1>& A,
		      Matrix<double, Prop2, RowMajor, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& alpha_real,
		      Vector<double, VectFull, Allocator4>& alpha_imag,
		      Vector<double, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N');
    int lwork = 8*n+16; Vector<double> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);
    dggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha_real.GetData(), alpha_imag.GetData(), beta.GetData(),
	   A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop1, RowMajor, Allocator1>& A,
				  Matrix<double, Prop2, RowMajor, Allocator2>& B,
				  Vector<double, VectFull, Allocator3>& alpha_real,
				  Vector<double, VectFull, Allocator4>& alpha_imag,
				  Vector<double, VectFull, Allocator5>& beta,
				  Matrix<double, Prop3, RowMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 8*n+16; Vector<double> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n, n);
    dggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha_real.GetData(), alpha_imag.GetData(), beta.GetData(),
	   V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    Transpose(V);
    // conjugate if necessary
    int i = 0;
    while (i < n)
      {
	if (i < (n-1))
	  if (alpha_real(i) == alpha_real(i+1))
	    {
	      for (int j = 0; j < n; j++)
		V(j,i+1) = -V(j,i+1);
	      
	      i++;
	    }
	
	i++;
      }
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<double>, Prop1, RowMajor, Allocator1>& A,
		      Matrix<complex<double>, Prop2, RowMajor, Allocator2>& B,
		      Vector<complex<double>, VectFull, Allocator4>& alpha,
		      Vector<complex<double>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N'); int lwork = 2*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha.GetData(), beta.GetData(), A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, RowMajor, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, RowMajor, Allocator2>& B,
				  Vector<complex<double>,
				  VectFull, Allocator4>& alpha,
				  Vector<complex<double>,
				  VectFull, Allocator5>& beta,
				  Matrix<complex<double>,
				  Prop3, RowMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 2*n; Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n, alpha.GetData(),
	   beta.GetData(), V.GetData(), &n, V.GetData(), &n, work.GetData(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
    TransposeConj(V);
  }
  
  
  /* ColMajor */
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3,
	   class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<float, Prop1, ColMajor, Allocator1>& A,
		      Matrix<float, Prop2, ColMajor, Allocator2>& B,
		      Vector<float, VectFull, Allocator3>& alpha_real,
		      Vector<float, VectFull, Allocator4>& alpha_imag,
		      Vector<float, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N');
    int lwork = 8*n+16; Vector<float> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);

    sggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha_real.GetData(), alpha_imag.GetData(), beta.GetData(),
	   A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<float, Prop1, ColMajor, Allocator1>& A,
				  Matrix<float, Prop2, ColMajor, Allocator2>& B,
				  Vector<float,
				  VectFull, Allocator3>& alpha_real,
				  Vector<float,
				  VectFull, Allocator4>& alpha_imag,
				  Vector<float, VectFull, Allocator5>& beta,
				  Matrix<float, Prop3, ColMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 8*n+16; Vector<float> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    sggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(),
	   &n, alpha_real.GetData(), alpha_imag.GetData(),
	   beta.GetData(), V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<float>, Prop1, ColMajor, Allocator1>& A,
		      Matrix<complex<float>, Prop2, ColMajor, Allocator2>& B,
		      Vector<complex<float>, VectFull, Allocator4>& alpha,
		      Vector<complex<float>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N'); int lwork = 2*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(),
	   &n, alpha.GetData(),
	   beta.GetData(), A.GetData(), &n, A.GetData(), &n, work.GetData(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<complex<float>,
				  Prop1, ColMajor, Allocator1>& A,
				  Matrix<complex<float>,
				  Prop2, ColMajor, Allocator2>& B,
				  Vector<complex<float>,
				  VectFull, Allocator4>& alpha,
				  Vector<complex<float>,
				  VectFull, Allocator5>& beta,
				  Matrix<complex<float>,
				  Prop3, ColMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('V'); int lwork = 2*n;
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    cggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(),
	   &n, alpha.GetData(),
	   beta.GetData(), V.GetData(), &n, V.GetData(), &n, work.GetData(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator3,
	   class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<double, Prop1, ColMajor, Allocator1>& A,
		      Matrix<double, Prop2, ColMajor, Allocator2>& B,
		      Vector<double, VectFull, Allocator3>& alpha_real,
		      Vector<double, VectFull, Allocator4>& alpha_imag,
		      Vector<double, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N');
    int lwork = 8*n+16; Vector<double> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);
    
    dggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(), &n,
	   alpha_real.GetData(), alpha_imag.GetData(), beta.GetData(),
	   A.GetData(), &n, A.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator3, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<double, Prop1, ColMajor, Allocator1>& A,
				  Matrix<double, Prop2, ColMajor, Allocator2>& B,
				  Vector<double,
				  VectFull, Allocator3>& alpha_real,
				  Vector<double,
				  VectFull, Allocator4>& alpha_imag,
				  Vector<double, VectFull, Allocator5>& beta,
				  Matrix<double, Prop3, ColMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('V'), jobvr('N');
    int lwork = 8*n+16; Vector<double> work(lwork);
    alpha_real.Reallocate(n);
    alpha_imag.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    
    dggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(),
	   &n, alpha_real.GetData(), alpha_imag.GetData(),
	   beta.GetData(), V.GetData(), &n, V.GetData(), &n,
	   work.GetData(), &lwork, &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif
    
  }
  
  
  template<class Prop1, class Prop2, class Allocator1,
	   class Allocator2, class Allocator4, class Allocator5>
  void GetEigenvalues(Matrix<complex<double>, Prop1, ColMajor, Allocator1>& A,
		      Matrix<complex<double>, Prop2, ColMajor, Allocator2>& B,
		      Vector<complex<double>, VectFull, Allocator4>& alpha,
		      Vector<complex<double>, VectFull, Allocator5>& beta,
		      LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('N'); int lwork = 2*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);

    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(),
	   &n, alpha.GetData(),
	   beta.GetData(), A.GetData(), &n, A.GetData(), &n, work.GetData(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  template<class Prop1, class Prop2, class Prop3, class Allocator1,
	   class Allocator2, class Allocator4,
	   class Allocator5, class Allocator6>
  void GetEigenvaluesEigenvectors(Matrix<complex<double>,
				  Prop1, ColMajor, Allocator1>& A,
				  Matrix<complex<double>,
				  Prop2, ColMajor, Allocator2>& B,
				  Vector<complex<double>,
				  VectFull, Allocator4>& alpha,
				  Vector<complex<double>,
				  VectFull, Allocator5>& beta,
				  Matrix<complex<double>,
				  Prop3, ColMajor, Allocator6>& V,
				  LapackInfo& info = lapack_info)
  {
    int n = A.GetM();
    char jobvl('N'), jobvr('V'); int lwork = 2*n;
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(8*n);
    alpha.Reallocate(n);
    beta.Reallocate(n);
    V.Reallocate(n,n);
    zggev_(&jobvl, &jobvr, &n, A.GetData(), &n, B.GetData(),
	   &n, alpha.GetData(),
	   beta.GetData(), V.GetData(), &n, V.GetData(), &n, work.GetData(),
	   &lwork, rwork.GetData(), &info.GetInfoRef());

#ifdef SELDON_LAPACK_CHECK_INFO
    if (info.GetInfo() != 0)
      throw LapackError(info.GetInfo(), "GetEigenvalues",
			"Failed to find eigenvalues ");
#endif

  }
  
  
  // GENERALIZED EIGENVALUE PROBLEM //
  ////////////////////////////////////
  
  
  //////////////////////////////////
  // SINGULAR VALUE DECOMPOSITION //
  
  
  /* RowMajor */
  
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<float, Prop1, RowMajor, Allocator1>& A,
	      Vector<float, VectFull, Allocator4>& lambda,
	      Matrix<float, General, RowMajor, Allocator2>& u,
	      Matrix<float, General, RowMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetM();
    char jobl('A'), jobr('A');
    int lwork = max(3*min(m,n)+max(m,n), 5*min(m,n));
    Vector<float> work(lwork);
    sgesvd_(&jobl, &jobr, &m, &n, A.GetData(), &n, lambda.GetData(),
	    u.GetData(), &n, v.GetData(), &n, work.GetData(),
	    &lwork, &info.GetInfoRef());
  }
  
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<complex<float>, Prop1, RowMajor, Allocator1>& A,
	      Vector<complex<float>, VectFull, Allocator4>& lambda,
	      Matrix<complex<float>, General, RowMajor, Allocator2>& u,
	      Matrix<complex<float>, General, RowMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetM();
    char jobl('A'), jobr('A');
    int lwork = 2*min(m,n)+max(m,n);
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(5*min(m,n));
    cgesvd_(&jobl, &jobr, &m, &n, A.GetDataVoid(), &n, lambda.GetDataVoid(),
	    u.GetDataVoid(), &n, v.GetDataVoid(), &n, work.GetDataVoid(),
	    &lwork, rwork.GetData(), &info.GetInfoRef());
  }
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<double, Prop1, RowMajor, Allocator1>& A,
	      Vector<double, VectFull, Allocator4>& lambda,
	      Matrix<double, General, RowMajor, Allocator2>& u,
	      Matrix<double, General, RowMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetM();
    char jobl('A'), jobr('A');
    int lwork = max(3*min(m,n)+max(m,n), 5*min(m,n));
    Vector<double> work(lwork);
    dgesvd_(&jobl, &jobr, &m, &n, A.GetData(), &n, lambda.GetData(),
	    u.GetData(), &n, v.GetData(), &n, work.GetData(),
	    &lwork, &info.GetInfoRef());
  }
  
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<complex<double>, Prop1, RowMajor, Allocator1>& A,
	      Vector<complex<double>, VectFull, Allocator4>& lambda,
	      Matrix<complex<double>, General, RowMajor, Allocator2>& u,
	      Matrix<complex<double>, General, RowMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetM();
    char jobl('A'), jobr('A');
    int lwork = 2*min(m,n)+max(m,n);
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(5*min(m,n));
    zgesvd_(&jobl, &jobr, &m, &n, A.GetDataVoid(), &n, lambda.GetDataVoid(),
	    u.GetDataVoid(), &n, v.GetDataVoid(), &n, work.GetDataVoid(),
	    &lwork, rwork.GetData(), &info.GetInfoRef());
  }
  
  
  /* ColMajor */
  
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<float, Prop1, ColMajor, Allocator1>& A,
	      Vector<float, VectFull, Allocator4>& lambda,
	      Matrix<float, General, ColMajor, Allocator2>& u,
	      Matrix<float, General, ColMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetN();
    char jobl('A'), jobr('A');
    int lwork = max(3*min(m,n)+max(m,n), 5*min(m,n));
    Vector<float> work(lwork);
    lambda.Reallocate(min(m, n)); lambda.Zero();
    u.Reallocate(m, m); u.Zero();
    v.Reallocate(n, n); v.Zero();
    sgesvd_(&jobl, &jobr, &m, &n, A.GetData(), &m, lambda.GetData(),
	    u.GetData(), &m, v.GetData(), &n, work.GetData(),
	    &lwork, &info.GetInfoRef());
  }
  
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<complex<float>, Prop1, ColMajor, Allocator1>& A,
	      Vector<complex<float>, VectFull, Allocator4>& lambda,
	      Matrix<complex<float>, General, ColMajor, Allocator2>& u,
	      Matrix<complex<float>, General, ColMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetN();
    char jobl('A'), jobr('A');
    int lwork = 2*min(m,n)+max(m,n);
    Vector<complex<float> > work(lwork);
    Vector<float> rwork(5*min(m,n));
    lambda.Reallocate(min(m, n)); lambda.Zero();
    u.Reallocate(m, m); u.Zero();
    v.Reallocate(n, n); v.Zero();
    cgesvd_(&jobl, &jobr, &m, &n, A.GetDataVoid(), &m, lambda.GetDataVoid(),
	    u.GetDataVoid(), &m, v.GetDataVoid(), &n, work.GetDataVoid(),
	    &lwork, rwork.GetData(), &info.GetInfoRef());
  }
  
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<double, Prop1, ColMajor, Allocator1>& A,
	      Vector<double, VectFull, Allocator4>& lambda,
	      Matrix<double, General, ColMajor, Allocator2>& u,
	      Matrix<double, General, ColMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetN();
    char jobl('A'), jobr('A');
    int lwork =10*max(m,n);
    Vector<double> work(lwork);
    lambda.Reallocate(min(m, n)); lambda.Zero();
    u.Reallocate(m, m); u.Zero();
    v.Reallocate(n, n); v.Zero();
    dgesvd_(&jobl, &jobr, &m, &n, A.GetData(), &m, lambda.GetData(),
	    u.GetData(), &m, v.GetData(), &n, work.GetData(),
	    &lwork, &info.GetInfoRef());
    
  }
  
  
  template<class Prop1, class Allocator1, class Allocator2,
	   class Allocator3, class Allocator4>
  void GetSVD(Matrix<complex<double>, Prop1, ColMajor, Allocator1>& A,
	      Vector<complex<double>, VectFull, Allocator4>& lambda,
	      Matrix<complex<double>, General, ColMajor, Allocator2>& u,
	      Matrix<complex<double>, General, ColMajor, Allocator3>& v,
	      LapackInfo& info = lapack_info)
  {
    int m = A.GetM(); int n = A.GetN();
    char jobl('A'), jobr('A');
    int lwork = 2*min(m,n)+max(m,n);
    Vector<complex<double> > work(lwork);
    Vector<double> rwork(5*min(m,n));
    lambda.Reallocate(min(m, n)); lambda.Zero();
    u.Reallocate(m, m); u.Zero();
    v.Reallocate(n, n); v.Zero();
    zgesvd_(&jobl, &jobr, &m, &n, A.GetDataVoid(), &m, lambda.GetDataVoid(),
	    u.GetDataVoid(), &m, v.GetDataVoid(), &n, work.GetDataVoid(),
	    &lwork, rwork.GetData(), &info.GetInfoRef());
  }
  
  
  // pseudo inverse
  template<class Prop1, class Allocator1>
  void GetPseudoInverse(Matrix<double, Prop1, ColMajor, Allocator1>& A,
			double epsilon, LapackInfo& info = lapack_info)
  {
    int m = A.GetM(), n = A.GetN();
    Vector<double, VectFull, Allocator1> lambda;
    Matrix<double, General, ColMajor, Allocator1> U;
    Matrix<double, General, ColMajor, Allocator1> V;
    
    GetSVD(A, lambda, U, V);
    
    A.Reallocate(n, m); A.Fill(0);
    // computation of A = V Sigma U^*
    for (int k = 0; k < min(m, n); k++)
      if (abs(lambda(k)) > epsilon)
	{
	  lambda(k) = 1.0/lambda(k);
	  for (int i = 0; i < n; i++)
	    for (int j = 0; j < m; j++)
	      A(i, j) += V(k, i)*lambda(k)*U(j, k);
	}
  }

  
  // SINGULAR VALUE DECOMPOSITION //
  //////////////////////////////////
  
  

  void GetHessian(Matrix<complex<double>, General, ColMajor>& A,
		  Matrix<complex<double>, General, ColMajor>& B,
		  Matrix<complex<double>, General, ColMajor>& Q,
		  Matrix<complex<double>, General, ColMajor>& Z)
  {
    char compq('V'), compz('I');
    int n = A.GetM(), ilo = 1, ihi = n, info, lwork = 4 * n;
    Vector<complex<double> > tau(n);
    Vector<complex<double> > work(lwork);
    zgeqrf_(&n, &n, B.GetDataVoid(), &n, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info);
    
    Q = B;
    zungqr_(&n, &n, &n, Q.GetDataVoid(), &n, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info);
    
    char side('L'), trans('C');
    zunmqr_(&side, &trans, &n, &n, &n, B.GetDataVoid(), &n, tau.GetDataVoid(),
	    A.GetDataVoid(), &n, work.GetData(), &lwork, &info);
    
    for (int i = 0; i < n; i++)
      for (int j = 0; j < i; j++)
	B(i, j) = 0;
    
    zgghrd_(&compq, &compz, &n, &ilo, &ihi, A.GetDataVoid(), &n,
	    B.GetDataVoid(), &n, Q.GetDataVoid(), &n, Z.GetDataVoid(),
	    &n, &info);
    
  }
  

  void GetQZ(Matrix<complex<double>, General, ColMajor>& A,
	     Matrix<complex<double>, General, ColMajor>& B,
	     Matrix<complex<double>, General, ColMajor>& Q,
	     Matrix<complex<double>, General, ColMajor>& Z)
  {
    char compq('V'), compz('I');
    int n = A.GetM(), ilo = 1, ihi = n, info, lwork = 4*n;
    Vector<complex<double> > tau(n);
    Vector<complex<double> > work(lwork);
    zgeqrf_(&n, &n, B.GetDataVoid(), &n, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info);
    
    Q = B;
    zungqr_(&n, &n, &n, Q.GetDataVoid(), &n, tau.GetDataVoid(),
	    work.GetDataVoid(), &lwork, &info);
    
    char side('L'), trans('C');
    zunmqr_(&side, &trans, &n, &n, &n, B.GetDataVoid(), &n, tau.GetDataVoid(),
	    A.GetDataVoid(), &n, work.GetData(), &lwork, &info);
    
    for (int i = 0; i < n; i++)
      for (int j = 0; j < i; j++)
	B(i,j) = 0;
    
    zgghrd_(&compq, &compz, &n, &ilo, &ihi, A.GetDataVoid(), &n,
	    B.GetDataVoid(), &n, Q.GetDataVoid(), &n, Z.GetDataVoid(),
	    &n, &info);
    
    char job('S');
    compq = 'V';
    compz = 'V';
    Vector<complex<double> > alpha(n), beta(n);
    Vector<double> rwork(lwork);
    zhgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, A.GetDataVoid(), &n,
	    B.GetDataVoid(), &n, alpha.GetDataVoid(), beta.GetDataVoid(),
	    Q.GetDataVoid(), &n, Z.GetDataVoid(), &n, work.GetDataVoid(),
	    &lwork, rwork.GetData(), &info);
  }
  
  
  void SolveSylvester(Matrix<complex<double>, General, ColMajor>& A,
		      Matrix<complex<double>, General, ColMajor>& B,
		      Matrix<complex<double>, General, ColMajor>& C,
		      Matrix<complex<double>, General, ColMajor>& D,
		      Matrix<complex<double>, General, ColMajor>& E)
  {
    complex<double> one(1), zero(0);
    int n = A.GetM();
    Matrix<complex<double>, General, ColMajor> Q1(n, n), Q2(n, n),
      Z1(n, n), Z2(n, n);
    Matrix<complex<double>, General, ColMajor> Y(n, n), F(n, n);
    
    GetHessian(A, C, Q1, Z1);
    GetQZ(D, B, Q2, Z2);
    
    Y.Zero();
    MltAdd(one, SeldonConjTrans, Q1, SeldonNoTrans, E, zero, Y);
    MltAdd(one, SeldonNoTrans, Y, SeldonNoTrans, Q2, zero, F);
    Y.Zero();
    
    Vector<complex<double> > ftemp(n), Yvec(n), ytmp(n), xtmp(n);
    Vector<int> pivot(n);
    for (int k = n-1; k >= 0; k--)
      {
	for (int j = 0; j < n; j++)
	  ftemp(j) = F(j,k);
	
	for (int j = k+1; j < n; j++)
	  {
	    for (int i = 0; i < n; i++)
	      Yvec(i) = Y(i, j);
	    
	    Mlt(A, Yvec, xtmp); Mlt(C, Yvec, ytmp);
	    for (int i = 0; i < n; i++)
	      ftemp(i) -= (conj(B(k, j)) * xtmp(i) + conj(D(k, j)) * ytmp(i));
	    
	  }
	
	for (int i = 0; i < n; i++)
	  for (int j = 0; j < n; j++)
	    E(i,j) = conj(B(k, k)) * A(i, j) + conj(D(k, k)) * C(i, j);
	
	GetLU(E, pivot);
	Yvec = ftemp;
	SolveLU(E, pivot, Yvec);
	
	for (int i = 0; i < n; i++)
	  Y(i, k) = Yvec(i);
      }
    
    MltAdd(one, SeldonNoTrans, Y, SeldonConjTrans, Z2, zero, F);
    MltAdd(one, Z1, F, zero, E);
  }

  
} // end namespace

#define SELDON_FILE_LAPACK_EIGENVALUES_CXX
#endif
