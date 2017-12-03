/**
 * This file includes linear algebra utilities used in solving SVD and
 * other problems.
 *
 * @author Ara Hayrapetyan
 */

#include "StdAfx.h"

/**
 * Computes al eigenvalues and eigenvectors of a real symmetric matrix
 * a[1..n][1.n]. On output, elements above the diagonal are destroyed.
 * Algorithm taken from Numerical Recipies pg#467.
 * The output of the algorithm is unsorted.
 *
 * @param a the input matrix.
 * @param n the dimension of the matrix.
 * @param d the eigenvalues are returned here.
 * @param v the eigenvectors are returned here.
 * @param nrot the number of Jacobi rotations required.
 */
void jacobi (double ** a, int n, double * d, double ** v, int * nrot)
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

  // Allocate.
  b = new double[n];
  z = new double[n];

  // Initialize to the identity matrix.
  for (ip = 0; ip < n; ip ++)
    {
      for (iq = 0; iq < n; iq ++)
	v[ip][iq] = 0.0;
      v[ip][ip] = 1.0;
    }

  // Initialize b & d to the diagonal of a.
  for (ip = 0; ip < n; ip ++)
    {
      b[ip] = d[ip] = a[ip][ip];
      z[ip] = 0.0;
    }

  // Initialize the number of rotations.
  *nrot = 0;

  // Main loop.
  for (i = 0; i < 50; i ++)
    {
      // Sum off diagonal elements.
      sm = 0;
      for (ip = 0; ip < (n - 1); ip ++)
	for (iq = (ip + 1); iq < n; iq ++)
	  sm += (double) fabs (a[ip][iq]);

      // The normal return, which relies on the quadratic convergence
      // to machine overflow.
      if (sm == 0.0)
	{
	  delete[] b;
	  delete[] z;
	  return;
	}

      // On the first 3 sweeps.
      if (i < 3)
	tresh = (double) 0.2*sm/(n*n);
      else
	tresh = 0.0;

      for (ip = 0; ip < (n - 1); ip ++)
	{
	  for (iq = (ip + 1); iq < n; iq ++)
	    {
	      g = (double)( 100 * fabs (a[ip][iq]) );

	      // After four sweeps, skip the rotation if the
	      // off-diagonal element is small.
	      if ((i > 3) && ((double)(fabs (d[ip] + g)) == (double)fabs (d[ip]))
		  && (double)(fabs (d[iq] + g) == (double)fabs(d[iq])))
		a[ip][iq] = 0;
	      else if (fabs(a[ip][iq]) > tresh)
		{
		  h = d[iq] - d[ip];
		  if ((double)(fabs(h) + g) == (double)fabs (h))
		    t = (a[ip][iq]) / h;
		  else
		    {
		      theta = (double)( 0.5 * h / (a[ip][iq]) );
		      t = (double)( 1.0 / (fabs (theta) + sqrt (1.0 + theta * theta)) );
		      if (theta < 0.0) t = -t;
		    }
		  c = (double) ( 1.0 / sqrt (1.0 + t * t) );
		  s = t * c;
		  tau = (double) (s / (1.0 + c) );
		  h= t * a[ip][iq];
		  z[ip] -= h;
		  z[iq] += h;
		  d[ip] -= h;
		  d[iq] += h;
		  a[ip][iq] = 0;

		  // Case of rotations 0<=j<p.
		  for (j = 0; j <= (ip - 1); j ++)
		    {
		      ROTATE(a,j,ip,j,iq);
		    }
		  // Case of rotations p<j<q.
		  for (j = (ip + 1); j <= (iq - 1); j ++)
		    {
		      ROTATE(a,ip,j,j,iq);
		    }
		  // Case of rotations q<j<n.
		  for (j = (iq + 1); j < n; j ++)
		    {
		      ROTATE(a,ip,j,iq,j);
		    }
		  
		  for (j = 0; j < n; j ++)
		    {
		      ROTATE(v,j,ip,j,iq);
		    }

		  ++ (*nrot);
		}
	    }
	}

      // Update d and reinitialize z.
      for (ip = 0; ip < n; ip ++)
	{
	  b[ip] += z[ip];
	  d[ip] = b[ip];
	  z[ip] = 0.0;
	}
    }
  
  // Deallocate.
  delete[] b;
  delete[] z;
  
  printf("Too many iterations in the jacobi routine.\n");
  exit (1);
}

/**
 * Given the eigenvalues d and eigenvectors v as output from Jacobi,
 * this routine sorts the eigenvalues into descending order, and
 * rearranges the columns of V correspondingly. The method is straight
 * insertion. Notice that O(N^2) is affordable since Jacobi is cubic.
 * Algorithm taken from Numerical Recipies pg#468.
 *
 * @param d the eigenvalues.
 * @param v the eigenvalues.
 * @param n the dimension of the matrix.
 */
void eigsrt (double * d, double ** v, int n)
{
  int i, j, k;
  double p;

  for (i = 0; i < n; i ++)
    {
      p = d[k = i];
      for (j = (i + 1); j < n; j ++)
	if (d[j] >= p)
	  p = d[k = j];

      if (i != k)
	// The current is not the maximum.
	{
	  d[k] = d[i];
	  d[i] = p;
	  for (j = 0; j < n; j ++)
	    {
	      p = v[j][i];
	      v[j][i] = v[j][k];
	      v[j][k] = p;
	    }
	}
    }
}

/**
 * Calculates the eigenvalues and eigenvectors of a matrix.
 * Assumes that the matrix is symmetric.
 *
 * @param A the matrix.
 * @param n the dimension of the matrix.
 * @param u the eigenvectors.
 * @param w the eigenvalues.
 */
void eigen (double ** A, int n, double ** u, double * w)
{
  int i, j;
  double ** a;
  int nrot;
  
  // Copy a, since it is to be damaged by the jacobi calculation
  // funciton.
  a = new double *[n];
  for (i = 0; i < n; i ++)
    a[i] = new double[n];

  for (i = 0; i < n; i ++)
    for (j = 0; j < n; j ++)
      a[i][j] = A[i][j];

  // Call the jacobi routine.
  jacobi (a, n, w, u, &nrot);
  
  // Sort the eigenvectors & values to descending order.
  eigsrt (w, u, n);

  // Deallocate.
  for (i = 0; i < n; i ++)
    delete[] a[i];
  delete[] a;
  
  // Task complete, return.
  return;
}

/**
 * Computes the Euclidean distance between to vectors.
 *
 * @param v1 the first vector.
 * @param v2 the second vector.
 */
double distance (double * v1, double * v2)
{
  double * v = new double[3];
  v[0] = v1[0] - v2[0];
  v[1] = v1[1] - v2[1];
  v[2] = v1[2] - v2[2];

  double result = sqrt (  (v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2])  );
  delete[](v);
  return result;
}

/**
 * Normalizes the vector.
 *
 * @param vector the vector.
 */
void normalize (double * v)
{
  double length = (double) (sqrt (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) );

  if(length < .000001){
	  v[0] = 0.0;
	  v[1] = 0.0;
	  v[2] = 0.0;
	  return;
  }

  v[0] = v[0] / length;
  v[1] = v[1] / length;
  v[2] = v[2] / length;
}

/*
 * Performs matrix-vector multiplication.
 * Specialized for the 3D case.
 * This function is a helper for the RMSD calculation.
 *
 * @param a the matrix.
 * @param v the vector.
 * @param vr the resulting vector.
 */
void matrix_vector_mult3D(double ** a, double * v, double * vr)
{
  int i;
  for (i = 0; i < 3; i++)
    vr[i] = v[0] * a[0][i] + v[1] * a[1][i] + v[2] * a[2][i];
}

/*
 * Performs matrix-matrix multiplications.
 * Specialized for the 3D case.
 * This function is a helper for the RMSD calculation.
 *
 * @param a the first matrix.
 * @param b the second matrix.
 * @param c the result matrix.
 */
void matrix_matrix_mult3D(double ** a, double ** b, double ** c)
{
  int i, j, k;
  
  // Clear c.
  for (i = 0; i < 4; i ++)
    for (j = 0; j < 4; j ++)
      c[i][j] = 0;

  // Complete the multiplication.
  for (i = 0; i < 4; i ++)
    for (j = 0; j < 4; j ++)
      for (k = 0; k < 4; k ++)
	c[i][j] += a[i][k] * b[k][j];
}

/**
 * Calculates the minimized RMSD between two molecular conformations.
 *
 * @param ref_orig the reference conformation.
 * @param fit_orig the fit conformation.
 * @param atoms the number of atoms.
 * @return the fit rmsd.
 *
 *###I HAVE MODIFIED THIS CODE TO RETURN AN ARRAY OF FLOATS###
 *###the first float is the RMSD, the rest are the atomic distances at optimal alignment
 */
double * min_rmsd (double ** ref_orig, double ** fit_orig, int atoms, double * transform)
{
//	bool testing = false;
//	if(transform[0] == 1){
//		testing = true;
//		printf("TESTING\n");
//	}

	int i, j, k;
	double * cmass_ref = new double[3];
	double * cmass_fit = new double[3];
	double ** ref;
	double ** fit;

	double * results = new double[atoms + 1];			///This is the result I have added

	// Copy the reference and fit to preserve the originals.
	ref = new double *[atoms];
	fit = new double *[atoms];
	for (i = 0; i < atoms; i ++){
		ref[i] = new double[3];
		fit[i] = new double[3];
	}
	for (i = 0; i < atoms; i ++){
		for (j = 0; j < 3; j ++){
			ref[i][j] = ref_orig[i][j];
			fit[i][j] = fit_orig[i][j];
		}
	}

	// Compute the center of mass of the reference and fit.
	cmass_ref[0] = cmass_ref[1] = cmass_ref[2] = 0;
	cmass_fit[0] = cmass_fit[1] = cmass_fit[2] = 0;

	for (i = 0; i < atoms; i ++){
		cmass_ref[0] += ref[i][0];
		cmass_ref[1] += ref[i][1];
		cmass_ref[2] += ref[i][2];
		cmass_fit[0] += fit[i][0];
		cmass_fit[1] += fit[i][1];
		cmass_fit[2] += fit[i][2];
    }

	cmass_ref[0] = cmass_ref[0] / ((double)atoms);
	cmass_ref[1] = cmass_ref[1] / ((double)atoms);
	cmass_ref[2] = cmass_ref[2] / ((double)atoms);
	cmass_fit[0] = cmass_fit[0] / ((double)atoms);
	cmass_fit[1] = cmass_fit[1] / ((double)atoms);
	cmass_fit[2] = cmass_fit[2] / ((double)atoms);


	// Compute R. R(i,j) is the dot product of i'th column of the
	// normalized reference and j'th column of the normalized fit.
	// But normalized we mean subtracting the appropriate Centers of
	// Masses. 
	double ** R = new double *[4];
	for (i = 0; i < 4; i ++){
		R[i] = new double[4];
	}
  
	for (i = 0; i < 4; i ++){
		for (j = 0; j < 4; j ++){
			R[i][j] = 0;
		}
	}
	R[3][3] = 1;

	for (i = 0; i < 3; i ++){
		for (j = 0; j < 3; j ++){
			for (k = 0; k < atoms; k ++){
				R[i][j] += (ref[k][i] - cmass_ref[i]) * (fit[k][j] - cmass_fit[j]);
			}
		}
	}

  // Transpose R.
  double ** Rt = new double *[4];
  for (i = 0; i < 4; i ++){
    Rt[i] = new double[4];
  }

  for (i = 0; i < 4; i ++){
	  for (j = 0; j < 4; j ++){
		   Rt[i][j] = R[j][i];
	  }
  }

  // Compute Rt * R.
  double ** RtR = new double *[4];
  for (i = 0; i < 4; i ++){
    RtR[i] = new double[4];
  }
  matrix_matrix_mult3D(Rt, R, RtR);

  // Compute the eigenvalues/eigenvectors of RtR.
  // The eigenvectors will come sorted.
  double * evalue = new double[3];
  double ** evector = new double*[3];
  for (i = 0; i < 3; i ++){
    evector[i] = new double[3];
  }
  eigen (RtR, 3, evector, evalue);

  // Transpose so that eigenvectors become rows.
  double temp_float;

  temp_float = evector[0][1];
  evector[0][1] = evector[1][0];
  evector[1][0] = temp_float;  
  temp_float = evector[0][2];
  evector[0][2] = evector[2][0];
  evector[2][0] = temp_float;
  temp_float = evector[1][2];
  evector[1][2] = evector[2][1];
  evector[2][1] = temp_float;

  // Determine b(i) = R*a(i).
  double ** b = new double *[3];
  for (i = 0; i < 3; i ++)
    b[i] = new double[3];

  matrix_vector_mult3D (Rt, evector[0], b[0]);
  matrix_vector_mult3D (Rt, evector[1], b[1]);
  matrix_vector_mult3D (Rt, evector[2], b[2]);

  normalize (b[0]);
  normalize (b[1]);
  normalize (b[2]);

  // Compute U = u(i,j) = sum(b(k,i) * a(k,j)) */
  double * *U = new double*[4];
  for (i = 0; i < 4; i ++){
	  U[i] = new double[4];
	  for (j = 0; j < 4; j ++){
		  U[i][j] = 0;
	  }
  }
  U[3][3] = 1;
  
	for (i = 0; i < 3; i ++){
		for (j = 0; j < 3; j ++) {
			for (k = 0; k < 3; k ++) {
				U[j][i] += (b[k][i] * evector[k][j]);
			}
		}
	}

/////inserted - cull numerical underflow
  for(i = 0; i<4; i++){
	  for(j = 0; j<4; j++){
		  if( fabs(U[i][j]) < .000001){
			U[i][j] = 0.0;
		  }
	  }
  }
/////inserted

/////inserted - catch degenerate rotations
  if( (U[0][0] == 0) && (U[0][1] == 0) && (U[0][2] == 0) && (U[1][0] == 0) && (U[2][0] == 0)){
	U[0][0] = 1.0;
  }
  if( (U[1][1] == 0) && (U[0][1] == 0) && (U[2][1] == 0) && (U[1][0] == 0) && (U[1][2] == 0)){
	U[1][1] = 1.0;
  }
  if( (U[2][2] == 0) && (U[2][1] == 0) && (U[2][0] == 0) && (U[1][2] == 0) && (U[0][2] == 0)){
	U[2][2] = 1.0;
  }
/////inserted


  // Update the rotational part of transform
  if(transform != NULL){
	  transform[0] = U[0][0];	transform[4] = U[1][0];		transform[8] = U[2][0];		transform[12] = U[3][0];
	  transform[1] = U[0][1];	transform[5] = U[1][1];		transform[9] = U[2][1];		transform[13] = U[3][1];
	  transform[2] = U[0][2];	transform[6] = U[1][2];		transform[10] = U[2][2];	transform[14] = U[3][2];
	  transform[3] = 0;			transform[7] = 0;			transform[11] = 0;			transform[15] = 1.0;
  }

  // Update coordinates.
  double x, y, z;
  for (i = 0; i < atoms; i ++)
    {
      x = fit[i][0] * U[0][0] + fit[i][1] * U[1][0]	+ fit[i][2] * U[2][0] + U[3][0];
      
      y = fit[i][0] * U[0][1] + fit[i][1] * U[1][1]	+ fit[i][2] * U[2][1] + U[3][1];
      
      z = fit[i][0] * U[0][2] + fit[i][1] * U[1][2]	+ fit[i][2] * U[2][2] + U[3][2];
      
      fit[i][0] = x;
      fit[i][1] = y;
      fit[i][2] = z;
    }

  // Compute the new center of mass.
  cmass_fit[0] = cmass_fit[1] = cmass_fit[2] = 0;
  for (i = 0; i < atoms; i ++)
    {
      cmass_fit[0] += fit[i][0];
      cmass_fit[1] += fit[i][1];
      cmass_fit[2] += fit[i][2];
    }
  cmass_fit[0] = cmass_fit[0] / ((double)atoms);
  cmass_fit[1] = cmass_fit[1] / ((double)atoms);
  cmass_fit[2] = cmass_fit[2] / ((double)atoms);

  double tx = cmass_ref[0] - cmass_fit[0];
  double ty = cmass_ref[1] - cmass_fit[1];
  double tz = cmass_ref[2] - cmass_fit[2];

/////Inserted
  if( fabs(tx) < 0.000001 ){ tx = 0.0;}
  if( fabs(ty) < 0.000001 ){ ty = 0.0;}
  if( fabs(tz) < 0.000001 ){ tz = 0.0;}
/////Inserted

  // Translate the fit.
  for (i = 0; i < atoms; i ++)
    {
      fit[i][0] += tx;
      fit[i][1] += ty;
      fit[i][2] += tz;
    }

  // Compound this onto the translation matrix, and set the scale to 1.
	if(transform != NULL){
		transform[12] += tx;
		transform[13] += ty;
		transform[14] += tz;
		transform[15] = 1.0;
	}

  // Return the RMSD of the fit and the reference.
  // Fill the Results Array with distance data
  double result = 0.0;
  for (i = 0; i < atoms; i++)
  {
	results[i+1] = distance (ref[i], fit[i]);
    result += (distance (ref[i], fit[i])) * (distance (ref[i], fit[i]));
  }
  result /= atoms;
  result = (double) (sqrt (result) );

  ///fill first results array position with RMSD
  results[0] = result;


#ifdef RMSD_CLAMP
  /////////////////////////////RMSD CLAMP/////////////////////////////////////////////////
  ////This hack fixes retarded RMSD issues assuming that geometry that can
  ////be very very closely aligned is actually identical.
  ////This is to catch the situation with RMSD = 1.234e-15 gives bad
  ////Transformations far from identity.
  /////////////////////////////RMSD CLAMP/////////////////////////////////////////////////
  if(results[0] < RMSD_CLAMP){
	  transform[0] = 1;	transform[4] = 0;		transform[8] = 0;		transform[12] = 0;
	  transform[1] = 0;	transform[5] = 1;		transform[9] = 0;		transform[13] = 0;
	  transform[2] = 0;	transform[6] = 0;		transform[10] = 1;		transform[14] = 0;
	  transform[3] = 0;	transform[7] = 0;		transform[11] = 0;		transform[15] = 1;
  }
  /////////////////////////////RMSD CLAMP/////////////////////////////////////////////////
  ////The problem is that while an excellent alignment can be calculated
  ////(i.e. with RMSD = epsilon) the transformation is not accurately 
  ////calculated and thus is not well defined numerically.  This hopefully
  ////resolves this issue.
  /////////////////////////////RMSD CLAMP/////////////////////////////////////////////////
#endif RMSD_CLAMP
  

  // Deallocate the allocated memory.
  for (i = 0; i < atoms; i ++)
    {
      delete[](ref[i]);
      delete[](fit[i]);
    }
  delete[](ref);
  delete[](fit);
  for (i = 0; i <  4; i ++)
    {
      delete[](R[i]);
      delete[](Rt[i]);
      delete[](RtR[i]);
    }
  delete[](R);
  delete[](Rt);
  delete[](RtR);
  for (i = 0; i < 3; i ++)
    {
      delete[](evector[i]);
      delete[](b[i]);
    }
  for (i = 0; i < 4; i ++)
    {
      delete[](U[i]);
    }
  delete[](U);
  delete[](evector);
  delete[](evalue);
  delete[](b);
  delete[](cmass_ref);
  delete[](cmass_fit);
  
  // Return.
  return results;
}


