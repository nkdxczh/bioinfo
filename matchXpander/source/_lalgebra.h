/**
 * This file includes linear algebra utilities used in solving SVD and
 * other problems.
 *
 * @author Ara Hayrapetyan
 */

#ifndef _LALGEBRA_H
#define _LALGEBRA_H

/** Rotation function - taken from Numerical Recipies. */
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);

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
void jacobi (double ** a, int n, double * d, double ** v, int * nrot);

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
void eigsrt (double * d, double ** v, int n);

/**
 * Calculates the eigenvalues and eigenvectors of a matrix.
 * Assumes that the matrix is symmetric.
 *
 * @param A the matrix.
 * @param n the dimension of the matrix.
 * @param u the eigenvectors.
 * @param w the eigenvalues.
 */
void eigen (double ** A, int n, double ** u, double * w);

/**
 * Computes the Euclidean distance between to vectors.
 *
 * @param v1 the first vector.
 * @param v2 the second vector.
 */
double distance (double * v1, double * v2);

/**
 * Normalizes the vector.
 *
 * @param vector the vector.
 */
void normalize (double * v);

/*
 * Performs matrix-vector multiplication.
 * Specialized for the 3D case.
 * This function is a helper for the RMSD calculation.
 *
 * @param a the matrix.
 * @param v the vector.
 * @param vr the resulting vector.
 */
void matrix_vector_mult3D(double ** a, double * v, double * vr);

/*
 * Performs matrix-matrix multiplications.
 * Specialized for the 3D case.
 * This function is a helper for the RMSD calculation.
 *
 * @param a the first matrix.
 * @param b the second matrix.
 * @param c the result matrix.
 */
void matrix_matrix_mult3D(double ** a, double ** b, double ** c);

/**
 * Calculates the minimized RMSD between two molecular conformations.
 *
 * @param ref_orig the reference conformation.
 * @param fit_orig the fit conformation.
 * @param atoms the number of atoms.
 * @return the fit rmsd.
 *
 *transform is to be a 16 float array to be filled with a 4x4 matrix in the following order
 *0  4  8 12 
 *1  5  9 13
 *2  6 10 14
 *3  7 11 15
 *
 */
double * min_rmsd (double ** ref_orig, double ** fit_orig, int atoms, double* transform);

#endif



