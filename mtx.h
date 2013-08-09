#ifndef MTX_H_
#define MTX_H_
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>


class mat
{
public:


    
	//constructors
	mat();//NULL constructor
	mat (int i,int j,double* in);
	mat(int i,int j);
	mat(const char *);

	//deconstructor
	~mat ();
	int write(const char *);


	double get(int i, int j);
	void set(int i, int j,double value);

	mat mult(mat B);
	//mat operator*(const mat B);
	void print();

	mat minus(mat sub);
	//mat operator-(const mat sub);
	mat plus(mat sub);
	//mat operator+(const mat plus);

	void replaceColumn(int col, mat matr);
	void replaceRow(int row, mat matr);
	mat transpose();
	/**
	Normalizes the vector using the frobenius norm
	**/
	mat normalize();
	/**
	Returns the frobenius norm of the vector
	**/
	double norm();
	/**
	Returns a deep copy of the current matrix
	**/
	mat deepCopy();
	/*
	Swaps row r1 and row r2
	*/
	void swapRows(int r1,int r2);

	//h and w are the dimensions of m
	//m is the input matrix of linear equations and the destination for the solution
	//
	mat solver(double);
	mat lstSqrs(mat y);
	//mat operator/(mat y);
	mat append(mat appende);
	mat subMat(int sRow,int pRow,int sCol,int pCol);
	double sumsqrs();
    double * getData();
    int getM();
    int getN();
    
private:
    double* data;
	int m;
	int n;

};






#endif /* MTX_H_ */
