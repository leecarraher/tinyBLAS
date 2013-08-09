/*Copyright 2010 Lee Carraher. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY Lee Carraher ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Lee Carraher OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those of the
authors and should not be interpreted as representing official policies, either expressed
or implied, of Lee Carraher.
*/



#include "mtx.h"


mat::mat (int i,int j,double* in) {
     m =i;
     n=j;
     data=in;
}

mat::~mat(){
	//both of these cause errors?
	//delete data;
    printf("freeing: %x\n",data);
    free(data);
	data=NULL;
}


mat::mat(int i,int j){
     data = (double*) malloc(sizeof(double)*i*j);
     m = i;
     n = j;
     j=i*j;
     for(i=0;i<j;i++)data[i]=0.0;
}

mat::mat(const char * infile){
     int length = 0;
     double d;
     FILE    *fptr;
     /* Open the file */
     if ((fptr = fopen(infile,"r")) == NULL) {
          fprintf(stderr,"Unable to open data file\n");
          printf("%s",infile);
          exit(0);
     }

     int row;
     int col;
     fscanf(fptr,"%i",&row);
     fscanf(fptr,"%i",&col);

     //printf("%i\n",row);
     //printf("%i\n\n\n",col);


     data = (double*) malloc(sizeof(double)*row*col);

     /* Read as many points as we can */
     while (fscanf(fptr,"%lf",&d) == 1) {
          if(length>row*col){
               fprintf(stderr,"Malformed Data File\n");
               exit(0);
          }
          data[length] = d;
          length++;
     }

     fclose(fptr);
     m = row;
     n = col;

}

mat::mat(){

}


double mat::get(int i, int j){return data[i*n+j];}
void mat::set(int i, int j,double value){data[i*n+j]=value;}

int mat::write(const char * outfile){

     int length = 0;
     double d;
     FILE    *fptr;
     /* Open the file */
     if ((fptr = fopen(outfile,"wt")) == NULL) {
          fprintf(stderr,"Unable to open data file\n");
          printf("%s",outfile);
	  return 1;
     }

	fprintf (fptr, "%i\n",m);
	fprintf (fptr, "%i\n",n);
	int i;
	for(i=0;i<m*n;i++)fprintf (fptr, "%f\n",data[i]);
	fclose (fptr);
	return 0;



}


//naive multiply, our matrices aren't big enough to exploit strassen's algo.
//{ai x aj},{bi x bj} are the dimensions of A and B respectively, C is ai x bj
//A and B are the input matrices, and C is a result matrix
mat mat::mult(mat B)
{
     if(n!=B.m){
          printf("Incompatible Dimensions %ix%i * %ix%i",m,n,B.m,B.n);
          return mat(0,0,NULL);
     }

     mat* ret = new mat(m,B.n);

     double * newdata = ret->data;//(double*) malloc(sizeof(double)*m*B.n);
     //mat rtr = mat(m,B.n,new double[m*B.n]);

     register int i,j,k;
     double sum;
     for(i = 0; i<m;i++)
     {
          for(j = 0; j<B.n;j++)
          {
               sum = 0.0;
               for(k = 0; k<n;k++){
                    sum+=data[i*n+k]*B.data[k*B.n+j];
               }
               newdata[i*B.n+j]=sum;
          }
     }
     return *ret;//mat(m,B.n,newdata);
}


//mat mat::operator*( mat B){
//     return mult(B);
//}

void mat::print()
{
     register int i,j;
     for(i=0;i<m;i++){
          for(j=0;j<n;j++)printf("%f ",data[i*n+j]);
          printf("\n");
     }
     printf("\n");
}


mat mat::minus(mat sub){


     
     if(m!=sub.m || n!=sub.n)
     {
          printf("incompatible dimensions");
          return mat(0,0);
     }
     mat* ret = new  mat(m,n);
     double * newData = ret->data;//double * newData = (double*) malloc(sizeof(double)*m*n);
     register int i;
     for(i=0;i<m*n;i++)newData[i] = data[i]-sub.data[i];

     return *ret;//mat(m,n,newData);

}
//mat mat::operator-( mat B){
//     return minus(B);
//}


mat mat::plus(mat sub){
     //double *newData;
     if(m!=sub.m || n!=sub.n)
     {
          printf("incompatible dimensions");
          return mat(0,0);
     }
     register int i;
     mat* ret =  new mat(m,n);
     double * newData = ret->data;//double * newData = (double*) malloc(sizeof(double)*m*n);
     //newData= (double*) malloc(sizeof(double)*m*n);
     for(i=0;i<m*n;i++)newData[i] = data[i]+sub.data[i];

     return *ret;//mat(m,n,newData);

}


//mat mat::operator+(const mat B){
//     return plus(B);
//}

void mat::replaceColumn(int col, mat matr)
{
     if(matr.m>m){
          printf("error imcompatible dimensions");
          return;
     }
     int i;
     for(i=0;i<matr.m;i++)data[i*n+col]=matr.data[i];
}

void mat::replaceRow(int row, mat matr)
{
     if(matr.n>n){
          printf("error imcompatible dimensions");
          return;
     }
     //int i;
     //double * mDat = matr.data;
     //for(i=0;i<(matr.n);i++)data[row*n+i]=mDat[i];

     memcpy(&data[row*n],matr.data,m*sizeof(double));
}


mat mat::transpose()
{    register int i,j;
    

     mat *ret = new mat(m,n);
     double * newData = ret->data;//double * newdata = (double*) malloc(sizeof(double)*m*n);
     //newData= (double*) malloc(sizeof(double)*m*n);
     //memcpy(newData,data,n*m*sizeof(double));
     for(i=0;i<m;i++)
          for(j=0;j<n;j++)
               newData[j*m+i] = data[i*n+j];

     //m^=n;
     //n^=m;
     //m^=n;
     
     return *ret;//mat(n,m,newdata);
}

/**
Normalizes the vector using the frobenius norm
**/
mat mat::normalize()
{
     register int i,j;
     
     double c= norm();
     
     mat *ret = new mat(m,n);
     double * newData = ret->data;//double * newData = (double*) malloc(sizeof(double)*m*n);
     c = 1.0/c; //mults are faster

     
     for(i = 0; i < m*n;i++)newData[i]=data[i]*c;
     return *ret;//mat(m,n,newData);
}

/**
Returns the frobenius norm of the vector
**/
double mat::norm(){
     register int i,j;
     double c=0.0;
     for(i=0;i<m*n;i++)c+=data[i]*data[i];
     return sqrt(c);
}

/**
Returns a deep copy of the current matrix
**/
mat mat::deepCopy()
{
     mat *ret = new mat(m,n);
     double * newData = ret->data;//newData= (double*) malloc(sizeof(double)*m*n);
     memcpy(newData,data,n*m*sizeof(double));
     return mat(m,n,newData);

}

/*
Swaps row r1 and row r2
*/
void mat::swapRows(int r1,int r2)
{
          double* temp = (double*) malloc(sizeof(double)*n);
          memcpy(temp,&data[r1*n],n*sizeof(double));
          memcpy(&data[r1*n],&data[r2*n],n*sizeof(double));
          memcpy(&data[r2*n],temp,n*sizeof(double));
          free(temp);
          temp=NULL;
}

//h and w are the dimensions of m
//m is the input matrix of linear equations and the destination for the solution
//
mat mat::solver(double eps=10e-32)
{
     if(n-m <=0)
     {
          printf("overdetermined system! try using least squares method (lstSqrs()).");
          return mat(0,0);
     }
     int x,y, y2,maxrow;
     double c;

     for(y=0;y<m;y++)
     {
         maxrow = y;

         for(y2=y+1;y2< m;y2++)// find pivot
          if(abs( data[y2*n+y]) > abs(data[maxrow*n+y]))maxrow = y2;

         swapRows(y,maxrow);


         if(abs(data[maxrow*n+y]) < eps){//check singularity
          printf("singular matrix: mat[%i,%i]=%f\n",maxrow,y,data[y*n+y]);
           return mat(0,0);
         }
         //eliminate column y
         for(y2=y + 1;y2< m;y2++)
         {
          c = data[y2*n+y] / data[y*n+y];
             for(x=y;x<n;x++)data[y2*n+x]-= data[y*n+x] * c;
         }

      }

      for(y=m-1;y>-1;y--)
      {
         c = data[y*n+y];
         for(y2=0;y2<y;y2++)
               for(x=n - 1;x> y - 1; x--) data[y2*n+x] -=  data[y*n+x] *  data[y2*n+y] / c;
         data[y*n+y] /= c;
         for(x=m;x< n;x++ ) data[y*n+x] /= c;
      }
      return subMat(0,m,n-(n-m),n);//return the last column of solutions
}

mat mat::lstSqrs(mat y){
//solve overdetermined systems
//we will use A^T*A*x=A^T*b
//A = A.subMat(0,A.m,1,A.n);
     mat me = deepCopy();
     mat ata=transpose().mult(me);
     mat atb=transpose().mult(y);
     mat aug = ata.append(atb);
     return aug.solver();
}


//mat mat::operator/(const mat B){
//     return lstSqrs(B);
//}

mat mat::append(mat appende)
{

     if(appende.m!=m){
          printf("error imcompatible dimensions: %ix%i + %ix%i \n",m,n,appende.m,appende.n);
          return mat(0,0);
     }
     int nn = n+appende.n;
     int i,j;
     mat *ret = new mat(m,nn);
     double * newData = ret->data;//double * newData = (double*) malloc(sizeof(double)*m*nn);
     
     //fill with old data
     for(i=0;i<m;i++)for(j=0;j<n;j++)newData[i*nn+j]=data[i*n+j];
     //fill with new data starting at m
     for(i=0;i<m;i++)for(j=n;j<nn;j++)newData[i*nn+j]=appende.data[i*appende.n+(j-n)];

     return *ret;//mat(m,nn,newData);

}

int mat::getM(){
    return m;
}
int mat::getN(){
    return n;
}

double * mat::getData(){
    return data;
}

mat mat::subMat(int sRow,int pRow,int sCol,int pCol)
{


     int nm = (pRow-sRow);
     int nn = (pCol-sCol);
     mat *ret = new mat(nm,nn);
     double * newData = ret->data;//double * newData = (double*) malloc(sizeof(double)*nm*nn);
     int i,j,ni,nj;
     ni=0;
     for(i = sRow;i<pRow;){
          nj = 0;
          for(j = sCol;j<pCol;){
               newData[ni*nn+nj]=data[i*n+j];
               nj++;
               j++;
          }
          ni++;
          i++;
     }
     return *ret;//mat(nm,nn,newData);
}



double mat::sumsqrs(){
     int i;
     double sum = 0.0;
     for(i=0;i<m*n;i++)sum+=data[i]*data[i];
     return sum;
}
