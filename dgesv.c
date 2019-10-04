#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "mkl_lapacke.h"

double *generate_matrix(int size,int seed)
{
    int i;
    double *matrix = (double *)malloc(sizeof(double) * size * size);
    srand(seed);

    for (i = 0; i < size * size; i++)
    {
        matrix[i] = rand() % 100;
    }

    return matrix;
}

//create an id matrix
void IdMatrix(int size, double** IdMat) {
//fill in the matrix
for (int i=0;i<size;i++) {
    for (int j=0;j<size;j++) {
      if (i==j) {
        IdMat[i][j] =  1;
      } else {
        IdMat[i][j] = 0;
      }
    }
  }
}

void print_matrix(const char *name, double *matrix, int size)
{
    int i, j;
    printf("matrix: %s \n", name);

    for (i = 0; i < size; i++)
    {
            for (j = 0; j < size; j++)
            {
                printf("%f ", matrix[i * size + j]);
            }
            printf("\n");
    }
}


int check_result(double *bref, double *b, int size) {
    int i;
    for(i=0;i<size*size;i++) {
        if (abs(bref[i]-b[i])>0.000005) return 0;
    }
    return 1;
}

//convert vector to Matrix
void vectorToMatrix(double* vect,double** matRes, int size) {
  for (int i=0;i<size;i++) {
    for (int j=0;j<size;j++) {
      matRes[i][j] = vect[i * size + j];
    }
  }
}



//matrix dilation
void dilatation (double **double_ttMat, int int_line, int int_N, double double_k) {
  int j;   //variable for column loop
    for (j=0;j<int_N;j++) {
        double_ttMat[int_line][j] = double_k * double_ttMat[int_line][j];
      } 
    }



//matrix transvection
void transvection(double **double_ttMat, int int_N, double double_k,int int_line1,int int_line2) {  
  int j;    
    for (j=0;j<int_N;j++) {
	      double_ttMat[int_line1][j] = double_ttMat[int_line1][j] + (double_k * double_ttMat[int_line2][j]);   
      } 
    }

//method of gauss jordan
//add zeros below the pivot
void underPivot(double **double_ttMat,double**double_ttMatId, int int_N, int int_indicePivot) {
  int i;
  double double_k;  //dilation coefficient
  double double_valeur;    //transvection coefficient

  //if the pivot is equal to 0 we quit
  if (double_ttMat[int_indicePivot][int_indicePivot] == 0) {
    exit(-1);
    //otherwise we define it as 1
  } else {
    double_k = 1 / (double_ttMat[int_indicePivot][int_indicePivot]);
  if(int_indicePivot == int_N - 1) {   //if the pivot is the first or last cell of the matrix, we use the dilation to define it as 1
    //computations on the matrix and the id matrix
    dilatation(double_ttMat,int_indicePivot,int_N,double_k);
    dilatation(double_ttMatId,int_indicePivot,int_N,double_k);  
  } else {  //otherwise we dilate and we write zeros below
    double_k = 1 / (double_ttMat[int_indicePivot][int_indicePivot]);  //to define the pivot as 1
    dilatation(double_ttMat,int_indicePivot,int_N,double_k);
    dilatation(double_ttMatId,int_indicePivot,int_N,double_k);
    //write zeros below the pivot
      for(i=(int_indicePivot + 1);i<int_N;i++) {  //we start at indicePivot + 1 to be below the diagonal
        double_valeur = - (double_ttMat[i][int_indicePivot]);
        transvection(double_ttMat,int_N,double_valeur,i,int_indicePivot);   
        transvection(double_ttMatId,int_N,double_valeur,i,int_indicePivot);   
      }
    }
  }
  
}

//add zeros above the pivot
void abovePivot(double **double_ttMat, double**double_ttMatId, int int_N, int int_indicePivot) {
  int i;
  double double_valeur;    //transvection coefficient
  //write zeros above the pivot
  for(i=int_indicePivot;i>0;i--) {  //we start by the pivot at the end of the diagonal
    double_valeur = - (double_ttMat[i-1][int_indicePivot]);  //we decrease by one to be above the diagonal
    transvection(double_ttMat,int_N,double_valeur,i-1,int_indicePivot);
    transvection(double_ttMatId,int_N,double_valeur,i-1,int_indicePivot);
    //we use the result matrix of the precedent loop
  }
}

 
//compute the reverse of the ttMat matrix, the result is stocked in ttMatResId 
void inverseG(double **double_ttMat,double **double_ttMatId, int int_N) {
    
  int i;

  for (i=0;i<int_N;i++) {
    underPivot(double_ttMat, double_ttMatId,int_N,i);  //lower triangularization     
  }
 
  for(i=(int_N)-1;i>0;i--) {
    abovePivot(double_ttMat,double_ttMatId,int_N,i);  //upper triangularization
  }
 
}

 

//convert matrix to vector
void MatrixToVector(double** mat,double* vectRes, int size) {
  for (int i=0;i<size;i++) {
    for (int j=0;j<size;j++) {
      vectRes[i * size + j] = mat[i][j];
    }
  }
}

//matrix product sotcked in matRes
void MatrixProduct (double**mat1, double**mat2, double**matRes, int size) {
  int int_tmp; //to temporarily stock the value of the scalar product of each cell
  for (int i=0;i<size;i++) {
    for (int j=0;j<size;j++) {
      for (int k=0;k<size;k++) {  //we compute the scalar product between a line of the first matrix and a column of the second matrix
    	matRes[i][j] += mat1[i][k] * mat2[k][j];
      }
    }
  }
}



    void main(int argc, char *argv[])
    {

        int size = atoi(argv[1]);

        double *a, *aref;
        double *b, *bref;

        //initialization of an empty matrix to stock the result of vectorToMatrix for matrix a
        double** matA = malloc(size * sizeof(double *));
        for (int i=0; i<size;i++) {
            matA[i] = malloc(size * sizeof(double));
        }
        //initialization of an empty matrix to stock the result of vectorToMatrix for matrix b
        double** matB = malloc(size * sizeof(double *));
        for (int i=0; i<size;i++) {
            matB[i] = malloc(size * sizeof(double));
        }
        //intialization of an empty matrix to stock the id matrix
        double** matId = malloc(size * sizeof(double *));
        for (int i=0; i<size;i++) {
            matId[i] = malloc(size * sizeof(double));
        }
        //to stock the result of the system 
        double ** matRes = malloc(size * sizeof(double *));
        for (int i=0; i<size;i++) {
            matRes[i] = malloc(size * sizeof(double));
        }

        //to stock X
        double*X=malloc(size * size * sizeof(double));

        a = generate_matrix(size,2);
        aref = generate_matrix(size,2);        
        b = generate_matrix(size,5);
        bref = generate_matrix(size,5);
        

        //print_matrix("A", a, size);
        //print_matrix("B", b, size);

        // Using MKL to solve the system
        MKL_INT n = size, nrhs = size, lda = size, ldb = size, info;
        MKL_INT *ipiv = (MKL_INT *)malloc(sizeof(MKL_INT)*size);

        clock_t tStart = clock();
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, n, nrhs, aref, lda, ipiv, bref, ldb);
        printf("Time taken by MKL: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);

        //beginning of the sequencial code
        
        //first, I convert vectors into matrices
        vectorToMatrix(a,matA,size);
        vectorToMatrix(b,matB,size);
        //I create the id matrix
        IdMatrix(size,matId);
        //I compute the reverse matrix of A
        inverseG(matA,matId,size); //the reverse matrix is stocked in the variable matA

        tStart = clock();    
        //I complete a matrix product between the reverse matrix of A and the b matrice to find X
        MatrixProduct(matId,matB,matRes,size);
        //I convert X into a vector
        MatrixToVector(matRes,X,size);
        printf("Time taken by my implementation: %.2fs\n", (double)(clock() - tStart) / CLOCKS_PER_SEC);
        
        if (check_result(bref,X,size)==1)
            printf("Result is ok!\n");
        else    
            printf("Result is wrong!\n");
        
        //print_matrix("X", X, size);
        //print_matrix("Xref", bref, size);

}
    
