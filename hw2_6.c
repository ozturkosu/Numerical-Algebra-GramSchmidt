/* This is the template file for HW2_6. You need to complete the code by 
 * implementing the parts marked by FIXME. 
 *
 * Submit your completed code to the TA by email. Also submit the plots 
 * and your conclusions of the comparative study.
 *
 * To compile, use the UNIX command 
 *      make
 * in the source directory, and it would then invoke the compilation command 
 * with the included makefile.
 *
 * This will generate an executable hw2_6. The -g option is optional and 
 * is needed only if you need to debug the program in a debugger such as ddd.
 * The -Wall option would enable compiler's warning messages.
 *
 * Run the program with command
 *      ./hw2_6
 * which would generate a M-file, which you can use to generate the plots by
 * issueing the UNIX command 
 *      make plot
 * in the source directory. It will then generate two PDF files, which you
 * should submit to the TA along with your source code and your conclusions.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Define types for matrix and vector. 
 * Access the (i,j)-entries of a matrix A using A.vals[i][j], and 
 * access the ith entries of a vector x using x.vals[i], where 
 * i and j start from 0.
 */

typedef struct {
    float** vals;
    int m;
    int n;
} Matrix;

typedef struct {
    float*  vals;
    int m;
} Vector;

/* Prototypes for local functions. You need to implement these functions. */
float check_orthogonality(const Matrix Q);
void  dgs(Matrix Q, Matrix R, const Matrix A);
void  cgs(Matrix Q, Matrix R, const Matrix A);
void  mgs(Matrix Q, Matrix R, const Matrix A);

/* Define prototypes of helper functions */
Matrix allocate_matrix(int m, int n);
Vector allocate_vector( int m);
void   deallocate_matrix( Matrix);
void   deallocate_vector( Vector);

void  gen_hilbert(Matrix H);
void  write_vec(FILE *fid, char *name, const float *vec, int n);
float tictoc (int n);

#define N   12

/*************************************************************
 *
 * FUNCTION: hw2_driver
 *
 * This is the driver function
 *************************************************************/

int main(int argc, char **argv) {
    Matrix R, H, Q, CheckH;
    float digits_cgs[N], digits_mgs[N], digits_dgs[N], times_cgs[N], times_mgs[N];
    float times_dgs[N];
    int i, j;
    FILE *fid;
    //int niter = 10000;

    H = allocate_matrix(N, N);
    Q = allocate_matrix(N, N);
    R = allocate_matrix(N, N);

    CheckH = allocate_matrix(N,N) ;
   
    for (i=2; i<=N; ++i) {

        /* Use first i rows and first i columns of H, Q, and R */
        H.m = H.n = Q.m = Q.n = R.m = R.n = i;

        /* Generate Hilbert matrix */
        gen_hilbert(H);
        /*%%%%%%%%%%%%%%%%%% Run classical Gram-Schmidt */

        //Print H matrix

        /* print the matrix H  */
        printf("H = \n");
        for(int k = 0; k< i; k++) {
            for(int l = 0; l < i; l++) {

                printf("%9.6g ", H.vals[k][l]);
            }
            printf("\n");
        }
        printf("\n");


        tictoc(1);

    
        cgs(Q, R, H);

         /* print the matrix H  */
        printf("Q = \n");
        for(int k = 0; k< i; k++) {
            for(int l = 0; l < i; l++) {

                printf("%9.6g ", Q.vals[k][l]);
            }
            printf("\n");
        }
        printf("\n");


         /* print the matrix H  */
        printf("R = \n");
        for(int k = 0; k< i; k++) {
            for(int l = 0; l < i; l++) {

                printf("%9.6g ", R.vals[k][l]);
            }
            printf("\n");
        }
        printf("\n");
        
        //Lets Check if Q*R = H ????
        float sumM=0;

        for ( int k = 0; k < i; k++){
            for (int l = 0; l < i; l++){
                sumM =0 ;
                for ( int m = 0; m < i; m++){
                    sumM = sumM + Q.vals[k][m]*R.vals[m][l] ;
                }
                CheckH.vals[k][l] = sumM;
            }
        }

        /* print the matrix CheckH  */
        printf("CheckH = \n");
        for(int k = 0; k< i; k++) {
            for(int l = 0; l < i; l++) {

                printf("%9.6g ", CheckH.vals[k][l]);
            }
            printf("\n");
        }
        printf("\n");


        times_cgs[i - 2] =  tictoc(2) ;
        digits_cgs[i - 2] = check_orthogonality(Q);
        /*%%%%%%%%%%%%%%%%%% Run modified Gram-Schmidt */

        tictoc(1);

        /* Repeat for many times for more accurate timing */
        mgs(Q, R, H);
        
	   
        times_mgs[i - 2] =  tictoc(2) ;
        digits_mgs[i - 2] = check_orthogonality(Q);
        /*%%%%%%%%%%%%%%%%%% Run float Gram-Schmidt */

        tictoc(1);
       

        dgs(Q, R, H);
        
        times_dgs[i - 2] =  tictoc(2) ;
        digits_dgs[i - 2] = check_orthogonality(Q);
    }

    deallocate_matrix(H);
    deallocate_matrix(Q);
    deallocate_matrix(R);
    deallocate_matrix(CheckH) ;

    /* Write out results into a Matlab file */
    fid = fopen("results.m", "w");
    write_vec(fid, "digits_cgs ", digits_cgs, N-1);
    write_vec(fid, "digits_mgs ", digits_mgs, N-1);
    write_vec(fid, "digits_dgs ", digits_dgs, N-1);
    write_vec(fid, "times_cgs", times_cgs, N-1);
    write_vec(fid, "times_mgs", times_mgs, N-1);
    write_vec(fid, "times_dgs", times_dgs, N-1);
    fclose(fid);
    return 0;
}

/*************************************************************
 *
 * FUNCTION: write_vec
 *
 * Write out a vector into file
 *************************************************************/

void  write_vec(FILE *fid, char *name, const float *vec, int n) {
    int i;

    fprintf(fid, "%s=[", name);
    for (i=0; i<n; ++i) {
        fprintf(fid, "%g; ", vec[i]);
    }
    fprintf(fid, "];\n");
}

/*************************************************************
 *
 * FUNCTION: check_orthogonality
 *
 * Check the orthogonality of a given unitary matrix.
 *************************************************************/
//Calculate digit accuracy -log10(|| I - Q'*Q||F)

float check_orthogonality(const Matrix Q) {
    float sqnrm=0.0;
    int n=Q.m, i,j,k;

    Matrix resultMatrix = allocate_matrix(n,n);

    //Calculate Q'*Q

    float sum=0;
    for ( i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            
            sum=0;
            for ( k = 0; k < n; k++){
               sum = sum + Q.vals[k][i] * Q.vals[k][j] ;
            }
            resultMatrix.vals[i][j] = sum;

        }
    }

    /* Compute the Frobenius norm of I-Q'*Q */
    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j) {
            
            /* FIXME */
        	if ( i == j ){
        		sqnrm = sqnrm + (1 - resultMatrix.vals[i][j]) * (1- resultMatrix.vals[i][j]) ;
        	}
        	else{
        		sqnrm = sqnrm + resultMatrix.vals[i][j] * resultMatrix.vals[i][j] ;
        	}	
        }
    }
	
	deallocate_matrix(resultMatrix);
    return (-log10(sqrt(sqnrm)));
}

/*************************************************************
 *
 * FUNCTION: dgs
 *
 * float Classical Gram-Schdmit Algorithm
 *************************************************************/

void  dgs(Matrix Q, Matrix R, const Matrix A) {
    int n=A.m, i,j,k;

    Matrix buf1 = allocate_matrix(n,n);
    Matrix buf2 = allocate_matrix(n,n);

    cgs(Q, buf1, A);
    cgs(Q, buf2, Q);

    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            for (k = 0; k < n; ++k)
                R.vals[i][j] = R.vals[i][j] + buf1.vals[i][k] * buf2.vals[k][j];

    deallocate_matrix(buf1);
    deallocate_matrix(buf2);
}

/*************************************************************
 *
 * FUNCTION: cgs
 *
 * Classical Gram-Schdmit Algorithm
 *************************************************************/

void  cgs(Matrix Q, Matrix R, const Matrix A) {

    int n=A.m, i,j,k;
    float v[N];

    float result=0;

    for (j=0; j<n; ++j) {
        /* v = A(:,j); */
        for (i = 0; i < n; ++i)
            v[i] = A.vals[i][j];



        for (i = 0; i < j; ++i) {
            
            /* R(i,j) = Q(:, i)' * A(:, j); */

        	for ( k = 0; k < n; ++k){
        		result = result + Q.vals[k][i]*A.vals[k][j] ;
        	}

        	R.vals[i][j] = result ;
            result=0;

            /* v = v - R(i, j) * Q(:,i); */

        	for (int k = 0; k < n; ++k){
        		v[k] = v[k] - R.vals[i][j]*Q.vals[k][i];
        	}


        }

        /* R(j, j) = norm(v, 2); */
       
        result = 0 ;
        for ( k = 0; k < n; ++k) {
        	result = result + v[k]*v[k] ;
        }

        R.vals[j][j] = sqrt(result) ;


        /* Q(:,j) = v / R(j, j); */
        for (i=0; i<n; ++i)
            Q.vals[i][j] = v[i] / R.vals[j][j];
    }
}

/*************************************************************
 *
 * FUNCTION: mgs
 *
 * Modified Gram-Schdmit Algorithm
 *************************************************************/

void  mgs(Matrix Q, Matrix R, Matrix A) {
    int n=A.m, i,j,k;
    float result=0;

    /* Q = A; */
    for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            Q.vals[i][j] = A.vals[i][j];

    for (i=0; i<n; ++i) {
        /* R(i,i) = norm (Q(:,i), 2); */
        /* FIXME */

        for (k = 0; k < n; ++k)
        {
        	/* code */
        	result = result + Q.vals[k][i]*	Q.vals[k][i];
        }

        R.vals[i][i] = sqrt(result) ;

        /* Q(:,i) = Q(:,i) ./ R(i,i); */
        for (k = 0; k < n; ++k)
            Q.vals[k][i] = Q.vals[k][i] / R.vals[i][i];
		
		result=0;
        for (j=i+1; j<n; ++j) {
            /* R(i,j) = Q(:,i)' * Q(:,j); */
            /* FIXME */

        	for (k = 0; k < n; ++k)
        	{
        		/* code */
        		result = result + Q.vals[k][i]*Q.vals[k][j];
        	}
        	R.vals[i][j] = result ;

            /* Q(:,j) = Q(:,j) - R(i,j) * Q(:,i); */
            /* FIXME */

            for (k = 0; k < n; ++k)
            {
            	/* code */
            	Q.vals[k][j] = Q.vals[k][j] - R.vals[i][j]*Q.vals[k][i] ;
            }
        }
    }
}

/*************************************************************
 *
 * FUNCTION: gen_hilbert
 *
 * Generating nxn Hilbert matrix.
 *************************************************************/

void  gen_hilbert(Matrix H) {
    int n = H.m, i, j;

    for (i=0; i<n; ++i) {
        for (j=0; j<n; ++j) {
            H.vals[i][j] = 1.0 / (1 + j + i);
        }
    }
}

/*-----------------------------------------------------------------------------
 * FUNCTION: allocate_matrix - Allocate memory for a given matrix.
 * --------------------------------------------------------------------------*/
/* Disclaimer: The approach used here is not the most efficient way for
 * implementing matrices, but it is adopted here for convenience.
 */
Matrix allocate_matrix(int m, int n) {
    Matrix A;
    float *ptmp;
    int i;

    A.m = m; A.n = n;
    ptmp = (float *)malloc( sizeof(float) * m * n);
    A.vals = (float **)malloc( sizeof(float **) * m);

    for (i=0; i<m; ++i) {
        A.vals[i] = ptmp + i*n;
    }

    return A;
}

/*-----------------------------------------------------------------------------
 * FUNCTION: deallocate_matrix - De-allocate memory for a matrix.
 * --------------------------------------------------------------------------*/
void deallocate_matrix( Matrix A) {
    free( A.vals[0]);
    free( A.vals); A.vals=NULL;
}


/*-----------------------------------------------------------------------------
 * FUNCTION: allocate_vector - Allocate memory for an m vector.
 * --------------------------------------------------------------------------*/
Vector allocate_vector( int m) {
    Vector x;

    x.m = m;
    x.vals = (float *)malloc( sizeof(float) * m);
    return x;
}

/*-----------------------------------------------------------------------------
 * FUNCTION: deallocate_vector - De-allocate the memory for an m vector. 
 * --------------------------------------------------------------------------*/
void deallocate_vector( Vector x) {
    free( x.vals); x.vals=NULL;
}

#include <time.h>
#include <sys/time.h>

/*-----------------------------------------------------------------------------
 * FUNCTION: tictoc - Initialize (1) or get (2) the elapsed time since in seconds
 * --------------------------------------------------------------------------*/
float tictoc (int n)
{
    float    y = 0.0;
    static    struct timeval start_time, end_time; 

    if (n == 1) { 
        gettimeofday(&start_time, NULL);
    }
    else if (n == 2) { 
        gettimeofday(&end_time, NULL);

        y = (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec-start_time.tv_usec)*1.e-6;
    }
    return (y); 
}
