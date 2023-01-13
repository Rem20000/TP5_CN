#include<cblas.h>
#include<stdio.h>

//## Programme main##

int main(){

	double a[]={1.0,2.0,3.0,4.0};
	double  b[]={5.0,2.0,2.0,1.0};
	double c[4];

// ## Multiplication matricielle en utilisant dgemm##

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2,2,2,1.0,a,2,b, 2, 0.0,c,2);

	printf("Le résultat obtenu est:");
	for(int i=0;i<4;i++){
		printf("%lf",c[i]);
	}

	return 0;

}
