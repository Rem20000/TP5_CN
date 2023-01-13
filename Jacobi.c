#include <stdio.h>
#include <math.h>
#include<cblas.h>
#include<lapack.h>

void Jacobi(double A[][10], double b[], double X[], int n) {
    double r = 0, s = 0, eps = 1e-15;              // r est le residus pour savoir l'itération "stop"
    int i, j, k = 0;
    for (i = 0; i < n; i++)
        X[i] = rand();            // valeurs aléatoires initialisées au début pour évaluer le vecteur solution x
  do {
        for (i = 0; i < n; i++) {
            s = 0;
            for (j = 0; j < n; j++) {
                if (j != i)
                    s += A[i][j] * X[j];
            }
            X[i] += (b[i] - s) / A[i][i];
        }
        r = 0;
        for (i = 0; i < n; i++)
            r += (b[i] - A[i][i] * X[i]) * (b[i] - A[i][i] * X[i]);
        r = sqrt(r) / norm(b);     // le residus
    k++;
    } while (r > eps);
    printf("k = %d\nr = %lf\n", k, r);     // les résultats obtenus
}

int main() {
    double A[10][10], b[10], X[10];
    int n, i, j;
    printf("entrer la taille de A : ");
    scanf("%d", &n);
    printf("Enter les éléments de A : \n");
    for (i = 0; i < n; i++) {for (j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }
    printf("Enter les elements du second membre b: \n");
    for (i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }
    Jacobi(A, b, X, n);
    printf("X = [");            // afficher la solution approchée après l'ensemble des itérations 
    for (i = 0; i < n; i++) {
        printf("%lf ", X[i]);
    }
   
    return 0;
}




