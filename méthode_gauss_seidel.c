#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<cblas.h>
#include<lapack.h>


void Seidel(double A[][10], double b[], double X[], int n) {   //résolution par une méthode itérative d'un système de n inconnues AX=b
    double r = 0, s1 = 0, s2 = 0, eps = 1e-15;     //eps est la tolérence de la méthode, et r est le residus
    int i, j, k = 0;
    for (i = 0; i < n; i++)
        X[i] = rand();                 //On initialise le vecteur solution initial par des valeurs aléatoires choisies par le système
    do {
        for (i = 0; i < n; i++) {                  //boucle pour calculer la solution
            s1 = 0;
            s2 = 0;
            for (j = 0; j < i; j++) {
                if (j != i - 1)
                    s1 += A[i][j] * X[j];
            }
            for (j = i + 1; j < n; j++) {
                if (j != i)
                    s2 += A[i][j] * X[j];
              }
            X[i] += (b[i] - s1 - s2) / A[i][i];
        }
        r = 0;
        for (i = 0; i < n; i++)
            r += (b[i] - A[i][i] * X[i]) * (b[i] - A[i][i] * X[i]);
        r = sqrt(r) / norm(b);
        k++;
    } while (r > eps);
    printf("k = %d\nr = %lf\n", k, r);
}
int main() {
    double A[10][10], b[10], X[10];
    int n, i, j;
    printf("Entrer la taille de A: ");
    scanf("%d", &n);
    printf("Entrer les elements de A: \n");
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            scanf("%lf", &A[i][j]);
        }
    }
  printf("Entrer b: \n");
    for (i = 0; i < n; i++) {
        scanf("%lf", &b[i]);
    }
  
  //Application de la méthode
    gauss_seidel(A, b, X, n);
    printf(" la solution finale est X = [");
    for (i = 0; i < n; i++) {
        printf("%lf ", X[i]);
    }
    
    return 0;
}
