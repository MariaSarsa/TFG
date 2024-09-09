#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>

#define N 10000
#define EPS 0.7
#define BETA 0.5
#define PASOS (2 * N)
#define SIMUS 50



//Función para generar un número random entre 0 y 1
double rand_double() {
    return rand() / (double) RAND_MAX;
}


int main() {
    // Semilla para la comprobación de resultados
    srand(123);
    //srand((unsigned int) time(NULL));

    // Parámetros
    const double TOLERANCIA = 1e-9;
    char filename3[100];

    int N_Su0 = N, N_Sv0 = N;
    int j;

    double eps = 0.7;
    double eps_delta = 0.9;

    double **Time = (double **)malloc(SIMUS * sizeof(double *));
    int **N_Su = (int **)malloc(SIMUS * sizeof(int *));
    int **N_Sv = (int **)malloc(SIMUS * sizeof(int *));

    FILE *file1 = fopen("simulaciones.txt", "w");
    FILE *file2 = fopen("promedios.txt", "w");
    sprintf(filename3, "individuos_Leaky_N%i_S%i.txt", N, SIMUS);
    FILE *file3 = fopen(filename3, "w");

    fprintf(file3, "Eps Simu Individuo Vacuna_Status Eficacia Tiempo Infeccion_Status\n");

    if (file1 == NULL || file2 == NULL || file3 == NULL) {
        printf("Error al abrir el archivo para escritura.\n");
        return 1;
    }

    for (int i = 0; i < SIMUS; i++) {
        Time[i] = (double *)malloc((PASOS + 1) * sizeof(double));
        N_Su[i] = (int *)malloc((PASOS + 1) * sizeof(int));
        N_Sv[i] = (int *)malloc((PASOS + 1) * sizeof(int));
    }

    double *Time_ave = (double *)malloc((PASOS + 1) * sizeof(double));
    double *N_Su_ave = (double *)malloc((PASOS + 1) * sizeof(double));
    double *N_Sv_ave = (double *)malloc((PASOS + 1) * sizeof(double));

    //Bucle barrido eficacia promedio
    while (eps <= 1 + TOLERANCIA){

        // Condiciones iniciales
        for (int i = 0; i < SIMUS; i++) {
            Time[i][0] = 0.0;
            N_Su[i][0] = N_Su0;
            N_Sv[i][0] = N_Sv0;
        }
        for (int i = 0; i <= PASOS; i++) {
            Time_ave[i] = 0.0;   //Promedio tiempo
            N_Su_ave[i] = 0.0;   //Promedio CONTROL
            N_Sv_ave[i] = 0.0;
        }

    /** ALGORITMO DE GILLESPIE **/
        for (int i = 0; i < SIMUS; i++) {
            j=0;
            while(j != PASOS){

                // Suma de propensidades total
                double P = BETA * N_Su[i][j] + BETA * (1 - eps) * N_Sv[i][j];

                // Tiempo al siguiente evento
                double r1 = rand_double();
                double tau = (1 / P) * log(1 / r1);

                if (P < TOLERANCIA){ /** A este if solo se entra cuando todos los de CONTROL se han infectado y eps = 1**/
                    Time[i][j + 1] = Time[i][j];
                    N_Su[i][j + 1] = N_Su[i][j];
                    N_Sv[i][j + 1] = N_Sv[i][j];

                    fprintf(file3, "%f %i %i %i %f %f %i\n",eps, i+1, j + 1, 1, eps, Time[i][j + 1], 0);

                    j++;
                }

                else if ( isinf(tau) == 0) {

                    // Guardamos tiempo
                    Time[i][j + 1] = Time[i][j] + tau;
                    //printf("Tau: %f Tiempo %f:\n", tau, Time[i][j + 1]);

                    // Seleccionamos el próximo evento
                    double Pu = (BETA * N_Su[i][j]) / P;
                    double r2 = rand_double();

                    //Infectado de la cohorte de control
                    if (r2 < Pu) {
                        N_Su[i][j + 1] = N_Su[i][j] - 1;
                        N_Sv[i][j + 1] = N_Sv[i][j];

                        fprintf(file3, "%f %i %i %i %f %f %i\n", eps,i+1, j + 1, 0, 0.0, Time[i][j + 1], 1);

                    } else {//Infectado de la cohorte vacunada
                        N_Su[i][j + 1] = N_Su[i][j];
                        N_Sv[i][j + 1] = N_Sv[i][j] - 1;

                        fprintf(file3, "%f %i %i %i %f %f %i\n",eps,i+1, j + 1, 1, eps, Time[i][j + 1], 1);
                    }
                    j++;
                }
            }
        } /** FIN DE 1 BARRIDO DE eps **/

        eps += eps_delta;

    }

    //Liberación de memoria
    for (int i = 0; i < SIMUS; i++) {
            free(Time[i]);
            free(N_Su[i]);
            free(N_Sv[i]);
        }

        free(Time);
        free(N_Su);
        free(N_Sv);
        free(Time_ave);
        free(N_Su_ave);
        free(N_Sv_ave);


    fclose(file1);
    fclose(file2);
    fclose(file3);

    printf("Simulaciones guardadas.\n");


    return 0;
}
