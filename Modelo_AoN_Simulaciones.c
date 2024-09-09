/*
 * Programa de Simulación de Infección en una Población con Protección Heterogénea Máxima
 *
 * Descripción:
 * Este programa realiza simulaciones de la dinámica de infección en una población de individuos,
 * considerando una cohorte de control y una cohorte vacunada AoN. Utiliza el algoritmo de Gillespie
 * para simular la evolución de la infección a lo largo del tiempo. Los resultados de las simulaciones
 * se guardan en archivos de texto para su análisis posterior.
 *
 * Parámetros:
 * - N: Número total de individuos en la población.
 * - BETA: Tasa de infección.
 * - PASOS: Número total de pasos en la simulación (2*N).
 * - SIMUS: Número de simulaciones a realizar.
 * - eps: Eficacia de la vacuna.
 * - eps_delta: Incremento en la eficacia de la vacuna para el barrido.
 *
 * Archivos:
 * - individuos_AoN_N{N}_S{SIMUS}.txt: Archivo para guardar detalles de los individuos y su estado
 *   a lo largo del tiempo en cada simulación.
 *
 * Funciones:
 * - rand_double: Genera un número aleatorio de tipo double entre 0 y 1.
 *
 */




#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N 1000
#define BETA 0.5
#define PASOS 2*N
#define SIMUS 100

double rand_double() {
    return (double) rand() / RAND_MAX;
}

int main() {
    // Semilla para comprobación de resultados
    srand(123);
    //srand((unsigned int) time(NULL));

    const double TOLERANCIA = 1e-9;
    char filename3[100];

    // Parámetros
    int N_Su0 = N, N_Sv0 = N;
    int j;
    int pasos;

    double P = 0;
    double eps = 0.0;
    double eps_delta = 0.01;

    //Ficheros
    FILE *file1 = fopen("simulaciones.txt", "w");
    FILE *file2 = fopen("promedios.txt", "w");
    sprintf(filename3, "individuos_AoN_N%i_S%i.txt", N, SIMUS);
    FILE *file3 = fopen(filename3, "w");

    fprintf(file3, "Eps Simu Individuo Vacuna_Status Eficacia Tiempo Infeccion_Status\n");

    if (file1 == NULL || file2 == NULL || file3 == NULL) {
        printf("Error al abrir el archivo para escritura.\n");
        return 1;
    }


    // Arrays para almacenar la información
    double **Time = (double **)malloc(SIMUS * sizeof(double *));
    int **N_Su = (int **)malloc(SIMUS * sizeof(int *));
    int **N_Sv = (int **)malloc(SIMUS * sizeof(int *));
    int **N_Sv_np = (int **)malloc(SIMUS * sizeof(int *));

    for (int i = 0; i < SIMUS; i++) {
        Time[i] = (double *)malloc((PASOS + 1) * sizeof(double));
        N_Su[i] = (int *)malloc((PASOS + 1) * sizeof(int));
        N_Sv[i] = (int *)malloc((PASOS + 1) * sizeof(int));
        N_Sv_np[i] = (int *)malloc((PASOS + 1) * sizeof(int));
    }

    double *Time_ave = (double *)malloc((PASOS + 1) * sizeof(double));
    double *N_Su_ave = (double *)malloc((PASOS + 1) * sizeof(double));
    double *N_Sv_ave = (double *)malloc((PASOS + 1) * sizeof(double));


    //Bucle barrido eficacia promedio
    while (eps <= 1.0 + TOLERANCIA){

      pasos = (int)round((N + N * (1 - eps)));

        // Condiciones iniciales individuos
        for (int i = 0; i < SIMUS; i++) {
            N_Su[i][0] = N_Su0;
            N_Sv[i][0] = N_Sv0;
            N_Sv_np[i][0] = (int) round((N_Sv0 * (1 - eps)));
        }

        /** ALGORITMO DE GILLESPIE **/
        for (int i = 0; i < SIMUS; i++) {
            j=0;

            //Bucle infección de individuos por completo
            while(j != PASOS){

                // Suma de propensidades total
                P = BETA * (double)N_Su[i][j] + BETA * (double)N_Sv_np[i][j];

                // Tiempo al siguiente evento
                double r1 = rand_double();
                double tau = 1 / P * log(1 / r1);


                if (P < TOLERANCIA){/** A este if solo se entra cuando todos los de CONTROL se han infectado y eps = 1**/

                    Time[i][j + 1] = Time[i][j];
                    N_Su[i][j + 1] = N_Su[i][j];
                    N_Sv[i][j + 1] = N_Sv[i][j];
                    N_Sv_np[i][j + 1] = N_Sv_np[i][j];

                    fprintf(file3, "%f %i %i %i %i %f %i\n",eps, i+1, j + 1, 1, 1, Time[i][j + 1], 0);

                    j++;
                }

                else if ( isinf(tau) == 0){

                    // Guardar tiempo
                    Time[i][j + 1] = Time[i][j] + tau;

                    // Seleccionar el próximo evento
                    double Pu = BETA * N_Su[i][j] / P;
                    double r2 = rand_double();

                    //Infectado de la cohorte de control
                    if (r2 < Pu) {

                        N_Su[i][j + 1] = N_Su[i][j] - 1;
                        N_Sv[i][j + 1] = N_Sv[i][j];
                        N_Sv_np[i][j + 1] = N_Sv_np[i][j];

                        fprintf(file3, "%f %i %i %i %i %f %i\n",eps, i+1,  j + 1, 0, 0, Time[i][j + 1], 1);

                    } else { //Infectado de la cohorte vacunada
                        N_Su[i][j + 1] = N_Su[i][j];
                        N_Sv[i][j + 1] = N_Sv[i][j] - 1;
                        N_Sv_np[i][j + 1] = N_Sv_np[i][j] - 1;

                        fprintf(file3, "%f %i %i %i %i %f %i\n",eps, i+1, j + 1, 1, 0, Time[i][j + 1], 1);

                    }

                    j++;
                }
            }
        }/** FIN DE 1 BARRIDO de eps**/

     eps += eps_delta;
    }

    //Liberación memoria
    for (int i = 0; i < SIMUS; i++) {
        free(Time[i]);
        free(N_Su[i]);
        free(N_Sv[i]);
        free(N_Sv_np[i]);
    }

    free(Time);
    free(N_Su);
    free(N_Sv);
    free(N_Sv_np);
    free(Time_ave);
    free(N_Su_ave);
    free(N_Sv_ave);

    fclose(file1);
    fclose(file2);
    fclose(file3);

    printf("Simulaciones guardadas.\n");


    return 0;
}
