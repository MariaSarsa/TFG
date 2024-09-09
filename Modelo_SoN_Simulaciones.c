/*
 * Programa de Simulación de Infección en una Población con Protección Heterogénea Intermedia
 *
 * Descripción:
 * Este programa realiza simulaciones de la dinámica de infección en una población de individuos,
 * considerando una cohorte de control y una cohorte vacunada SoN. Utiliza el algoritmo de Gillespie
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
 * - individuos_SoN_N{N}_S{SIMUS}.txt: Archivo para guardar detalles de los individuos y su estado
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
#include <float.h>

#define N 1000
#define BETA 0.5
#define PASOS (2 * N)
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

    double eps = 0.0;
    double eps_delta = 0.05;
    double var = 0.0;
    double var_delta = 0.01;
    double p;
    double mu;

    // Arrays para almacenar la información
    double **Time = (double **)malloc(SIMUS * sizeof(double *));
    int **N_Su = (int **)malloc(SIMUS * sizeof(int *));
    int **N_Sv_p = (int **)malloc(SIMUS * sizeof(int *));
    int **N_Sv_np = (int **)malloc(SIMUS * sizeof(int *));

    if (Time == NULL || N_Su == NULL || N_Sv_p == NULL || N_Sv_np == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    FILE *file1 = fopen("simulaciones.txt", "w");
    FILE *file2 = fopen("promedios.txt", "w");
    sprintf(filename3, "individuos_SoN_N%i_S%i_eps02_var02.txt", N, SIMUS);
    FILE *file3 = fopen(filename3, "w");

    fprintf(file3, "Eps Var Simu Individuo Vacuna_Status Eficacia Tiempo Infeccion_Status\n");

    if (file1 == NULL || file2 == NULL || file3 == NULL) {
        printf("Error al abrir el archivo para escritura.\n");
        return 1;
    }

    for (int i = 0; i < SIMUS; i++) {
        Time[i] = (double *)malloc((PASOS + 1) * sizeof(double));
        N_Su[i] = (int *)malloc((PASOS + 1) * sizeof(int));
        N_Sv_p[i] = (int *)malloc((PASOS + 1) * sizeof(int));
        N_Sv_np[i] = (int *)malloc((PASOS + 1) * sizeof(int));

        if (Time[i] == NULL || N_Su[i] == NULL || N_Sv_p[i] == NULL || N_Sv_np[i] == NULL) {
            printf("Memory allocation failed.\n");
            return 1;
        }
    }

    double *Time_ave = (double *)malloc((PASOS + 1) * sizeof(double));
    double *N_Su_ave = (double *)malloc((PASOS + 1) * sizeof(double));
    double *N_Sv_ave = (double *)malloc((PASOS + 1) * sizeof(double));

    if (Time_ave == NULL || N_Su_ave == NULL || N_Sv_ave == NULL) {
        printf("Memory allocation failed.\n");
        return 1;
    }

    //Bucle barrido eficacia promedio
    while (eps <= (double)(1)+ TOLERANCIA){

        var = 0;

        //Bucle barrido varianza
        while (var  <=  (eps * (1 - eps)+ TOLERANCIA)){

        //Caso especial donde tanto eps como var son 0
        if (eps == 0 && var == 0){

            p = 0;
            mu = 0;

        }else{

            p = eps * eps /(var + eps * eps);
            mu = eps / p;

        }

        // Condiciones iniciales
        for (int i = 0; i < SIMUS; i++) {
            //printf("hola\n");
            N_Su[i][0] = N_Su0;
            N_Sv_p[i][0] = (int) round((N_Sv0 * p));
            N_Sv_np[i][0] = N_Su0 - (int) round((N_Sv0 * p));

        }

        //Bucle de simulación de ensayos
        for (int i = 0; i < SIMUS; i++) {

            j =0;

            while(j != PASOS ){

                // Suma de propensidades total
                double P;
                P = BETA * N_Su[i][j] + BETA * (1 - mu) * N_Sv_p[i][j] + BETA * N_Sv_np[i][j];

                // Tiempo al siguiente evento
                double r1 = rand_double();
                double tau = 1 / P * log(1 / r1);

                if (P < TOLERANCIA){

                            Time[i][j + 1] = Time[i][j];
                            N_Su[i][j + 1] = N_Su[i][j];
                            N_Sv_p[i][j + 1] = N_Sv_p[i][j];
                            N_Sv_np[i][j + 1] = N_Sv_np[i][j];

                            fprintf(file3, "%f %f %i %i %i %f %f %i\n",eps, var, i+1, j +1, 1, 1.0, Time[i][j + 1], 0);

                            j++;

                }else if ( isinf(tau) == 0) {

                    // Guardar tiempo
                    Time[i][j + 1] = Time[i][j] + tau;

                    // Seleccionar el próximo evento
                    double Pu = BETA * N_Su[i][j] / P;
                    double Pv_p = BETA * (1 - mu) * N_Sv_p[i][j] / P;
                    double r2 = rand_double();

                    //Se infecta CONTROL
                    if (r2 < Pu) {

                        N_Su[i][j + 1] = N_Su[i][j] - 1;
                        N_Sv_p[i][j + 1] = N_Sv_p[i][j];
                        N_Sv_np[i][j + 1] = N_Sv_np[i][j];

                        fprintf(file3, "%f %f %i %i %i %f %f %i\n",eps, var, i+1, j +1, 0, 0.0, Time[i][j + 1], 1);

                    } else if(r2 < (Pu + Pv_p) || (N_Sv_p[i][j] != 0 && N_Sv_np[i][j]==0)){//Se infecta alguien de los SOMETHING

                        N_Su[i][j + 1] = N_Su[i][j];
                        N_Sv_p[i][j + 1] = N_Sv_p[i][j] - 1;
                        N_Sv_np[i][j + 1] = N_Sv_np[i][j];

                        fprintf(file3, "%f %f %i %i %i %f %f %i\n",eps, var, i+1, j+1, 1, mu, Time[i][j + 1], 1);

                      } else{ //Se infecta alguien de los NOTHING

                        N_Su[i][j + 1] = N_Su[i][j];
                        N_Sv_p[i][j + 1] = N_Sv_p[i][j];
                        N_Sv_np[i][j + 1] = N_Sv_np[i][j] - 1;

                        fprintf(file3, "%f %f %i %i %i %f %f %i\n",eps, var, i+1, j+1, 1, 0.0, Time[i][j + 1], 1);

                        }
                        j++;
                    }//Fin else if

                }//Fin while
            }//Fin simus---> Fin 1 valor de var

            var += var_delta;

        }/** FIN BARRIDO var**/

        eps += eps_delta;

    }/** FIN BARRIDO eps**/

        fclose(file1);
        fclose(file2);
        fclose(file3);

        //Liberación de memoria
      for (int i = 0; i < SIMUS; i++) {
        free(Time[i]);
        free(N_Su[i]);
        free(N_Sv_p[i]);
        free(N_Sv_np[i]);
    }

        free(Time);
        free(N_Su);
        free(N_Sv_p);
        free(N_Sv_np);
        free(Time_ave);
        free(N_Su_ave);
        free(N_Sv_ave);



    printf("Simulaciones guardadas.\n");

    return 0;



}
































