#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "Interpolazione.h"


void interp_alpha( double AoA, double alpha_des, double *Cp, double Cs){
    //printf("\n Interpolazione of %lf and %lf gives:\n",*Cp,Cs);
    *Cp=((Cs-*Cp)/(0.2))*(alpha_des-AoA)+*Cp;
    //printf("%lf\n", *Cp);
    // NB Salvo i coefficienti interpolati che userò nella simulazione nella posizione dello stesso vettore indicata da Alpha in modo da usare variabili già in memoria e non allocare ulteriormente altre variabili nella stessa
}

void interpol_coeffs(int*index, int nAoA, double*ALPHA , double alpha_tr, double css[][6],double cxd[][7], double cyd[][6], double czd[][7],
                     double cld[][6], double cmd[][7],double cnd[][6],double cuf[][6],double cum[][6],double crd[][6] ){

    const double pi=4*atan(1);
    double AoA;
    FILE*fin;

    fin=fopen("Validation/Interpolated_coefficients.txt","w"); // File di Validazione

    if (fin==NULL){
        printf("\nError in opening the alpha interp. validation file\n");
        exit(6);
    }

    alpha_tr=alpha_tr*180/pi;
    for (int i = 0; i < nAoA-1; i++) {
        //printf("\nalpha is: %lf, index is: %d\n",ALPHA[i],i);
        if (alpha_tr>=ALPHA[i] && alpha_tr<=ALPHA[i+1]) {
            *index=i; AoA=ALPHA[i];
            break; // If Flag is 0 then alpha_des coincides with ne of the given values in the dba
        }
    }// It returns the index of closest smaller alpha

    //fprintf(fin,"Validation File for Alpha interpolation wrt the alpha trim value\n");
    //fprintf(fin,"SS\t\tX\t\tY\t\tZ\t\tL\t\tM\t\tN\t\tFu\t\tMu\t\tRd\n");

    for(int i=0; i<7; i++){
        if(i<6) {
            interp_alpha(AoA, alpha_tr, &css[*index][i], css[*index+1][i]);
            interp_alpha(AoA, alpha_tr, &cyd[*index][i], cyd[*index+1][i]);
            interp_alpha(AoA, alpha_tr, &cld[*index][i], cld[*index+1][i]);
            interp_alpha(AoA, alpha_tr, &cnd[*index][i], cnd[*index+1][i]);
            interp_alpha(AoA, alpha_tr, &cuf[*index][i], cuf[*index+1][i]);
            interp_alpha(AoA, alpha_tr, &cum[*index][i], cum[*index+1][i]);
            interp_alpha(AoA, alpha_tr, &crd[*index][i], crd[*index+1][i]);
        }
        interp_alpha(AoA, alpha_tr, &cmd[*index][i], cmd[*index+1][i]);
        interp_alpha(AoA, alpha_tr, &cxd[*index][i], cxd[*index+1][i]);
        interp_alpha(AoA, alpha_tr, &czd[*index][i], czd[*index+1][i]);
        if(i<6) {
            fprintf(fin, "%.6lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\t\t%.6lf\n",
                    css[*index][i], cxd[*index][i], cyd[*index][i], czd[*index][i], cld[*index][i], cmd[*index][i],
                    cnd[*index][i], cuf[*index][i], cum[*index][i], crd[*index][i]);
        }else{
            fprintf(fin, "0.000000\t\t%.6lf\t\t0.000000\t\t%.6lf\t\t0.000000\t\t%.6lf\t\t0.000000\t\t0.000000\t\t0.000000\t\t0.000000\n",cxd[*index][i], czd[*index][i], cmd[*index][i]);
        }

    } // Saving the needed coefficient in the i-th position

    fclose(fin);
}

double **MatrixCoeff(int r,int c) {

    double** C = (double**)malloc(r*sizeof(double*));
    double C_tmp[10][10];
    FILE *fp_tmp;

    fp_tmp = fopen("Validation/Interpolated_coefficients.txt","r");
    if (fp_tmp==NULL){
        printf("[!] File temporaneo per i coefficienti aerodinamici non aperto correttamente");
        exit(4);
    }

    // Salvataggio coefficienti interpolati in una matrice 10x10
    for (int row = 0; row <10; row++) {
        C[row]=(double*)malloc(r*sizeof(double*));
        for (int col = 0; col < 10; col++) {
            if (row > 6){
                C[row][col] = 0;
            }else {
                fscanf(fp_tmp, "%lf", &C[row][col]);
            }
            C_tmp[row][col] = C[row][col];
        }
    }

    //Trasposizione della matrice 10x10 (fatto per convenienza)
    for (int row = 0; row < 10; row++){
        for (int col = 0; col < 10; col++){
            C[col][row] = C_tmp[row][col];
        }
    }

    fclose(fp_tmp);

    return C;
};


