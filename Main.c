#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "Functions/Atmosfera/Atmosfera.h"
#include "Functions/ManageFile/ManageFile.h"
#include "Functions/Interpolazione/Interpolazione.h"
#include "Functions/Condizioni_trim/Alpha_trim.h"
#include "Functions/Condizioni_trim/Condizioni_trim.h"
#include "Functions/Calcolo_spinta_RPM/Calcolo_spinta_RPM.h"
#include "Functions/Integrazione/Integrazione.h"
#include "Functions/Routh/Routh.h"
#include "Functions/Modi/Longitudinali.h"

int main() {

    double m = 0.973, S = 0.247;
    double Jx = 1.088315E-2, Jy = 1.1922E-2, Jz = 2.232462E-2, Jxz = 4.5905E-5, b = 0.872, c = 0.296;

    double press0 = 101325, temp0 = 15, rho0 = 1.225, vsuono0 = 340;
    double press_h, temp_h, rho_h, vsuono_h, Pmax_h = 0.16; //[KW]

    double n_min, n_max, v0, a0_batt, a_stall, q_stall; // File engine.txt
    double C_max, eta_s; // File battery.txt
    double datipropeller[11], CSI[125], RD[125], CH_AD[125], BA[125]; // File propeller.txt

    double CSS[110][6], CXD[110][7], CYD[110][6], CZD[110][7], CLD[110][6], CMD[110][7], CND[110][6],
            CUF[110][6], CUM[110][6], CRD[110][6];
    double alpha[110];

    double CI[3] = {0,0,0}; // Condizioni iniziali
    double Vmin;

    char FileEngine[40] = "Drone/engine.ini";
    char FileBattery[40] = "Drone/battery.ini";
    char FilePropeller[40] = "Drone/propeller.ini";

    int flagatm, flagtrim=0, inp=0, alpha_index, nAoA=110;

    double alpha_trim, delta_e_trim;

    //Cambio della working directory
    char projectPath[256];
    strcpy(projectPath, __FILE__);
    WorkingDirectory(projectPath);

    //Inizio del simulatore
    do {
        if (flagtrim!=0) {
            system("read -n 1 -s -p \"Premi invio per riconfigurare il nuovo set di condizioni iniziali...\""); //MACOS
            system("pause"); //WINDOWS
        }
        system("cls");

        loadCI(CI,&inp);

        // Lettura file salvataggio in variabili
        DbaRead(CI[1],alpha,CSS,CXD,CYD,CZD,CLD,CMD,CND,CUF,CUM,CRD);

        // Dati atmosferici e condizioni iniziali
        AtmosphereChoice(&press0, &temp0, &rho0, &vsuono0, &press_h, &temp_h, &rho_h, &vsuono_h, CI,
                         &flagatm); //fa scegliere all'utente le opzioni atmosferiche
        AtmosphereCalc(CI, &Pmax_h, &press0, &temp0, &press_h, &temp_h, &rho_h, &flagatm);

        // Calcolo condizioni di TRIM
        // La funzione per il calcolo di alpha_trim funziona per il volo longitudinale
        flagtrim = calcolo_alpha_trim(alpha,m,S,rho_h,&CI[0],&CI[2],CSS,CZD,CUF,CMD,CUM,&alpha_trim,&delta_e_trim,&Vmin);

    }while(flagtrim!=0);

    // Lettura dei file statici
    ReadAndSavePropeller(FilePropeller, datipropeller, CSI, RD, CH_AD, BA);
    ReadAndSaveStatic(FileEngine, FileBattery, &n_min, &n_max, &v0, &a0_batt, &a_stall, &q_stall, &C_max, &eta_s);

    // Interpolazione dei coefficienti nel dba
    interpol_coeffs(&alpha_index,nAoA,alpha,alpha_trim,CSS,CXD,CYD,CZD,CLD,CMD,CND,CUF,CUM,CRD);

    double trim_cond[11];
    double trim_control[3];

    Condizioni_trim(trim_cond, trim_control, alpha_trim, CI[0], CI[2], CI[1], delta_e_trim,
                    m, S, rho_h, n_min, n_max, datipropeller, BA);

    //Calcolo della stabilita'
    Routh(m,S,rho_h,alpha_trim,CI[0]);

    //Calcolo dei modi longitudinali
    long_modi(m,rho_h,S,CI[0],alpha_trim);

    //Integrazione numerica
    integration(trim_cond, trim_control, alpha_trim, rho_h, S, b, c, m, Jx, Jy, Jz, Jxz, CI[0], CI[1],Vmin,inp);

    printf("--------------------------------------------------------------------------------");
    printf("\nFINE: Simulazione terminata correttamente. Ora verra' automaticamente aperto il "
           "\nfile Plot.m per l'esecuzione della fase finale di post-processing.\n\n");

    // Apertura automatica del file Matlab
    char filename[256] = "/Plot.m";
    char command[256];
    sprintf(command, "open \"%s%s\"", projectPath,filename);  // Genera il comando per aprire il file
    system(command);  // Esegue il comando di sistema per aprire il file

    return 0;
}