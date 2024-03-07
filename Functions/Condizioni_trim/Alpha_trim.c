#include <stdio.h>
#include <math.h> // Per le funzioni trigonometriche
#include <stdlib.h> // Per poter usare il comando exit

int calcolo_alpha_trim(double alpha[], double m, double S, double rho, double* V, double *gamma, double css[][6], double czd[][7],
                        double cuf[][6],double cmd[][7], double cum[][6], double* alpha_trim, double* delta_e_trim, double *Vmin)  {

    // Costanti
    const double g = 9.81; // accelerazione di gravi
    const double pi = atan(1) * 4; // pi greco
    const double alpha_min = -5 * pi / 180, alpha_max = 16.8 * pi / 180; // [rad] incidenza minima e massima
    const double delta_e_min = -20 * pi / 180, delta_e_max = 20 * pi / 180; // [rad] angolo equilibratore minimo e massimo
    const double delta_alpha_dba = 0.2 * pi / 180; // [rad] intervallo tra due alpha consecutivi nel file dba

    // Parametri dell'algoritmo (scelti dall'utente)
    const double delta_alpha = 0.001 * pi / 180;
    const float res_Z = 0.01; // [N] massimo residuo equilibrio forza lungo z (MODIFICABILE)

    // Variabili
    double Cz_tot; // coefficiente forza lungo Z totale
    int flag = 0; // flag per uscire dal ciclo while
    int index; // Indice del vettore alpha corrispondente all'angolo immediatamente inferiore o uguale a alpha_trim
    double diff_alpha; // differenza tra alpha_trim e alpha[index]
    double k; // coefficiente per il calcolo dei coefficienti interpolati
    double Cz_alpha, Cz_delta_e, Cm_ss, Cm_alpha, Cm_delta_e, Cz_ss, Cx_ss; // Coefficienti (presi dal DBA)

    *gamma = *gamma * pi / 180; // Converto gamma da gradi a radianti

    /* NOTA SU INIZIALIZZAZIONE ANGOLI : gli angoli alpha_trim e delta_e_trim sono inizializzati al loro valore minimo meno l'intervallo
    di cui si aggiornano dopo ogni ciclo. Questa scelta   dovuta al fatto che in ogni ciclo while i valori degli angoli vengono
    incrementati prima dello svolgimento dei calcoli. In questo modo, una volta che gli angoli di trim vengono trovati, questi non
    saranno modificati prima di uscire dai cicli while. */


    *alpha_trim = alpha_min - delta_alpha; // Inizializzo alpha trim


    while ( *alpha_trim <= alpha_max && flag == 0) {

        *alpha_trim = *alpha_trim + delta_alpha; // incremento alpha_trim

        // Calcolo le due variabili per l'interpolazione
        index = (int)((*alpha_trim - alpha_min) / delta_alpha_dba);
        diff_alpha = *alpha_trim - (alpha_min + index * delta_alpha_dba);

        // Calcolo i coefficienti interpolati
        k = (diff_alpha * 180 / pi) / (alpha[index + 1] - alpha[index]);
        Cz_alpha = czd[index][0] + (czd[index + 1][0] - czd[index][0]) * k;
        Cz_delta_e = cuf[index][2] + (cuf[index + 1][2] - cuf[index][2]) * k;
        Cm_ss = css[index][4] + (css[index + 1][4] - css[index][4]) * k;
        Cm_alpha = cmd[index][0] + (cmd[index + 1][0] - cmd[index][0]) * k;
        Cm_delta_e = cum[index][2] + (cum[index + 1][2] - cum[index][2]) * k;
        Cz_ss = css[index][2] + (css[index + 1][2] - css[index][2]) * k;
        Cx_ss = css[index][0] + (css[index + 1][0] - css[index][0]) * k;
        // NOTA: la dimensione dei vettori dei coefficienti   maggiore del numero di angoli alpha nel dba.
        // Pertanto, non vi   il rischio che [index+1] sia maggiore della dimensione dei vettori dei coefficienti.

        *delta_e_trim = -(Cm_ss + Cm_alpha*(*alpha_trim))/(Cm_delta_e);

        Cz_tot = Cz_ss + Cz_alpha * (*alpha_trim) + Cz_delta_e * (*delta_e_trim);

        if (fabs(m * g * cos(*alpha_trim + *gamma) + 0.5 * rho * (*V) * (*V) * S * Cz_tot) < res_Z) { // Condizione di equilibrio lungo z
            flag = 1; // => se le condizioni sono soddisfatte si esce dal while
        }
    }

    printf("--------------------------------------------------------------------------------\n");
    printf("CONDIZIONI DI TRIM:\n");

    // Calcolo della velocità minima
    double Cz_ss_max = css[109][2];
    double Cx_ss_max = css[109][0];
    double CLmax = -(-Cx_ss_max*sin(alpha_max)+Cz_ss_max*cos(alpha_max));
    *Vmin = sqrt(2*m*g/(S*rho*CLmax));

    if (*V <= *Vmin){
        printf("\nERRORE[!] La velocita' è minore o uguale della velocità di stallo pari a %lf. Reinserire le nuove condizioni iniziali.\n", Vmin);
        return 1;
    }

    if (flag == 0) {
        printf("\nERRORE[!] Condizioni di trim non trovate: nessuna deflessione dell'elevone soddisfa l'equilibrio longitudinale. \n");
        printf("\nInserire nuovamente quota, velocita' e angolo di rampa.\n\n");
        return 1;
    }
    else {
        if (*delta_e_trim > delta_e_max || *delta_e_trim < delta_e_min) {
            printf("\nERRORE[!] Deflessione dell'elevone non ragguingibile.\nInserire nuovamente quota, velocita' e angolo di rampa.\n\n");
            return 1;
        }
        else {
            // printf("\nIl residuo dell'equilibrio lungo z e' %.3f N\n", fabs(m * g * cos(*alpha_trim + *gamma) + 0.5 * rho * *V * *V * S * Cz_tot));
            // printf("Il residuo dell'equilibrio attorno a y e' %.6f\n\n", fabs(Cm_ss + Cm_alpha * (*alpha_trim) + Cm_delta_e * (*delta_e_trim)));
             printf("\nalpha_trim:\t%.3f [deg]\t(%.3f [rad]) \n", *alpha_trim * 180 / pi, *alpha_trim);
             printf("delta_e_trim:\t%.3f [deg]\t(%.3f [rad]) \n\n", *delta_e_trim * 180 / pi, *delta_e_trim);
             return 0;
        }
    }

}