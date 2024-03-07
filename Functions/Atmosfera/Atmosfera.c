#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "Atmosfera.h"

// il codice mostra i dati atmosferici di default a video e permette di modificarli.

int AtmosphereChoice (double *press0,double *temp0,double *rho0,double *vsuono0,
                      double *press_h,double *temp_h,double *rho_h,double *vsuono_h, double *CI, int *flagatm)
{
    int input;

    printf("--------------------------------------------------------------------------------");
    printf("\nATMOSFERA ISA:\n");
    printf("\nLa simulazione fa riferimento al modello atmosferico ISA, avente i segeunti valori per quota h=0m (sea level):\n\n");
    printf("\tPressione: \t\tP = %g Pa\n",*press0);
    printf("\tTemperatura: \t\tT = %g C\n",*temp0);
    printf("\tDensita': \t\trho = %g kg/m^3\n",*rho0);
    printf("\tVelocita' del Suono: \ta = %g m/s\n",*vsuono0);

    printf("\nSe non si desidera procedere con i suddetti parametri e' possibile modificarli, reinserendoli manualmente o scegliendo una quota differente.\n");
    do {
        printf("\nPremere:\n");
        printf("\t[1] se si desidera procedere con i valori iniziali precedentemente indicati\n"); // mantiene i valori iniziali di default
        printf("\t[2] se si desidera modificare i valori iniziali a h=0\n"); //modifica i valori iniziali di default
        printf("\t[3] se di desidera immettere i valori manualmente a una quota specifica (richiesta in seguito)\n");// immette i valori manualmente a una data quota senza richiederla poi

        scanf("%d", &input);
        if(input!=1 && input!=2 && input!=3)
        {
            scanf("%*[^\n]"); //svuota il buffer di scanf se immetto un carattere invece che un numero
            printf("[!]WARNING immettere un numero da 1 a 3\n");
        }
        switch(input)
        {
            case 1: // mantiene i valori iniziali di default
                *flagatm = 1;
                break;
            case 2: //modifica i valori iniziali di default
                printf("Inserire un valore di pressione (h=0) in Pa:");
                scanf("%lf", press0);
                if(*press0<0) {
                    do { //controllo che il numero immesso sia maggiore di zero
                        printf("[!]WARNING immettere un numero positivo\n");
                        printf("Inserire un valore di pressione (h=0) in Pa:");
                        scanf("%lf", press0);
                    } while (*press0<0);
                }

                printf("Inserire un valore di temperatura (h=0) in C:");
                scanf("%lf", temp0);

                printf("Inserire un valore di densita' (h=0) in kg/m^3:");
                scanf("%lf", rho0);
                if(*rho0<0) {
                    do { //controllo che il numero immesso sia maggiore di zero
                        printf("[!]WARNING immettere un numero positivo\n");
                        printf("Inserire un valore di densita' (h=0) in kg/m^3:");
                        scanf("%lf", rho0);
                    } while (*rho0 < 0);
                }

                printf("Inserire un valore di velocita' (h=0) del suono in m/s:");
                scanf("%lf", vsuono0);
                if(*vsuono0<0) {
                    do { //controllo che il numero immesso sia maggiore di zero o effettivamente un numero
                        printf("[!]WARNING immettere un numero positivo\n");
                        printf("Inserire un valore di velocita' (h=0) del suono in m/s:");
                        scanf("%lf", vsuono0);
                    } while (*vsuono0 < 0);
                }

                *flagatm = 1;
                break;
            case 3: // immette i valori manualmente a una data quota senza richiederla poi
                printf("Inserire un valore di pressione in Pa:");
                scanf("%lf", press_h);
                if(*press_h<0) {
                    do { //controllo che il numero immesso sia maggiore di zero
                        printf("[!]WARNING immettere un numero positivo\n");
                        printf("Inserire un valore di pressione in Pa:");
                        scanf("%lf", press_h);
                    } while (*press_h < 0);
                }

                printf("Inserire un valore di temperatura in C:");
                scanf("%lf", temp_h);

                printf("Inserire un valore di densita' in kg/m^3:");
                scanf("%lf", rho_h);
                if(*rho_h<0) {
                    do { //controllo che il numero immesso sia maggiore di zero
                        printf("[!]WARNING immettere un numero positivo\n");
                        printf("Inserire un valore di densita' in kg/m^3:");
                        scanf("%lf", rho_h);
                    } while (*rho_h < 0);
                }

                printf("Inserire un valore di velocita' del suono in m/s:");
                scanf("%lf", vsuono_h);
                if(*vsuono_h<0) {
                    do { //controllo che il numero immesso sia maggiore di zero o effettivamente un numero
                        printf("[!]WARNING immettere un numero positivo\n");
                        printf("Inserire un valore di velocita' del suono in m/s:");
                        scanf("%lf", vsuono_h);
                    } while (*vsuono_h < 0);

                }
                /*if(CI[1]<0)
                {
                    printf("[!]ERRORE: I valori inseriti corrispondono ad una quota minore di zero\nSIMULAZIONE TERMINATA!");
                    system("PAUSE");
                    exit(14);
                }
                else if(CI[1]>4116)
                {
                    printf("[!]ERRORE: I valori inseriti corrispondono ad una quota maggiore di quella di tangenza\nSIMULAZIONE TERMINATA!\n");
                    system("PAUSE");
                    exit(14);
                }*/
                *flagatm = 2;
                break;
        }
    }
    while(input!=1 && input!=2 && input!=3);

    return 0;
}

int AtmosphereCalc (double *CI, double *Pmax_h,double *press0,double *temp0,
                    double *press_h,double *temp_h,double *rho_h, int *flagatm)
{

    double R=287.05;
    switch(*flagatm)
    {
        case 1: // calcola le condizioni atmosferiche e di pot per i valori a SL di default o inseriti dall'utente
            *temp0= *temp0 + 273.15;
            *temp_h=*temp0 - 0.0065*CI[1];
            *press_h=*press0*(pow((*temp_h/(*temp0)),5.2561));
            *rho_h=*press_h/(R*(*temp_h));

            printf("\nI dati atmosferici e di potenza per la quota di %f m sono:\n\n",CI[1]);
            printf("Temperatura: \t\t%f [K]\n",*temp_h);
            printf("Pressione: \t\t%f [Pa]\n",*press_h);
            printf("Desita': \t\t%f [kg/m^3]\n",*rho_h);
            printf("Potenza: \t\t%f [kW]\n",*Pmax_h);
            break;
        case 2:// usa i valori scelti dall'utente
            *temp_h=*temp_h+273.15;

            printf("\nI dati atmosferici e di potenza per la quota di %f m sono:\n\n",CI[1]);
            printf("Temperatura: \t\t%f [K]\n",*temp_h);
            printf("Pressione: \t\t%f [Pa]\n",*press_h);
            printf("Densita': \t\t%f [kg/m^3]\n",*rho_h);
            printf("Potenza: \t\t%f [kW]\n",*Pmax_h);
            break;
    }
    return 0;
}

void loadCI(double *CI, int* inp){

    double V = 0; int h = 0; double gamma = 0;
    double V_max = 25;  // Limite di velocità [m/s]
    double h_t = 2000; // Quota di tangenza [m]

    printf("\nSIMULATORE DI VOLO: \n");

    while(*inp!=1 && *inp!=2) {
        printf("\nSelezionare una tra le seguenti scelte:\n"
               "[1] Simulazione del mantenimento del trim.\n"
               "[2] Simulazione di un percorso predefinito.\n");
        scanf("%d", inp);
        printf("--------------------------------------------------------------------------------");
    }

    if (*inp == 2) {
        V = 13.5;
        h = 100;
        gamma = 0;
        printf("\nCONDIZIONI INIZIALI: \n\n");
        printf("Le condizioni iniziali sono state automaticamente settate a:\n"
               "V = 13.5 [m/s]\n"
               "h = 100 [m]\n"
               "gamma = 0 [deg]\n");

    }
    else if (*inp == 1){
        printf("\nCONDIZIONI INIZIALI: \n\n");
        printf("Inserisci la quota alla quale stai volando in [m]: ");
        scanf("%d", &h);

        if (h<0) {
            //printf("\n[!] Inserisci un valore positivo della quota in [m]: ");
            //scanf("%d", &h);
            printf("Inserito un valore di quota negativo.\n");
            exit(0);
        }

        if (h == h_t) {
            printf("\n[!] La quota inserita e' quella di tangenza.\n");
        }
        if (h > h_t) {
            printf("\n[!] La quota inserita e' superiore rispetto a quella di tangenza.\n");
            exit(1);
        }

        printf("Inserisci il valore della velocita' in [m/s]: ");
        scanf("%lf", &V);

        do {
            if (V > V_max) {
                printf("\nLa velocità inserita e' al di sopra del limite di velocità consentito. Inserire una velocità inferiore a %lf [m/s]: ", V_max);
                scanf("%lf", &V);
            }
        } while (V > V_max);

        printf("Inserire l'angolo di rampa in [deg]: ");
        scanf("%lf", &gamma);

    }else{
        printf("\nIl numero inserito non e' corretto, riprovare.\n\n");
    }



    *CI = V; CI++;
    *CI = h; CI++;
    *CI = gamma;
}