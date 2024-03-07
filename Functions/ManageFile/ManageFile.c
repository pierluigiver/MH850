#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct coeff{
    double CSS[110][6], CXD[110][7], CYD[110][6], CZD[110][7], CLD[110][6], CMD[110][7], CND[110][6], CUF[110][6], CUM[110][6], CRD[110][6];
}COEFF;

void WorkingDirectory(const char* projectPath){
    printf("\nSETTAGGIO DELLA WORKING DIRECTORY:\n");

    char WorkingDirectory[100];

    char* lastSeparator = strrchr(projectPath, '/');
    if (lastSeparator == NULL) {
        lastSeparator = strrchr(projectPath, '\\');
    }

    if (lastSeparator != NULL) {
        *lastSeparator = '\0';
        printf("\nPercorso della cartella del progetto: %s\n", projectPath);
        sprintf(projectPath,"%s/cmake-build-debug",projectPath);
        chdir(projectPath);
        printf("--------------------------------------------------------------------------------");
    } else {
        printf("\nImpossibile ottenere il percorso della cartella del progetto.\n");
        printf("Aprire la cartella 'Simulatore' e trascinare qui la cartella 'cmake-build-debug': ");
        scanf("%s",WorkingDirectory);
        chdir(WorkingDirectory);
        printf("\n--------------------------------------------------------------------------------");
    }

}

void ReadAndSaveStatic (char *file_engine, char *file_battery, double *n_min, double *n_max, double *v0,
                        double *a0, double *a_stall, double *q_stall, double *C_max, double *eta_s){

    char line[100];
    double vec_var[10]; // Vettore dove vado a salvare i valori delle variabili
    int i=0 , j=0;

    FILE *fp[2];//Creo un vettore di due puntatori fp[0] e fp[1]
    fp[0] = fopen(file_engine, "r");
    fp[1] = fopen(file_battery, "r");
    
    FILE *fp_val[2];
    fp_val[0] = fopen("Validation/engine.ini","w");
    fp_val[1] = fopen("Validation/battery.ini","w");
   
    if (fp[0] == NULL){
        printf("[!]ERRORE nella lettura del file '%s'\n", file_engine);
        exit(1);
    }

    if (fp[1] == NULL){
        printf("[!]ERRORE nella lettura del file '%s'\n", file_battery);
        exit(2);
    }

    if (fp_val[0] == NULL) {
        printf("[!]ERRORE nella lettura del file di validazione '%s'\n", file_engine);
        exit(3);
    }

    if (fp_val[1] == NULL) {
        printf("[!]ERRORE nella lettura del file di validazione '%s'\n", file_battery);
        exit(4);
    }

    // Lettura dei file e salvataggio dei dati
    for (j=0; j<2; j++){ //Ciclo for aggiorna il puntatore

        while (fgets(line,sizeof(line), fp[j]) != NULL){

            if (line[0] == '*') {
                continue;  //Se trovi un '*' vai alla riga successiva
            }
            else{
                sscanf(line, "%lf", &vec_var[i]);  // Prendo il il numero e lo salvo nel vettore
                i++;
            }
        } //Ciclo while salva le varibili in un vettore
        fclose(fp[j]);
    }

    // Assegnazioni dei valori alle variabili e stampa di questi ultimi a video
    *n_min = vec_var[0];
    *n_max = vec_var[1];
    *v0 = vec_var[2];
    *a0 = vec_var[3];
    *a_stall = vec_var[4];
    *q_stall = vec_var[5];
    *C_max = vec_var[6];
    *eta_s = vec_var[7];

    fprintf(fp_val[0],"I dati letti nel file '%s' sono i seguenti : \n\n", file_engine);
    fprintf(fp_val[0],"NUMERO DI GIRI MINIMO DEL MOTORE  [rpm] = %.2f\n", *n_min);
    fprintf(fp_val[0],"NUMERO DI GIRI MASSIMO DEL MOTORE [rpm] SENZA CARICO = %.2f\n", *n_max);
    fprintf(fp_val[0],"NOMINAL VOLTAGE = %.2f\n", *v0);
    fprintf(fp_val[0],"NO LOAD CURRENT = %.2f\n", *a0);
    fprintf(fp_val[0],"STALL CURRENT = %.2f\n", *a_stall);
    fprintf(fp_val[0],"STALL TORQUE = %.2f\n\n\n", *q_stall);
    
    fclose(fp_val[0]);

    fprintf(fp_val[1],"I dati letti nel file '%s' sono i seguenti: \n\n", file_battery);
    fprintf(fp_val[1],"CAPACITA' NOMINALE \t[mAh] = %.2f\n", *C_max);
    fprintf(fp_val[1],"RENDIMENTO DI SCARICA = %.2f\n\n", *eta_s);
    
    fclose(fp_val[1]);

}

double dbainterp(double Cp, double Cs, double *h){

    int quote[4] = {100,1000,2000,3000};
    double C, k;
    int i, hp, hs;

    for (i=0;i<(sizeof(quote)/sizeof(quote[0]))-1;i++) {
        if (*h>=quote[i] && *h<=quote[i+1]){
            hp = quote[i];
            hs = quote[i+1];
        }
    }

    k = (Cs - Cp) / (hs-hp);
    C = Cp + k * (*h - hp);

    return C;
}

/*
void ReadAndSaveDba (double *h, double *alpha, double *CX, double *CY, double *CZ, double *Cl, double *Cm, double *Cn,
              double *CXA, double *CXAP, double *CXU, double *CXQ, double *CXB, double *CXP, double *CXR,
              double *CYB, double *CYBP, double *CYP, double *CYR, double *CYA, double *CYQ,
              double *CZALPHA, double *CZAP, double *CZU, double *CZQ, double *CZB, double *CZP, double *CZR,
              double *ClB, double *ClBP, double *ClP, double *ClR, double *ClA, double *ClQ,
              double *CmA, double *CmAP, double *CmU, double *CmQ, double *CmB, double *CmP, double *CmR,
              double *CnB, double *CnBP, double *CnP, double *CnR, double *CnA, double *CnQ,
              double *CXde, double *CXdle, double *CZde, double *CZdle, double *CYda, double *CYdr,
              double *Clda, double *Cldr, double *Cmde, double *Cmdle, double *Cnda, double *Cndr,
              double *CXom, double *CYom, double *CZom, double *Clom, double *Cmom, double *Cnom){

    struct coefficients{
        double CX[125], CY[125], CZ[125], Cl[125], Cm[125], Cn[125];
        double CXA[125], CXAP[125], CXU[125], CXQ[125], CXB[125], CXP[125], CXR[125];
        double CYB[125], CYBP[125], CYP[125], CYR[125], CYA[125], CYQ[125];
        double CZALPHA[125], CZAP[125], CZU[125], CZQ[125], CZB[125], CZP[125], CZR[125];
        double ClB[125], ClBP[125], ClP[125], ClR[125], ClA[125], ClQ[125];
        double CmA[125], CmAP[125], CmU[125], CmQ[125], CmB[125], CmP[125], CmR[125];
        double CnB[125], CnBP[125], CnP[125], CnR[125], CnA[125], CnQ[125];
        double CXde[125], CXdle[125], CZde[125], CZdle[125], CYda[125], CYdr[125];
        double Clda[125], Cldr[125], Cmde[125], Cmdle[125], Cnda[125], Cndr[125];
        double CXom[125], CYom[125], CZom[125], Clom[125], Cmom[125], Cnom[125];
    };

    struct coefficients p;
    struct coefficients s;

    int line_number = 1;
    int vec_in[10],vec_fin[10],i=0,j=0;
    double a[130];
    int tmp = 0;
    char line[200];
    char file[20];
    FILE *fp1 = NULL;
    FILE *fp2 = NULL;
    FILE *fp_val;


    int quote[4] = {100,1000,2000,3000};

    for (i = 0; i<sizeof(quote)/sizeof(quote[0]); i++){

        if ((quote[i]<=*h)&&(quote[i+1]>=*h)){

            sprintf(file,"Drone/dba_%d.txt",quote[i]);
            fp1 = fopen(file,"r");
            sprintf(file,"Drone/dba_%d.txt",quote[i+1]);
            fp2 = fopen(file,"r");

            if (fp1 == NULL){
                printf("[!]ERRORE nella lettura del file txt.\n");
                exit(3);
            }
            if (fp2 == NULL){
                printf("[!]ERRORE nella lettura del file txt.\n");
                exit(3);
            }
            break;
        }
    }


    fp_val = fopen("Validation/dba.ini","w");

    if (fp_val == NULL){
        printf("\n[!]Problemi con il file dba di validazione.\n");
        exit(4);
    }

    // Calcolo delle posizioni nel file in cui si trovano i dati che devono essere salvati nei vettori
    while (fgets(line, sizeof(line), fp1)) {
        if (strstr(line, "ALPHA") != NULL) {
            vec_in[i] = line_number + 1;
            tmp = 1;
        }
        if (tmp == 1) {
            //fprintf(fp_tmp, "%s", line); // Salvataggio dati su un nuovo file di testo
            if (line[1] == '*') {
                vec_fin[i] = line_number - 1;
                i++;
                tmp = 0;
            }
        }
        line_number++;
    }
    vec_fin[i] = line_number-1;

    // Salvataggio variabili nelle apposite memorie allocate, salvataggio in file temporaneo e stampa su terminale
    fprintf(fp_val,"Valori presenti nel file dba:\n\n");

    // Salvataggio STEADY STATE COEFFICIENT [0]
    fprintf(fp_val, "STEADY STATE COEFFICIENT: \n");
    fseek(fp1, 0, SEEK_SET);
    fseek(fp2, 0, SEEK_SET);
    for (i = 1; i < vec_in[0]; i++) {
        fgets(line, sizeof(line), fp1);
        fgets(line, sizeof(line), fp2);
    }
    for (i=vec_in[0]; i<vec_fin[0]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%lf %lf %lf %lf %lf %lf %lf", alpha,&p.CX[i],&p.CY[i],&p.CZ[i],&p.Cl[i],&p.Cm[i],&p.Cn[i]);
        fgets(line, sizeof(line), fp2);
        sscanf(line,"%lf %lf %lf %lf %lf %lf %lf", alpha,&s.CX[i],&s.CY[i],&s.CZ[i],&s.Cl[i],&s.Cm[i],&s.Cn[i]);
        a[j] = *alpha;
        *CX = dbainterp(p.CX[i],s.CX[i],h);
        *CY = dbainterp(p.CY[i],s.CY[i],h);
        *CZ = dbainterp(p.CZ[i],s.CZ[i],h);
        *Cl = dbainterp(p.Cl[i],s.Cl[i],h);
        *Cm = dbainterp(p.Cm[i],s.Cm[i],h);
        *Cn = dbainterp(p.Cn[i],s.Cn[i],h);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n", *alpha,*CX,*CY,*CZ,*Cl,*Cm,*Cn); // Salvataggio dati su un nuovo file di testo
        alpha++;CX++;CY++;CZ++;Cl++;Cm++;Cn++;j++;
    }

    // Salvataggio X FORCE DERIVATIVES [1]
    fprintf(fp_val, "\nX FORCE DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[1]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[1]; i<vec_fin[1]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf %lf",CXA,CXAP,CXU,CXQ,CXB,CXP,CXR);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*CXA,*CXAP,*CXU,*CXQ,*CXB,*CXP,*CXR); // Salvataggio dati su un nuovo file di testo
        CXA++;CXAP++;CXU++;CXQ++;CXB++;CXP++;CXR++;j++;
    }

    // Salvataggio Y FORCE DERIVATIVES [2]
    fprintf(fp_val, "\nY FORCE DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[2]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[2]; i<vec_fin[2]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf",CYB,CYBP,CYP,CYR,CYA,CYQ);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*CYB,*CYBP,*CYP,*CYR,*CYA,*CYQ); // Salvataggio dati su un nuovo file di testo
        CYB++;CYBP++;CYP++;CYR++;CYA++;CYQ++;j++;
    }

    // Salvataggio Z FORCE DERIVATIVES [3]
    fprintf(fp_val, "\nZ FORCE DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[3]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[3]; i<vec_fin[3]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf %lf",CZALPHA,CZAP,CZU,CZQ,CZB,CZP,CZR);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*CZALPHA,*CZAP,*CZU,*CZQ,*CZB,*CZP,*CZR); // Salvataggio dati su un nuovo file di testo
        CZALPHA++;CZAP++;CZU++;CZQ++;CZB++;CZP++;CZR++;j++;
    }

    // Salvataggio ROLLING MOMENT DERIVATIVES [4]
    fprintf(fp_val, "\nROLLING MOMENT DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[4]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[4]; i<vec_fin[4]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf",ClB,ClBP,ClP,ClR,ClA,ClQ);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*ClB,*ClBP,*ClP,*ClR,*ClA,*ClQ); // Salvataggio dati su un nuovo file di testo
        ClB++;ClBP++;ClP++;ClR++;ClA++;ClQ++;j++;
    }

    // Salvataggio PITCHING MOMENT DERIVATIVES [5]
    fprintf(fp_val, "\nPITCHING MOMENT DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[5]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[5]; i<vec_fin[5]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf %lf",CmA,CmAP,CmU,CmQ,CmB,CmP,CmR);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*CmA,*CmAP,*CmU,*CmQ,*CmB,*CmP,*CmR); // Salvataggio dati su un nuovo file di testo
        CmA++;CmAP++;CmU++;CmQ++;CmB++;CmP++;CmR++;j++;
    }

    // Salvataggio YAWING MOMENT DERIVATIVES [6]
    fprintf(fp_val, "\nYAWING MOMENT DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[6]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[6]; i<vec_fin[6]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf",CnB,CnBP,CnP,CnR,CnA,CnQ);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*CnB,*CnBP,*CnP,*CnR,*CnA,*CnQ); // Salvataggio dati su un nuovo file di testo
        CnB++;CnBP++;CnP++;CnR++;CnA++;CnQ++;j++;
    }

    // Salvataggio CONTROL FORCE DERVATIVES [7]
    fprintf(fp_val, "\nCONTROL FORCE DERVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[7]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[7]; i<vec_fin[7]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf",CXde,CXdle,CZde,CZdle,CYda,CYdr);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*CXde,*CXdle,*CZde,*CZdle,*CYda,*CYdr); // Salvataggio dati su un nuovo file di testo
        CXde++;CXdle++;CZde++;CZdle++;CYda++;CYdr++;j++;
    }

    // Salvataggio CONTROL MOMENT DERIVATIVES [8]
    fprintf(fp_val, "\nCONTROL MOMENT DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[8]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[8]; i<vec_fin[8]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf",Clda,Cldr,Cmde,Cmdle,Cnda,Cndr);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*Clda,*Cldr,*Cmde,*Cmdle,*Cnda,*Cndr); // Salvataggio dati su un nuovo file di testo
        Clda++;Cldr++;Cmde++;Cmdle++;Cnda++;Cndr++;j++;
    }

    // Salvataggio ROTARY DERIVATIVES [9]
    fprintf(fp_val, "\nROTARY DERIVATIVES: \n");
    fseek(fp1, 0, SEEK_SET);
    j=0;
    for (i = 1; i < vec_in[9]; i++) {
        fgets(line, sizeof(line), fp1);
    }
    for (i=vec_in[9]; i<vec_fin[9]+1; i++) {
        fgets(line, sizeof(line), fp1);
        sscanf(line,"%*lf %lf %lf %lf %lf %lf %lf",CXom,CYom,CZom,Clom,Cmom,Cnom);
        fprintf(fp_val,"%f\t %f\t %f\t %f\t %f\t %f\t %f\n",a[j],*CXom,*CYom,*CZom,*Clom,*Cmom,*Cnom); // Salvataggio dati su un nuovo file di testo
        CXom++;CYom++;CZom++;Clom++;Cmom++;Cnom++;j++;
    }

    fclose(fp1);
    fclose(fp_val);

    printf("\n");
}
*/

double C_h_interp(double Cp, double Cs, double hp, double hs, double h){
    double C;
    C=Cp+(Cs-Cp)/(hs-hp)*(h-hp);
    return C;
};

void DbaRead(double h, double*ALPHA, double css[][6],double cxd[][7], double cyd[][6], double czd[][7],
              double cld[][6], double cmd[][7],double cnd[][6],double cuf[][6],double cum[][6],double crd[][6] ){


    char s1[255], s2[255], FN[20];
    int quota[4]={100, 1000, 2000, 3000}; // Possible Dbas
    int qp,qs;
    int i=-1, flag=0;
    COEFF h1, h2; // One for height 1 and the other for height 2 to interpolate
    FILE *fp1, *fp2; //File pointer

//File selection based on height

while( i<3 &&  flag==0) {
        i++;
        if (h >= quota[i] && h <= quota[i + 1]) {
            qp=quota[i]; qs=quota[i+1];
            sprintf(FN,"Drone/dba_%d.txt", qp);
            fp1 = fopen(FN,"r");
            sprintf(FN,"Drone/dba_%d.txt", qs);
            fp2 = fopen(FN,"r");
            flag = 1;
        }
    }

//Check correctness in opening
    if (fp1==NULL){
        printf("\n[!]ERRORE nell'apertura del file dba realtivo alla quota %d.\n\n", quota[i]);
        exit(3);
    }
    if(fp2==NULL){
        printf("\n[!]ERRORE nell'apertura del file dba realtivo alla quota %d.\n\n", quota[i+1]);
        exit(4);
    }

//READING COEFFICIENTS THAT ARE FUNCTIONS OF ALPHA
    int j=0; //Table number
    while(fgets(s1,sizeof(s1),fp1) && fgets(s2,sizeof(s2),fp2)){

        if (strstr(s1, "ALPHA")!=NULL && strstr(s2, "ALPHA")!=NULL){// This function searches for the substring ALPHA in the whole string s

            for ( i=0; i<110; i++) { //110 or 106=number of AoA
                fgets(s1, sizeof(s1), fp1);
                fgets(s2,sizeof(s2),fp2);
                //Copying variables in a file
                if (s1[0]!='*' && s2[0]!='*') {
                    //fputs(s, fp_tmp);
                    switch (j) {
                        case 0: // Steady state coeffs
                            sscanf(s1, "%lf %lf %lf %lf %lf %lf %lf",&ALPHA[i], &h1.CSS[i][0], &h1.CSS[i][1], &h1.CSS[i][2], &h1.CSS[i][3], &h1.CSS[i][4], &h1.CSS[i][5]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf",  &h2.CSS[i][0], &h2.CSS[i][1], &h2.CSS[i][2], &h2.CSS[i][3], &h2.CSS[i][4], &h2.CSS[i][5]);
                            break;
                        case 1: //X force derivatives
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf %lf", &h1.CXD[i][0], &h1.CXD[i][1], &h1.CXD[i][2], &h1.CXD[i][3], &h1.CXD[i][4], &h1.CXD[i][5], &h1.CXD[i][6]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf %lf", &h2.CXD[i][0], &h2.CXD[i][1], &h2.CXD[i][2], &h2.CXD[i][3], &h2.CXD[i][4], &h2.CXD[i][5],&h2.CXD[i][6]);
                            break;
                        case 2: //Y Force der
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf",&h1.CYD[i][0], &h1.CYD[i][1], &h1.CYD[i][2], &h1.CYD[i][3], &h1.CYD[i][4], &h1.CYD[i][5]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf", &h2.CYD[i][0], &h2.CYD[i][1], &h2.CYD[i][2], &h2.CYD[i][3], &h2.CYD[i][4], &h2.CYD[i][5]);
                            break;
                        case 3: //Z F der
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf %lf", &h1.CZD[i][0], &h1.CZD[i][1], &h1.CZD[i][2], &h1.CZD[i][3], &h1.CZD[i][4], &h1.CZD[i][5], &h1.CZD[i][6]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf %lf", &h2.CZD[i][0], &h2.CZD[i][1], &h2.CZD[i][2], &h2.CZD[i][3], &h2.CZD[i][4], &h2.CZD[i][5],&h2.CZD[i][6]);
                            break;
                        case 4: //Roll Mom der
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf", &h1.CLD[i][0], &h1.CLD[i][1], &h1.CLD[i][2], &h1.CLD[i][3], &h1.CLD[i][4], &h1.CLD[i][5]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf", &h2.CLD[i][0], &h2.CLD[i][1], &h2.CLD[i][2], &h2.CLD[i][3], &h2.CLD[i][4], &h2.CLD[i][5]);
                            break;
                        case 5: //Pitch Mom der
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf %lf", &h1.CMD[i][0], &h1.CMD[i][1], &h1.CMD[i][2], &h1.CMD[i][3], &h1.CMD[i][4], &h1.CMD[i][5], &h1.CMD[i][6]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf %lf", &h2.CMD[i][0], &h2.CMD[i][1], &h2.CMD[i][2], &h2.CMD[i][3], &h2.CMD[i][4], &h2.CMD[i][5], &h2.CMD[i][6]);
                            break;
                        case 6: //Yaw Mom der
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf", &h1.CND[i][0], &h1.CND[i][1], &h1.CND[i][2], &h1.CND[i][3], &h1.CND[i][4], &h1.CND[i][5]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf",  &h2.CND[i][0], &h2.CND[i][1], &h2.CND[i][2], &h2.CND[i][3], &h2.CND[i][4], &h2.CND[i][5]);
                            break;
                        case 7: //Control force der.
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf", &h1.CUF[i][0], &h1.CUF[i][1], &h1.CUF[i][2], &h1.CUF[i][3], &h1.CUF[i][4], &h1.CUF[i][5]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf",  &h2.CUF[i][0], &h2.CUF[i][1], &h2.CUF[i][2], &h2.CUF[i][3], &h2.CUF[i][4], &h2.CUF[i][5]);
                            break;
                        case 8: //Control mom. der.
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf", &h1.CUM[i][0], &h1.CUM[i][1], &h1.CUM[i][2], &h1.CUM[i][3], &h1.CUM[i][4], &h1.CUM[i][5]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf", &h2.CUM[i][0], &h2.CUM[i][1], &h2.CUM[i][2], &h2.CUM[i][3], &h2.CUM[i][4], &h2.CUM[i][5]);
                            break;
                        case 9: //Rotary mom. der.
                            sscanf(s1, "%*lf %lf %lf %lf %lf %lf %lf", &h1.CRD[i][0], &h1.CRD[i][1], &h1.CRD[i][2], &h1.CRD[i][3], &h1.CRD[i][4], &h1.CRD[i][5]);
                            sscanf(s2, "%*lf %lf %lf %lf %lf %lf %lf", &h2.CRD[i][0], &h2.CRD[i][1], &h2.CRD[i][2], &h2.CRD[i][3], &h2.CRD[i][4], &h2.CRD[i][5]);
                            break;
                        default:
                            break;
                    }
                }
            }
            j++; // Look up table number increment
            //fprintf(fp_tmp, "\n\n"); //Separation between tables
        }
    }

    fclose(fp1);
    fclose(fp2);

    // Height interpolation
    for(i=0;i<110;i++){
        for(j=0; j<7; j++){ // Sovrascrivo dati interpolati su h1
            if(j<6) {
                h1.CSS[i][j] = C_h_interp(h1.CSS[i][j], h2.CSS[i][j], qp, qs, h);
                h1.CYD[i][j] = C_h_interp(h1.CYD[i][j], h2.CYD[i][j], qp,qs, h);
                h1.CLD[i][j] = C_h_interp(h1.CLD[i][j], h2.CLD[i][j], qp,qs,h);
                h1.CND[i][j] = C_h_interp(h1.CND[i][j], h2.CND[i][j], qp,qs, h);
                h1.CUF[i][j] = C_h_interp(h1.CUF[i][j], h2.CUF[i][j], qp,qs, h);
                h1.CUM[i][j] = C_h_interp(h1.CUM[i][j], h2.CUM[i][j], qp,qs, h);
                h1.CRD[i][j] = C_h_interp(h1.CRD[i][j], h2.CRD[i][j], qp,qs, h);
                // Assigning variables
                css[i][j]=h1.CSS[i][j];
                cyd[i][j]=h1.CYD[i][j];
                cld[i][j]=h1.CLD[i][j];
                cnd[i][j]=h1.CND[i][j];
                cuf[i][j]=h1.CUF[i][j];
                cum[i][j]=h1.CUM[i][j];
                crd[i][j]=h1.CRD[i][j];

            }
            h1.CXD[i][j]=C_h_interp(h1.CXD[i][j], h2.CXD[i][j], qp,qs, h);
            h1.CZD[i][j]=C_h_interp(h1.CZD[i][j], h2.CZD[i][j], qp,qs, h);
            h1.CMD[i][j]=C_h_interp(h1.CMD[i][j], h2.CMD[i][j], qp,qs, h);
            // Assigning variables
            cxd[i][j]=h1.CXD[i][j];
            czd[i][j]=h1.CZD[i][j];
            cmd[i][j]=h1.CMD[i][j];
        }
    }
};

void ReadAndSavePropeller(char *file_propeller, double *datipropeller, double *CSI, double *RD, double *CH_AD, double *BA){

    char line[100];
    FILE *fp,*fp_val;

    fp = fopen(file_propeller,"r");
    fp_val = fopen("Validation/propeller.ini","w");

    if (fp == NULL){
        printf("[!]ERRORE nella lettura del file '%s'\n", file_propeller);
        exit(5);
    }

    if (fp_val == NULL){
        printf("\n[!]ERRORE Nella creazione/apertura file dba di validazione.\n");
        exit(6);
    }

    while (fgets(line, sizeof(line), fp)) {

        if (line[0] == '*') {
            continue;  // Se trovo '*' vado alla riga successiva
        }

        if (strstr(line, "CSI") != NULL) { //Cerco l'inizio e la fine dei vettori nel file in modo da puntarli direttamente dopo

            while (fgets(line, sizeof(line), fp)) {
                if (line[0] != '*') {
                    sscanf(line, "%lf %lf %lf %lf", CSI, RD, CH_AD, BA);
                    CSI++;RD++;CH_AD++;BA++;
                }
            }
        }

        fprintf(fp_val, "%s", line); // Salvo tutto nel file di verifica "propeller.ini"
        sscanf(line, "%lf", datipropeller); // Salvo i valori presenti nel file in delle variabili
        datipropeller++;
    }

    fclose(fp);
    fclose(fp_val);
};

