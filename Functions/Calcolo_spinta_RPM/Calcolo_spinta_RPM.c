#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

//#include "propel.h"

#define pi 3.14159265

double propel_trim(double T_trim, double alpha_trim, double deltae_trim, double Vel, double rho1){

    double ch_ad[101] = {0.1055,0.1131,0.1207,0.1283,0.1360,0.1436,0.1512,0.1588,0.1664,0.1740,0.1816,0.1892,0.1968,0.2006,
                         0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,0.2006,
                         0.2006,0.2006,0.2004,0.1998,0.1993,0.1987,0.1981,0.1976,0.1970,0.1964,0.1958,0.1953,0.1947,0.1941,0.1936,0.1930,0.1924,0.1918,
                         0.1913,0.1907,0.1901,0.1896,0.1890,0.1884,0.1879,0.1873,0.1854,0.1835,0.1816,0.1797,0.1778,0.1759,0.1740,0.1721,0.1702,0.1683,
                         0.1664,0.1645,0.1626,0.1607,0.1588,0.1569,0.1550,0.1531,0.1512,0.1493,0.1474,0.1445,0.1416,0.1386,0.1356,0.1326,0.1296,0.1267,
                         0.1237,0.1207,0.1177,0.1148,0.1118,0.1088,0.1058,0.1028,0.0999,0.0969,0.0939,0.0750,0.0545,0.0339,0.0134};

    double BA[101] = {39.1010,38.0058,36.9617,35.9658,35.0154,34.1080,33.2412,32.4126,31.6200,30.8616,30.1354,29.4395,28.7723,
                      28.1323,27.5180,26.9279,26.3608,25.8155,25.2909,24.7858,24.2993,23.8305,23.3784,22.9422,22.5211,22.1145,21.7216,21.3418,20.9744,20.6190,20.2748,19.9415,19.6186,19.3056,19.0020,
                      18.7075,18.4217,18.1442,17.8746,17.6127,17.3582,17.1106,16.8699,16.6356,16.4075,16.1855,15.9692,15.7585,15.5532,15.3530,15.1577,
                      14.9673,14.7815,14.6001,14.4231,14.2502,14.0813,13.9163,13.7550,13.5974,13.4432,13.2925,13.1450,13.0008,12.8596,12.7214,12.5861,12.4535,
                      12.3237,12.1966,12.0720,11.9499,11.8302,11.7128,11.5977,11.4848,11.3741,11.2654,11.1588,11.0542,10.9515,10.8506,10.7516,10.6543,10.5588,
                      10.4649,10.3727,10.2821,10.1930,10.1054,10.0194,9.9347,9.8515,9.7696,9.6891,9.6098,9.5319,9.4552,9.3797,9.3053,9.2322};

    double CS[101] = {0.2000,0.2080,0.2160,0.2240,0.2320,0.2400,0.2480,0.2560,0.2640,
                      0.2720,0.2800,0.2880,0.2960,0.3040,0.3120,0.3200,0.3280,0.3360,0.3440,0.3520,0.3600,0.3680,
                      0.3760,0.3840,0.3920,0.4000,0.4080,0.4160,0.4240,0.4320,0.4400,0.4480,0.4560,0.4640,0.4720,0.4800,
                      0.4880,0.4960,0.5040,0.5120,0.5200,0.5280,0.5360,0.5440,0.5520,0.5600,0.5680,0.5760,0.5840,0.5920,
                      0.6000,0.6080,0.6160,0.6240,0.6320,0.6400,0.6480,0.6560,0.6640,0.6720,0.6800,0.6880,0.6960,0.7040,
                      0.7120,0.7200,0.7280,0.7360,0.7440,0.7520,0.7600,0.7680,0.7760,0.7840,0.7920,0.8000,0.8080,0.8160,
                      0.8240,0.8320,0.8400,0.8480,0.8560,0.8640,0.8720,0.8800,0.8880,0.8960,0.9040,0.9120,0.9200,0.9280,
                      0.9360,0.9440,0.9520,0.9600,0.9680,0.9760,0.9840,0.9920,1.0000};

    double diam = 0.125; //diametro elica
    double Raggio=diam/2.0; //raggio elica

    double cl_alpha = 3.984900;
    double a_0 = -0.09616302;
    double RPM_ref=3600; //RPM di riferimento
    double rpm_trim = 0;


    while(RPM_ref <= 30000){

        double n=RPM_ref/60.0; //round-per-second
        double omega=n*2.0*pi; //velocià angolare [rad/s]

        //double coef1=(tip-hub)/(xt-xs); //coefficiente #1 di supporto al calcolo dell'angolo di svergolamento theta
        //double coef2=hub-coef1*xs+pitch; //coefficiente #2 di supporto al calcolo dell'angolo di svergolamento theta
        double r1[101]; //creazione vettore delle 60 stazioni meno le stazioni dell'ogiva (-> CSI in propeller.txt)

        double theta1;
        double alpha1;
        double t2[101], a2[101], b2[101];
        double th, phi1, eff, DtDr,DqDr,cl,cd,CT,CQ, tem1,tem2;
        double a, anew;
        double b, bnew;
        //double chord;

        //double coef1_new;
        //double coef2_new;
        double rstep_new;

        int j=0;
        int finished=0;
        /*Inizializzazione r1*/
        for(j=0;j<101;j++){
            //r1[j]=xs+j*rstep;               //    !! NB !! : r1 non è esattamente = a CSI
            //printf("\n r1=%f",r1[j]);   //      prova --- rimuovere riga
            r1[j] = CS[j]*Raggio;
            //r1[j] = rd[j];//CS[j]*Raggio;
        }

        double rad;

        double Vlocal,V0,V2;
        double thrust=0.0; //inizializzazione vettore spinta
        double torque=0.0;//inizializzazione vettore coppia che abbiamo già definito
        //*P_alone=0.0; // Potenza all'albero

        for(j=1; j<101; j++)
        {

            double chord =ch_ad[j]*Raggio;

            double xt_new = CS[j]*Raggio;
            double xs_new = CS[j-1]*Raggio;
            //double xt_new = rd[j]*Raggio;//CS[j]*Raggio;
            //double xs_new = rd[j-1]*Raggio;//CS[j-1]*Raggio;
            rstep_new = xt_new-xs_new;
            double tip_new = BA[j];
            double hub_new = BA[j-1];
            //coef1_new = (tip_new-hub_new)/(xt_new-xs_new); //coefficiente #1 di supporto al calcolo dell'angolo di svergolamento theta
            //coef2_new = hub_new - coef1_new*xs_new+pitch; //coefficiente #2 di supporto al calcolo dell'angolo di svergolamento theta

            //double chord=0.198;
            rad=r1[j]; //coordinata j-esima stazione (-> CSI in propeller.txt)
            //double rad_new = xt_new;
            //rad = rad_new;
            theta1=BA[j] + 10;
            //double theta1_bis=coef1*rad+coef2; //calcolo angolo di svergolamento della j-esima stazione       ---->  ERRATO
            //double theta1_new=coef1_new*rad_new+coef2_new; //calcolo angolo di svergolamento della j-esima stazione       ---->  ERRATO

            //sleep(1.0);
            //printf("theta1 da BA is %lf\n",theta1);
            //printf("theta1 da coef1 is %lf\n",theta1_bis);
            //printf("theta1_new is %lf\n\n",theta1_new);

            //printf("rad da CS is %lf\n",rad);
            //printf("rad_new is %lf\n\n",rad_new);

            t2[j]=theta1; //angolo di svergolamento della j-esima stazione (-> BA su propeller.txt)
            th=theta1/180.0*pi; //angolo di svergolamento [rad]
            a=0.1; //inizializzazione axial inflow factor (vedi pag.4 PROPEL.pdf)
            b=0.01; //inizializzazione angular inflow (swirl) factor (vedi pag.4 PROPEL.pdf)
            finished=0; //inizializzione flag
            int sum=1; //inizializzione variabile di supporto
//      if(j==0){
//          printf("\n\n\nCHORD_END=%lf\n\n\n",chord);
//          printf("\n\n\nBA_END=%lf\n\n\n",BA[12]);}
//      if(j==47)
//          printf("\n\n\nCSI=%f\nRD=%f\nCH=%f\nBA=%f\n",CSI[j+12],RD[j+12],CH[12+j],BA[12+j]);
//      if(j==46)

            //printf("chord is %lf\n",chord);
            //printf("xt_new is %lf\n",xt_new);
            //printf("xs_new is %lf\n",xs_new);
            //printf("rstep_new is %lf\n",rstep_new);
            //printf("coef1_new is %lf\n",coef1_new);
            //printf("coef2_new is %lf\n",coef2_new);

            while (finished==0)
            {
                V0=Vel*(1+a); //componente del flusso all'incirca uguale alla velocità di avanzamento del velivolo (Vinf), aumentata tramite l'axial inflow factor
                V2=omega*rad*(1-b); //componente del flusso all'incirca uguale alla velocità angolare della sezione della pala (omega*rad), ridotta tramite l'angular inflow factor
                phi1=atan2(V0,V2); //angolo tra le due componenti del flusso V0 e V2
                alpha1=th - phi1 - a_0; //angolo di attacco raltivo alla j-esima sezione della pala

                //sleep(1.0);
                //printf("V0 %lf\n",V0);
                //printf("V2 %lf\n",V2);
                //printf("phi1 %lf\n",phi1);
                //printf("alpha1 %lf\n",alpha1);
                //printf("a %lf\n\n",a);

                // NB: i valori di Cl e CD sarebbero da parametrizzare in fnz dei dati ottenuti dal file di testo propeller.txt

                cl=0.3832+ cl_alpha *alpha1; //L coefficiente di portanza
                cd=0.0235+0.154*alpha1+1.0476*alpha1*alpha1; // CD coefficiente di resistenza CD = CD0+CD1*CL+CD2*CL^2 (NB nel nostro caso, CD = CD0+CD_alpha*alpha+CD_alpha2*alpha^2 -> slide lezione 2)

                Vlocal=sqrt(V0*V0+V2*V2); // velocità locale del flusso
                CT = cl*cos(phi1)-cd*sin(phi1); //CT coefficiente di spinta adimensionale
                DtDr=0.5*rho1*Vlocal*Vlocal*2.0*chord*CT; //contributo di spinta della j-esima sezione
                CQ = cd*cos(phi1)+cl*sin(phi1); //CQ coefficiente di coppia adimensionale
                DqDr=0.5*rho1*Vlocal*Vlocal*2.0*chord*rad*CQ; //contributo di coppia della j-esima sezione
                tem1=DtDr/(4.0*pi*rad*rho1*Vel*Vel*(1+a)); //fattore correttivo del coefficiente "a"
                tem2=DqDr/(4.0*pi*rad*rad*rad*rho1*Vel*(1+a)*omega); //fattore correttivo del coefficiente "b"
                anew=0.5*(a+tem1); //nuovo valore coefficiente "a"
                bnew=0.5*(b+tem2); //nuovo valore coefficiente "b"
                //processo iterativo per arrivare a convergenza
                if (fabs(anew-a)<1/100000){
                    if (fabs(bnew-b)<1/100000){
                        finished=1;
                    }
                }
                a=anew; //definizione valore finale coefficiente "a"
                b=bnew; //definizione valore finale coefficiente "b"
                sum=sum+1;
                if (sum>500){
                    finished=1;
                }
            }
            a2[j]=a; //definizione valore finale coefficiente "a" per la j-esima stazione
            b2[j]=b; //definizione valore finale coefficiente "b" per la j-esima stazione
            thrust=thrust+DtDr*rstep_new; //sommatoria dei contributi di spinta dalla stazione 1 alla stazione j
            torque=torque+DqDr*rstep_new; //sommatoria dei contributi di coppia dalla stazione 1 alla stazione j

            //printf("thrust temp is %lf\n",thrust);
        }

        //printf("**********************\n");
        //printf("** %f RPM ** \n",RPM_ref);
        //printf("Thrust: %f N\n",thrust);
        //printf("Torque: %f Nm\n",torque);
        //printf("coef1: %f Nm\n",coef1_new);
        //printf("coef2: %f Nm\n",coef2_new);
        //printf("alpha: %f Nm\n",alpha1);
        //printf("Vel: %f Nm\n",Vel);
        //printf("rstep: %f Nm\n",rstep_new);
        //printf("**********************\n\n\n");

        double t = thrust/(rho1*n*n*diam*diam*diam*diam); //coefficiente di spinta adimensionale
        double q_torque= torque/(rho1*n*n*diam*diam*diam*diam*diam); //coefficiente di coppia adimensionale
        double J = Vel/(n*diam); //rapporto di avanzamento;
        if (t < 0){
            eff = 0.0; //efficienza elica
        }else{
            eff = t / q_torque*J/(2.0*M_PI); //efficienza elica
        }

        if(fabs(T_trim-thrust) < 0.01){
            rpm_trim = RPM_ref;
        }

        RPM_ref = RPM_ref+100;
    }



    printf("\n");
    printf("Spinta di trim: %lf [N]\n", T_trim);
    printf("RPM di trim: %lf \n", rpm_trim);

    // Calcolo Potenza all'albero
    //*P_alone=torque*omega/1000;




    //printf("**********************\n");
    //printf("Thrust: %f N \n",thrust);
    //printf("Torque: %f Nm \n",torque);
    //printf("**********************\n");

    //double t=thrust/(rho1*n*n*diam*diam*diam*diam); //coefficiente di spinta adimensionale
    //double q=torque/(rho1*n*n*diam*diam*diam*diam*diam); //coefficiente di coppia adimensionale
    //double J=Vel/(n*diam); //rapporto di avanzamento
    //if (t<0){
    //    eff=0.0; //efficienza elica
    //}else{
    //    eff=t/q*J/(2.0*pi); //efficienza elica
    //}


    return rpm_trim;
}

double propel(double n, double Vel, double rho1) {

    double ch_ad[101] = {0.1055, 0.1131, 0.1207, 0.1283, 0.1360, 0.1436, 0.1512, 0.1588, 0.1664, 0.1740, 0.1816, 0.1892,
                         0.1968, 0.2006,
                         0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006,
                         0.2006, 0.2006, 0.2006, 0.2006, 0.2006, 0.2006,
                         0.2006, 0.2006, 0.2004, 0.1998, 0.1993, 0.1987, 0.1981, 0.1976, 0.1970, 0.1964, 0.1958, 0.1953,
                         0.1947, 0.1941, 0.1936, 0.1930, 0.1924, 0.1918,
                         0.1913, 0.1907, 0.1901, 0.1896, 0.1890, 0.1884, 0.1879, 0.1873, 0.1854, 0.1835, 0.1816, 0.1797,
                         0.1778, 0.1759, 0.1740, 0.1721, 0.1702, 0.1683,
                         0.1664, 0.1645, 0.1626, 0.1607, 0.1588, 0.1569, 0.1550, 0.1531, 0.1512, 0.1493, 0.1474, 0.1445,
                         0.1416, 0.1386, 0.1356, 0.1326, 0.1296, 0.1267,
                         0.1237, 0.1207, 0.1177, 0.1148, 0.1118, 0.1088, 0.1058, 0.1028, 0.0999, 0.0969, 0.0939, 0.0750,
                         0.0545, 0.0339, 0.0134};

    double BA[101] = {39.1010, 38.0058, 36.9617, 35.9658, 35.0154, 34.1080, 33.2412, 32.4126, 31.6200, 30.8616, 30.1354,
                      29.4395, 28.7723,
                      28.1323, 27.5180, 26.9279, 26.3608, 25.8155, 25.2909, 24.7858, 24.2993, 23.8305, 23.3784, 22.9422,
                      22.5211, 22.1145, 21.7216, 21.3418, 20.9744, 20.6190, 20.2748, 19.9415, 19.6186, 19.3056, 19.0020,
                      18.7075, 18.4217, 18.1442, 17.8746, 17.6127, 17.3582, 17.1106, 16.8699, 16.6356, 16.4075, 16.1855,
                      15.9692, 15.7585, 15.5532, 15.3530, 15.1577,
                      14.9673, 14.7815, 14.6001, 14.4231, 14.2502, 14.0813, 13.9163, 13.7550, 13.5974, 13.4432, 13.2925,
                      13.1450, 13.0008, 12.8596, 12.7214, 12.5861, 12.4535,
                      12.3237, 12.1966, 12.0720, 11.9499, 11.8302, 11.7128, 11.5977, 11.4848, 11.3741, 11.2654, 11.1588,
                      11.0542, 10.9515, 10.8506, 10.7516, 10.6543, 10.5588,
                      10.4649, 10.3727, 10.2821, 10.1930, 10.1054, 10.0194, 9.9347, 9.8515, 9.7696, 9.6891, 9.6098,
                      9.5319, 9.4552, 9.3797, 9.3053, 9.2322};

    double CS[101] = {0.2000, 0.2080, 0.2160, 0.2240, 0.2320, 0.2400, 0.2480, 0.2560, 0.2640,
                      0.2720, 0.2800, 0.2880, 0.2960, 0.3040, 0.3120, 0.3200, 0.3280, 0.3360, 0.3440, 0.3520, 0.3600,
                      0.3680,
                      0.3760, 0.3840, 0.3920, 0.4000, 0.4080, 0.4160, 0.4240, 0.4320, 0.4400, 0.4480, 0.4560, 0.4640,
                      0.4720, 0.4800,
                      0.4880, 0.4960, 0.5040, 0.5120, 0.5200, 0.5280, 0.5360, 0.5440, 0.5520, 0.5600, 0.5680, 0.5760,
                      0.5840, 0.5920,
                      0.6000, 0.6080, 0.6160, 0.6240, 0.6320, 0.6400, 0.6480, 0.6560, 0.6640, 0.6720, 0.6800, 0.6880,
                      0.6960, 0.7040,
                      0.7120, 0.7200, 0.7280, 0.7360, 0.7440, 0.7520, 0.7600, 0.7680, 0.7760, 0.7840, 0.7920, 0.8000,
                      0.8080, 0.8160,
                      0.8240, 0.8320, 0.8400, 0.8480, 0.8560, 0.8640, 0.8720, 0.8800, 0.8880, 0.8960, 0.9040, 0.9120,
                      0.9200, 0.9280,
                      0.9360, 0.9440, 0.9520, 0.9600, 0.9680, 0.9760, 0.9840, 0.9920, 1.0000};

    double rd[101] = {0.0119, 0.0124, 0.0129, 0.0134, 0.0139, 0.0143, 0.0148, 0.0153, 0.0158, 0.0162, 0.0167,
                      0.0172, 0.0177, 0.0182, 0.0186, 0.0191, 0.0196, 0.0201, 0.0205, 0.0210, 0.0215, 0.0220, 0.0224,
                      0.0229, 0.0234,
                      0.0239, 0.0244, 0.0248, 0.0253, 0.0258, 0.0263, 0.0267, 0.0272, 0.0277, 0.0282, 0.0287, 0.0291,
                      0.0296, 0.0301,
                      0.0306, 0.0310, 0.0315, 0.0320, 0.0325, 0.0330, 0.0334, 0.0339, 0.0344, 0.0349, 0.0353, 0.0358,
                      0.0363, 0.0368,
                      0.0373, 0.0377, 0.0382, 0.0387, 0.0392, 0.0396, 0.0401, 0.0406, 0.0411, 0.0415, 0.0420, 0.0425,
                      0.0430, 0.0435,
                      0.0439, 0.0444, 0.0449, 0.0454, 0.0458, 0.0463, 0.0468, 0.0473, 0.0478, 0.0482, 0.0487, 0.0492,
                      0.0497, 0.0501,
                      0.0506, 0.0511, 0.0516, 0.0521, 0.0525, 0.0530, 0.0535, 0.0540, 0.0544, 0.0549, 0.0554, 0.0559,
                      0.0564, 0.0568,
                      0.0573, 0.0578, 0.0583, 0.0587, 0.0592, 0.0595};

    double diam = 0.125; //diametro elica
    double Raggio = diam / 2.0; //raggio elica
    double cl_alpha = 3.984900;
    double a_0 = -0.09616302;

        double omega = n * 2.0 * pi/60; //velocià angolare [rad/s]
        double r1[101]; //creazione vettore delle 60 stazioni meno le stazioni dell'ogiva (-> CSI in propeller.txt)

        double theta1;
        double alpha1;
        double t2[101], a2[101], b2[101];
        double th, phi1, eff, DtDr, DqDr, cl, cd, CT, CQ, tem1, tem2;
        double a, anew;
        double b, bnew;
        //double chord;

        //double coef1_new;
        //double coef2_new;
        double rstep_new;

        int j = 0;
        int finished = 0;
        /*Inizializzazione r1*/
        for (j = 0; j < 101; j++) {
            //r1[j]=xs+j*rstep;               //    !! NB !! : r1 non è esattamente = a CSI
            //printf("\n r1=%f",r1[j]);   //      prova --- rimuovere riga
            r1[j] = CS[j] * Raggio;
            //r1[j] = rd[j];//CS[j]*Raggio;
        }

        double rad;

        double Vlocal, V0, V2;
        double thrust = 0.0; //inizializzazione vettore spinta
        double torque = 0.0;//inizializzazione vettore coppia che abbiamo già definito
        //*P_alone=0.0; // Potenza all'albero

        for (j = 1; j < 101; j++) {

            double chord = ch_ad[j] * Raggio;

            double xt_new = CS[j] * Raggio;
            double xs_new = CS[j - 1] * Raggio;
            rstep_new = xt_new - xs_new;
            double tip_new = BA[j];
            double hub_new = BA[j - 1];

            rad = r1[j]; //coordinata j-esima stazione (-> CSI in propeller.txt)

            theta1 = BA[j] + 10;

            t2[j] = theta1; //angolo di svergolamento della j-esima stazione (-> BA su propeller.txt)
            th = theta1 / 180.0 * pi; //angolo di svergolamento [rad]
            a = 0.1; //inizializzazione axial inflow factor (vedi pag.4 PROPEL.pdf)
            b = 0.01; //inizializzazione angular inflow (swirl) factor (vedi pag.4 PROPEL.pdf)
            finished = 0; //inizializzione flag
            int sum = 1; //inizializzione variabile di supporto


            while (finished == 0) {
                V0 = Vel * (1 + a); //componente del flusso all'incirca uguale alla velocità di avanzamento del velivolo (Vinf), aumentata tramite l'axial inflow factor
                V2 = omega * rad * (1 - b); //componente del flusso all'incirca uguale alla velocità angolare della sezione della pala (omega*rad), ridotta tramite l'angular inflow factor
                phi1 = atan2(V0, V2); //angolo tra le due componenti del flusso V0 e V2
                alpha1 = th - phi1 - a_0; //angolo di attacco raltivo alla j-esima sezione della pala

                // NB: i valori di Cl e CD sarebbero da parametrizzare in fnz dei dati ottenuti dal file di testo propeller.txt

                cl = 0.3832 + cl_alpha * alpha1; //L coefficiente di portanza
                cd = 0.0235 + 0.154 * alpha1 + 1.0476 * alpha1 *
                                               alpha1; // CD coefficiente di resistenza CD = CD0+CD1*CL+CD2*CL^2 (NB nel nostro caso, CD = CD0+CD_alpha*alpha+CD_alpha2*alpha^2 -> slide lezione 2)

                Vlocal = sqrt(V0 * V0 + V2 * V2); // velocità locale del flusso
                CT = cl * cos(phi1) - cd * sin(phi1); //CT coefficiente di spinta adimensionale
                DtDr = 0.5 * rho1 * Vlocal * Vlocal * 2.0 * chord * CT; //contributo di spinta della j-esima sezione
                CQ = cd * cos(phi1) + cl * sin(phi1); //CQ coefficiente di coppia adimensionale
                DqDr = 0.5 * rho1 * Vlocal * Vlocal * 2.0 * chord * rad * CQ; //contributo di coppia della j-esima sezione
                tem1 = DtDr / (4.0 * pi * rad * rho1 * Vel * Vel * (1 + a)); //fattore correttivo del coefficiente "a"
                tem2 = DqDr / (4.0 * pi * rad * rad * rad * rho1 * Vel * (1 + a) * omega); //fattore correttivo del coefficiente "b"
                anew = 0.5 * (a + tem1); //nuovo valore coefficiente "a"
                bnew = 0.5 * (b + tem2); //nuovo valore coefficiente "b"
                //processo iterativo per arrivare a convergenza
                if (fabs(anew - a) < 1 / 100000) {
                    if (fabs(bnew - b) < 1 / 100000) {
                        finished = 1;
                    }
                }
                a = anew; //definizione valore finale coefficiente "a"
                b = bnew; //definizione valore finale coefficiente "b"
                sum = sum + 1;
                if (sum > 500) {
                    finished = 1;
                }
            }
            a2[j] = a; //definizione valore finale coefficiente "a" per la j-esima stazione
            b2[j] = b; //definizione valore finale coefficiente "b" per la j-esima stazione
            thrust = thrust + DtDr * rstep_new; //sommatoria dei contributi di spinta dalla stazione 1 alla stazione j
            torque = torque + DqDr * rstep_new; //sommatoria dei contributi di coppia dalla stazione 1 alla stazione j

        }

    return thrust;
}