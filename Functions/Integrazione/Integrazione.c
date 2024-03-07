#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "Integrazione.h"
#include "../Interpolazione/Interpolazione.h"
#include "../Calcolo_spinta_RPM/Calcolo_spinta_RPM.h"


void integration(double trim_cond[], double trim_control[],double alpha_trim,double rho, double S, double b,
                 double c, double m, double Jx, double Jy, double Jz, double Jxz, double V_ref, double h_ref, double Vmin, int inp){

    const double pi=4*atan(1);
    const double g=9.81; // [m/s^2]
    const double d2r=pi/180;
    double tf;
    double dt=0.01; //Passo di integrazione fissato
    int n; //Numero di passi
    double u,v,w,p,q,r,phi,theta,psi,h,x,y,alpha,beta,Da,De,Dr=0,Dth; // Initialized from trim conditions
    double du,dv,dw,dp,dq,dr,dphi,dtheta,dpsi,dh,dx,dy, drp=0,dpp=0;
    double X,Y,Z,L,M,N,T; //Forze e Momenti
    double de_max=15*d2r, de_min=-15*d2r;//[rad]
    double de_dot_max=4, de_dot_min=-4; //[deg/sec]
    double dth_max=0.55, dth_min=-0.45;
    double dth_dot_max=4, dth_dot_min=-4;
    double h_max=2000; //Quota di tangenza
    double RPM;
    double Time=0;
    double **C=MatrixCoeff(10,10);

    u=trim_cond[0];v=trim_cond[1];w=trim_cond[2];p=trim_cond[3];q=trim_cond[4];r=trim_cond[5];phi=trim_cond[6];theta=trim_cond[7]*d2r;psi=trim_cond[8];
    h=trim_cond[9];x=trim_cond[10];y=trim_cond[11];
    De=trim_control[0];Da=trim_control[1];Dth=trim_control[2];
    double psi_ref = psi;
    double psi_vec[35000];
    double V = sqrt(u*u+v*v+w*w);

//Aerodynamic coefficients definition
    double CX=C[0][0], CY=C[0][1], CZ=C[0][2], CL=C[0][3], CM=C[0][4], CN=C[0][5];
    double CXA=C[1][0], CXB=C[1][4], CXP=C[1][5], CXQ=C[1][3], CXR=C[1][6], CXDE=C[7][0];
    double CYA=C[2][4], CYB=C[2][0], CYP=C[2][2], CYQ=C[2][5], CYR=C[2][3], CYDA=C[7][4];
    double CZA=C[3][0], CZB=C[3][4], CZP=C[3][5], CZQ=C[3][3], CZR=C[3][6], CZDE=C[7][2];
    double CLA=C[4][4],CLB=C[4][0],CLP=C[4][2],CLQ=C[4][5],CLR=C[4][3],CLDA=C[8][0],CLDR=C[8][1];
    double CMA=C[5][0], CMB=C[5][4], CMP=C[5][5], CMQ=C[5][3], CMR=C[5][6], CMDE=C[8][2];
    double CNA=C[6][4],CNB=C[6][0],CNP=C[6][2],CNQ=C[6][5],CNR=C[6][3],CNDA=C[8][4],CNDR=C[8][5];

//File to write results down
    FILE *f_int,*f_int2;
    f_int=fopen("Integrazione.txt","w");
    f_int2=fopen("Control_vector.txt","w");

    if(f_int==NULL){
        printf("\n[X]ERRORE: apertura file di salvataggio variabili di stato non riuscita\n");
        exit(6);
    }
    if(f_int2==NULL){
        printf("\n[X]ERRORE: apertura file di salvataggio variabili di controllo non riuscita\n");
        exit(7);
    }

//Mettere errori per passo e tempo di integrazione
    // Chiedere estremi di integrazione a utente
    if (inp == 1) {
        printf("\nINTEGRAZIONE: \n");
        printf("\nInserire tempo integrazione in [s]: ");
        scanf("%lf", &tf);
    }

    if (inp == 2){
        //Variabili algoritmo guida (Lettura del vettore di psi)
        ReadPsi(psi_vec,&tf,&x,&y);
    }

    n=(int)(tf/dt+1);

    fprintf(f_int,"Time[s]   u[m/s]  v[m/s]  w[m/s]  p[rad/s]  q[rad/s]  r[rad/s]  phi[rad]  theta[rad]  psi[rad]  h[m]  x[m]  y[m]  V[m/s]  alpha[rad]  beta[rad]\n\n");
    fprintf(f_int,"%0.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", Time,u,v,w,p,q,r,phi,theta,psi,h,x,y,V_ref,alpha_trim,asin(v/V_ref));

    fprintf(f_int2,"Times[s]\t\tDe[°]\t\tDa[°]\t\tDth[%]\n");
    fprintf(f_int2,"%0.5lf\t%0.5lf\t%0.5lf\t%0.5lf\n",Time,De*180/pi,Da*180/pi,Dth);

    //Variabili necessarie per le funzioni del controllore
    double I_psi=0, I_phi=0, D_psi=0, D_phi=0, e_phip=0, e_psip=0,Da_c=0; // Variabili controllore latero_direzionale
    double I_V=0, I_h=0, I_th=0, D_th=0, D_V=0, D_h=0, e_Vp=0, e_thp=0, e_hp=0,De_c=0,Dth_c=0; // Variabili controllore longitudinale

    //Variabili attuatore Elevon/Manetta + Flags di saturazione
    double u_a,x1a=Da,x2a=0,z1h=Dth,z2h=0,u_h;
    double u_e,x1e=De,x2e=0;

    int fl_de=0, fl_de_p=0, fl_da=0, fl_da_p=0;
    int fl_dth=0, fl_dth_p=0;
    int first_time_V=0, first_time_De=0, first_time_Da=0, first_time_Dth=0;


    // Ricorda di inizializzare X,U,alpha,beta

    for(int i=0; i<n; i++) {

        if(h>=h_max){
            printf("\n[X]ERRORE: La quota non rientra nei limiti fisici [Superata la quota di tangenza]\n");
            exit(7);
        }else if (h<=0){
            printf("\n[X]ERRORE: Hai raggiunto il suolo! \n");
            exit(7);
        }

        double Da_max=0.1745;
        double De_max=0.1745;
        //Autopilots inputs to actuators
        if(inp==2) {
            //Psi comandato
            psi_ref = psi_vec[i];
            //Autopilot
            long_controller(i, dt, V_ref, h_ref, theta, V, h, &I_V, &D_V, &e_Vp, &I_th, &D_th, &e_thp, &I_h, &D_h,
                            &e_hp, &Dth_c, &De_c);
            lat_dir_controller(i, dt, psi_ref, psi, phi, &I_psi, &I_phi, &D_psi, &D_phi, &e_phip, &e_psip,
                               &Da_c); //Psi_ref from guidance
            Da = Da_c+trim_control[1];
            De = De_c+trim_control[0];
            Dth = Dth_c+trim_control[2];
        }
        saturation(10*d2r,-10*d2r,&Da,&fl_da,&fl_da_p);
        saturation(10*d2r,-10*d2r,&De,&fl_de,&fl_de_p);
        saturation(1,0.1,&Dth,&fl_dth,&fl_dth_p);

        if(fl_da==1 && fl_da_p==0 && first_time_Da==0){
            printf("\n[!]WARNING: Deflessione di elevone dedicata all'alettone satura durante la manovra.\n");
            first_time_Da = 1;
        }
        if(fl_de==1 && fl_de_p==0 && first_time_De==0){
            printf("\n[!]WARNING: Deflessione di elevone dedicata all'equilibratore satura durante la manovra.\n");
            first_time_De = 1;
        }
        if(fl_dth==1 && fl_dth_p==0 && first_time_Dth==0){
            printf("\n[!]WARNING: Manetta satura durante la manovra.\n");
            first_time_Dth = 1;
        }

        //Definizione Forze N.B. Variabili adimensionali come p_,q_ da dichiarare sopra!!!
        X=0.5*rho*S*(V*V)*(CX+CXA*alpha+CXB*beta+CXP*(p*b/(2*V))+CXQ*(q*c/(2*V))+CXR*(r*b/(2*V))+CXDE*De); //CXDLE Non sappiamo cosa sia
        Y=0.5*rho*S*(V*V)*(CY+CYA*alpha+CYB*beta+CYP*(p*b/(2*V))+CYQ*(q*c/(2*V))+CYR*(r*b/(2*V))+CYDA*Da);
        Z=0.5*rho*S*(V*V)*(CZ+CZA*alpha+CZB*beta+CZP*(p*b/(2*V))+CZQ*(q*c/(2*V))+CZR*(r*b/(2*V))+CZDE*De);
        L=0.5*rho*S*b*(V*V)*(CL+CLA*alpha+CLB*beta+CLP*(p*b/(2*V))+CLQ*(q*c/(2*V))+CLR*(r*b/(2*V))+CLDA*Da);
        M=0.5*rho*S*c*(V*V)*(CM+CMA*alpha+CMB*beta+CMP*(p*b/(2*V))+CMQ*(q*c/(2*V))+CMR*(r*b/(2*V))+CMDE*De); //cmdle? ETCC.. VERIFICARE TUTTI I COEFFNON USATI
        N=0.5*rho*S*b*(V*V)*(CN+CNA*alpha+CNB*beta+CNP*(p*b/(2*V))+CNQ*(q*c/(2*V))+CNR*(r*b/(2*V))+CNDA*Da);

        // saturation(1,0.1,&Dth,&fl_th,&fl_rate_th);
        //RPM = 3600+(30000-3600)*(Dth+0.5)/(0.55+0.5);
        //RPM = 3600+(30000-3600)*(trim_control[2]-0.1)/(1-0.1);
        //Dth = fabs(Dth);
        RPM = 3600+(30000-3600)*(Dth-0.1)/(1-0.1);

        T=propel(RPM,V,rho);
        //printf("%lf\t%lf\t %lf\n",T,RPM,Dth);

        //N.B. per il primo ciclo xp (stato precedente) viene dalle trim_conditions
        // N.B. dr e dp a dx di eq. da condizioni iniziali
        //Jxz=0;
        du=(r*v-q*w)-g*sin(theta)+X/m+T/m;
        dv=(p*w-r*u)+g*sin(phi)*cos(theta)+Y/m;
        dw=(q*u-p*v)+g*cos(phi)*cos(theta)+Z/m;
        dp=-(Jz-Jy)*q*r/Jx+(p*q+drp)*Jxz/Jx+L/Jx;
        dq=-(Jx-Jz)*p*r/Jy-(p*p-r*r)*Jxz/Jy+M/Jy;
        dr=-(Jy-Jx)*p*q/Jz-(q*r-dpp)*Jxz/Jz+N/Jz;
        dphi=p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta);
        dtheta=q*cos(phi)-r*sin(phi);
        dpsi=q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta);
        dh=(-u*sin(theta)+v*cos(theta)*sin(phi)+w*cos(theta)*cos(phi));
        dx=u*cos(psi)*cos(theta)+v*(cos(psi)*sin(theta)*sin(phi)-sin(psi)*cos(phi))+w*(cos(psi)*sin(theta)*cos(psi)+sin(psi)*sin(phi));
        dy=u*sin(psi)*cos(theta)+v*(sin(psi)*sin(theta)*sin(phi)+cos(psi)*cos(phi))+w*(sin(psi)*sin(theta)*cos(phi)-cos(psi)*sin(phi));

        // Saving precedent r_dot and p_dot for next iteration
        drp=dr;
        dpp=dp;
        //Updating state variables
        u=u+dt*du;
        v=v+dt*dv;
        w=w+dt*dw;
        p=p+dt*dp;
        q=q+dt*dq;
        r=r+dt*dr;
        phi=phi+dt*dphi;
        theta=theta+dt*dtheta;
        psi=psi+dt*dpsi;
        h=h+dt*dh;
        x=x+dx*dt;
        y=y+dy*dt;

        V=sqrt(u*u+v*v+w*w);
        if (V<Vmin){
            if (first_time_V==0) {
                printf("\n[!]WARNING La velocita' è minore di quella di stallo e quindi è stata settata pari a V_min = %lf\n",
                       Vmin);
            }
            V = Vmin;
            first_time_V = 1;
        }

        alpha=atan2(w,u);
        beta=asin(v/V);

        Time=Time+dt;

        fprintf(f_int,"%0.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\t%.5lf\n", Time,u,v,w,p,q,r,phi,theta,psi,h,x,y,V,alpha,beta);
        fprintf(f_int2,"%0.5lf\t%0.5lf\t%0.5lf\t%0.5lf\n",Time,De*180/pi,Da*180/pi,Dth);

    }


    fclose(f_int);
}

void lat_dir_controller(int i, double dt, double psi_ref, double psi, double phi,double*I_psi, double*I_phi, double*D_psi, double*D_phi,
                        double*e_phip, double*e_psip, double*da_c){
    // Controllo latero direzionale
    double pi = 4*atan(1);
    double Kp_psi=1.5, Kp_phi=0.12, Ki_psi=0.005, Ki_phi=0.0005, Kd_phi=0.001, Kd_psi=0.01, N_phi=2*pi*10, N_psi=2*pi*10;
    double phi_ref, e_phi, e_psi;
    double PID_psi, PID_phi;

    //Attenzione a gradi e radianti
    // State space controllore
    e_psi=psi_ref-psi;

    if (e_psi<-pi)
    {
        e_psi=e_psi+2*pi;
    }
    else if (e_psi>pi)
    {
        e_psi=e_psi-2*pi;
    }


    PID(i,Kp_psi ,Ki_psi,Kd_psi, N_psi, dt, e_psi,I_psi,D_psi,e_psip,&PID_psi);
    phi_ref=PID_psi; //=PI_psi

    e_phi=phi_ref-phi;
    PID(i,Kp_phi ,Ki_phi,Kd_phi, N_phi, dt, e_phi,I_phi,D_phi,e_phip,&PID_phi);
    *da_c=PID_phi;//Input attuatore
}

void long_controller(int i,double dt, double V_ref, double h_ref, double theta, double V, double h, double*I_V, double*D_V, double*e_Vp,
                     double*I_th, double*D_th, double*e_thp,double*I_h, double*D_h, double*e_hp, double*dth_c, double*de_c){

    // Controllo longitudinale
    const double pi=4*atan(1);
    double Kp_V=-0.0021 ,Ki_V=-0.00087, Kd_V=-0.0015, tau_V=1/(2*pi*10), N_V=1/tau_V;
    double Kp_th=-0.3, Ki_th=-3.25, Kd_th=-0.01, tau_th=1/(2*pi*10), N_th=1/tau_th;
    double Kp_h=0.019, Ki_h=0.0002, Kd_h=0.01, tau_h=1/(2*pi*1), N_h=1/tau_h;
    double e_V, theta_ref, e_theta, u_V, e_h;
    double PID_th, PID_V, PID_h;

//V channel
    e_V=V_ref-V;
    PID(i,Kp_V,Ki_V,Kd_V,N_V,dt,e_V,I_V,D_V,e_Vp,&PID_V); //I puntatori devono essere definiti nella funzione principale di integrazione
    theta_ref=PID_V;

    e_theta=theta_ref-theta;
    PID(i,Kp_th,Ki_th,Kd_th,N_th,dt,e_theta,I_th,D_th,e_thp,&PID_th);
    u_V=PID_th;
    *de_c=u_V; //Commanded delta_elevator

//h channel
    e_h=-(h_ref-h);
    PID(i,Kp_h,Ki_h,Kd_h,N_h,dt,e_h,I_h,D_h,e_hp,&PID_h);
    *dth_c=PID_h; //Commanded throttle
    // printf("%lf  %lf\n", *e_hp, e_h);
}

void PID(int i, double Kp, double Ki, double Kd, double N, double dt, double e, double *I, double *D, double *ep, double *PID){
    double P;
    P=Kp*e;
    (*I)=(*I)+Ki*(*ep)*dt;
    (*D)=(*D-Kd*N*(*ep))-N*dt*(*D)+Kd*N*e;
    *ep=e;
    if (i==0) {
        *I = 0;
        *D = 0;
    }
    *PID=P+(*I)+(*D);
}

void saturation(double umax, double umin, double*u, int*sat_flag,int*flagp){
    *flagp=*sat_flag;
    if(*u>umax){
        *u=umax;
        *sat_flag=1;
    }else if(*u<umin){
        *u=umin;
        *sat_flag=-1;
    }else{
        *u=*u;
        *sat_flag=0;
    }
}

/*
double search_T(double dth, double*RPM, double*Torque, double*SPINTA){
    const double pi=4*atan(1);
    int N_RPM=260, i;
    double P_max=160; //[W]
    double PP,PS;
    double T;

    for(i=0;i<N_RPM;i++) {
        PS = Torque[i] * RPM[i] * (2 * pi) / 60;
        if(fabs(PS)>fabs(dth*P_max)){
            if (i==0){
                T = SPINTA[i];
            }else {
                PP = Torque[i - 1] * RPM[i - 1] * (2 * pi) / 60;
                T = SPINTA[i - 1] + ((SPINTA[i] - SPINTA[i - 1]) / (PS - PP)) * (fabs(dth) * P_max - PP);
            }
            return T;
        }
    }
}
*/

void ReadPsi(double psi[], double *tf, double *x, double *y) {

    FILE *fp;
    char line[50];
    int i = 0, flag;
    int p;

    printf("--------------------------------------------------------------------------------\n");
    printf("SCELTA DEL PERCORSO:\n");

    do {
        printf("\nScegliere un percorso tra: "
               "\n[1] Square"
               "\n[2] Butterfly"
               "\n[3] Snake\n");
        scanf("%d", &p);

        flag = 0;
        switch (p) {
            case 1:
                fp = fopen("Drone/PERCORSI/Square.txt", "r");
                *tf = 129.5;
                *y = *y - 300;
                break;
            case 2:
                fp = fopen("Drone/PERCORSI/Butterfly.txt", "r");
                *tf = 160.99;
                *x = *x + 250;
                *y = *y - 250;
                break;
            case 3:
                fp = fopen("Drone/PERCORSI/Snake.txt", "r");
                *tf = 325.99;
                *x = *x + 250;
                *y = *y - 250;
                break;
            default:
                printf("\nNumero inserito non corretto.");
                flag = 1;
        }
    } while (flag == 1);


    if (fp == NULL) {
        printf("Errore apertura file PsiReference.txt");
    }

    while (fgets(line, 50, fp) != NULL) {
        sscanf(line, "%*lf %lf", &psi[i]);
        i++;
    }

    fclose(fp);
}