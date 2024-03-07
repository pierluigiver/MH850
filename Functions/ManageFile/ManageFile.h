#ifndef FLIGHTSIMULATOR_MANAGEFILE_H
#define FLIGHTSIMULATOR_MANAGEFILE_H

void WorkingDirectory(const char* projectPath);

void ReadAndSaveStatic (char *file_engine, char *file_battery, double *n_min, double *n_max, double *v0,
                        double *a0, double *a_stall, double *q_stall, double *C_max, double *eta_s);

void ReadAndSaveDba (double *h, double *alpha, double *CX, double *CY, double *CZ, double *Cl, double *Cm, double *Cn,
              double *CXA, double *CXAP, double *CXU, double *CXQ, double *CXB, double *CXP, double *CXR,
              double *CYB, double *CYBP, double *CYP, double *CYR, double *CYA, double *CYQ,
              double *CZALPHA, double *CZAP, double *CZU, double *CZQ, double *CZB, double *CZP, double *CZR,
              double *ClB, double *ClBP, double *ClP, double *ClR, double *ClA, double *ClQ,
              double *CmA, double *CmAP, double *CmU, double *CmQ, double *CmB, double *CmP, double *CmR,
              double *CnB, double *CnBP, double *CnP, double *CnR, double *CnA, double *CnQ,
              double *CXde, double *CXdle, double *CZde, double *CZdle, double *CYda, double *CYdr,
              double *Clda, double *Cldr, double *Cmde, double *Cmdle, double *Cnda, double *Cndr,
              double *CXom, double *CYom, double *CZom, double *Clom, double *Cmom, double *Cnom);

void ReadAndSavePropeller(char *file_propeller, double *datipropeller, double *CSI, double *RD, double *CH_AD, double *BA);

double dbainterp(double Cp, double Cs, double *h);

void DbaRead(double h, double*ALPHA, double css[][6],double cxd[][7], double cyd[][6], double czd[][7],
             double cld[][6], double cmd[][7],double cnd[][6],double cuf[][6],double cum[][6],double crd[][6] );

double C_h_interp(double Cp, double Cs, double hp, double hs, double h);

#endif //FLIGHTSIMULATOR_MANAGEFILE_H
