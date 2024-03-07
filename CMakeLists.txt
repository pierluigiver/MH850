cmake_minimum_required(VERSION 3.17)
project(FlightSimulator C)

set(CMAKE_C_STANDARD 99)

add_executable(FlightSimulator
        Functions/Atmosfera/Atmosfera.c
        Functions/Atmosfera/Atmosfera.h
        Main.c
        Functions/ManageFile/ManageFile.h
        Functions/ManageFile/ManageFile.c
        "Functions/Condizioni_trim/Alpha_trim.c"
        "Functions/Condizioni_trim/Alpha_trim.h"
        "Functions/Calcolo_spinta_RPM/Calcolo_spinta_RPM.c"
        "Functions/Calcolo_spinta_RPM/Calcolo_spinta_RPM.h"
        "Functions/Condizioni_trim/Condizioni_trim.c"
        "Functions/Condizioni_trim/Condizioni_trim.h"
        "Functions/Routh/Routh.c"
        "Functions/Routh/Routh.h"
        Functions/Modi/Longitudinali.h
        Functions/Modi/Longitudinali.c
        Functions/Integrazione/Integrazione.c
        Functions/Integrazione/Integrazione.h
        Functions/Interpolazione/Interpolazione.c
        Functions/Interpolazione/Interpolazione.h)
