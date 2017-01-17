#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>

#include <survey.h>


#include "util.h"
// Read ASCII file containing geophone information, i.e., ID, X, Y, Z
bool readGeophones(const char* fileName, int num, vector<double> &rx, vector<double> &ry, vector<double> &rz);


// Read ASCII file containing velocity information
bool readVelModel(const char* fileName, vector<double> &depthTreatWell, vector<double> &vpTreatWell, vector<double> &vsTreatWell, vector<double> &depthMonitorWell, vector<double> &vpMonitorWell, vector<double> &vsMonitorWell);


// Read ASCII file containing monitoring well trajectory information
bool readWellDev(const char* fileName, vector<WellDevPara*> &devData);


#endif // IO_H
