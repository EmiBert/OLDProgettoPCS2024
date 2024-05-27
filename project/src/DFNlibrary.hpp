#pragma once

#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace DFNlibrary {

    //int NumFracture = 0; ///< number of fracture
struct Fracture
{
    int id = 0; ///< numero identificativo della frattura
    int NumVertices = 0; ///< NumVertices = numero di vertici
    vector<double> coordx = {}; ///< coordx[i] = coordinata x del vertice i
    vector<double> coordy = {}; ///< coordy[i] = coordinata y del vertice i
    vector<double> coordz = {}; ///< coordz[i] = coordinata z del vertice i
};




struct Traces
{   int id = 0; ///< numero identificativo della traccia
    int FractureID1 = {}; ///< numero identificativo della frattura 1
    int FractureID2 = {}; ///< numero identificativo della frattura 2
    Vector3d P1 = {}; ///< P1 = estremo 1
    Vector3d P2 = {}; ///< P2 = estremo 2
    bool Tips1 = {}; ///< Tips1 falso se la traccia è passante rispetto alla frattura ID1, vero se è non-passante
    bool Tips2 = {}; ///< Tips2 falso se la traccia è passante rispetto alla frattura ID2, vero se è non-passante
    double Length = {}; ///< lunghezza della traccia
};




}
