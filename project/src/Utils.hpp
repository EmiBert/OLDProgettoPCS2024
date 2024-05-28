#pragma once

#include <iostream>
#include "DFNlibrary.hpp"

using namespace std;
namespace DFNlibrary{

bool ImportaFratture(const string& filepath,
                     vector<Fracture>& fract,
                     int& numF);



Vector3d CalcoloNormale(const Fracture F);



bool CercaIntersezioni(Vector3d P,
                       Vector3d t,
                       Fracture F,
                       double& c1,
                       double& c2);



void InserisciTraccia(double alpha,
                      double beta,
                      double gamma,
                      double delta,
                      vector<Traces>& tracesContainer,
                      Vector3d P,
                      Vector3d t,
                      int Fid1,
                      int Fid2);




void CercaTracce(const Fracture F1,
                 const Fracture F2,
                 vector<Traces>& tracesContainer);



void StampaTracce(vector<Traces> tracesContainer,
                  int numF);


void Scambia(vector<Vector2d>& A,
             int i,
             int j);



int Distribuzione(vector<Vector2d>& A,
                  int sinistra,
                  int destra);


void QuickSort(vector<Vector2d>& A,
               int sinistra,
               int destra);



void QuickSort(vector<Vector2d>& A);


void StampaTracceOrdinate(vector<Traces> tracesContainer,
                          int numFracture);



}
