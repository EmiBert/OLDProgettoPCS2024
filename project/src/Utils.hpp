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



template<typename T>
void Scambia(vector<T>& A, int i, int j){
    T temp = A[j];
    A[j] = A[i];
    A[i] = temp;
}


template<typename T>
int Distribuzione(vector<T>& A, int sinistra, int destra){
    T x = A[destra];
    int i = sinistra-1;
    for(int j=sinistra; j<destra; j++){
        if(A[j] >= x){
            i++;
            Scambia(A,i,j);
        }
    }
    Scambia(A,i+1,destra);
    return i+1;
}



//    2 pre: 0≤sinistra,destra≤n−1
template<typename T>
void QuickSort(vector<T>& A, int sinistra, int destra){
    if (sinistra < destra){
        // il pivot è l'ultimo indice, "destra"
        int rango = Distribuzione(A, sinistra, destra);
        QuickSort(A, sinistra, rango-1);
        QuickSort(A, rango+1, destra);
    }
}


template<typename T>
void QuickSort(vector<T>& A){
    int sinistra = 0;
    int destra = A.size()-1;
    QuickSort(A, sinistra, destra);
}


void StampaTracceOrdinate(vector<Traces> tracesContainer,
                          int numFracture,
                          map<int, vector<int>>& sortedPassanti,
                          map<int, vector<int>>& sortedNonPassanti);


void IntersezioneEdges(edges L1,
                       edges L2,
                       double& alpha,
                       double& beta);


void TagliaFratture(Fracture F,
                    vector<Traces> contenitoreTracce,
                    map<int, vector<int>> tPO,
                    map<int, vector<int>> tNPO,
                    vector<edges>& latiBordo,
                    vector<edges>& latiInterni);


void CercaEstremo(Vector3d P,
                  int& id,
                  vector<Vector3d> CoordinateNodi,
                  vector<int> idNodi,
                  bool& flag);


void CercaEstremo(Vector3d P,
                  int& id,
                  vector<Vector3d> CoordinateNodi,
                  vector<int> idNodi);


void CaricamentoCell0e1D(vector<edges>& latiBordo,
                         vector<edges>& latiInterni,
                         PolygonalMesh& mesh,
                         vector<int>& idBordo,
                         vector<int>& idInterno);



Vector3d CalcoloNormaleMesh(PolygonalMesh mesh);

void LatoSuccessivo(Vector3d& CurrentEdgeTan,
                    int& CurrentNode,
                    int& CurrentEdgeId,
                    vector<int>& nodiPoly,
                    vector<int>& latiPoly,
                    int& inverti,
                    bool& chiuso,
                    PolygonalMesh mesh,
                    Vector3d N);


void CaricamentoCell2D (PolygonalMesh& mesh,
                       vector<int>& idBordo,
                       vector<int>& idInterno);



}
