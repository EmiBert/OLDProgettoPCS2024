#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

namespace DFNlibrary {
const double tau = 1e-6;

bool ImportaFratture(const string& filepath,
                     vector<Fracture>& fract,
                     int& numF){

    ifstream file;
    file.open(filepath);

    if(file.fail())
        return false;

    vector<string> inputLines;
    string line;
    while (getline(file, line))
        inputLines.push_back(line);

    file.close();


    string s = inputLines[1];

    numF = stoi(s);

    if (numF == 0)
    {
        cerr << "non ci sono fratture" << endl;
        return false;
    }

    fract.reserve(numF);

    int id;
    int numV;

    for (int i=0; i<numF; i++){

        Fracture F;
        int pos = i*6+3;
        istringstream converter(inputLines[pos]);
        getline(converter, s, ';');
        id = stoi(s);
        F.id = id;

        getline(converter, s);
        numV = stoi(s);
        F.NumVertices = numV;

        F.coordx.reserve(numV);
        F.coordy.reserve(numV);
        F.coordz.reserve(numV);

        pos+=2;
        istringstream converterX(inputLines[pos]);

        for(int k =0; k<numV; k++){
            if(k == numV-1){
                getline(converterX, s);
                F.coordx.push_back(stod(s));
            }
            else{
                getline(converterX, s, ';');
                F.coordx.push_back(stod(s));
            }
        }

        pos++;
        istringstream converterY(inputLines[pos]);

        for(int k =0; k<numV; k++){
            if(k == numV-1){
                getline(converterY, s);
                F.coordy.push_back(stod(s));
            }
            else{
                getline(converterY, s, ';');
                F.coordy.push_back(stod(s));
            }
        }

        pos++;
        istringstream converterZ(inputLines[pos]);

        for(int k =0; k<numV; k++){
            if(k == numV-1){
                getline(converterZ, s);
                F.coordz.push_back(stod(s));
            }
            else{
                getline(converterZ, s, ';');
                F.coordz.push_back(stod(s));
            }
        }

        fract.push_back(F);

    }

    file.close();
    return true;
}

Vector3d CalcoloNormale(const Fracture F){

    // calcolo del vettore ortonormale N al piano contenente F

    Vector3d v; // vettore (F vertice 1 - F vertice 0)
    v[0] = F.coordx[1] - F.coordx[0]; // coordinata x di v
    v[1] = F.coordy[1] - F.coordy[0]; // coordinata y di v
    v[2] = F.coordz[1] - F.coordz[0]; // coordinata z di v

    Vector3d u; // vettore (F vertice 2 - F vertice 0)
    u[0] = F.coordx[2] - F.coordx[0]; // coordinata x di u
    u[1] = F.coordy[2] - F.coordy[0]; // coordinata y di u
    u[2] = F.coordz[2] - F.coordz[0]; // coordinata z di u

    Vector3d N; // v x u (prodotto scalare)
    N[0] = v[1]*u[2] - v[2]*u[1];
    N[1] = v[2]*u[0] - v[0]*u[2];
    N[2] = v[0]*u[1] - u[0]*v[1];

    // normalizzazione
    N = N/ N.norm();

    return N;
}




bool CercaIntersezioni(Vector3d P, Vector3d t, Fracture F, double& c1, double& c2){
    bool interseca = false;
    for(int i = 0; i< F.NumVertices; i++){
        Vector3d k = {F.coordx[(i+1)%F.NumVertices]-F.coordx[i], F.coordy[(i+1)%F.NumVertices]-F.coordy[i], F.coordz[(i+1)%F.NumVertices]-F.coordz[i]};
        // il modulo serve ad associare all'ultima iterata i vertici con indice pari a F.NumVertices-1 e 0 (F.NumVertices % F.NumVertices)
        // ovvero per un quadrilatero il vertice 3 e il vertice 0
        if((k[0]/t[0] == k[1]/t[1]) && (k[2]/t[2] == k[1]/t[1])) // k e t sono paralleli, non puo' esserci intersezione
        {
            continue;
        }

        Vector3d b = {F.coordx[i]-P[0], F.coordy[i]-P[1], F.coordz[i]-P[2]};
        MatrixXd A(3, 2);
        A.col(0) = t;
        A.col(1) = k;
        Vector2d sol = A.fullPivLu().solve(b);
        if((0< -sol(1)) && (-sol(1)<1)){
            if (interseca){
                c2 = sol(0);
                break;
            }
            else{
                interseca = true;
                c1 = sol(0);
            }
        }
    }

    return interseca;
}



void InserisciTraccia(double alpha, double beta, double gamma, double delta,
                      vector<Traces>& tracesContainer, Vector3d P, Vector3d t, int Fid1, int Fid2){
    double max1 = max(alpha,beta);
    double min1 = min(alpha,beta);
    double max2 = max(gamma,delta);
    double min2 = min(gamma,delta);

    if( (min1>max2) || (min2>max1) ){
        return; // le due fratture non si intersecano
    }
    Vector3d estremo1 = {};
    Vector3d estremo2 = {};
    bool T1;
    bool T2;

    if( (max1>=max2) && (min1<=min2)){
        estremo1 = P + max2*t;
        estremo2 = P + min2*t;
        T1 = true;
        T2 = false;
        if( (max1 == max2) && (min1 == min2) )
            T2 = false;
    }

    else if( (max2>=max1) && (min2<=min1)){
        estremo1 = P + max1*t;
        estremo2 = P + min1*t;
        T1 = false;
        T2 = true;
        if( (max1 == max2) && (min1 == min2) )
            T2 = false;
    }

    else if( max1>max2 ){
        estremo1 = P + max2*t;
        estremo2 = P + min1*t;
        T1 = true;
        T2 = true;
    }

    else{
        estremo1 = P + max1*t;
        estremo2 = P + min2*t;
        T1 = true;
        T2 = true;
    }

    Traces T = {};
    T.id = tracesContainer.size();
    T.FractureID1 = Fid1;
    T.FractureID2 = Fid2;
    T.P1 = estremo1;
    T.P2 = estremo2;
    T.Tips1 = T1;
    T.Tips2 = T2;
    T.Length = (estremo1-estremo2).norm();
    tracesContainer.push_back(T);
}














void CercaTracce(const Fracture F1, const Fracture F2, vector<Traces>& tracesContainer){

    // calcolo dei vettori ortonormali ai piani F1, F2
    Vector3d N1 = CalcoloNormale(F1);
    Vector3d N2 = CalcoloNormale(F2);

    // calcolo direzione t della retta intersezione dei due piani,
    // prodotto vettoriale tra N1 e N2
    Vector3d t;
    t[0] = N1[1]*N2[2] - N1[2]*N2[1];
    t[1] = N1[2]*N2[0] - N1[0]*N2[2];
    t[2] = N1[0]*N2[1] - N2[0]*N1[1];

    if (t.norm() < tau) // condizione: t == 0, in algebra finita con tolleranza tau
    {
        return; // i due piani sono paralleli, non possono esserci tracce
    }

    // prodotto scalare tra N1 e F1 vertice 0
    double d1 = N1[0]*F1.coordx[0] + N1[1]*F1.coordy[0] + N1[2]*F1.coordz[0];

    // prodotto scalare tra N2 e F2 vertice 0
    double d2 = N2[0]*F2.coordx[0] + N2[1]*F2.coordy[0] + N2[2]*F2.coordz[0];

    // creazione matrice A = {N1, N2, t}
    Matrix3d A = {};
    A.row(0) = N1;
    A.row(1) = N2;
    A.row(2) = t;

    // creazione vettore termine noto b = {d1, d2, 0}'
    Vector3d b = {d1, d2, 0};

    // risoluzione del sistema tramite la decomposizione PALU,
    // P punto appartente alla retta di intersezione dei due piani
    Vector3d P = A.fullPivLu().solve(b);

    double alpha = 0;
    double beta = 0;

    if(!CercaIntersezioni(P,t,F1,alpha,beta))
        return;


    double gamma = 0;
    double delta = 0;

    if(!CercaIntersezioni(P,t,F2,gamma,delta))
        return;


    InserisciTraccia(alpha, beta, gamma, delta, tracesContainer, P, t, F1.id, F2.id);

}





}
