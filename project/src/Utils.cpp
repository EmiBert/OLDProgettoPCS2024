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
        double valCond = abs(t[0]*k[0] + t[1]*k[1] + t[2]*k[2]) - t.norm()*k.norm();
        if(abs(valCond) < tau) // k e t sono paralleli, non puo' esserci intersezione
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


void StampaTracce(vector<Traces> tracesContainer, int numF){
    // fileName = "traces_FRX_data.csv", X = numero fratture nel file considerato
    string fileName = "traces_FR"+to_string(numF)+"_data.csv";
    cout<<fileName<<endl;
    ofstream outputTraces (fileName);
    outputTraces << "# Number of Traces"<<endl;    //intestazione
    outputTraces << tracesContainer.size()<<endl;
    outputTraces << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2"<<endl;

    for(auto& t: tracesContainer){
        outputTraces<<t.id<<";"<<t.FractureID1<<";"<<t.FractureID2<<";";
        outputTraces<<t.P1[0]<<";"<<t.P1[1]<<";"<<t.P1[2]<<";"<<t.P2[0]<<";"<<t.P2[1]<<";"<<t.P2[2]<<endl;
    }

    outputTraces.close();
}




void Scambia(vector<Vector2d>& A, int i, int j){
    Vector2d temp = A[j];
    A[j] = A[i];
    A[i] = temp;
}


int Distribuzione(vector<Vector2d>& A, int sinistra, int destra){
    Vector2d x = A[destra];
    int i = sinistra-1;
    for(int j=sinistra; j<destra; j++){
        if(A[j][1] >=x[1]){
            i++;
            Scambia(A,i,j);
        }
    }
    Scambia(A,i+1,destra);
    return i+1;
}


//    2 pre: 0≤sinistra,destra≤n−1
void QuickSort(vector<Vector2d>& A, int sinistra, int destra){
    if (sinistra < destra){
        // il pivot è l'ultimo indice, "destra"
        int rango = Distribuzione(A, sinistra, destra);
        QuickSort(A, sinistra, rango-1);
        QuickSort(A, rango+1, destra);
    }
}


void QuickSort(vector<Vector2d>& A){
    int sinistra = 0;
    int destra = A.size()-1;
    QuickSort(A, sinistra, destra);
}




void StampaTracceOrdinate(vector<Traces> tracesContainer, int numFracture,
                          map<int, vector<int>>& sortedPassanti,
                          map<int, vector<int>>& sortedNonPassanti){

    // fileName = "sorted_traces_FRX_data.csv", X = numero fratture nel file considerato
    string fileName = "sorted_traces_FR"+to_string(numFracture)+"_data.csv";
    cout<<fileName<<endl;
    ofstream outputSortedTraces (fileName);

    for(int i=0; i<numFracture; i++){
        // contenitori di tracce in forma breve, Vector2d: [TraceId, Length]
        vector<Vector2d> tPassanti = {};
        vector<Vector2d> tNonPassanti = {};

        for(auto& t: tracesContainer){

            if(t.FractureID1 == i){
                Vector2d shortT = {t.id, t.Length};
                if(t.Tips1 == true){
                    tNonPassanti.push_back(shortT);
                }
                else{
                    tPassanti.push_back(shortT);
                }
            }

            if(t.FractureID2 == i){
                Vector2d shortT = {t.id, t.Length};
                if(t.Tips2 == true){
                    tNonPassanti.push_back(shortT);
                }
                else{
                    tPassanti.push_back(shortT);
                }
            }
        }

        if(tPassanti.size() != 0){
            QuickSort(tPassanti);
        }

        if(tNonPassanti.size() != 0){
            QuickSort(tNonPassanti);
        }

        vector<int> idTP = {}; // id tracce passanti
        vector<int> idTNP = {}; // id tracce non passanti


        if(tPassanti.size() + tNonPassanti.size()!= 0){
            outputSortedTraces << "# FractureId; NumTraces"<<endl;
            outputSortedTraces << i<<";"<< tPassanti.size() + tNonPassanti.size()<<endl;
            outputSortedTraces << "# TraceId; Tips; Length"<<endl;

            if(tNonPassanti.size() !=0){
                for(auto& t: tNonPassanti){
                    outputSortedTraces<<t[0]<<";1;"<<t[1]<<endl;
                    idTNP.push_back(t[0]);
                }
                sortedNonPassanti.insert({i, idTNP});
            }

            if(tPassanti.size() !=0){
                for(auto& t: tPassanti){
                    outputSortedTraces<<t[0]<<";0;"<<t[1]<<endl;
                    idTP.push_back(t[0]);
                }
                sortedPassanti.insert({i, idTP});
            }
        }
    }

    outputSortedTraces.close();
}



void IntersezioneEdges(edges L1, edges L2, double& alpha, double& beta){

    double valCond = abs(L1.t[0]*L2.t[0] + L1.t[1]*L2.t[1] + L1.t[2]*L2.t[2]) - L1.t.norm()*L2.t.norm();
    if(abs(valCond)<1e-6) // k e t sono paralleli, non puo' esserci intersezione
    {
        return;
    }
    Vector3d b = L2.P - L1.P;
    MatrixXd A(3, 2);
    A.col(0) = L1.t;
    A.col(1) = L2.t;
    Vector2d sol = A.fullPivLu().solve(b);
    alpha = sol[0];
    beta = -sol[1];
}


void TagliaFratture(Fracture F, vector<Traces> contenitoreTracce,
                    map<int, vector<int>> tPO, map<int, vector<int>> tNPO,
                    vector<edges>& latiBordo, vector<edges>& latiInterni){
    // tPO = tracce Passanti Ordinate, tNPO = tracce Non-Passanti Ordinate
    for(int i=0; i<F.NumVertices; i++){
        Vector3d P = {F.coordx[i], F.coordy[i], F.coordz[i]};
        Vector3d t = {F.coordx[(i+1)%F.NumVertices]-F.coordx[i], F.coordy[(i+1)%F.NumVertices]-F.coordy[i], F.coordz[(i+1)%F.NumVertices]-F.coordz[i]};
        edges e = {};
        e.P = P;
        e.t = t;
        latiBordo.push_back(e);
    }

    if(tPO.find(F.id) != tPO.end()){
        for(int i=0; i<tPO[F.id].size(); i++){
            edges l = {};
            int idTraccia = tPO[F.id][i];
            for(auto& T: contenitoreTracce){
                if(T.id == idTraccia){
                    l.P = T.P1;
                    l.t = T.P2 - T.P1;
                    break;
                }
            }

            bool flagStart = false;
            bool flagEnd = false;

            for(auto& b: latiBordo){
                double alpha = 0;
                double beta = 0;
                IntersezioneEdges(b,l,alpha,beta);
                if(flagStart && flagEnd){
                    break;
                }
                if(abs(beta)<tau){
                    b.intersection.push_back(alpha);
                    flagStart = true;
                }
                else if(abs(beta-1)<tau){
                    b.intersection.push_back(alpha);
                    flagEnd = true;
                }
            }


            for(auto& i: latiInterni){
                double alpha = 0;
                double beta = 0;
                IntersezioneEdges(i,l,alpha,beta);
                if(0<alpha && alpha<1 && 0<beta && beta<1){
                    i.intersection.push_back(alpha);
                    l.intersection.push_back(beta);
                }
            }

            latiInterni.push_back(l);
        }
    }


}



}
