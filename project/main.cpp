#include <iostream>
#include "Utils.hpp"
#include "DFNlibrary.hpp"
#include "Eigen/Eigen"


using namespace std;
using namespace Eigen;
using namespace DFNlibrary;

int main()
{
    vector<Fracture> contenitoreFratture;
    int numeroFratture = 0;
    string filepath = "DFN/FR3_data.txt";

    if(!ImportaFratture(filepath, contenitoreFratture, numeroFratture)){
        return 1;
    }

    vector<Traces> contenitoreTracce = {};

    for(int i =0; i<numeroFratture-1; i++){
        for(int j =i+1; j<numeroFratture; j++){
            CercaTracce(contenitoreFratture[i],contenitoreFratture[j],contenitoreTracce);
        }
    }


    cout<<" # Number of Fractures"<<endl;
    cout<<numeroFratture<<endl;

/*    for(auto f: contenitoreFratture){
        cout<<"# FractureId; NumVertices"<<endl;
        cout<<f.id<<"; "<<f.NumVertices<<endl;
        for(auto x: f.coordx){
            cout<<x<<"; ";
        }
        cout<<endl;
        for(auto y: f.coordy){
            cout<<y<<"; ";
        }
        cout<<endl;
        for(auto z: f.coordz){
            cout<<z<<"; ";
        }
        cout<<endl;
    }




    cout<<" # TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2; Tips1; Tips2; length"<<endl;
    for(auto t: contenitoreTracce){
        cout<< t.id<<";";
        cout<< t.FractureID1<<";";
        cout<< t.FractureID2<<";";
        cout<< t.P1[0]<<";";
        cout<< t.P1[1]<<";";
        cout<< t.P1[2]<<";";
        cout<< t.P2[0]<<";";
        cout<< t.P2[1]<<";";
        cout<< t.P2[2]<<";";
        cout<< t.Tips1<<";";
        cout<< t.Tips2<<";";
        cout<< t.Length<<";"<<endl;
    }

*/

    StampaTracce(contenitoreTracce,filepath);

return 0;
}
