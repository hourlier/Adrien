// compile with g++ GenerateChimeraTargetFile.cc -o newTargetFile
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <string>
#include <time.h>
#include <random>


int main(int argc, char **argv){
    std::string filename = "Target.csv";
    int event = 0;
    int subrun = 1;
    int run = 1;


    int Nevt  = 10;
    int Npart = 2;

    if(argc >= 2){Nevt  = atoi(argv[1]);}
    if(argc >= 3){Npart = atoi(argv[2]);}

    std::pair<int,int> Xdet[3];
    Xdet[0].first = 0;
    Xdet[0].second = 250;
    Xdet[1].first = -100;
    Xdet[1].second = 100;
    Xdet[2].first = 0;
    Xdet[2].second = 1000;
    double X0[3]; // vertex position
    std::vector<double> L(Npart); // length of the various particles.
    std::vector<double> theta(Npart);
    std::vector<double> phi(Npart);
    int Run(run), SubRun(SubRun),Event(event);
    int eventnumberrange = 1e4*Nevt;


    double pi = acos(-1.);
    srand(time(NULL));
    std::default_random_engine generator;
    std::normal_distribution<double> distribLength(10.0,5.0);
    std::normal_distribution<double> distribTheta(pi/4.,pi/8.);
    std::normal_distribution<double> distribPhi(pi/2.,pi/4.);

    std::ofstream targetfile(filename);
    if(!targetfile){std::cout << "ERROR, could not open the output file" << std::endl;return 0;}
    if(Event!=1)targetfile << std::endl;
    for(int ievt = 0;ievt<Nevt;ievt++){
        Event++;
        if(Event == 10000){
            SubRun++;
            Event = 0;
        }
        if(SubRun == 100){
            Run++;
            SubRun=0;
            Event=0;
        }
        targetfile << Run << ",  " << SubRun << ",  " << Event;
        for(int dim=0;dim<3;dim++){
            X0[dim] = rand()%(Xdet[dim].second-Xdet[dim].first)-Xdet[dim].first;
            targetfile <<",  "<< X0[dim];
        }
        for(int ipart = 0;ipart<Npart;ipart++){
            L[ipart] = distribLength(generator);
            while (L[ipart]<6) {L[ipart] = distribLength(generator);}

            theta[ipart] = distribTheta(generator);
            while(theta[ipart] > pi/2.){theta[ipart]-=pi/2.;}
            if(theta[ipart]<0)theta[ipart]*=-1.;

            phi[ipart] = distribPhi(generator);
            if(phi[ipart]<0)phi[ipart]*=-1.;
            while(phi[ipart] > pi){phi[ipart]-=pi;}
            phi[ipart]-=pi/2.;

            targetfile << ",  " << L[ipart] << ",  " << theta[ipart] << ",  " << phi[ipart];
        }
        targetfile << " z";
        if(ievt < Nevt-1)targetfile << std::endl;
    }

    targetfile.close();
    return 0;
}
