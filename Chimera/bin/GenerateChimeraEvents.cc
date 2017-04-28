#include <iostream>
#include <string>
#include <vector>

#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/storage_manager.h"

#include "../ChimeraTrackEvaluator.h"
#include "../ChimeraFinding.h"
#include "../ChimeraPatching.h"

#include "TVector3.h"

std::vector<std::vector<double> > ReadTargetFile(std::string _TargetFile, int Nevt);

int main(int nargs, char** argv){
    int Nevt = -1;
    if(nargs > 1){
        Nevt = atoi(argv[1]);
    }

    std::string basedir = "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE";

    std::string protonfilename = basedir+"/DeepLearning/DLwork/DataFiles/larlite_pandoraNu_20170211_205524_331124.root";
    std::string muonfilename   = basedir+"/mylarlite/UserDev/Chimera_Adrien/ChimeraMuonTracks/data/larlite_kalmanhit_0000.root";
    //std::string muonfilename = basedir+"/DeepLearning/DLwork/DataFiles/larlite_pandoraNu_20170211_205524_331124.root";

    std::string muonTrackFile  = basedir+"/mylarlite/UserDev/Chimera_Adrien/ChimeraMuonTracks/data/muon_track_list.csv";
    std::string protonTrackFile= basedir+"/DeepLearning/DLwork/DataFiles/passedGBDT_extBNB_AnalysisTrees_cosmic_trained_only_on_mc_score_0.99.csv";

    std::string targetfilename = basedir+"/mylarlite/UserDev/Adrien/Chimera/bin/Target.csv";

    int NparticlesInEvt = 2;

    std::vector<std::vector<double> > Scores;
    std::vector<std::vector<std::vector<larlite::hit> > > BestHitClusters;
    std::vector<std::vector<larlite::track> > BestTracks;

    std::vector< larlite::storage_manager > storage(2);
    std::vector< ChimeraFinding > findTrack(2);

    std::vector<std::string> larliteInputFile(2);
    larliteInputFile[0] = protonfilename;
    larliteInputFile[1] = muonfilename;

    std::vector< std::string > TrackFile(2);
    TrackFile[0] = protonTrackFile;
    TrackFile[1] = muonTrackFile;

    std::vector< std::string > particleType(2);
    particleType[0] = "proton";
    particleType[1] = "muon";

    std::vector< std::string > trackGenerator(2);
    trackGenerator[0] = "pandoraNu";
    trackGenerator[1] = "trackkalmanhit";
    //trackGenerator[1] = "pandoraNu";

    if(NparticlesInEvt > larliteInputFile.size() || NparticlesInEvt > TrackFile.size()){
        std::cout << "ERROR : not enough input files" << std::endl;
        return 0;
    }

    std::vector<std::vector<double> > parSigmas(2);
    std::vector<double> parSigmasProton(3); parSigmasProton[0] = 10;parSigmasProton[1] = 0.1; parSigmasProton[2] = 10*TMath::Pi()/180.; //σ_R, σ_L,  σ_angles
    std::vector<double> parSigmasMuon(3);   parSigmasMuon[0] =100;   parSigmasMuon[1] = 2.;    parSigmasMuon[2] = 90*TMath::Pi()/180.; //σ_R, σ_L,  σ_angles
    parSigmas[0] = parSigmasProton;
    parSigmas[1] = parSigmasMuon;

    std::vector<std::vector<double> > FullTargetParameters = ReadTargetFile(targetfilename, Nevt);
    if(FullTargetParameters.size() == 0){std::cout << "ERROR, empty target vector" << std::endl;return 0;}

    //=====================
    //=====================
    // Read Data File =====
    //=====================
    //=====================

    for(int ipart = 0;ipart<NparticlesInEvt;ipart++){

        std::cout << "Read " << particleType[ipart] << " file" << std::endl;
        storage[ipart].set_verbosity(larlite::msg::kNORMAL);
        storage[ipart].set_io_mode(larlite::storage_manager::kREAD);
        storage[ipart].add_in_filename(larliteInputFile[ipart]);
        storage[ipart].open();

        if(!storage[ipart].is_open()){
            std::cerr << "File open failed!" << std::endl;
            return 0;
        }
        if(!storage[ipart].is_ready_io()){
            std::cerr << "I/O preparation failed!" << std::endl;
        }

        findTrack[ipart].SetTargetFile(targetfilename);
        findTrack[ipart].SetTargetVector(FullTargetParameters);
        findTrack[ipart].SetTrackFile(TrackFile[ipart]);
        findTrack[ipart].SetTrackGenerator(trackGenerator[ipart]);
        findTrack[ipart].SetParticleType(particleType[ipart]);
        findTrack[ipart].initialize();
        double sigmaParameters[3] = {parSigmas[ipart][0],parSigmas[ipart][1],parSigmas[ipart][2]};
        findTrack[ipart].SetSigmaEval(sigmaParameters);
        int ievent = 0;
        int Ntot = storage[ipart].get_entries();
        while(storage[ipart].next_event()){
            if((int)(ievent*2000./Ntot)%(200) == 0){std::cout << ievent*100/Ntot << "\%" << std::endl;}
            ievent++;
            findTrack[ipart].analyze(storage[ipart]);
        } // Loop over events in the file

        Scores.push_back(findTrack[ipart].GetScores());
        BestHitClusters.push_back(findTrack[ipart].GetBestHitClusters());
        BestTracks.push_back(findTrack[ipart].GetBestTracks());
        
        findTrack[ipart].finalize();
    }

    //=====================
    //=====================
    // Print Chimera Evt ==
    //=====================
    //=====================
    if(BestHitClusters.size() == 0){std::cout << "ERROR : not Hit cluster found" << std::endl;return 0;}
    if(BestTracks.size() == 0){std::cout << "ERROR : not Tracks found" << std::endl;return 0;}

    ChimeraPatching patch;
    int Npart = Scores.size();
    patch.Initialize();
    patch.SetNpartperevent(NparticlesInEvt);
    patch.SetScoreLimit(1e-8);
    std::string evtID;
    for(int ipart = 0;ipart<Npart-1;ipart++){
        if(BestHitClusters[ipart].size()!=BestHitClusters[ipart+1].size()
           || BestHitClusters[ipart].size()!=BestTracks[ipart].size()
           || BestHitClusters[ipart].size()!=BestTracks[ipart+1].size()){
            std::cout << "ERROR : not the same number of tracks and hit clusters" << std::endl;
            return 0;
        }
    }

    for(int iEvt = 0;iEvt<BestHitClusters[0].size();iEvt++){
        //std::cout << "Drawing event #" << iEvt << std::endl;
        //std::cout << "X0 set" << std::endl;
        evtID = Form("1_0_%03d",iEvt);
        patch.NewEvent(evtID,FullTargetParameters.at(iEvt));
        for(int ipart=0;ipart<Npart;ipart++){
            patch.AddTrack(BestTracks[ipart][iEvt],BestHitClusters[ipart][iEvt],Scores[ipart][iEvt]);
        }
        patch.TranslateClusters();
        patch.DrawEvent();
    }
    return 0;
}

std::vector<std::vector<double> > ReadTargetFile(std::string _TargetFile, int Nevt){
    if(Nevt == -1)Nevt = 1e9;
    std::vector< std::vector<double> >FullTargetParameters;
    std::ifstream file(_TargetFile);
    if(!file){std::cout << "ERROR, could not open file" << std::endl;return FullTargetParameters;}
    int Run, SubRun, Event;
    double X0,Y0,Z0,L_p,theta_p,phi_p,L_u,theta_u,phi_u;
    char coma;
    std::vector<double> eventTarget(9);
    bool goOn = true;
    int nevtfound = 0;
    while(goOn && nevtfound < Nevt){
        file >> Run >> coma >> SubRun >> coma >> Event >> coma >> X0 >> coma >> Y0 >> coma >> Z0 >> coma >> L_p >> coma >> theta_p >> coma >> phi_p >> coma >> L_u >> coma >> theta_u >> coma >> phi_u;
        eventTarget[0] = X0;
        eventTarget[1] = Y0;
        eventTarget[2] = Z0;
        eventTarget[3] = L_p;
        eventTarget[4] = theta_p;
        eventTarget[5] = phi_p;
        eventTarget[6] = L_u;
        eventTarget[7] = theta_u;
        eventTarget[8] = phi_u;
        FullTargetParameters.push_back(eventTarget);
        nevtfound++;
        /*std::cout << Run << "\t" << SubRun << "\t" << Event << "\t";
         for(auto ipar:eventTarget){std::cout << ipar << "\t";}
         std::cout << std::endl;*/
        if(file.eof()){goOn=false;break;}
    }
    file.close();
    return FullTargetParameters;
}
