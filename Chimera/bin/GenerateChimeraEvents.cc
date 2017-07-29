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
void EndOfCode();

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

    std::vector<std::vector<double> > Scores;
    std::vector<std::vector<std::vector<larlite::hit> > > BestHitClusters;
    std::vector<std::vector<larlite::track> > BestTracks;


    int NparticlesInEvt = 2;
    if(nargs > 2){
        NparticlesInEvt = nargs - 2;
    }

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

    std::vector< larlite::storage_manager > storage(NparticlesInEvt);
    std::vector< ChimeraFinding > findTrack(NparticlesInEvt);
    std::vector< std::vector<std::string> > TypeParticle(NparticlesInEvt);
    std::vector<std::vector<double> > parSigmas(NparticlesInEvt);

    std::vector<double> parSigmasProton(3); parSigmasProton[0] = 10;parSigmasProton[1] = 0.1; parSigmasProton[2] = 10*TMath::Pi()/180.; //σ_R, σ_L,  σ_angles
    std::vector<double> parSigmasMuon(3);   parSigmasMuon[0] =1000;   parSigmasMuon[1] = 2.;    parSigmasMuon[2] = 10*TMath::Pi()/180.; //σ_R, σ_L,  σ_angles

    for(int ipart = 0;ipart<NparticlesInEvt;ipart++){
        int partIndex;
        if((std::string)argv[ipart+2] == "proton"){   partIndex = 0; parSigmas[ipart] = parSigmasProton;}
        else if((std::string)argv[ipart+2] == "muon"){partIndex = 1; parSigmas[ipart] = parSigmasMuon;}
        else{std::cout << "ERROR, not recognized particle type : " << argv[ipart+2] << std::endl;return 1;}
        TypeParticle[ipart].resize(4);
        TypeParticle[ipart][0] = argv[ipart+2];               // particle type
        TypeParticle[ipart][1] = trackGenerator[partIndex];   // track generator
        TypeParticle[ipart][2] = larliteInputFile[partIndex]; // input file
        TypeParticle[ipart][3] = TrackFile[partIndex];        // track selection file
    }
    std::cout << Nevt << " events" << std::endl;
    for(int ipart = 0;ipart<NparticlesInEvt;ipart++){
        std::cout << TypeParticle[ipart][0] << std::endl;
    }

    //if(NparticlesInEvt > larliteInputFile.size()){std::cout << "ERROR : not enough input files" << std::endl;return 2;}
    //if(NparticlesInEvt > TrackFile.size()){std::cout << "ERROR : not enough input files" << std::endl;return 2;}

    //=====================
    //=====================
    // Read Track File ====
    //=====================
    //=====================

    std::vector<std::vector<double> > FullTargetParameters = ReadTargetFile(targetfilename, Nevt);
    if(FullTargetParameters.size() == 0){std::cout << "ERROR, empty target vector" << std::endl;return 3;}
    if(FullTargetParameters[0].size() < 3+3*NparticlesInEvt){std::cout << "ERROR : Target file cannot handle " << NparticlesInEvt << " particle per events for now" << std::endl; return 4;}

    //=====================
    //=====================
    // Read Data File =====
    //=====================
    //=====================

    for(int ipart = 0;ipart<NparticlesInEvt;ipart++){

        std::cout << "Read " << TypeParticle[ipart][0] << " file" << std::endl;
        storage[ipart].set_verbosity(larlite::msg::kNORMAL);
        storage[ipart].set_io_mode(larlite::storage_manager::kREAD);

        storage[ipart].add_in_filename(TypeParticle[ipart][2]);
        storage[ipart].open();

        if(!storage[ipart].is_open()){
            std::cerr << "File open failed!" << std::endl;
            return 5;
        }
        if(!storage[ipart].is_ready_io()){
            std::cerr << "I/O preparation failed!" << std::endl;
        }

        findTrack[ipart].SetTargetFile(targetfilename);
        findTrack[ipart].SetTargetVector(FullTargetParameters);
        findTrack[ipart].SetTrackFile(TypeParticle[ipart][3]);
        findTrack[ipart].SetTrackGenerator(TypeParticle[ipart][1]);
        findTrack[ipart].SetParticleType(TypeParticle[ipart][0]);
        findTrack[ipart].SetPartIndex(ipart);
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
        storage[ipart].close();
    }

    std::cout << "Read all input files, found tracks, ready to plot" << std::endl;

    //=====================
    //=====================
    // Print Chimera Evt ==
    //=====================
    //=====================
    if(BestHitClusters.size() == 0){std::cout << "ERROR : not Hit cluster found" << std::endl;return 6;}
    if(BestTracks.size() == 0){std::cout << "ERROR : not Tracks found" << std::endl;return 7;}

    ChimeraPatching patch;
    int Npart = Scores.size();
    patch.Initialize();
    patch.SetNpartperevent(NparticlesInEvt);
    std::vector<std::string> partypes(NparticlesInEvt);
    for(int ipart = 0;ipart<NparticlesInEvt;ipart++){
        partypes[ipart] = TypeParticle[ipart][0];
    }
    patch.SetPartTypes(partypes);
    patch.SetScoreLimit(1e-12);
    std::string evtID;
    std::cout << "about to set tracks to plot" << std::endl;
    for(int ipart = 0;ipart<Npart-1;ipart++){
        if(BestHitClusters[ipart].size()!=BestHitClusters[ipart+1].size()
           || BestHitClusters[ipart].size()!=BestTracks[ipart].size()
           || BestHitClusters[ipart].size()!=BestTracks[ipart+1].size()){
            std::cout << "ERROR : not the same number of tracks and hit clusters" << std::endl;
            return 8;
        }
    }

    for(int iEvt = 0;iEvt<BestHitClusters[0].size();iEvt++){
        evtID = Form("1_0_%03d",iEvt);
        patch.NewEvent(evtID,FullTargetParameters.at(iEvt));
        for(int ipart=0;ipart<Npart;ipart++){
            patch.AddTrack(BestTracks[ipart][iEvt],BestHitClusters[ipart][iEvt],Scores[ipart][iEvt]);
        }
        patch.TranslateClusters();
        patch.DrawEvent();
    }

    //=====================
    //=====================
    // End of main    =====
    //=====================
    //=====================
    EndOfCode();
    return 0;
}

std::vector<std::vector<double> > ReadTargetFile(std::string _TargetFile, int Nevt){
    if(Nevt == -1)Nevt = 1e9;
    std::vector< std::vector<double> >FullTargetParameters;
    std::ifstream file(_TargetFile);
    if(!file){std::cout << "ERROR, could not open file" << std::endl;return FullTargetParameters;}
    int Run, SubRun, Event;
    double X0,Y0,Z0,L,theta,phi;
    char coma;
    std::vector<double> eventTarget;
    bool goOn = true;
    int nevtfound = 0;
    while(goOn && nevtfound < Nevt){
        eventTarget.clear();
        file >> Run >> coma >> SubRun >> coma >> Event >> coma >> X0 >> coma >> Y0 >> coma >> Z0 >> coma;
        eventTarget.push_back(X0);
        eventTarget.push_back(Y0);
        eventTarget.push_back(Z0);
        while(coma != 'z'){
            file >> L >> coma >> theta >> coma >> phi >> coma;
            eventTarget.push_back(L);
            eventTarget.push_back(theta);
            eventTarget.push_back(phi) ;
        }
        FullTargetParameters.push_back(eventTarget);
        nevtfound++;

        if(file.eof()){goOn=false;break;}
    }
    file.close();
    return FullTargetParameters;
}

void EndOfCode(){
    std::cout << "=================" << std::endl << std::endl;
    std::cout << "Known issues to be addressed : " << std::endl;
    std::cout << "\t 1) : Track reconstruction sometimes clearly doesn't fit the hit distribution, can I do something about it?" << std::endl;
    std::cout << "\t 2) : Charge deposition at vertex may not be the exact sum of ionization by each individual particles" << std::endl;
    std::cout << "\t 3) : How to save the event?" << std::endl;
    std::cout << "\t 4) : Need to use wire info instead of hits" << std::endl;
    std::cout << "\t 5) : Implement the union of the dead wires for all tracks" << std::endl;
    std::cout << std::endl << "=================" << std::endl;
    std::cout << "Gracefuly stopped" << std::endl;
    std::cout << "=================" << std::endl;
}
