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


int main(int nargs, char** argv){
    std::string protonfilename = "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/DLwork/DataFiles/larlite_pandoraNu_20170211_205524_331124.root";
    std::string muonfilename   = "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/mylarlite/UserDev/Chimera_Adrien/ChimeraMuonTracks/data/larlite_kalmanhit_0000.root";
    std::string targetfilename = "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/mylarlite/UserDev/Adrien/Chimera/Targets.csv";


    //=====================
    //=====================
    // Read Proton File ===
    //=====================
    //=====================

    std::cout << "READ PROTON FILE" << std::endl;
    larlite::storage_manager p_storage;
    p_storage.set_verbosity(larlite::msg::kNORMAL);
    p_storage.set_io_mode(larlite::storage_manager::kREAD);
    p_storage.add_in_filename(protonfilename);
    p_storage.open();

    if(!p_storage.is_open()){
        std::cerr << "File open failed!" << std::endl;
        return 0;
    }
    if(!p_storage.is_ready_io()){
        std::cerr << "I/O preparation failed!" << std::endl;
    }

    ChimeraFinding findProton;
    findProton.SetTargetFile(targetfilename);
    findProton.SetTrackGenerator("pandoraNu");
    findProton.SetParticleType("proton");
    findProton.initialize();

    while(p_storage.next_event()){
        findProton.analyze(p_storage);
    }

    std::vector<std::vector<larlite::hit> > BestHitClusters_p = findProton.GetBestHitClusters();
    std::vector<larlite::track> BestTracks_p = findProton.GetBestTracks();

    findProton.finalize();

    //=====================
    //=====================
    // Read Muon File =====
    //=====================
    //=====================

    std::cout << "READ MUON FILE" << std::endl;
    larlite::storage_manager u_storage;
    u_storage.set_verbosity(larlite::msg::kNORMAL);
    u_storage.set_io_mode(larlite::storage_manager::kREAD);
    u_storage.add_in_filename(muonfilename);
    u_storage.open();

    if(!u_storage.is_open()){
        std::cerr << "File open failed!" << std::endl;
        return 0;
    }
    if(!u_storage.is_ready_io()){
        std::cerr << "I/O preparation failed!" << std::endl;
    }

    ChimeraFinding findMuon;
    findMuon.SetTargetFile(targetfilename);
    findMuon.SetTrackGenerator("trackkalmanhit");
    findMuon.SetParticleType("muon");
    findMuon.initialize();
    double parSigmasMuon[3] =  {50,1.,90*TMath::Pi()/180.}; //σ_R, σ_L,  σ_angles
    findMuon.SetSigmaEval(parSigmasMuon);

    while(u_storage.next_event()){
        findMuon.analyze(u_storage);
    }

    std::vector<std::vector<larlite::hit> > BestHitClusters_u = findMuon.GetBestHitClusters();
    std::vector<larlite::track> BestTracks_u = findMuon.GetBestTracks();

    findMuon.finalize();

    //=====================
    //=====================
    // Print Chimera Evt ==
    //=====================
    //=====================

    ChimeraPatching patch;
    std::vector<std::vector<double> > ParTarget = findMuon.GetParTarget();
    patch.Initialize();
    std::string evtID;
    if(BestHitClusters_p.size()!=BestHitClusters_u.size()
       || BestHitClusters_p.size()!=BestTracks_p.size()
       || BestHitClusters_p.size()!=BestTracks_u.size()){
        std::cout << "ERROR : not the same number of muon and proton tracks and hit clusters" << std::endl;
        return 0;
    }

    for(int iEvt = 0;iEvt<BestHitClusters_p.size();iEvt++){
        TVector3 X0(ParTarget.at(iEvt)[0],ParTarget.at(iEvt)[1],ParTarget.at(iEvt)[2]);
        evtID = Form("1_0_%03d",iEvt);
        patch.NewEvent(evtID, X0);
        patch.AddTrack(BestTracks_p[iEvt],BestHitClusters_p[iEvt]);
        patch.AddTrack(BestTracks_u[iEvt],BestHitClusters_u[iEvt]);
        patch.TranslateClusters();
        patch.DrawEvent();
    }
    return 0;
}
