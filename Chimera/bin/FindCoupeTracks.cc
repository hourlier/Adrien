#include <iostream>
#include <string>
#include <vector>

#include "DataFormat/opflash.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/storage_manager.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"

#include "../ChimeraPatching.h"
#include "../ChimeraTrackEvaluator.h"
#include "utils.cpp"

#include "TVector3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"


int main( int nargs, char** argv){

    std::string protonfilename = "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/DLwork/DataFiles/larlite_pandoraNu_20170211_205524_331124.root";
    std::string muonfilename   = "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/mylarlite/UserDev/Chimera_Adrien/ChimeraMuonTracks/data/larlite_kalmanhit_0000.root";
    std::string targetfilename = "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/mylarlite/UserDev/Adrien/Chimera/Targets.csv";
    std::string evtID;

    //=====================
    //=====================
    // Read .csv File =====
    //=====================
    //=====================
    std::cout << "reading .csv file" << std::endl;
    std::vector<std::vector<double> > FullTargets = ReadTargetFile(targetfilename);
    std::cout << ".csv file read OK " << FullTargets.size() << " events requested" << std::endl;

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

    larlite::ChimeraTrackEvaluator evaluator;

    double parSigmas[3] =  {10,           0.1,   10*TMath::Pi()/180.}; //σ_R,       σ_L,  σ_angles
    //double parTargets[6] = {100,75, 365,  20,   0.3,0.67};            // X0,Y0,Z0, L0,   theta,phi

    evaluator.SetSigmas(parSigmas);
    double max_p_Score = 0;
    int best_p_Trackindex(0),best_p_Run(0),best_p_SubRun(0),best_p_Event(0);
    larlite::track best_p_track;
    std::vector<larlite::hit> best_p_HitCluster;

    for(size_t iEvt = 0;iEvt<FullTargets.size();iEvt++){
        evtID = Form("%03zu",iEvt);
        double parTargets[6] = {
            FullTargets[iEvt][0],
            FullTargets[iEvt][1],
            FullTargets[iEvt][2],
            FullTargets[iEvt][3],
            FullTargets[iEvt][4],
            FullTargets[iEvt][5]};
        evaluator.SetTargets(parTargets);
        while(p_storage.next_event()){
            auto ev_track = p_storage.get_data<larlite::event_track>("pandoraNu"); // tracks
            if(!ev_track)  throw larlite::DataFormatException("Could not locate event_track data product!");

            int run =     p_storage.run_id();
            int subrun =  p_storage.subrun_id();
            int event =   p_storage.event_id();
            larlite::event_hit* ev_hit=nullptr;

            auto const& track_to_hit = p_storage.find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
            if(!ev_hit) throw larlite::DataFormatException("Could not find associated hit data product!");

            double score = 0;

            for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
                auto p_track = track_index;
                auto const& track = (*ev_track)[track_index];
                //
                // Evaluate track score
                score = evaluator.EvalTrack(track);
                if(score > max_p_Score){
                    max_p_Score=score;
                    best_p_Run        = run;
                    best_p_SubRun     = subrun;
                    best_p_Event      = event;
                    best_p_Trackindex = p_track;
                    best_p_track      = track;
                    best_p_HitCluster.clear();
                    for(auto const& hit_index : track_to_hit[track_index]) {
                        best_p_HitCluster.push_back( (*ev_hit)[hit_index] );
                    }
                }
            }
        }

        std::cout << "Max score for proton track : " << max_p_Score << std::endl;

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
        if(!u_storage.is_open()){std::cerr << "File open failed!" << std::endl;return 0;}
        if(!u_storage.is_ready_io()){std::cerr << "I/O preparation failed!" << std::endl;}

        double parSigmasMuon[3] =  {20,0.5,20*TMath::Pi()/180.}; //σ_R, σ_L,  σ_angles
        double parTargetsMuon[6] = {
            FullTargets[iEvt][0],
            FullTargets[iEvt][1],
            FullTargets[iEvt][2],
            FullTargets[iEvt][6],
            FullTargets[iEvt][7],
            FullTargets[iEvt][8]};
        evaluator.SetTargets(parTargetsMuon);
        evaluator.SetSigmas(parSigmasMuon);

        double max_u_Score = 0;
        int best_u_Trackindex(0),best_u_Run(0),best_u_SubRun(0),best_u_Event(0);
        larlite::track best_u_track;
        std::vector<larlite::hit> best_u_HitCluster;

        while(u_storage.next_event()){
            auto ev_track = u_storage.get_data<larlite::event_track>("trackkalmanhit"); // tracks
            if(!ev_track)  throw larlite::DataFormatException("Could not locate event_track data product!");

            int run =     u_storage.run_id();
            int subrun =  u_storage.subrun_id();
            int event =   u_storage.event_id();
            // Get associated hits + association info
            larlite::event_hit* ev_hit=nullptr;

            auto const& track_to_hit = u_storage.find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
            if(!ev_hit) throw larlite::DataFormatException("Could not find associated hit data product!");

            double score = 0;

            for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
                auto u_track = track_index;
                auto const& track = (*ev_track)[track_index];
                //
                // Evaluate track score
                score = evaluator.EvalTrack(track);
                if(score > max_u_Score){
                    max_u_Score=score;
                    best_u_Run        = run;
                    best_u_SubRun     = subrun;
                    best_u_Event      = event;
                    best_u_Trackindex = u_track;
                    best_u_track      = track;
                    best_u_HitCluster.clear();
                    for(auto const& hit_index : track_to_hit[track_index]) {
                        best_u_HitCluster.push_back( (*ev_hit)[hit_index] );
                    }
                }
            }
        }

        std::cout << std::endl << "Max score for muon track : " << max_u_Score << std::endl;

        //=====================
        //=====================
        // Merge tracks =======
        //=====================
        //=====================
        std::vector<std::vector<larlite::hit> > GlobalHitCluster;
        std::vector<TVector3> points;
        TVector3 X0(parTargets[0],parTargets[1],parTargets[2]);
        TVector3 X0proton = best_p_track.Vertex();
        TVector3 X0muon   = best_u_track.Vertex();
        points.push_back(X0);
        //points.push_back(X0proton);
        //points.push_back(X0muon);

        std::vector<larlite::hit> newProtonHits = TranslateHit(best_p_HitCluster, X0, best_p_track.Vertex());
        std::vector<larlite::hit> newMuonHits   = TranslateHit(best_u_HitCluster, X0, best_u_track.Vertex());

        //GlobalHitCluster.push_back(best_p_HitCluster);
        //GlobalHitCluster.push_back(best_u_HitCluster);
        GlobalHitCluster.push_back(newProtonHits);
        GlobalHitCluster.push_back(newMuonHits);

        DrawTrack(GlobalHitCluster,points, evtID);
    }
    
    
    return 0;
}
