#ifndef LARLITE_CHIMERAPATCHING_CXX
#define LARLITE_CHIMERAPATCHING_CXX

#include "ChimeraPatching.h"

namespace larlite {

    bool ChimeraPatching::initialize() {
        hX0     = new TH1D("hX0",    "hX0",    100,0,250);
        hY0     = new TH1D("hY0",    "hY0",    100,-120,120);
        hZ0     = new TH1D("hY0",    "hY0",    1100,0,1100);
        hL0     = new TH1D("hL0",    "hL0",    100,0,1000);
        hTheta0 = new TH1D("hTheta0","hTheta0",90,0,0.5*3.141597);
        hPhi0   = new TH1D("hPhi0",  "hPhi0",  180,-0.5*3.141597,0.5*3.141597);

        hX     = new TH1D("hX",    "hX",    100,-200,200);
        hY     = new TH1D("hY",    "hY",    100,-200,200);
        hZ     = new TH1D("hY",    "hY",    1000,0,10000);
        hL     = new TH1D("hL",    "hL",    100,0,1000);
        hTheta = new TH1D("hTheta","hTheta",90,0,0.5*3.141597);
        hPhi   = new TH1D("hPhi",  "hPhi",  180,-0.5*3.141597,0.5*3.141597);

        hScore = new TH1D("hScore","hScore",100,0,1);

        double parSigmas[3] = {20,100,10*TMath::Pi()/180.}; //σ_R,σ_L,σ_angles
        double parTargets[6] = {100,75,365,20,0.3,0.67}; // X0,Y0,Z0,L0,theta,phi
        _evaluator.SetSigmas(parSigmas);
        _evaluator.SetTargets(parTargets);
        return true;
    }

    bool ChimeraPatching::analyze(storage_manager* storage) {
        auto ev_track    = storage->get_data<event_track>(_track_producer); // tracks
        auto ev_gaushit  = storage->get_data<event_hit>(_hit_producer);     // gaushits

        if(!ev_track)  throw DataFormatException("Could not locate event_track data product!");
        //if(!ev_gaushit)throw DataFormatException("Could not locate event_gaushit data product!");

        // Get associated hits + association info
        larlite::event_hit* ev_hit=nullptr;
        auto const& track_to_hit = storage->find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
        if(!ev_hit) throw DataFormatException("Could not find associated hit data product!");

        for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
            auto const& track = (*ev_track)[track_index];
            hX0->Fill(track.Vertex().X());
            hY0->Fill(track.Vertex().Y());
            hZ0->Fill(track.Vertex().Z());
            hL0->Fill(track.Length(0));
            hTheta0->Fill(track.Vertex().Theta());
            hPhi0->Fill(track.Vertex().Phi());
            double score = _evaluator.EvalTrack(track);
            std::cout << "Track #" << track_index << " :\t" << score << std::endl;
            hX->Fill(track.Vertex().X());
            hY->Fill(track.Vertex().Y());
            hZ->Fill(track.Vertex().Z());
            hL->Fill(track.Length(0));
            hTheta->Fill(track.Vertex().Theta());
            hPhi->Fill(track.Vertex().Phi());
            hScore->Fill(score);
        }

        return true;
    }

    bool ChimeraPatching::finalize() {
        _fout->cd();
        _fout->Write();
        return true;
    }
    
}
#endif
