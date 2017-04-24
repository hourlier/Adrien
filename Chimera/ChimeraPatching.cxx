#ifndef LARLITE_CHIMERAPATCHING_CXX
#define LARLITE_CHIMERAPATCHING_CXX

#include "ChimeraPatching.h"
#include "TH2D.h"
#include "TCanvas.h"

namespace larlite {

    bool ChimeraPatching::initialize() {
        hX0     = new TH1D("hX0",    "hX0",    100,0,250);
        hY0     = new TH1D("hY0",    "hY0",    100,-120,120);
        hZ0     = new TH1D("hZ0",    "hZ0",    1100,0,1100);
        hL0     = new TH1D("hL0",    "hL0",    100,0,1000);
        hTheta0 = new TH1D("hTheta0","hTheta0",90,0,0.5*3.141597);
        hPhi0   = new TH1D("hPhi0",  "hPhi0",  180,-0.5*3.141597,0.5*3.141597);

        hX     = new TH1D("hX",    "hX",    100,0,250);
        hY     = new TH1D("hY",    "hY",    100,-120,120);
        hZ     = new TH1D("hZ",    "hZ",    1100,0,1100);
        hL     = new TH1D("hL",    "hL",    100,0,1000);
        hTheta = new TH1D("hTheta","hTheta",90,0,0.5*3.141597);
        hPhi   = new TH1D("hPhi",  "hPhi",  180,-0.5*3.141597,0.5*3.141597);

        hScore = new TH1D("hScore","hScore",100,0,1);

        double parSigmas[3] =  {10,           0.1,   10*TMath::Pi()/180.}; //σ_R,       σ_L,  σ_angles
        double parTargets[6] = {100,75, 365,  20,   0.3,0.67};            // X0,Y0,Z0, L0,   theta,phi
        _evaluator.SetSigmas(parSigmas);
        _evaluator.SetTargets(parTargets);

        maxScore = 0;
        bestTrackindex=0;

        return true;
    }
    //____________________________________________________________
    bool ChimeraPatching::analyze(storage_manager* storage) {
        auto ev_track    = storage->get_data<event_track>(_track_producer); // tracks
        auto ev_gaushit  = storage->get_data<event_hit>(_hit_producer);     // gaushits

        _run = storage->run_id();
        _subrun = storage->subrun_id();
        _event = storage->event_id();

        if(!ev_track)  throw DataFormatException("Could not locate event_track data product!");
        if(!ev_gaushit)throw DataFormatException("Could not locate event_gaushit data product!");

        // Get associated hits + association info
        larlite::event_hit* ev_hit=nullptr;

        auto const& track_to_hit = storage->find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
        if(!ev_hit) throw DataFormatException("Could not find associated hit data product!");

        double score = 0;

        for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
            _track = track_index;
            auto const& track = (*ev_track)[track_index];
            //
            // Evaluate track score
            score = _evaluator.EvalTrack(track);

            hX0->Fill(track.Vertex().X());
            hY0->Fill(track.Vertex().Y());
            hZ0->Fill(track.Vertex().Z());
            hL0->Fill(track.Length(0));
            hTheta0->Fill(track.Vertex().Theta());
            hPhi0->Fill(track.Vertex().Phi());
            hScore->Fill(score);

            //if(score < 0.05)continue;
            if(score > maxScore){
                maxScore=score;
                bestRun        = _run;
                bestSubRun     = _subrun;
                bestEvent      = _event;
                bestTrackindex = _track;
                bestTrack      = track;
                bestHitCluster.clear();
                for(auto const& hit_index : track_to_hit[track_index]) {
                    bestHitCluster.push_back( (*ev_hit)[hit_index] );
                }
            }

            if(score > 0.05){
                hX->Fill(track.Vertex().X());
                hY->Fill(track.Vertex().Y());
                hZ->Fill(track.Vertex().Z());
                hL->Fill(track.Length(0));
                hTheta->Fill(track.Vertex().Theta());
                hPhi->Fill(track.Vertex().Phi());
            }
        }
        return true;
    }
    //____________________________________________________________
    bool ChimeraPatching::finalize() {
        DrawTrack(bestHitCluster);
        std::vector<double>targets = _evaluator.GetTargets();
        if(maxScore < 0.05){std::cout << "ERROR, could not find suitable track" << std::endl;}
        std::cout << "best candidate is track : " << bestRun << "\t" << bestSubRun << "\t" << bestEvent << "\t" << bestTrackindex << " with a score of " << maxScore << std::endl;
        std::cout << "X0    = " << bestTrack.Vertex().X()     << " : " << targets[0] << std::endl;
        std::cout << "Y0    = " << bestTrack.Vertex().Y()     << " : " << targets[1] << std::endl;
        std::cout << "Z0    = " << bestTrack.Vertex().Z()     << " : " << targets[2] << std::endl;
        std::cout << "L     = " << bestTrack.Length(0)        << " : " << targets[3] << std::endl;
        std::cout << "theta = " << bestTrack.Vertex().Theta() << " : " << targets[4] << std::endl;
        std::cout << "phi   = " << bestTrack.Vertex().Phi()   << " : " << targets[5] << std::endl;
        _fout->cd();
        _fout->Write();
        return true;
    }
    //____________________________________________________________
    void ChimeraPatching::DrawTrack(const std::vector<larlite::hit> HitCluster){
        TCanvas *cTrack = new TCanvas("cTrack","cTrack",1200,300);
        cTrack->Divide(3,1);
        TH2D *hTrack[3];
        int timebounds[2] = {(int)(1e9),0};
        int wirebounds[3][2];
        for(size_t iPlane=0;iPlane<3;iPlane++){
            wirebounds[iPlane][0] = (int)(1e9) ;
            wirebounds[iPlane][1] = 0;

            for(larlite::hit h:HitCluster){
                if(h.WireID().Plane!=iPlane)continue;
                if(h.PeakTime()    < timebounds[0]){timebounds[0] = h.PeakTime();}
                if(h.PeakTime()    > timebounds[1]){timebounds[1] = h.PeakTime();}
                if(h.WireID().Wire < wirebounds[iPlane][0]){wirebounds[iPlane][0] = h.WireID().Wire;}
                if(h.WireID().Wire > wirebounds[iPlane][1]){wirebounds[iPlane][1] = h.WireID().Wire;}
            }
            int marginwire = (int)(0.5*(wirebounds[iPlane][1]-wirebounds[iPlane][0]));
            wirebounds[iPlane][0] -= marginwire;
            wirebounds[iPlane][1] += marginwire;
        }
        int margintime = (int)(0.5*(timebounds[1]-timebounds[0]));
        timebounds[0]-=margintime;
        timebounds[1]+=margintime;

        for(size_t iPlane = 0;iPlane<3;iPlane++){
            hTrack[iPlane] = new TH2D(Form("hTrack_%zu",iPlane),Form("hTrack_%zu;wire;tick",iPlane),wirebounds[iPlane][1]-wirebounds[iPlane][0],wirebounds[iPlane][0],wirebounds[iPlane][1],timebounds[1]-timebounds[0],timebounds[0],timebounds[1]);
            for(larlite::hit h:HitCluster){
                if(h.WireID().Plane!=iPlane)continue;
                hTrack[iPlane]->SetBinContent(h.WireID().Wire+1-wirebounds[iPlane][0],h.PeakTime()+1-timebounds[0],h.Integral());
            }
            cTrack->cd(iPlane+1);
            hTrack[iPlane]->Draw("colz");
        }
        _fout->cd();
        cTrack->Write();
    }
}
#endif
