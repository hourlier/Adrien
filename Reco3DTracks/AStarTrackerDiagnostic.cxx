#ifndef LARLITE_AStarTrackerDiagnostic_CXX
#define LARLITE_AStarTrackerDiagnostic_CXX

#include "AStarTrackerDiagnostic.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/mctrack.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/SCE/SpaceChargeMicroBooNE.h"

namespace larlite {

    bool AStarTrackerDiagnostic::initialize() {
        Window = new TCanvas("Window","Window",1200,400);
        Window->Divide(3,1);
        dist2MC_all = new TH1D("dist2MC_all","dist2MC_all",100,0,20);
        return true;
    }

    bool AStarTrackerDiagnostic::analyze(storage_manager* storage) {
        auto ev_track    = storage->get_data<event_track>(_track_producer);       // tracks
        auto ev_chstatus = storage->get_data<event_chstatus>(_chstatus_producer); // chstatus
        auto ev_mct      = storage->get_data<event_mctrack>(_mctrack_producer);   // mctracks

        _run = storage->run_id();
        _subrun = storage->subrun_id();
        _event = storage->event_id();

        std::cout << _run << "\t" << _subrun << "\t" << _event << std::endl;

        //
        // Check validity
        if(!ev_track){      throw DataFormatException("Could not locate event_track data product!");}
        if(!ev_chstatus){   throw DataFormatException("Could not locate event_chstatus data product!");}
        if(ev_track->empty()){std::cout << "No tracks in ev_track" << std::endl; return true;}
        if(!ev_mct){        throw DataFormatException("Could not locate event_mctrack data product!");}

        //
        // Get associated hits + association info
        larlite::event_hit* ev_hit=nullptr;
        auto const& track_to_hit = storage->find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
        if(!ev_hit) throw DataFormatException("Could not find associated hit data product!");

        //
        // Loop over tracks
        for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
            thisTrack = ev_track->at(track_index);
            _trackID = Form("%05d_%05d_%05d_%d",_run,_subrun,_event,thisTrack.ID());
            if(thisTrack.NumberTrajectoryPoints() == 0)continue;
            auto const& start_pt = thisTrack.Vertex();
            auto const& end_pt   = thisTrack.End();
            std::cout << "\t\t" << thisTrack.ID() << "\t" << thisTrack.Length(0) << std::endl;

            //
            // Find corresponding MC track
            for(auto const& mct : *ev_mct) {
                double DRstart = sqrt(pow(start_pt.X()-mct.Start().X(),2)+pow(start_pt.Y()-mct.Start().Y(),2)+pow(start_pt.Z()-mct.Start().Z(),2));
                double DRend = sqrt(pow(end_pt.X()-mct.End().X(),2)+pow(end_pt.Y()-mct.End().Y(),2)+pow(end_pt.Z()-mct.End().Z(),2));
                if(DRstart > 20 || DRend > 20)continue;
                thisTrueTrack = mct;
            }
            if(thisTrueTrack.size() == 0){std::cout << "Could not find MC track corresponding to provided start and end points => stop" << std::endl;continue;}

            //
            // Get Hit vector
            thisHit_v.clear();
            for(auto const& hit_index : track_to_hit[track_index]) {
                auto const& h = (*ev_hit)[hit_index];
                thisHit_v.push_back(h);
            }
            //
            // Correct SCE
            thisCorrectedTrack=CorrectSCE();
            //
            // Create track histogram
            VisualizeTrack();
            //
            // Compare corrected track to MC in 3D
            CompareRecoCorr2MC();
        }

        return true;
    }

    bool AStarTrackerDiagnostic::finalize() {
        _fout->cd();
        //_fout->Write();
        for(TCanvas *c:Window_v)c->Write();
        for(TCanvas *c:dist2MC_v)c->Write();
        dist2MC_all->Write();
        _fout->Close();
        return true;
    }

    void AStarTrackerDiagnostic::VisualizeTrack(){
        //
        // Get Histogram limits
        std::cout << "VisualizeTrack " << _trackID << std::endl;
        auto const& start_pt = thisTrack.Vertex();
        auto const& end_pt   = thisTrack.End();
        std::vector<std::pair<double, double> > timeLimit(_numPlanes);
        std::vector<std::pair<double, double> > wireLimit(_numPlanes);
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            timeLimit[iPlane].first  = 1e9;
            timeLimit[iPlane].second = 0;
            wireLimit[iPlane].first  = 1e9;
            wireLimit[iPlane].second = 0;
            double wireProjStartPt = larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w/0.3;
            double timeProjStartPt = X2Tick(start_pt.X(),iPlane);
            if(wireProjStartPt<wireLimit[iPlane].first ){wireLimit[iPlane].first  = wireProjStartPt;}
            if(wireProjStartPt>wireLimit[iPlane].second){wireLimit[iPlane].second = wireProjStartPt;}
            if(timeProjStartPt<timeLimit[iPlane].first ){timeLimit[iPlane].first  = timeProjStartPt;}
            if(timeProjStartPt>timeLimit[iPlane].second){timeLimit[iPlane].second = timeProjStartPt;}

            double wireProjEndPt = larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w/0.3;
            double timeProjEndPt = X2Tick(end_pt.X(),iPlane);
            if(wireProjEndPt<wireLimit[iPlane].first ){wireLimit[iPlane].first  = wireProjEndPt;}
            if(wireProjEndPt>wireLimit[iPlane].second){wireLimit[iPlane].second = wireProjEndPt;}
            if(timeProjEndPt<timeLimit[iPlane].first ){timeLimit[iPlane].first  = timeProjEndPt;}
            if(timeProjEndPt>timeLimit[iPlane].second){timeLimit[iPlane].second = timeProjEndPt;}

            for(larlite::hit h:thisHit_v){
                double wireProjHit = h.WireID().Wire;
                double timeProjHit = h.PeakTime();
                if(h.WireID().Plane != iPlane)continue;
                if(wireProjHit<wireLimit[iPlane].first ){wireLimit[iPlane].first  = wireProjHit;}
                if(wireProjHit>wireLimit[iPlane].second){wireLimit[iPlane].second = wireProjHit;}
                if(timeProjHit<timeLimit[iPlane].first ){timeLimit[iPlane].first  = timeProjHit;}
                if(timeProjHit>timeLimit[iPlane].second){timeLimit[iPlane].second = timeProjHit;}

            }
            for(size_t iNode=0;iNode<thisTrueTrack.size();iNode++){
                TVector3 node(thisTrueTrack[iNode].X(),thisTrueTrack[iNode].Y(),thisTrueTrack[iNode].Z());
                double wireProjmc = larutil::GeometryHelper::GetME()->Point_3Dto2D(node,iPlane).w/0.3;
                double timeProjmc = X2Tick(node.X(),iPlane);
                if(wireProjmc<wireLimit[iPlane].first ){wireLimit[iPlane].first  = wireProjmc;}
                if(wireProjmc>wireLimit[iPlane].second){wireLimit[iPlane].second = wireProjmc;}
                if(timeProjmc<timeLimit[iPlane].first ){timeLimit[iPlane].first  = timeProjmc;}
                if(timeProjmc>timeLimit[iPlane].second){timeLimit[iPlane].second = timeProjmc;}
            }
        }
        //
        // Equalize time margins
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            for(size_t jPlane=0;jPlane<_numPlanes;jPlane++){
                if(timeLimit[jPlane].first <timeLimit[iPlane].first ){timeLimit[iPlane].first =timeLimit[jPlane].first;}
                if(timeLimit[jPlane].second>timeLimit[iPlane].second){timeLimit[iPlane].second=timeLimit[jPlane].second;}
            }
        }
        //
        // Add margin to the track, declare histograms ans graphs
        for(size_t iPlane = 0;iPlane<_numPlanes;iPlane++){
            timeLimit[iPlane].first  -=20;
            timeLimit[iPlane].second +=20;
            wireLimit[iPlane].first  -=20;
            wireLimit[iPlane].second +=20;
            if(timeLimit[iPlane].first<0)timeLimit[iPlane].first=0;
            if(wireLimit[iPlane].first<0)wireLimit[iPlane].first=0;
            hHitImage2D[iPlane] = new TH2D(Form("hHitImage_%s_%zu",_trackID.c_str(),iPlane),Form("hHitImage_%s_%zu;wire;time",_trackID.c_str(),iPlane),wireLimit[iPlane].second-wireLimit[iPlane].first+1,wireLimit[iPlane].first,wireLimit[iPlane].second,timeLimit[iPlane].second-timeLimit[iPlane].first+1,timeLimit[iPlane].first,timeLimit[iPlane].second);
            gRecoedTrack[iPlane] = new TGraph();
            gRecoedTrack[iPlane]->SetNameTitle(Form("gRecoedTrack_%s_%zu",_trackID.c_str(),iPlane),Form("gRecoedTrack_%s_%zu",_trackID.c_str(),iPlane));
            gTrueTrack[iPlane] = new TGraph();
            gTrueTrack[iPlane]->SetNameTitle(Form("gTrueTrack_%s_%zu",_trackID.c_str(),iPlane),Form("gTrueTrack_%s_%zu",_trackID.c_str(),iPlane));
            gCorrectedTrack[iPlane] = new TGraph();
            gCorrectedTrack[iPlane]->SetNameTitle(Form("gCorrectedTrack_%s_%zu",_trackID.c_str(),iPlane),Form("gCorrectedTrack_%s_%zu",_trackID.c_str(),iPlane));

        }
        //
        // Fill histograms
        for(larlite::hit h:thisHit_v){
            size_t iPlane = h.WireID().Plane;
            if(!(wireLimit[iPlane].first < h.WireID().Wire)){std::cout << iPlane << "\t" << wireLimit[iPlane].first << " !< " << h.WireID().Wire << std::endl;continue;}
            if(!(h.WireID().Wire <wireLimit[iPlane].second)){std::cout << iPlane << "\t" << h.WireID().Wire << " !< " << wireLimit[iPlane].second << std::endl;continue;}
            if(!(h.PeakTime() > timeLimit[iPlane].first)   ){std::cout << iPlane << "\t" << timeLimit[iPlane].first << " !< " << h.PeakTime() << std::endl;continue;}
            if(!(h.PeakTime() < timeLimit[iPlane].second)  ){std::cout << iPlane << "\t" << h.PeakTime() << " !< " << timeLimit[iPlane].second << std::endl;continue;}
            hHitImage2D[iPlane]->SetBinContent(h.WireID().Wire+1-wireLimit[iPlane].first,h.PeakTime()+1-timeLimit[iPlane].first,h.SummedADC());
        }
        //
        // Fill graphs
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            for(size_t iNode=0;iNode<thisTrack.NumberTrajectoryPoints();iNode++){
                gRecoedTrack[iPlane]->SetPoint(iNode,larutil::GeometryHelper::GetME()->Point_3Dto2D(thisTrack.LocationAtPoint(iNode),iPlane).w/0.3,X2Tick(thisTrack.LocationAtPoint(iNode).X(),iPlane));
            }
            for(size_t iNode=0;iNode<thisCorrectedTrack.NumberTrajectoryPoints();iNode++){
                gCorrectedTrack[iPlane]->SetPoint(iNode,larutil::GeometryHelper::GetME()->Point_3Dto2D(thisCorrectedTrack.LocationAtPoint(iNode),iPlane).w/0.3,X2Tick(thisCorrectedTrack.LocationAtPoint(iNode).X(),iPlane));
            }
            for(size_t iNode=0;iNode<thisTrueTrack.size();iNode++){
                TVector3 node(thisTrueTrack[iNode].X(),thisTrueTrack[iNode].Y(),thisTrueTrack[iNode].Z());
                gTrueTrack[iPlane]->SetPoint(iNode,larutil::GeometryHelper::GetME()->Point_3Dto2D(node,iPlane).w/0.3,X2Tick(thisTrueTrack[iNode].X(),iPlane));
            }
        }
        //
        // Fill Window
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            Window->cd(iPlane+1);
            hHitImage2D[iPlane]->Draw("colz");
            gRecoedTrack[iPlane]->Draw("same LP");
            gTrueTrack[iPlane]->SetLineWidth(2);
            gTrueTrack[iPlane]->SetLineColor(2);
            gTrueTrack[iPlane]->Draw("same LP");
            gCorrectedTrack[iPlane]->SetLineColor(4);
            gCorrectedTrack[iPlane]->Draw("same LP");
        }
        //
        //Load canvas in container
        TCanvas *thisWindow = (TCanvas*)Window->Clone(Form("thisWindow_%s",_trackID.c_str()));
        thisWindow->SetName(Form("thisWindow_%s",_trackID.c_str()));
        thisWindow->SetTitle(Form("thisWindow_%s",_trackID.c_str()));
        Window_v.push_back(thisWindow);
    }

    double AStarTrackerDiagnostic::X2Tick(double x, size_t plane) const {

        auto ts = larutil::TimeService::GetME();
        auto larp = larutil::LArProperties::GetME();

        // (X [cm] / Drift Velocity [cm/us] - TPC waveform tick 0 offset) ... [us]
        double tick = (x / larp->DriftVelocity() - ts->TriggerOffsetTPC() - _speedOffset);
        // 1st plane's tick
        if(plane==0) return tick * 2;// 2 ticks/us
        // 2nd plane needs more drift time for plane0=>plane1 gap (0.3cm) difference
        tick -= 0.3 / larp->DriftVelocity(larp->Efield(1));
        // 2nd plane's tick
        if(plane==1) return tick * 2;
        // 3rd plane needs more drift time for plane1=>plane2 gap (0.3cm) difference
        tick -= 0.3 / larp->DriftVelocity(larp->Efield(2));
        return tick * 2;
    }

    larlite::track AStarTrackerDiagnostic::CorrectSCE(){
        std::cout << "correct SCE" << std::endl;
        // larlite::track => std::vector<TVector3>
        std::vector<TVector3> newPath;
        TVector3 node;
        for(size_t iNode = 0;iNode<thisTrack.NumberTrajectoryPoints();iNode++){
            node = thisTrack.LocationAtPoint(iNode);
            newPath.push_back(node);
        }

        // correct for SCE
        std::vector<TVector3> newNodes = CorrectSCE(newPath);
        //std::vector<TVector3> newNodes = newPath;

        // std::vector <TVector3> => larlite::track
        larlite::track recoTrack;
        for(size_t iNode=0;iNode<newNodes.size();iNode++){
            recoTrack.add_vertex(newNodes[iNode]);
            if(iNode<newNodes.size()-1){
                double norm = sqrt(pow(newNodes[iNode+1][0]-newNodes[iNode][0],2)
                                   +pow(newNodes[iNode+1][1]-newNodes[iNode][1],2)
                                   +pow(newNodes[iNode+1][2]-newNodes[iNode][2],2));
                TVector3 direction((newNodes[iNode+1][0]-newNodes[iNode][0])/norm,
                                   (newNodes[iNode+1][1]-newNodes[iNode][1])/norm,
                                   (newNodes[iNode+1][2]-newNodes[iNode][2])/norm);
                recoTrack.add_direction(direction);
            }
            else{
                double norm = sqrt(pow(newNodes[iNode][0]-newNodes[0][0],2)
                                   +pow(newNodes[iNode][1]-newNodes[0][1],2)
                                   +pow(newNodes[iNode][2]-newNodes[0][2],2));
                TVector3 direction((newNodes[iNode][0]-newNodes[0][0])/norm,
                                   (newNodes[iNode][1]-newNodes[0][1])/norm,
                                   (newNodes[iNode][2]-newNodes[0][2])/norm);
                recoTrack.add_direction(direction);
            }
        }
        return recoTrack;
    }
    std::vector<TVector3> AStarTrackerDiagnostic::CorrectSCE(std::vector<TVector3> originpath){
        std::vector<TVector3> newPath;
        larlitecv::SpaceChargeMicroBooNE *mySCEcorr = new larlitecv::SpaceChargeMicroBooNE();
        for(size_t iNode=0;iNode<originpath.size();iNode++){
            std::vector<double> PosCorr = mySCEcorr->GetPosOffsets(originpath.at(iNode).X(),originpath.at(iNode).Y(),originpath.at(iNode).Z());
            TVector3 newNode(originpath.at(iNode).X()+(PosCorr[0]-0.7), originpath.at(iNode).Y()-PosCorr[1], originpath.at(iNode).Z()-PosCorr[2]  );
            newPath.push_back(newNode);
        }
        return newPath;
    }

    void AStarTrackerDiagnostic::CompareRecoCorr2MC(){

        larlite::track recoTrack = thisCorrectedTrack;
        larlite::mctrack trueTrack = thisTrueTrack;
        TH1D *hDistance2MC = new TH1D(Form("hDistance2MC_%s",_trackID.c_str()),Form("hDistance2MC_%s",_trackID.c_str()),100,0,20);
        std::cout << "comparing Reco and MC tracks" << std::endl;
        double trueLength=0;
        for(size_t iNode = 0;iNode<trueTrack.size()-1;iNode++){
            trueLength+=sqrt(pow(trueTrack[iNode+1].X()-trueTrack[iNode].X(),2)+pow(trueTrack[iNode+1].Y()-trueTrack[iNode].Y(),2)+pow(trueTrack[iNode+1].Z()-trueTrack[iNode].Z(),2));
        }
        std::cout << "reconstructed length : " << recoTrack.Length(0) << " [cm]" << std::endl;
        std::cout << "true length          : " << trueLength          << " [cm]" << std::endl;

        std::cout << "Reco'ed nodes" << std::endl;
        for(size_t iNode=0; iNode<recoTrack.NumberTrajectoryPoints();iNode++){
            std::cout << Form("   node (X,Y,Z) : %3.2f \t %3.2f \t %3.2f",recoTrack.LocationAtPoint(iNode).X(),recoTrack.LocationAtPoint(iNode).Y(),recoTrack.LocationAtPoint(iNode).Z()) << std::endl;
        }
        std::cout << std::endl;
        std::cout << "MC nodes" << std::endl;
        for(auto const& step : trueTrack) {
            std::cout << Form("MC node (X,Y,Z) : %3.2f \t %3.2f \t %3.2f",step.X(),step.Y(),step.Z()) << std::endl;
        }

        for(size_t iNode = 0; iNode<recoTrack.NumberTrajectoryPoints();iNode++){
            double dPointSegment = 1e9;
            double dPointSegmentX = 1e9;
            double dPointSegmentY = 1e9;
            double dPointSegmentZ = 1e9;
            for(size_t i=0;i<trueTrack.size()-1;i++){
                double X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0;
                X1 = trueTrack[i].X();
                Y1 = trueTrack[i].Y();
                Z1 = trueTrack[i].Z();
                X2 = trueTrack[i+1].X();
                Y2 = trueTrack[i+1].Y();
                Z2 = trueTrack[i+1].Z();
                X0 = recoTrack.LocationAtPoint(iNode).X();
                Y0 = recoTrack.LocationAtPoint(iNode).Y();
                Z0 = recoTrack.LocationAtPoint(iNode).Z();

                double distPoint2thisSegment =  sqrt(pow((Y0-Y1)*(Z0-Z2)-(Z0-Z1)*(Y0-Y2),2)
                                                     +pow((Z0-Z1)*(X0-X2)-(X0-X1)*(Z0-Z2),2)
                                                     +pow((X0-X1)*(Y0-Y2)-(Y0-Y1)*(X0-X2),2))/sqrt(pow(X2-X1,2)+pow(Y2-Y1,2)+pow(Z2-Z1,2));
                if(distPoint2thisSegment < dPointSegment){
                    dPointSegment=distPoint2thisSegment;
                    dPointSegmentX=((Y0-Y1)*(Z0-Z2)-(Z0-Z1)*(Y0-Y2))/sqrt(pow(X2-X1,2)+pow(Y2-Y1,2)+pow(Z2-Z1,2));
                    dPointSegmentY=((Z0-Z1)*(X0-X2)-(X0-X1)*(Z0-Z2))/sqrt(pow(X2-X1,2)+pow(Y2-Y1,2)+pow(Z2-Z1,2));
                    dPointSegmentZ=((X0-X1)*(Y0-Y2)-(Y0-Y1)*(X0-X2))/sqrt(pow(X2-X1,2)+pow(Y2-Y1,2)+pow(Z2-Z1,2));
                }
            }
            hDistance2MC-> Fill(dPointSegment);
            dist2MC_all->Fill(dPointSegment);
            //hDistance2MCX->Fill(dPointSegmentX);
            //hDistance2MCY->Fill(dPointSegmentY);
            //hDistance2MCZ->Fill(dPointSegmentZ);
            if(dPointSegment > 20){std::cout << "STRANGE, TOO FAR FROM MC" << std::endl;}
        }
        TCanvas *c = new TCanvas(Form("dist2MC_%s",_trackID.c_str()),Form("dist2MC_%s",_trackID.c_str()),800,600);
        hDistance2MC->Draw();
        dist2MC_v.push_back(c);
    }

}
#endif
