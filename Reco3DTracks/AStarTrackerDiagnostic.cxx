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
    //------------------------------------------------------
    bool AStarTrackerDiagnostic::initialize() {
        Window = new TCanvas("Window","Window",1200,400);
        Window->Divide(3,1);
        hLengthComparison = new TH2D("hLengthComparison","hLengthComparison;L_{true} [cm];(L_{reco}-L_{true})/L_{true}",100,0,20,100,-1,1);
        dist2MC_all = new TH1D("dist2MC_all","dist2MC_all;dist MC-reco(SCE corr.) [cm]",100,0,20);
        hdQdx       = new TH1D("hdQdx","hdQdx",25,0,25);
        hdQdxEntries = (TH1D*)hdQdx->Clone("hdQdxEntries");
        return true;
    }
    //------------------------------------------------------
    bool AStarTrackerDiagnostic::analyze(storage_manager* storage) {
        auto ev_track    = storage->get_data<event_track>(_track_producer);       // tracks
        auto ev_chstatus = storage->get_data<event_chstatus>(_chstatus_producer); // chstatus
        auto ev_mct      = storage->get_data<event_mctrack>(_mctrack_producer);   // mctracks
        auto ev_gaushit  = storage->get_data<event_hit>("gaushit");               // gaushits

        _run = storage->run_id();
        _subrun = storage->subrun_id();
        _event = storage->event_id();

        std::cout << _run << "\t" << _subrun << "\t" << _event << std::endl;
        //
        // Check validity
        if(!ev_track){      throw DataFormatException("Could not locate event_track data product!");}
        //if(!ev_chstatus){   throw DataFormatException("Could not locate event_chstatus data product!");}
        if(ev_track->empty()){std::cout << "No tracks in ev_track" << std::endl; return true;}
        //if(!ev_mct){        throw DataFormatException("Could not locate event_mctrack data product!");}
        //
        // Get associated hits + association info
        larlite::event_hit* ev_hit=nullptr;
        auto const& track_to_hit = storage->find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
        if(!ev_hit) throw DataFormatException("Could not find associated hit data product!");
        //
        // Get empty channels
        if(ev_chstatus){
            for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
                for(size_t ch=0;ch<3456;ch++){chstat[iPlane][ch]=-1;}
                for(auto const& chstatus_v : *ev_chstatus) {
                    if( chstatus_v.plane().Plane != iPlane) continue;
                    for(size_t ch=0; ch<chstatus_v.status().size(); ++ch) {
                        auto const& status = chstatus_v.status()[ch];
                        chstat[iPlane][ch]=1;
                        if(status == larlite::larch::kGOOD){chstat[iPlane][ch]=0;}
                    }
                }
            }
        }
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
            if(ev_mct){
                for(auto const& mct : *ev_mct) {
                    double DRstart = sqrt(pow(start_pt.X()-mct.Start().X(),2)+pow(start_pt.Y()-mct.Start().Y(),2)+pow(start_pt.Z()-mct.Start().Z(),2));
                    double DRend = sqrt(pow(end_pt.X()-mct.End().X(),2)+pow(end_pt.Y()-mct.End().Y(),2)+pow(end_pt.Z()-mct.End().Z(),2));
                    if(DRstart > 20 || DRend > 20)continue;
                    thisTrueTrack = mct;
                }
                if(thisTrueTrack.size() == 0){std::cout << "Could not find MC track corresponding to provided start and end points => stop" << std::endl;continue;}
            }

            //
            // Get Hit vectors
            thisHit_v.clear();
            thisGausHit_v.clear();
            for(auto const& h:*ev_gaushit){thisGausHit_v.push_back(h);}
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
            // Compare corrected track to MC and hits
            CompareReco2Hits();
            //CompareRecoCorr2MC();
            //CompareLengths();
            DrawdQdX();
        }

        return true;
    }
    //------------------------------------------------------
    bool AStarTrackerDiagnostic::finalize() {
        DrawPlots();
        _fout->cd();
        //_fout->Write();
        for(TCanvas *c:Window_v)c->Write();
        for(TCanvas *c:dist2MC_v)c->Write();
        dist2MC_all->Write();
        _fout->Close();
        return true;
    }
    //------------------------------------------------------
    void AStarTrackerDiagnostic::VisualizeTrack(){
        //
        // Get Histogram limits
        auto const& start_pt = thisTrack.Vertex();
        auto const& end_pt   = thisTrack.End();
        std::vector<std::pair<double, double> > timeLimit(_numPlanes);
        std::vector<std::pair<double, double> > wireLimit(_numPlanes);
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            timeLimit[iPlane].first  = 1e9;
            timeLimit[iPlane].second = 0;
            wireLimit[iPlane].first  = 1e9;
            wireLimit[iPlane].second = 0;
            //
            // assure that the start and end points of the reconstructed tracks are within bounds
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
            //
            // assure that all hits are within bounds
            /*for(larlite::hit h:thisHit_v){
                double wireProjHit = h.WireID().Wire;
                double timeProjHit = h.PeakTime();
                if(h.WireID().Plane != iPlane)continue;
                if(wireProjHit<wireLimit[iPlane].first ){wireLimit[iPlane].first  = wireProjHit;}
                if(wireProjHit>wireLimit[iPlane].second){wireLimit[iPlane].second = wireProjHit;}
                if(timeProjHit<timeLimit[iPlane].first ){timeLimit[iPlane].first  = timeProjHit;}
                if(timeProjHit>timeLimit[iPlane].second){timeLimit[iPlane].second = timeProjHit;}

            }*/
            //
            // assure that all true nodes are within bounds
            if(thisTrueTrack.size() != 0){
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
            double marginTime = 1*(timeLimit[iPlane].second - timeLimit[iPlane].first);
            double marginWire = 1*(wireLimit[iPlane].second - wireLimit[iPlane].first);
            timeLimit[iPlane].first  -=marginTime;
            timeLimit[iPlane].second +=marginTime;
            wireLimit[iPlane].first  -=marginWire;
            wireLimit[iPlane].second +=marginWire;
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
        for(larlite::hit h:thisHit_v){ //to fill with all gaushits, just change thisHit_v by thisGausHit_v
            size_t iPlane = h.WireID().Plane;
            if(!(wireLimit[iPlane].first < h.WireID().Wire))continue;
            if(!(h.WireID().Wire <wireLimit[iPlane].second))continue;
            if(!(h.PeakTime() > timeLimit[iPlane].first)   )continue;
            if(!(h.PeakTime() < timeLimit[iPlane].second)  )continue;
            hHitImage2D[iPlane]->SetBinContent(h.WireID().Wire+1-wireLimit[iPlane].first,h.PeakTime()+1-timeLimit[iPlane].first,h.SummedADC());
        }
        //
        // Fill graphs
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            double x,y;
            for(size_t iNode=0;iNode<thisTrack.NumberTrajectoryPoints();iNode++){
                x = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisTrack.LocationAtPoint(iNode),iPlane).w/0.3;
                y = X2Tick(thisTrack.LocationAtPoint(iNode).X(),iPlane);
                gRecoedTrack[iPlane]->SetPoint(iNode,x,y);
            }
            for(size_t iNode=0;iNode<thisCorrectedTrack.NumberTrajectoryPoints();iNode++){
                x = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisCorrectedTrack.LocationAtPoint(iNode),iPlane).w/0.3;
                y = X2Tick(thisCorrectedTrack.LocationAtPoint(iNode).X(),iPlane);
                gCorrectedTrack[iPlane]->SetPoint(iNode,x,y);
            }
            if(thisTrueTrack.size()!=0){
                for(size_t iNode=0;iNode<thisTrueTrack.size();iNode++){
                    TVector3 node(thisTrueTrack[iNode].X(),thisTrueTrack[iNode].Y(),thisTrueTrack[iNode].Z());
                    gTrueTrack[iPlane]->SetPoint(iNode,larutil::GeometryHelper::GetME()->Point_3Dto2D(node,iPlane).w/0.3,X2Tick(thisTrueTrack[iNode].X(),iPlane));
                }
            }
        }
        //
        // Fill TLines for dead channels
        std::vector<TLine> deadch[3];
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            for(size_t ch=0;ch<3456;ch++){
                if(ch < wireLimit[iPlane].first || ch > wireLimit[iPlane].second)continue;
                if(chstat[iPlane][ch] == -1)continue;
                if(chstat[iPlane][ch] == 0) continue;
                TLine line(ch,timeLimit[iPlane].first,ch, timeLimit[iPlane].second);
                line.SetLineColor(kGray);
                line.SetLineWidth(2);
                deadch[iPlane].push_back(line);
            }
        }
        //
        // Draw dQdX integration area
        //DrawdQdxSurface();
        //
        // Fill Window
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            Window->cd(iPlane+1);
            hHitImage2D[iPlane]->Draw("colz");
            for(size_t deadchannel = 0;deadchannel<deadch[iPlane].size();deadchannel++){deadch[iPlane][deadchannel].Draw("same");}
            gRecoedTrack[iPlane]->Draw("same LP");
            if(thisTrueTrack.size()!=0){
                gTrueTrack[iPlane]->SetLineWidth(2);
                gTrueTrack[iPlane]->SetLineColor(2);
                gTrueTrack[iPlane]->SetMarkerStyle(4);
                gTrueTrack[iPlane]->SetMarkerColor(2);
                gTrueTrack[iPlane]->Draw("same LP");
            }
            gCorrectedTrack[iPlane]->SetLineColor(4);
            gCorrectedTrack[iPlane]->SetLineWidth(2);
            gCorrectedTrack[iPlane]->SetMarkerStyle(4);
            gCorrectedTrack[iPlane]->SetMarkerColor(4);
            gCorrectedTrack[iPlane]->Draw("same LP");
        }
        //
        //Load canvas in container
        TCanvas *thisWindow = (TCanvas*)Window->Clone(Form("thisWindow_%s",_trackID.c_str()));
        thisWindow->SetName(Form("thisWindow_%s",_trackID.c_str()));
        thisWindow->SetTitle(Form("thisWindow_%s",_trackID.c_str()));
        Window_v.push_back(thisWindow);
    }
    //------------------------------------------------------
    void AStarTrackerDiagnostic::DrawPlots(){
        for(TCanvas *c:Window_v){c->Draw();c->SaveAs(Form("plot/%s.pdf",c->GetName()));}
        for(TCanvas *c:dist2MC_v){c->Draw();c->SaveAs(Form("plot/%s.pdf",c->GetName()));}
        TCanvas *c = new TCanvas("dist2MC_all","dist2MC_all",800,600);
        dist2MC_all->Draw();
        c->SaveAs(Form("plot/%s.pdf",c->GetName()));
        hdQdx->Divide(hdQdxEntries);
        hdQdx->SetMarkerStyle(20);
        hdQdx->Draw();
        c->SaveAs("plot/dQdXprofile.pdf");
        hLengthComparison->Draw("colz");
        c->SaveAs("plot/LengthComparison.pdf");
    }
    //------------------------------------------------------
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
    //------------------------------------------------------
    larlite::track AStarTrackerDiagnostic::CorrectSCE(){
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
    //------------------------------------------------------
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
    //------------------------------------------------------
    void AStarTrackerDiagnostic::CompareRecoCorr2MC(){

        larlite::track recoTrack = thisCorrectedTrack;
        larlite::mctrack trueTrack = thisTrueTrack;
        TH1D *hDistance2MC = new TH1D(Form("hDistance2MC_%s",_trackID.c_str()),Form("hDistance2MC_%s",_trackID.c_str()),100,0,20);
        double trueLength=0;
        for(size_t iNode = 0;iNode<trueTrack.size()-1;iNode++){
            trueLength+=sqrt(pow(trueTrack[iNode+1].X()-trueTrack[iNode].X(),2)+pow(trueTrack[iNode+1].Y()-trueTrack[iNode].Y(),2)+pow(trueTrack[iNode+1].Z()-trueTrack[iNode].Z(),2));
        }

        for(size_t iNode = 0; iNode<recoTrack.NumberTrajectoryPoints();iNode++){
            double distance2track  = 1e9;
            double distance2point1 = 1e9;
            double distance2point2 = 1e9;
            double alpha;
            double ux,uy,uz;
            for(size_t i=0;i<trueTrack.size()-1;i++){
                double X1,Y1,Z1,X2,Y2,Z2,X0,Y0,Z0,Xp,Yp,Zp;
                X1 = trueTrack[i].X();
                Y1 = trueTrack[i].Y();
                Z1 = trueTrack[i].Z();
                X2 = trueTrack[i+1].X();
                Y2 = trueTrack[i+1].Y();
                Z2 = trueTrack[i+1].Z();
                X0 = recoTrack.LocationAtPoint(iNode).X();
                Y0 = recoTrack.LocationAtPoint(iNode).Y();
                Z0 = recoTrack.LocationAtPoint(iNode).Z();
                // 1) check first if the point projects inside the segment.
                // 2) if so, the relevant distance is the distance to the segment
                // 3) if not, the relevant distance is the distance to one of the edges fo the segment

                // 1) check first if the point projects inside the segment.
                ux=(X2-X1); uy=(Y2-Y1); uz=(Z2-Z1);
                alpha=( ux*(X0-X1)+uy*(Y0-Y1)+uz*(Z0-Z1) )/( pow(ux,2)+pow(uy,2)+pow(uz,2) );
                // 2) if so, the relevant distance is the distance to the segment
                if(alpha >= 0 && alpha <= 1){
                    Xp=X1+alpha*ux;
                    Yp=Y1+alpha*uy;
                    Zp=Z1+alpha*uz;
                    double distPoint2thisSegment = sqrt(pow(Xp-X0,2)+pow(Yp-Y0,2)+pow(Zp-Z0,2));
                    if(distPoint2thisSegment < distance2track)distance2track=distPoint2thisSegment;
                }
                // 3) if not, the relevant distance is the distance to one of the edges fo the segment
                else{
                    distance2point1 = sqrt(pow(X0-X1,2)+pow(Y0-Y1,2)+pow(Z0-Z1,2));
                    distance2point2 = sqrt(pow(X0-X2,2)+pow(Y0-Y2,2)+pow(Z0-Z2,2));
                    if(distance2point1<distance2track)distance2track=distance2point1;
                    if(distance2point2<distance2track)distance2track=distance2point2;
                }
            }
            hDistance2MC-> Fill(distance2track);
            dist2MC_all->Fill(distance2track);
        }
        TCanvas *c = new TCanvas(Form("dist2MC_%s",_trackID.c_str()),Form("dist2MC_%s",_trackID.c_str()),800,600);
        hDistance2MC->Draw();
        dist2MC_v.push_back(c);
    }
    //------------------------------------------------------
    void AStarTrackerDiagnostic::CompareLengths(){
        double trueLength=0;
        for(size_t iNode=0;iNode<thisTrueTrack.size()-1;iNode++){
            trueLength+=sqrt(pow(thisTrueTrack[iNode+1].X()-thisTrueTrack[iNode].X(),2)
                             +pow(thisTrueTrack[iNode+1].Y()-thisTrueTrack[iNode].Y(),2)
                             +pow(thisTrueTrack[iNode+1].Z()-thisTrueTrack[iNode].Z(),2));
        }
        double recoLength=thisTrack.Length(0);
        double correctedLength=thisCorrectedTrack.Length(0);
        std::cout << "True Length = " << Form("%.2f",trueLength) << " [cm]" << std::endl;
        std::cout << "Reco Length = " << Form("%.2f",recoLength) << " [cm] \t\t=>\t" << Form("%.2f",(recoLength-trueLength)*100./trueLength) << " \% from MC" << std::endl;
        std::cout << "Corr Length = " << Form("%.2f",correctedLength) << " [cm] \t\t=>\t" << Form("%.2f",(correctedLength-trueLength)*100./trueLength) << " \% from MC" << std::endl;
        hLengthComparison->Fill(trueLength,(correctedLength-trueLength)/trueLength);
    }
    //------------------------------------------------------
    void AStarTrackerDiagnostic::CompareReco2Hits(){
        std::cout << "Write CompareReco2Hits" << std::endl;
    }
    //------------------------------------------------------
    void AStarTrackerDiagnostic::DrawdQdX(){
        for(size_t iNode = 0;iNode<thisTrack.NumberTrajectoryPoints()-1;iNode++){
            double dqdx = thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)+thisTrack.DQdxAtPoint(iNode,larlite::geo::kU)+thisTrack.DQdxAtPoint(iNode,larlite::geo::kV);
            /*
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)==0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU)!=0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)!=0)dqdx*=3./2.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)!=0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU)==0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)!=0)dqdx*=3./2.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)!=0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU)!=0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)==0)dqdx*=3./2.;

            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)==0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)!=0)dqdx*=3.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)!=0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)==0)dqdx*=3.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)==0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) != 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)==0)dqdx*=3.;
             */
            hdQdx->Fill(thisTrack.Length(iNode),dqdx);
            hdQdxEntries->Fill(thisTrack.Length(iNode));
        }
    }
    //------------------------------------------------------
    void AStarTrackerDiagnostic::DrawdQdxSurface(){
        for(size_t iPlane=0;iPlane<_numPlanes;iPlane++){
            double X0,Y0,X1,Y1,X2,Y2,xp,yp,alpha,theta0,theta1,dist;
            double Npx = 5;

            for(size_t iNode=0;iNode<thisTrack.NumberTrajectoryPoints()-1;iNode++){
                for(int row = 0;row<hHitImage2D[iPlane]->GetNbinsX();row++){

                    for(int col = 0;col<hHitImage2D[iPlane]->GetNbinsY();col++){

                        X1 = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisTrack.LocationAtPoint(iNode  ),iPlane).w/0.3;
                        X2 = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisTrack.LocationAtPoint(iNode+1),iPlane).w/0.3;

                        Y1 = X2Tick(thisTrack.LocationAtPoint(iNode  ).X(),iPlane);
                        Y2 = X2Tick(thisTrack.LocationAtPoint(iNode+1).X(),iPlane);

                        if(X1 == X2){X1 -= 0.001;}
                        if(Y1 == Y2){Y1 -= 0.001;}

                        X0=row+hHitImage2D[iPlane]->GetXaxis()->GetXmin();
                        Y0=col+hHitImage2D[iPlane]->GetYaxis()->GetXmin();

                        alpha = ( (X2-X1)*(X0-X1)+(Y2-Y1)*(Y0-Y1) )/( pow(X2-X1,2)+pow(Y2-Y1,2) );
                        if(alpha < -0.1 || alpha > 1.1) continue;
                        xp = X1+alpha*(X2-X1);
                        yp = Y1+alpha*(Y2-Y1);
                        dist = sqrt(pow(xp-X0,2)+pow(yp-Y0,2));
                        if(dist > Npx)continue;

                        theta0 = std::acos( ((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0)) / (sqrt(pow(X1-X0,2)+pow(Y1-Y0,2))*sqrt(pow(X2-X0,2)+pow(Y2-Y0,2))) );

                        bool bestCandidate = true;

                        for(size_t jNode=0;jNode<thisTrack.NumberTrajectoryPoints()-1;jNode++){

                            if(iNode==jNode)continue;
                            double x1,y1,x2,y2;
                            x1 = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisTrack.LocationAtPoint(jNode  ),iPlane).w/0.3;
                            x2 = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisTrack.LocationAtPoint(jNode+1),iPlane).w/0.3;

                            y1 = X2Tick(thisTrack.LocationAtPoint(jNode  ).X(),iPlane);
                            y2 = X2Tick(thisTrack.LocationAtPoint(jNode+1).X(),iPlane);

                            theta1 = std::acos( ((x1-X0)*(x2-X0)+(y1-Y0)*(y2-Y0)) / (sqrt(pow(x1-X0,2)+pow(y1-Y0,2))*sqrt(pow(x2-X0,2)+pow(y2-Y0,2))) );
                            if(theta1 > theta0){bestCandidate=false;}

                        }

                        if(!bestCandidate)continue;
                        hHitImage2D[iPlane]->SetBinContent(row+1,col+1,1);
                    }
                }
            }
        }
    }
    //------------------------------------------------------
}
#endif
