#ifndef LARLITE_AStarTracker_CXX
#define LARLITE_AStarTracker_CXX

#include "AStarTracker.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/mctrack.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/core/DataFormat/ChStatus.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgoProton.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/SCE/SpaceChargeMicroBooNE.h"

namespace larlite {

    bool AStarTracker::initialize() {
        ReadProtonTrackFile();
        _eventTreated = 0;
        _eventSuccess = 0;
        hDistance2Hit = new TH1D("hDistance2Hit","hDistance2Hit",100,0,50);
        hdQdx = new TH1D("hdQdx","hdQdx",100,0,100);
        hdQdxEntries = new TH1D("hdQdxEntries","hdQdxEntries",100,0,100);
        return true;
    }

    bool AStarTracker::analyze(storage_manager* storage) {
        gStyle->SetOptStat(0);
        //
        // Retrieve data products
        //
        auto ev_track    = storage->get_data<event_track>(_track_producer);       // tracks
        auto ev_chstatus = storage->get_data<event_chstatus>(_chstatus_producer); // chstatus
        auto ev_mct      = storage->get_data<event_mctrack>(_mctrack_producer);   // mctracks
        auto ev_gaushit  = storage->get_data<event_hit>("gaushit");               // gaushits

        // Check validity
        if(!ev_track)    throw DataFormatException("Could not locate event_track data product!");
        //if(!ev_chstatus) throw DataFormatException("Could not locate event_chstatus data product!");
        //if(!ev_mct)      throw DataFormatException("Could not locate event_mctrack data product!");
        if(!ev_gaushit)  throw DataFormatException("Could not locate event_gaushit data product!");
        if(ev_track->empty()){tellMe("ERROR, ev_track empty",0);return true;}

        _run = storage->run_id();
        _subrun = storage->subrun_id();
        _event = storage->event_id();

        tellMe(Form("%d\t%d\t%d",_run,_subrun,_event),1);
        //if(!(_run == 6146 && _subrun ==	2 && _event ==	114)) return true;

        // if running on the proton files from Erez____________________________________________________________//
        bool okevent = false;//________________________________________________________________________________//
        for(auto trackinfo:_SelectableTracks){//_______________________________________________________________//
            if(_run == trackinfo[0] && _subrun == trackinfo[1] && _event == trackinfo[2]){okevent = true; break;}//___//
        }//____________________________________________________________________________________________________//
        if(!okevent) return true;//____________________________________________________________________________//
        //
        //

        // Get associated hits + association info
        larlite::event_hit* ev_hit=nullptr;
        auto const& track_to_hit = storage->find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
        if(!ev_hit) throw DataFormatException("Could not find associated hit data product!");

        //
        // Loop over each proto-track and call reconstruction algorithm
        // 0) retrieve 3D start/end points
        // 1) construct hit image2d
        // 2) construct chstatus image2d
        // 3) call reconstruct function
        //


        for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
            _track = track_index;
            if( !(_SelectableTracks.size()!=0 && IsGoodTrack()) ) continue;
            tellMe(Form("%d\t%d\t%d\t%d",_run,_subrun,_event,_track),0);
            //if((_run == 5161 && _subrun ==	18 && _event ==	937  && _track == 1)) continue;

            if (track_to_hit[track_index].size() == 0){
                tellMe("no hits associated to this track, stop here for now",0);
                continue;
            }

            //
            // 0) retrieve 3D start/end points
            //

            auto const& track = (*ev_track)[track_index];
            start_pt = track.Vertex();
            end_pt   = track.End();
            int thisTrackID = track.ID();

            bool endpointOutOfRange = false;
            if(start_pt.X() < 0    || end_pt.X() < 0   ){endpointOutOfRange = true;}
            if(start_pt.X() > 260  || end_pt.X() > 260 ){endpointOutOfRange = true;}

            if(start_pt.Y() < -115 || end_pt.Y() < -115){endpointOutOfRange = true;}
            if(start_pt.Y() >  115 || end_pt.Y() >  115){endpointOutOfRange = true;}

            if(start_pt.Z() < 0    || end_pt.Z() < 0   ){endpointOutOfRange = true;}
            if(start_pt.Z() > 1036 || end_pt.Z() > 1036){endpointOutOfRange = true;}

            if(endpointOutOfRange){tellMe("Point out of range",0); continue;}

            //
            // 0-bis) find corresponding MC track
            if(ev_mct){
                larlite::mctrack trueTrack;
                for(auto const& mct : *ev_mct) {
                    double DRstart = sqrt(pow(start_pt.X()-mct.Start().X(),2)+pow(start_pt.Y()-mct.Start().Y(),2)+pow(start_pt.Z()-mct.Start().Z(),2));
                    double DRend = sqrt(pow(end_pt.X()-mct.End().X(),2)+pow(end_pt.Y()-mct.End().Y(),2)+pow(end_pt.Z()-mct.End().Z(),2));
                    if(DRstart > 10 || DRend > 10)continue;
                    trueTrack = mct;
                }
                //if(trueTrack.size() == 0){std::cout << "Could not find MC track corresponding to provided start and end points => stop" << std::endl;continue;}
            }
            tellMe("11111",1);


            //
            // 0-ter) get vector of hits and channel status to create image2D in next step
            //

            std::vector<larlite::hit> trackHit_v;
            for(auto const& hit_index : track_to_hit[track_index]) {
                auto const& h = (*ev_hit)[hit_index];
                trackHit_v.push_back(h);
            }

            tellMe("22222",1);

            std::vector< std::vector< std::pair<float,float> > > chstatus_vector;
            if(ev_chstatus){
                for(int iPlane = 0;iPlane<3;iPlane++){
                    std::vector< std::pair<float,float> > chstatus_vector_plane;
                    for(auto const& chstatus_v : *ev_chstatus) {
                        if( chstatus_v.plane().Plane != iPlane) continue;

                        for(size_t ch=0; ch<chstatus_v.status().size(); ++ch) {
                            std::pair<float, float> pairch;
                            auto const& status = chstatus_v.status()[ch];
                            pairch.first = ch;
                            if(status == larlite::larch::kGOOD){pairch.second = 0;}
                            else {pairch.second = 1;}
                            chstatus_vector_plane.push_back(pairch);
                        }
                        break;
                    }
                    chstatus_vector.push_back(chstatus_vector_plane);
                }
            }

            tellMe("33333",1);

            //
            // 1) construct hit image2d ... start from figuring out the range
            // 2) construct chstatus image2d ... do this in the same loop
            //

            bool imageOK = false;
            std::pair< std::vector<larcv::Image2D>, std::vector<larcv::Image2D> > imagesPair = CreateImages(trackHit_v,chstatus_vector,imageOK);
            if(imageOK  == false){tellMe("imageOK is false",1); continue;}

            std::vector<larcv::Image2D> hit_image_v      = imagesPair.first;
            std::vector<larcv::Image2D> chstatus_image_v = imagesPair.second;

            tellMe("44444",1);

            //
            // 3) call reconstruct! and set the return track ovject to override the proto-track
            //
            //std::cout << "about to call Reconstruct()" << std::endl;
            _eventTreated++ ;
            (*ev_track)[track_index] = Reconstruct(hit_image_v, chstatus_image_v);
            //std::cout << "Reconstruct called and executed" << std::endl;
            (*ev_track)[track_index].set_track_id(thisTrackID);

            tellMe("calling CompareReco2hits",1);
            CompareReco2hits((*ev_track)[track_index],hit_image_v);

            std::cout << std::endl << std::endl;
        }
        return true;
    }

    double AStarTracker::X2Tick(double x, size_t plane) const {

        auto ts = larutil::TimeService::GetME();
        auto larp = larutil::LArProperties::GetME();

        // (X [cm] / Drift Velocity [cm/us] - TPC waveform tick 0 offset) ... [us]
        double tick = (x / larp->DriftVelocity() - ts->TriggerOffsetTPC() - _speedOffset);
        // 1st plane's tick
        /*if(plane==0) return tick * 2;// 2 ticks/us
        // 2nd plane needs more drift time for plane0=>plane1 gap (0.3cm) difference
        tick -= 0.3 / larp->DriftVelocity(larp->Efield(1));
        // 2nd plane's tick
        if(plane==1) return tick * 2;
        // 3rd plane needs more drift time for plane1=>plane2 gap (0.3cm) difference
        tick -= 0.3 / larp->DriftVelocity(larp->Efield(2));*/
        return tick * 2;
    }

    double AStarTracker::Tick2X(double tick, size_t plane) const{
        //std::cout << "tick : " << tick;
        auto ts = larutil::TimeService::GetME();
        auto larp = larutil::LArProperties::GetME();
        // remove tick offset due to plane
        /*if(plane >= 1) {tick += 0.3/larp->DriftVelocity(larp->Efield(1));}
        if(plane >= 2) {
            tick += 0.3/larp->DriftVelocity(larp->Efield(2));
        }*/
        tick/=2.;
        double x = (tick+ts->TriggerOffsetTPC()+_speedOffset)*larp->DriftVelocity();
        //std::cout << "  =>  x : " << x << "  => new tick : " << X2Tick(x,plane) << std::endl;

        return x;
    }


    larlite::track AStarTracker::Reconstruct(const std::vector<larcv::Image2D>& hit_image_v, const std::vector<larcv::Image2D>& chstatus_image_v){
        tellMe("creating tag image",1);
        std::vector<larcv::Image2D> tag_image_v(3);
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            tag_image_v[iPlane] = hit_image_v[iPlane];
            tag_image_v[iPlane].paint(0);
        }

        RecoedPath.clear();
        int goal_reached = 0;

        //_______________________________________________
        // Get Start and end points
        //-----------------------------------------------
        tellMe("getting start and end points",1);
        std::vector<int> start_cols(3);
        std::vector<int> end_cols(3);
        int start_row, end_row;

        if(hit_image_v[1].meta().height()>0) {
            start_row = hit_image_v[1].meta().row(X2Tick(start_pt.X(),1));
            end_row = hit_image_v[1].meta().row(X2Tick(end_pt.X(),1));
        }
        else if(hit_image_v[0].meta().height()>0) {
            start_row = hit_image_v[0].meta().row(X2Tick(start_pt.X(),0));
            end_row = hit_image_v[0].meta().row(X2Tick(end_pt.X(),0));
        }
        else if(hit_image_v[2].meta().height()>0) {
            start_row = hit_image_v[2].meta().row(X2Tick(start_pt.X(),2));
            end_row = hit_image_v[2].meta().row(X2Tick(end_pt.X(),2));
        }
        else{tellMe("error: 0 rows",0); return larlite::track();}

        //tellMe("ok here");

        for(size_t iPlane=0;iPlane<3;iPlane++){

            if(hit_image_v[iPlane].meta().width()==0)continue;

            double wireProjStartPt = larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w / 0.3;
            double wireProjEndPt   = larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w / 0.3;

            //if(iPlane == 1 && wireProjStartPt >= 2400){wireProjStartPt = 2399;std::cout << "had to move the start point wire projection" << std::endl;}
            //if(iPlane == 1 && wireProjEndPt   >= 2400){wireProjEndPt   = 2399;std::cout << "had to move the end point wire projection"   << std::endl;}

            //std::cout << iPlane << "\t" << wireProjStartPt << "\t" << wireProjEndPt << std::endl;

            if(wireProjStartPt < 0) wireProjStartPt = 0;
            if(wireProjEndPt   < 0) wireProjEndPt   = 0;

            start_cols[iPlane] = hit_image_v[iPlane].meta().col(wireProjStartPt);
            end_cols[iPlane]   = hit_image_v[iPlane].meta().col(wireProjEndPt);
        }

        //std::cout << "start and end rows OK" << std::endl;

        //_______________________________________________
        // Configure A* algo
        //-----------------------------------------------
        larlitecv::AStar3DAlgoConfig config;
        config.accept_badch_nodes = true;
        config.astar_threshold.resize(3,0);
        config.astar_neighborhood.resize(3,5); //can jump over n empty pixels
        config.astar_start_padding = 10;// allowed region around the start point
        config.astar_end_padding = 10;  // allowed region around the end point
        config.lattice_padding = 5; // margin around the edges
        config.min_nplanes_w_hitpixel = 2;
        config.restrict_path = false; // do I want to restrict to a cylinder around the strainght line of radius defined bellow
        config.path_restriction_radius = 30.0;

        //_______________________________________________
        // Define A* algo
        //-----------------------------------------------
        larlitecv::AStar3DAlgoProton algo( config );
        algo.setVerbose(0);
        algo.setPixelValueEsitmation(true);

        tellMe("starting A*",1);
        RecoedPath = algo.findpath( hit_image_v, chstatus_image_v, tag_image_v, start_row, end_row, start_cols, end_cols, goal_reached );

        if(goal_reached == 1){
            _eventSuccess++;
        }

        //_______________________________________________
        // Make track out of 3D nodes
        //-----------------------------------------------
        tellMe("making newTrack",1);
        larlite::track newTrack = MakeTrack();

        //_______________________________________________
        // Compute dQdx for newTrack
        //-----------------------------------------------
        tellMe("calling ComputedQdX()",1);
        newTrack = ComputedQdX(newTrack,hit_image_v,chstatus_image_v);

        return newTrack;
    }

    void AStarTracker::DrawdQdX(larlite::track thisTrack){
        for(size_t iNode = 0;iNode<thisTrack.NumberTrajectoryPoints()-1;iNode++){
            double dqdx =   thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)
                            +thisTrack.DQdxAtPoint(iNode,larlite::geo::kU)
                            +thisTrack.DQdxAtPoint(iNode,larlite::geo::kV);

            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) != 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)!=0) dqdx*=3./2.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ) != 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)!=0) dqdx*=3./2.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ) != 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) != 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)==0) dqdx*=3./2.;

            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)!=0) dqdx*=3.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ) != 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)==0) dqdx*=3.;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ) == 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) != 0 && thisTrack.DQdxAtPoint(iNode,larlite::geo::kV)==0) dqdx*=3.;

            hdQdx->Fill(thisTrack.Length(iNode),dqdx);
            hdQdxEntries->Fill(thisTrack.Length(iNode));
        }
    }

    void AStarTracker::CompareReco2MC3D(const larlite::track recoTrackRaw, const larlite::mctrack trueTrack){
        //larlite::track recoTrack = CorrectSCE(recoTrackRaw);
        larlite::track recoTrack = recoTrackRaw;
        tellMe("comparing Reco and MC tracks",1);
        double trueLength=0;
        for(size_t iNode = 0;iNode<trueTrack.size()-1;iNode++){
            trueLength+=sqrt(pow(trueTrack[iNode+1].X()-trueTrack[iNode].X(),2)+pow(trueTrack[iNode+1].Y()-trueTrack[iNode].Y(),2)+pow(trueTrack[iNode+1].Z()-trueTrack[iNode].Z(),2));
        }
        tellMe(Form("reconstructed length : %.1f [cm]",recoTrack.Length(0)),1);
        tellMe(Form("true length          : %.1f [cm]",trueLength),1);

        tellMe("Reco'ed nodes",1);
        for(size_t iNode=0; iNode<recoTrack.NumberTrajectoryPoints();iNode++){
            tellMe(Form("   node (X,Y,Z) : %3.2f \t %3.2f \t %3.2f",recoTrack.LocationAtPoint(iNode).X(),recoTrack.LocationAtPoint(iNode).Y(),recoTrack.LocationAtPoint(iNode).Z()),1);
        }
        std::cout << std::endl;
        tellMe("MC nodes",1);
        for(auto const& step : trueTrack) {
            tellMe(Form("MC node (X,Y,Z) : %3.2f \t %3.2f \t %3.2f",step.X(),step.Y(),step.Z()),1);
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
            hDistance2MCX->Fill(dPointSegmentX);
            hDistance2MCY->Fill(dPointSegmentY);
            hDistance2MCZ->Fill(dPointSegmentZ);
            if(dPointSegment > 20){tellMe("STRANGE, TOO FAR FROM MC",1);}
        }

    }

    void AStarTracker::CompareReco2hits(const larlite::track recoTrack, const std::vector<larcv::Image2D>& hit_image_v){

        for(size_t iNode = 0; iNode<recoTrack.NumberTrajectoryPoints();iNode++){
            //std::cout << iNode << std::endl;
            //TVector3 this_point = recoTrack.LocationAtPoint(iNode);
            double dist=1e9;
            for(int iPlane=0;iPlane<3;iPlane++){
                //std::cout << iPlane << std::endl;
                double thiscol = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(recoTrack.LocationAtPoint(iNode),iPlane).w/0.3);
                double thisrow = hit_image_v[iPlane].meta().row(X2Tick(recoTrack.LocationAtPoint(iNode).X(),iPlane));
                //std::cout << thisrow << "\t" << thiscol << std::endl;
                //std::cout << "reading images for this point" << std::endl;
                for(size_t icol = 0;icol<hit_image_v[iPlane].meta().cols();icol++){
                    for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                        //std::cout << icol << "/" << hit_image_v[iPlane].meta().cols() << "\t" << irow << "/" << hit_image_v[iPlane].meta().rows() << std::endl;
                        if(hit_image_v[iPlane].pixel(irow,icol)<=1)continue;
                        //std::cout << "about to evaluate thisdist" << std::endl;
                        double thisdist = sqrt(pow(_compressionFactor_t*(thisrow-irow),2)+pow(_compressionFactor_w*(thiscol-icol),2));
                        //std::cout << "thisdist = "<< thisdist << std::endl;
                        if(dist > thisdist)dist = thisdist;
                    }
                }
            }
            //std::cout << "filling hDistance2Hit" << std::endl;
            hDistance2Hit->Fill(dist);
        }
    }

    bool AStarTracker::finalize() {
        tellMe("Finalizing",0);
        tellMe(Form("%d / %d = %.2f",_eventSuccess,_eventTreated,_eventSuccess*100./_eventTreated),0);
        TCanvas *cDist = new TCanvas();
        hDistance2Hit->Draw();
        cDist->SetLogy();
        //hDistance2Hit->Write();
        tellMe(Form("hDistance2Hit mean: %.2f , RMS : %.2f",hDistance2Hit->GetMean(),hDistance2Hit->GetRMS()),0);
        tellMe(Form("entries : %.1f",hDistance2Hit->GetEntries()),0);
        cDist->SaveAs("DistanceReco2Hit.pdf");
        return true;
    }

    larlite::track AStarTracker::MakeTrack(){
        larlite::track newTrack;
        double nodeX,nodeY,nodeZ,nodeDirX,nodeDirY,nodeDirZ,norm;

        for(size_t iNode=RecoedPath.size()-1;iNode < RecoedPath.size();iNode--){
            //std::cout << iNode << " / " << RecoedPath.size() << std::endl;

            nodeX = Tick2X(RecoedPath.at(iNode).tyz[0],0);
            nodeY = RecoedPath.at(iNode).tyz[1];
            nodeZ = RecoedPath.at(iNode).tyz[2];
            TVector3 node(nodeX,nodeY,nodeZ);

            if(iNode < RecoedPath.size()-1){
                nodeDirX =  Tick2X(RecoedPath.at(iNode+1).tyz[0],0)-Tick2X(RecoedPath.at(iNode).tyz[0],0);
                nodeDirY =  RecoedPath.at(iNode+1).tyz[1]-RecoedPath.at(iNode).tyz[1];
                nodeDirZ =  RecoedPath.at(iNode+1).tyz[2]-RecoedPath.at(iNode).tyz[2];
                norm = sqrt(nodeDirX*nodeDirX+nodeDirY*nodeDirY+nodeDirZ*nodeDirZ);
                TVector3 nodeDir(nodeDirX/norm,nodeDirY/norm,nodeDirZ/norm);
                newTrack.add_vertex(node);
                newTrack.add_direction(nodeDir);
            }
            else{
                nodeDirX =  Tick2X(RecoedPath.at(iNode).tyz[0],0)-Tick2X(RecoedPath.at(0).tyz[0],0);
                nodeDirY =  RecoedPath.at(iNode).tyz[1]-RecoedPath.at(0).tyz[1];
                nodeDirZ =  RecoedPath.at(iNode).tyz[2]-RecoedPath.at(0).tyz[2];
                norm = sqrt(nodeDirX*nodeDirX+nodeDirY*nodeDirY+nodeDirZ*nodeDirZ);
                TVector3 nodeDir(nodeDirX/norm,nodeDirY/norm,nodeDirZ/norm);
                newTrack.add_vertex(node);
                newTrack.add_direction(nodeDir);
            }
        }
        return newTrack;
    }

    larlite::track AStarTracker::ComputedQdX(larlite::track newTrack, const std::vector<larcv::Image2D>& hit_image_v, const std::vector<larcv::Image2D>& chstatus_image_v){
        double X0,Y0,X1,Y1,X2,Y2;
        std::vector<double> planedQdx;
        double Npx = 5;
        double dist,xp,yp,alpha;
        double localdQdx;
        double nodeDirX,nodeDirY,nodeDirZ,norm;
        TH2D *hImage[3];
        TH2D *hStatus[3];
        TH2D *hIntegral[3];
        TGraph *gTrack[3];
        TGraph *gStartNend[3];

        for(size_t iPlane=0;iPlane<3;iPlane++){
            //std::cout << "\t" << iPlane;
            gTrack[iPlane] = new TGraph();
            hImage[iPlane] = new TH2D(Form("hImage_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hImage_%d_%d_%d_%d_%zu;wire;time",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());
            hStatus[iPlane] = new TH2D(Form("hStatus_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hStatus_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());
            hIntegral[iPlane] = new TH2D(Form("hIntegral_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hIntegral_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());
            //std::cout << "\t histos for plane defined ";
            for(int icol=0;icol<hit_image_v[iPlane].meta().cols();icol++){
                for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                    hImage[iPlane]->SetBinContent(icol+1,irow+1,hit_image_v[iPlane].pixel(irow,icol));
                    hStatus[iPlane]->SetBinContent(icol+1,irow+1,chstatus_image_v[iPlane].pixel(irow,icol));
                }
            }
            //std::cout << " /  hImage filled";
            //std::cout << " / X0 : " << start_pt.X() << " ";
            //std::cout << " / Y0 : " << start_pt.Y() << " ";
            //std::cout << " / Z0 : " << start_pt.Z() << " ";
            //std::cout << " / Xf : " << end_pt.X() << " ";
            //std::cout << " / Yf : " << end_pt.Y() << " ";
            //std::cout << " / Zf : " << end_pt.Z() << " ";

            gStartNend[iPlane] = new TGraph();
            gStartNend[iPlane]->SetPoint(0,hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w/0.3),hit_image_v[iPlane].meta().row(X2Tick(start_pt.X(),iPlane)));
            //std::cout << " /  gStartNend filled first point ";

            double this_wire = larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w/0.3;
            //std::cout << " / wire : " << this_wire << " ";
            double this_end_col = hit_image_v[iPlane].meta().col(this_wire);
            //std::cout << " /  col : " << this_end_col << " ";
            double this_end_row = hit_image_v[iPlane].meta().row(X2Tick(end_pt.X(),iPlane));
            //std::cout << " /  row : " << this_end_row << " ";
            gStartNend[iPlane]->SetPoint(1,this_end_col,this_end_row);
            //std::cout << " /  gStartNend filled last point " << std::endl;
        }

        tellMe("computing dqdx",1);
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            //std::cout << "\t" << iPlane << std::endl;

            for(size_t iNode=0;iNode<newTrack.NumberTrajectoryPoints()-1;iNode++){
                localdQdx=0;
                X1 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(iNode  ),iPlane).w/0.3);
                Y1 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(iNode  ).X(),iPlane));

                X2 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(iNode+1),iPlane).w/0.3);
                Y2 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(iNode+1).X(),iPlane));

                gTrack[iPlane]->SetPoint(iNode,X1,Y1);
                gTrack[iPlane]->SetPoint(iNode+1,X2,Y2);

                if(X1 == X2){X1 -= 0.001;}
                if(Y1 == Y2){Y1 -= 0.001;}

                //std::cout << "\t \t" << X1 << "," << Y1 << "   " << X2 << "," << Y2 << std::endl;

                for(int icol=0;icol<hit_image_v[iPlane].meta().cols();icol++){
                    for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                        //if(hit_image_v[iPlane].pixel(irow,icol) == 0)continue;
                        //----------------------------
                        // first, check if the point can belong to this segment
                        // measure the distance to the line
                        // check the angles
                        //----------------------------
                        double theta0; // angle with which the segment [iNode, iNode+1] is ssen by the current point
                        double theta1; // angle with which the segment [jNode, jNode+1] is seen by the current point
                        X0 = icol+0.5;
                        Y0 = irow+0.5;
                        alpha = ( (X2-X1)*(X0-X1)+(Y2-Y1)*(Y0-Y1) )/( pow(X2-X1,2)+pow(Y2-Y1,2) );
                        if(alpha < -0.1 || alpha > 1.1) continue;
                        xp = X1+alpha*(X2-X1);
                        yp = Y1+alpha*(Y2-Y1);
                        dist = sqrt(pow(xp-X0,2)+pow(yp-Y0,2));
                        if(dist > Npx)continue;
                        theta0 = std::acos( ((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0)) / (sqrt(pow(X1-X0,2)+pow(Y1-Y0,2))*sqrt(pow(X2-X0,2)+pow(Y2-Y0,2))) );
                        // check all other nodes, if I find an angle bigger with a distance lower than Npx, reject the point for the current segment

                        bool bestCandidate = true;
                        for(size_t jNode=0;jNode<newTrack.NumberTrajectoryPoints()-1;jNode++){
                            if(jNode==iNode)continue;
                            double x1,y1,x2,y2;
                            x1 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(jNode  ),iPlane).w/0.3);
                            x2 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(jNode+1),iPlane).w/0.3);

                            y1 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(jNode  ).X(),iPlane));
                            y2 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(jNode+1).X(),iPlane));

                            theta1 = std::acos( ((x1-X0)*(x2-X0)+(y1-Y0)*(y2-Y0)) / (sqrt(pow(x1-X0,2)+pow(y1-Y0,2))*sqrt(pow(x2-X0,2)+pow(y2-Y0,2))) );
                            if(theta1 > theta0){bestCandidate=false;}
                        }
                        if(!bestCandidate)continue;


                        nodeDirX =  newTrack.LocationAtPoint(iNode+1).X()-newTrack.LocationAtPoint(iNode).X();//Tick2X(RecoedPath.at(iNode-1).tyz[0],0)-Tick2X(RecoedPath.at(iNode).tyz[0],0);
                        nodeDirY =  newTrack.LocationAtPoint(iNode+1).Y()-newTrack.LocationAtPoint(iNode).Y();;//RecoedPath.at(iNode-1).tyz[1]-RecoedPath.at(iNode).tyz[1];
                        nodeDirZ =  newTrack.LocationAtPoint(iNode+1).Z()-newTrack.LocationAtPoint(iNode).Z();;//RecoedPath.at(iNode-1).tyz[2]-RecoedPath.at(iNode).tyz[2];
                        norm = sqrt(nodeDirX*nodeDirX+nodeDirY*nodeDirY+nodeDirZ*nodeDirZ);
                        localdQdx+=hit_image_v[iPlane].pixel(irow,icol)/norm;
                        hIntegral[iPlane]->Fill(X0,Y0);
                    }
                }
                planedQdx.push_back(localdQdx);
            }
            newTrack.add_dqdx(planedQdx);
            planedQdx.clear();
        }
        TCanvas *c = new TCanvas(Form("c_%d_%d_%d_%d",_run,_subrun,_event,_track),Form("c_%d_%d_%d_%d",_run,_subrun,_event,_track),1200,300);
        c->Divide(3,1);
        for(int iPlane=0;iPlane<3;iPlane++){
            c->cd(iPlane+1);
            hImage[iPlane]->Draw("colz");
            //hStatus[iPlane]->Draw("same colz");
            hIntegral[iPlane]->Draw("same");
            hImage[iPlane]->Draw("same colz");
            gTrack[iPlane]->SetLineColor(2);
            gTrack[iPlane]->Draw("same LP");
            gStartNend[iPlane]->SetMarkerStyle(20);
            gStartNend[iPlane]->Draw("same P");
        }
        c->SaveAs(Form("%s.pdf",c->GetName()));
        for(int iPlane = 0;iPlane<3;iPlane++){
            hImage[iPlane]->Delete();
            hStatus[iPlane]->Delete();
            hIntegral[iPlane]->Delete();
            gTrack[iPlane]->Delete();
            gStartNend[iPlane]->Delete();
        }
        //c->Delete();
        return newTrack;
    }

    std::pair< std::vector<larcv::Image2D>, std::vector<larcv::Image2D> > AStarTracker::CreateImages(std::vector<larlite::hit> trackHit_v, std::vector< std::vector< std::pair<float,float> > > chstatus_vector, bool &imageOK){

        const size_t num_planes = 3;
        std::pair< std::vector<larcv::Image2D>, std::vector<larcv::Image2D> > imagesPair;
        std::vector<larcv::Image2D> hit_image_v;
        std::vector<larcv::Image2D> chstatus_image_v;
        imagesPair.first = hit_image_v;
        imagesPair.second = chstatus_image_v;

        std::vector<std::pair<double,double> > time_bounds(num_planes); // first=min, second=max
        std::vector<std::pair<double,double> > wire_bounds(num_planes); // first=min, second=max

        // Initialize range limits
        for(size_t plane=0; plane<num_planes; ++plane) {
            time_bounds[plane].first  = 1.e9;
            time_bounds[plane].second = 0.;
            wire_bounds[plane].first  = 1.e9;
            wire_bounds[plane].second = 0.;
        }

        // Figure out range limits based on hits
        /*for(auto const& h : trackHit_v) {
         // make sure current hit is within the range
         if(h.PeakTime() < time_bounds[h.WireID().Plane].first ) time_bounds[h.WireID().Plane].first  = h.PeakTime();
         if(h.PeakTime() > time_bounds[h.WireID().Plane].second) time_bounds[h.WireID().Plane].second = h.PeakTime();
         if(h.WireID().Wire < wire_bounds[h.WireID().Plane].first ) wire_bounds[h.WireID().Plane].first  = h.WireID().Wire;
         if(h.WireID().Wire > wire_bounds[h.WireID().Plane].second) wire_bounds[h.WireID().Plane].second = h.WireID().Wire;
         }*/

        for(size_t iPlane=0;iPlane<3;iPlane++){
            // make sure starting point is projected back in the range
            if(X2Tick(start_pt.X(),iPlane) < time_bounds[iPlane].first ) time_bounds[iPlane].first  = X2Tick(start_pt.X(),iPlane);
            if(X2Tick(start_pt.X(),iPlane) > time_bounds[iPlane].second) time_bounds[iPlane].second = X2Tick(start_pt.X(),iPlane);
            double wireProjStartPt = larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w / 0.3;
            //std::cout << "before any resizing : z_0(" << iPlane << ") : " << wireProjStartPt << std::endl;
            if(wireProjStartPt < wire_bounds[iPlane].first ) wire_bounds[iPlane].first  = wireProjStartPt;
            if(wireProjStartPt > wire_bounds[iPlane].second) wire_bounds[iPlane].second = wireProjStartPt;

            // make sure ending point is projected back in the range
            if(X2Tick(end_pt.X(),iPlane) < time_bounds[iPlane].first ) time_bounds[iPlane].first  = X2Tick(end_pt.X(),iPlane);
            if(X2Tick(end_pt.X(),iPlane) > time_bounds[iPlane].second) time_bounds[iPlane].second = X2Tick(end_pt.X(),iPlane);
            double wireProjEndPt   = larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w / 0.3;
            //std::cout << "before any resizing : z_1(" << iPlane << ") : " << wireProjEndPt << std::endl;
            if((larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w / 0.3) < wire_bounds[iPlane].first ) wire_bounds[iPlane].first  = wireProjEndPt;
            if((larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w / 0.3) > wire_bounds[iPlane].second) wire_bounds[iPlane].second = wireProjEndPt;
        }
        // equalize time ranges to intercept the same range on 3 views
        for(size_t iPlane=0;iPlane<3;iPlane++){
            for(size_t jPlane=0;jPlane<3;jPlane++){
                if(time_bounds[jPlane].first  <= time_bounds[iPlane].first) {time_bounds[iPlane].first  = time_bounds[jPlane].first;}
                if(time_bounds[jPlane].second >= time_bounds[iPlane].second){time_bounds[iPlane].second = time_bounds[jPlane].second;}
            }
        }
        // Update the range with margin
        double ImageMargin = 10.;
        for(size_t iPlane=0;iPlane<3;iPlane++){
            time_bounds[iPlane].first  -= std::max(ImageMargin,0.2*(time_bounds[iPlane].second-time_bounds[iPlane].first));
            time_bounds[iPlane].second += std::max(ImageMargin,0.2*(time_bounds[iPlane].second-time_bounds[iPlane].first));
            wire_bounds[iPlane].first  -= std::max(ImageMargin,0.2*(wire_bounds[iPlane].second-wire_bounds[iPlane].first));
            wire_bounds[iPlane].second += std::max(ImageMargin,0.2*(wire_bounds[iPlane].second-wire_bounds[iPlane].first));

            if(time_bounds[iPlane].first < 0) time_bounds[iPlane].first = 0;
            if(wire_bounds[iPlane].first < 0) wire_bounds[iPlane].first = 0;
            //if(iPlane < 2 && wire_bounds[iPlane].second >= 2400) wire_bounds[iPlane].second = 2399;
            //if(iPlane == 2 && wire_bounds[iPlane].second >= 3456) wire_bounds[iPlane].second = 3455;

            wire_bounds[iPlane].first  = (size_t)(wire_bounds[iPlane].first + 0.5);
            wire_bounds[iPlane].second = (size_t)(wire_bounds[iPlane].second + 0.5);
        }

        // make sure the number of rows and cols are divisible by _compressionFactor




        for(size_t iPlane=0;iPlane<3;iPlane++){
            //std::cout << "can the image be compressed in time ? => " << (size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first) << " % " << _compressionFactor_t << " = " << (size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first)%_compressionFactor_t << std::endl;
            while(!( (size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first)%_compressionFactor_t == 0)){
                //std::cout << "not compressible in time, adding " << (_compressionFactor_t-((size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first)%_compressionFactor_t)) << std::endl;
                time_bounds[iPlane].second+=(_compressionFactor_t-((size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first)%_compressionFactor_t));
            }
        }

        for(int iPlane = 0;iPlane<3;iPlane++){
            //std::cout << "can the image be compressed in wires ? => " << (size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first) << " % " << _compressionFactor_w << " = " << (size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first)%_compressionFactor_w << std::endl;

            while(!( (size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first)%_compressionFactor_w == 0)){
                //std::cout << "not compressible in wires, adding " << (_compressionFactor_w-((size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first)%_compressionFactor_w)) << std::endl;
                wire_bounds[iPlane].second+=(_compressionFactor_w-((size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first)%_compressionFactor_w));
            }
        }

        // check if any associated hit is within the limits
        bool HitsPresent = false;
        int nHits=0;
        for(auto const& h : trackHit_v){
            auto& time_bound = time_bounds[h.WireID().Plane];
            auto& wire_bound = wire_bounds[h.WireID().Plane];
            // make sure current hit is within the range
            if(h.PeakTime() < time_bound.first )continue;
            if(h.PeakTime() > time_bound.second)continue;
            if(h.WireID().Wire < wire_bound.first ) continue;
            if(h.WireID().Wire > wire_bound.second) continue;
            HitsPresent=true;
            nHits++;
        }
        if(HitsPresent == false){std::cout << "no hits in Image2D" << std::endl; imageOK=false; return imagesPair;}
        //else{std::cout << "found " << nHits << " hist in the image" << std::endl;}

        // Using the range, construct Image2D

        for(size_t plane=0; plane<num_planes; ++plane) {
            // Retrieve boundaries
            auto const& time_bound = time_bounds[plane];
            auto const& wire_bound = wire_bounds[plane];
            // If no hit on this plane, make an empty image
            if(wire_bound.second <= wire_bound.first ||
               time_bound.second <= time_bound.first ) {
                hit_image_v.push_back(larcv::Image2D());
                chstatus_image_v.push_back(larcv::Image2D());
                continue;
            }
            // Construct meta
            //std::cout << "construct meta" << std::endl;
            size_t image_cols = (size_t)(wire_bound.second - wire_bound.first/* + 1*/);
            size_t image_rows = (size_t)(time_bound.second - time_bound.first/* + 1*/);
            larcv::ImageMeta hit_meta((double)image_cols, (double)image_rows,
                                      image_rows, image_cols,
                                      (size_t) wire_bound.first,  // origin x = min x
                                      (size_t) time_bound.second, // origin y = max time
                                      plane);
            auto chstatus_meta = hit_meta;

            //std::cout << "hit_meta      size (row,col) : " << hit_meta.rows()      << "\t" << hit_meta.cols()      << std::endl;
            //std::cout << "chstatus_meta size (row,col) : " << chstatus_meta.rows() << "\t" << chstatus_meta.cols() << std::endl;


            // Prepare hit image data + fill
            //std::cout << "Prepare hit image data + fill" << std::endl;

            std::vector<float> image_data(hit_meta.rows() * hit_meta.cols(), 0.);
            size_t row,col;
            for(auto const& h : trackHit_v){
                if(h.WireID().Plane != plane) continue;
                if(h.PeakTime()    < time_bound.first || h.PeakTime()    >= time_bound.second) continue;
                if(h.WireID().Wire < wire_bound.first || h.WireID().Wire >= wire_bound.second) continue;
                row = (size_t)(time_bound.second - h.PeakTime() + 0.5);
                col = (size_t)(h.WireID().Wire - wire_bound.first + 0.5);
                image_data[hit_meta.rows() * col + row] = h.Integral();
            }
            larcv::Image2D hit_image(std::move(hit_meta),std::move(image_data));
            //std::cout << "after filling hit_image " << hit_image.meta().rows() << "\t" << hit_image.meta().cols() << std::endl;
            hit_image_v.emplace_back(std::move(hit_image));

            ///////////////////////
            std::vector<float> chstatus_data(chstatus_meta.rows() * chstatus_meta.cols(), 0.);
            if(chstatus_vector.size()!=0){
                for(int i=0;i<chstatus_vector[plane].size();i++) {
                    auto const& status = chstatus_vector[plane][i].second;
                    auto const& ch = chstatus_vector[plane][i].first;
                    if(ch < wire_bound.first || ch >= wire_bound.second) continue;
                    for(size_t row=0; row<chstatus_meta.rows(); ++row){
                        chstatus_data[(ch - wire_bound.first) * chstatus_meta.rows() + row] = status;
                    }
                    break;
                }
            }
            larcv::Image2D chstatus_image(std::move(chstatus_meta),std::move(chstatus_data));
            //std::cout << "after filling chstatus_image " << chstatus_image.meta().rows() << "\t" << chstatus_image.meta().cols() << std::endl;
            ///////////////////////

            //
            // when running with no chstatus
            //
            /*for(int icol=0;icol<hit_image_v[plane].meta().cols();icol++){
             bool isBadCh = true;
             for(int irow=0;irow<hit_image_v[plane].meta().rows();irow++){
             if(hit_image_v[plane].pixel(irow,icol)!=0){isBadCh=false;}
             }
             if(isBadCh==true){
             for(int irow=0;irow<hit_image_v[plane].meta().rows();irow++){
             chstatus_image.set_pixel(irow,icol,1);
             }
             }
             }*/
            //
            // end of running with no chstatus
            //
            chstatus_image_v.emplace_back(std::move(chstatus_image));
            //std::cout << "chstatus_image done" << std::endl;
        }

        // compress Images
        tellMe("compress Images",1);
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            tellMe(Form("%zu size : %zu \t %zu",iPlane,hit_image_v.at(iPlane).meta().rows(),hit_image_v.at(iPlane).meta().cols()),1);

            hit_image_v.at(iPlane).compress(hit_image_v.at(iPlane).meta().rows()/_compressionFactor_t,
                                            hit_image_v.at(iPlane).meta().cols()/_compressionFactor_w);

            chstatus_image_v.at(iPlane).compress(chstatus_image_v.at(iPlane).meta().rows()/_compressionFactor_t,
                                                 chstatus_image_v.at(iPlane).meta().cols()/_compressionFactor_w);
        }
        tellMe("image compressed",1);

        imagesPair.first = hit_image_v;
        imagesPair.second = chstatus_image_v;
        imageOK = true;
        return imagesPair;
    }

    //
    // Correct for Space-Charge Effects
    //
    larlite::track AStarTracker::CorrectSCE(larlite::track thisTrack){
        tellMe("correct SCE",1);
        // larlite::track => std::vector<TVector3>
        std::vector<TVector3> newPath;
        TVector3 node;
        for(size_t iNode = 0;iNode<thisTrack.NumberTrajectoryPoints();iNode++){
            node = thisTrack.LocationAtPoint(iNode);
            newPath.push_back(node);
        }

        // correct for SCE
        std::vector<TVector3> newNodes = CorrectSCE(newPath);

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
    std::vector<TVector3> AStarTracker::CorrectSCE(larlite::mctrack thisTrack){
        std::vector<TVector3> newPath;
        for(size_t iNode = 0;iNode<thisTrack.size();iNode++){
            TVector3 node(thisTrack[iNode].X(),thisTrack[iNode].Y(),thisTrack[iNode].Z());
            newPath.push_back(node);
        }
        return CorrectSCE(newPath);
    }
    std::vector<TVector3> AStarTracker::CorrectSCE(std::vector<TVector3> thisTrack){
        std::vector<TVector3> newPath;
        larlitecv::SpaceChargeMicroBooNE mySCEcorr;
        for(size_t iNode=0;iNode<thisTrack.size();iNode++){
            std::vector<double> PosCorr = mySCEcorr.GetPosOffsets(thisTrack.at(iNode).X(),thisTrack.at(iNode).Y(),thisTrack.at(iNode).Z());
            TVector3 newNode(thisTrack.at(iNode).X()-(PosCorr[0]+0.7), thisTrack.at(iNode).Y()+PosCorr[1], thisTrack.at(iNode).Z()+PosCorr[2]  );
            newPath.push_back(newNode);
        }
        return newPath;
    }

    //
    // When reading proton file
    //
    void AStarTracker::ReadProtonTrackFile(){
        _SelectableTracks.clear();
        std::vector<int> trackinfo(4);
        std::ifstream file("/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/DLwork/DataFiles/passedGBDT_extBNB_AnalysisTrees_cosmic_trained_only_on_mc_score_0.99.csv");
        if(!file){std::cout << "ERROR, could not open file of tracks to sort through" << std::endl;return;}
        std::string firstline;
        getline(file, firstline);
        bool goOn = true;
        int Run,SubRun,Event,TrackID,WireUMin,TimeUMin,WireUMax,TimeUMax,WireVMin,TimeVMin,WireVMax,TimeVMax,WireYMin,TimeYMin,WireYMax,TimeYMax;
        double BDTScore;
        while(goOn){
            file >> Run >> SubRun >> Event >> TrackID >> WireUMin >> TimeUMin >> WireUMax >> TimeUMax >> WireVMin >> TimeVMin >> WireVMax >> TimeVMax >> WireYMin >> TimeYMin >> WireYMax >> TimeYMax >> BDTScore;
            trackinfo[0] = Run;
            trackinfo[1] = SubRun;
            trackinfo[2] = Event;
            trackinfo[3] = TrackID;
            _SelectableTracks.push_back(trackinfo);
            if(file.eof()){goOn=false;break;}
        }
    }
    
    bool AStarTracker::IsGoodTrack(){
        for(auto trackinfo:_SelectableTracks){
            if(_run    != trackinfo[0])continue;
            if(_subrun != trackinfo[1])continue;
            if(_event  != trackinfo[2])continue;
            if(_track  != trackinfo[3])continue;
            return true;
            //if(_run == trackinfo[0] && _subrun == trackinfo[1] && _event == trackinfo[2] && _track == trackinfo[3]){return true;}
        }
        return false;
    }

    void AStarTracker::tellMe(std::string s, int verboseMin = 0){
        if(_verbose >= verboseMin) std::cout << s << std::endl;
    }

}
#endif
