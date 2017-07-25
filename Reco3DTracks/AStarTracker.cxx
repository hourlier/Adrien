#ifndef LARLITE_AStarTracker_CXX
#define LARLITE_AStarTracker_CXX

#include "AStarTracker.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/wire.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"
#include "TVector3.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/core/DataFormat/ChStatus.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgoProton.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/SCE/SpaceChargeMicroBooNE.h"

namespace larlite {
    
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
    //______________________________________________________
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
    //______________________________________________________
    void AStarTracker::tellMe(std::string s, int verboseMin = 0){
        if(_verbose >= verboseMin) std::cout << s << std::endl;
    }
    //______________________________________________________
    void AStarTracker::SetTimeAndWireBounds(){
        time_bounds.clear();
        wire_bounds.clear();
        time_bounds.reserve(3);
        wire_bounds.reserve(3);
        for(size_t iPlane=0; iPlane<3; ++iPlane) {
            time_bounds[iPlane].first  = 1.e9;
            time_bounds[iPlane].second = 0.;
            wire_bounds[iPlane].first  = 1.e9;
            wire_bounds[iPlane].second = 0.;
        }

        for(size_t iPlane=0;iPlane<3;iPlane++){
            // make sure starting point is projected back in the range
            if(X2Tick(start_pt.X(),iPlane) < time_bounds[iPlane].first ) time_bounds[iPlane].first  = X2Tick(start_pt.X(),iPlane);
            if(X2Tick(start_pt.X(),iPlane) > time_bounds[iPlane].second) time_bounds[iPlane].second = X2Tick(start_pt.X(),iPlane);
            double wireProjStartPt = larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w / 0.3;
            if(wireProjStartPt < wire_bounds[iPlane].first ) wire_bounds[iPlane].first  = wireProjStartPt;
            if(wireProjStartPt > wire_bounds[iPlane].second) wire_bounds[iPlane].second = wireProjStartPt;

            // make sure ending point is projected back in the range
            if(X2Tick(end_pt.X(),iPlane) < time_bounds[iPlane].first ) time_bounds[iPlane].first  = X2Tick(end_pt.X(),iPlane);
            if(X2Tick(end_pt.X(),iPlane) > time_bounds[iPlane].second) time_bounds[iPlane].second = X2Tick(end_pt.X(),iPlane);
            double wireProjEndPt   = larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w / 0.3;
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
        double ImageMargin = 50.;
        double fractionImageMargin = 0.25;
        for(size_t iPlane=0;iPlane<3;iPlane++){
            time_bounds[iPlane].first  -= std::max(ImageMargin,fractionImageMargin*(time_bounds[iPlane].second-time_bounds[iPlane].first));
            time_bounds[iPlane].second += std::max(ImageMargin,fractionImageMargin*(time_bounds[iPlane].second-time_bounds[iPlane].first));
            wire_bounds[iPlane].first  -= std::max(ImageMargin,fractionImageMargin*(wire_bounds[iPlane].second-wire_bounds[iPlane].first));
            wire_bounds[iPlane].second += std::max(ImageMargin,fractionImageMargin*(wire_bounds[iPlane].second-wire_bounds[iPlane].first));

            if(time_bounds[iPlane].first < 0) time_bounds[iPlane].first = 0;
            if(wire_bounds[iPlane].first < 0) wire_bounds[iPlane].first = 0;

            wire_bounds[iPlane].first  = (size_t)(wire_bounds[iPlane].first + 0.5);
            wire_bounds[iPlane].second = (size_t)(wire_bounds[iPlane].second + 0.5);
        }

        // make sure the number of rows and cols are divisible by _compressionFactor
        for(size_t iPlane=0;iPlane<3;iPlane++){
            while(!( (size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first)%_compressionFactor_t == 0)){
                time_bounds[iPlane].second++;//=(_compressionFactor_t-((size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first)%_compressionFactor_t));
            }
            tellMe(Form("%zu rows for %d compression factor",(size_t)(time_bounds[iPlane].second - time_bounds[iPlane].first),_compressionFactor_t),1);
        }

        for(int iPlane = 0;iPlane<3;iPlane++){

            while(!( (size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first)%_compressionFactor_w == 0)){
                wire_bounds[iPlane].second++;//=(_compressionFactor_w-((size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first)%_compressionFactor_w));
            }
            tellMe(Form("%zu cols for %d compression factor",(size_t)(wire_bounds[iPlane].second - wire_bounds[iPlane].first),_compressionFactor_w),1);
        }
    }
    //______________________________________________________
    bool AStarTracker::CheckEndPointsInVolume(){
        bool endpointInRange = true;
        if(start_pt.X() < 0    || end_pt.X() < 0   ){endpointInRange = false;}
        if(start_pt.X() > 260  || end_pt.X() > 260 ){endpointInRange = false;}

        if(start_pt.Y() < -115 || end_pt.Y() < -115){endpointInRange = false;}
        if(start_pt.Y() >  115 || end_pt.Y() >  115){endpointInRange = false;}

        if(start_pt.Z() < 1    || end_pt.Z() < 1   ){endpointInRange = false;}
        if(start_pt.Z() > 1036 || end_pt.Z() > 1036){endpointInRange = false;}

        return endpointInRange;
    }
    //______________________________________________________
    bool AStarTracker::initialize() {
        ReadProtonTrackFile();
        _eventTreated = 0;
        _eventSuccess = 0;
        time_bounds.reserve(3);
        wire_bounds.reserve(3);
        hit_image_v.reserve(3);
        hDistance2Hit = new TH1D("hDistance2Hit","hDistance2Hit",100,0,50);
        if(_DrawOutputs){
            hdQdx = new TH1D("hdQdx","hdQdx",100,0,100);
            hdQdxEntries = new TH1D("hdQdxEntries","hdQdxEntries",hdQdx->GetNbinsX(),hdQdx->GetXaxis()->GetXmin(),hdQdx->GetXaxis()->GetXmax());
            hdQdX2D = new TH2D("hdQdX2D","hdQdX2D",200,0,100,300,0,300);
        }
        return true;
    }
    //______________________________________________________
    bool AStarTracker::analyze(storage_manager* storage) {
        gStyle->SetOptStat(0);
        //
        // Retrieve data products
        //
        auto ev_track    = storage->get_data<event_track>   ( _track_producer    ); // tracks
        auto ev_mct      = storage->get_data<event_mctrack> ( _mctrack_producer  ); // mctracks
        auto ev_gaushit  = storage->get_data<event_hit>     ( _hit_producer      ); // gaushits
        auto ev_wire     = storage->get_data<event_wire>    ( _wire_producer     ); // caldata

        // Check validity
        if(!ev_track)    throw DataFormatException("Could not locate event_track data product!");
        //if(!ev_mct)      throw DataFormatException("Could not locate event_mctrack data product!");
        if(!ev_gaushit)  throw DataFormatException("Could not locate event_gaushit data product!");
        if(ev_track->empty()){tellMe("ERROR, ev_track empty",0);return true;}

        _run = storage->run_id();
        _subrun = storage->subrun_id();
        _event = storage->event_id();

        // Get associated hits + association info
        larlite::event_hit* ev_hit=nullptr;
        auto const& track_to_hit = storage->find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
        if(!ev_hit) throw DataFormatException("Could not find associated hit data product!");

        //
        //loop over tracks
        //
        for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
            _track = track_index;
            //
            // is this a track I should waste time on?
            //
            if( _track_producer!= "dl" && !(_SelectableTracks.size()!=0 && IsGoodTrack()) ) continue;
            tellMe(Form("%d\t%d\t%d\t%d",_run,_subrun,_event,_track),0);

            //
            // retrieve 3D start/end points
            //
            auto const& track = (*ev_track)[track_index];
            start_pt = track.Vertex();
            end_pt   = track.End();
            if(!CheckEndPointsInVolume()){tellMe("Point out of range",0); continue;}

            int thisTrackID = track.ID();

            //
            // find corresponding MC track
            //
            if(ev_mct && _track_producer == "dl"){
                larlite::mctrack trueTrack;
                for(auto const& mct : *ev_mct) {
                    double DRstart = sqrt(pow(start_pt.X()-mct.Start().X(),2)+pow(start_pt.Y()-mct.Start().Y(),2)+pow(start_pt.Z()-mct.Start().Z(),2));
                    double DRend = sqrt(pow(end_pt.X()-mct.End().X(),2)+pow(end_pt.Y()-mct.End().Y(),2)+pow(end_pt.Z()-mct.End().Z(),2));
                    if(DRstart > 10 || DRend > 10)continue;
                    trueTrack = mct;
                }
                if(trueTrack.size() == 0){tellMe("Could not find MC track corresponding to provided start and end points => stop",0);continue;}
            }


            //
            // 0-ter) get vector of hits and channel status to create image2D in next step
            //

            std::vector<larlite::hit> trackHit_v;
            for(auto const& hit_index : track_to_hit[track_index]) {
                auto const& h = (*ev_hit)[hit_index];
                trackHit_v.push_back(h);
            }
            tellMe("hits collected for this track",1);

            //
            // create Image 2D
            //

            SetTimeAndWireBounds();
            CreateDataImage(*ev_wire);
            tellMe("hit_image_v created", 1);
            if(_DrawOutputs){DrawROI();tellMe("ROI pre reco drawn",1);}

            //
            // Check the end points match the wire info
            //

            tellMe("About to check EndPoints",1);
            CheckEndPoints();

            //
            // 3) call reconstruct! and set the return track ovject to override the proto-track
            //
            //std::cout << "about to call Reconstruct()" << std::endl;
            _eventTreated++ ;
            (*ev_track)[track_index] = Reconstruct();
            (*ev_track)[track_index].set_track_id(thisTrackID);
            tellMe("Track reconstructed",1);

            //CompareReco2hits((*ev_track)[track_index]);
            if(_DrawOutputs){DrawTrack( (*ev_track)[track_index] );}

            std::cout << std::endl << std::endl;
        }
        return true;
    }
    //______________________________________________________
    larlite::track AStarTracker::Reconstruct(){

        tellMe("creating ch status and tag image",1);
        chstatus_image_v.clear();
        std::vector<larcv::Image2D> tag_image_v(3);
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            tag_image_v[iPlane] = hit_image_v[iPlane];
            tag_image_v[iPlane].paint(0);

            chstatus_image_v.push_back(hit_image_v[iPlane]);
            chstatus_image_v[iPlane].paint(0);
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
        config.astar_threshold.resize(3,1);    // min value for a pixel to be considered non 0
        config.astar_neighborhood.resize(3,5); //can jump over n empty pixels
        config.astar_start_padding = 20;       // allowed region around the start point
        config.astar_end_padding = 20;         // allowed region around the end point
        config.lattice_padding = 5;            // margin around the edges
        config.min_nplanes_w_hitpixel = 1;     // minimum allowed coincidences between non 0 pixels across planes
        config.restrict_path = false; // do I want to restrict to a cylinder around the strainght line of radius defined bellow ?
        config.path_restriction_radius = 30.0;

        //_______________________________________________
        // Define A* algo
        //-----------------------------------------------
        larlitecv::AStar3DAlgoProton algo( config );
        algo.setVerbose(0);
        algo.setPixelValueEsitmation(true);

        tellMe("starting A*",1);
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            tellMe(Form("plane %zu %zu rows and %zu cols before findpath", iPlane,hit_image_v[iPlane].meta().rows(),hit_image_v[iPlane].meta().cols()),1);
        }
        RecoedPath = algo.findpath( hit_image_v, chstatus_image_v, tag_image_v, start_row, end_row, start_cols, end_cols, goal_reached );
        tellMe("A* done",1);
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

        //-----------------------------------------------
        // Draw dQdX profile
        //-----------------------------------------------
        if(_DrawOutputs){
            tellMe("calling DrawdQdX",1);
            DrawdQdX(newTrack);
        }
        return newTrack;
    }
    //______________________________________________________
    void AStarTracker::DrawdQdX(larlite::track thisTrack){
        bool highdQdX = false;
        for(size_t iNode = 0;iNode<thisTrack.NumberTrajectoryPoints()-1;iNode++){
            double dqdx =   0;
            int Nplanes = 0;
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ) != 0){dqdx += thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ); Nplanes++;}
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kU) != 0){dqdx += thisTrack.DQdxAtPoint(iNode,larlite::geo::kU); Nplanes++;}
            if(thisTrack.DQdxAtPoint(iNode,larlite::geo::kV) != 0){dqdx += thisTrack.DQdxAtPoint(iNode,larlite::geo::kV); Nplanes++;}
            if(Nplanes == 0) continue;
            dqdx *= 3./Nplanes;

            if(dqdx > 300) highdQdX = true;

            hdQdx->Fill(thisTrack.Length(iNode),dqdx);
            hdQdxEntries->Fill(thisTrack.Length(iNode));
            hdQdX2D->Fill(thisTrack.Length(iNode),dqdx);
        }
        if(highdQdX){tellMe("\t \t high dQdX",1);}
    }
    //______________________________________________________
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
    //______________________________________________________
    void AStarTracker::CompareReco2hits(const larlite::track recoTrack){

        for(size_t iNode = 0; iNode<recoTrack.NumberTrajectoryPoints();iNode++){

            double dist=1e9;
            for(int iPlane=0;iPlane<3;iPlane++){

                double thiscol = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(recoTrack.LocationAtPoint(iNode),iPlane).w/0.3);
                double thisrow = hit_image_v[iPlane].meta().row(X2Tick(recoTrack.LocationAtPoint(iNode).X(),iPlane));

                for(size_t icol = 0;icol<hit_image_v[iPlane].meta().cols();icol++){
                    for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){

                        if(hit_image_v[iPlane].pixel(irow,icol)<=1)continue;

                        double thisdist = sqrt(pow(_compressionFactor_t*(thisrow-irow),2)+pow(_compressionFactor_w*(thiscol-icol),2));

                        if(dist > thisdist)dist = thisdist;
                    }
                }
            }
            hDistance2Hit->Fill(dist);
        }
    }
    //______________________________________________________
    bool AStarTracker::finalize() {
        tellMe("Finalizing",0);
        tellMe(Form("%d / %d = %.2f",_eventSuccess,_eventTreated,_eventSuccess*100./_eventTreated),0);

        if(hDistance2Hit->GetEntries() != 0){
            TCanvas *cDist = new TCanvas();
            hDistance2Hit->Draw();
            cDist->SetLogy();
            //hDistance2Hit->Write();
            tellMe(Form("hDistance2Hit mean: %.2f , RMS : %.2f",hDistance2Hit->GetMean(),hDistance2Hit->GetRMS()),0);
            tellMe(Form("entries : %.1f",hDistance2Hit->GetEntries()),0);
            cDist->SaveAs("DistanceReco2Hit.pdf");
        }

        if(_DrawOutputs){

            if(hdQdx->GetEntries() > 1){
                TCanvas *cdQdX = new TCanvas("cdQdX","cdQdX",800,600);
                for(int i = 0;i<hdQdx->GetNbinsX(); i++){
                    if(hdQdxEntries->GetBinContent(i+1) == 0) continue;
                    if(hdQdx->GetBinContent(i+1) == 0) continue;
                    double a = hdQdx->GetBinContent(i+1);
                    hdQdx->SetBinContent(i+1,a*1./hdQdxEntries->GetBinContent(i+1));
                    hdQdx->SetBinError(i+1,0.001*hdQdx->GetBinContent(i+1));

                }
                hdQdx->Draw();
                cdQdX->SaveAs("cdQdX.pdf");
            }
            if(hdQdX2D->GetEntries() > 1){
                TCanvas *cdQdX2D = new TCanvas("cdQdX2D","cdQdX2D",800,600);
                hdQdX2D->Draw("colz");
                cdQdX2D->SaveAs("cdQdX2D.pdf");
                cdQdX2D->SaveAs("cdQdX2D.root");
            }
        }
        return true;
    }
    //______________________________________________________
    bool AStarTracker::CheckEndPoints(){
        tellMe(Form("start point : (%f,%f,%f)",start_pt.X(),start_pt.Y(),start_pt.Z()),2);
        tellMe(Form("end   point : (%f,%f,%f)",end_pt.X(),end_pt.Y(),end_pt.Z()),2);
        TVector3 EndPoint[2] = {start_pt,end_pt};
        double dR = 0.1;
        for(size_t point = 0;point<2;point++){
            double minDist = EvalMinDist(EndPoint[point]);
            if(minDist <= 3){tellMe("Point OK",1);continue;}
            tellMe("Updating point",0);
            TVector3 newPoint = EndPoint[point];
            int iter = 0;
            while (minDist > 3 && iter <= 100) {
                iter++;
                std::vector<TVector3> openSet = GetOpenSet(newPoint,dR);
                double minDist = 1e9;
                for(auto neighbour : openSet){
                    double dist = EvalMinDist(neighbour);
                    if(dist < minDist){minDist=dist;newPoint = neighbour;}
                }
                tellMe(Form("iteration #%d",iter),1);
            }
            EndPoint[point] = newPoint;
            tellMe("Point updated",0);
        }
        return true;
    }
    //______________________________________________________
    bool AStarTracker::CheckEndPoints(std::vector< std::pair<int,int> > endPix){
        tellMe(Form("start point : (%f,%f,%f)",start_pt.X(),start_pt.Y(),start_pt.Z()),2);
        tellMe(Form("end   point : (%f,%f,%f)",end_pt.X(),end_pt.Y(),end_pt.Z()),2);
        TVector3 EndPoint[2] = {start_pt,end_pt};
        double dR = 0.1;

        for(size_t point = 0;point<2;point++){
            double minDist = EvalMinDist(EndPoint[point],endPix);
            if(minDist <= 3){tellMe("Point OK",1);continue;}
            tellMe("Updating point",0);
            TVector3 newPoint = EndPoint[point];
            int iter = 0;
            while (minDist > 3 && iter <= 100) {
                iter++;
                std::vector<TVector3> openSet = GetOpenSet(newPoint,dR);
                double minDist = 1e9;
                for(auto neighbour : openSet){
                    double dist = EvalMinDist(neighbour,endPix);
                    if(dist < minDist){minDist=dist;newPoint = neighbour;}
                }
                tellMe(Form("iteration #%d",iter),1);
            }
            EndPoint[point] = newPoint;
            tellMe("Point updated",0);
        }
        return true;
    }
    //______________________________________________________
    double AStarTracker::EvalMinDist(TVector3 point){
        tellMe(Form("EvalMinDist : (%f,%f,%f)",point.X(),point.Y(),point.Z()),2);
        double minDist = 1e9;
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            tellMe(Form("plane %zu :",iPlane),2);
            int pointCol = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(point,iPlane).w/0.3);
            int pointRow = hit_image_v[iPlane].meta().row(X2Tick(point.X(),iPlane));
            tellMe(Form("\t (%d,%d)",pointRow,pointCol),2);
            double dist = 1e9;
            for(int icol=0;icol<hit_image_v[iPlane].meta().cols();icol++){
                for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                    if(hit_image_v[iPlane].pixel(irow,icol) == 0)continue;
                    dist = sqrt(pow(pointCol-icol,2)+pow(pointRow-irow,2));
                    if(dist < minDist)minDist=dist;
                }
            }
        }
        return minDist;
    }
    //______________________________________________________
    double AStarTracker::EvalMinDist(TVector3 point, std::vector< std::pair<int,int> > endPix){
        tellMe(Form("EvalMinDist : (%f,%f,%f)",point.X(),point.Y(),point.Z()),2);
        double minDist = 1e9;
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            tellMe(Form("plane %zu :",iPlane),2);
            int pointCol = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(point,iPlane).w/0.3);
            int pointRow = hit_image_v[iPlane].meta().row(X2Tick(point.X(),iPlane));
            tellMe(Form("\t (%d,%d)",pointRow,pointCol),2);
            double dist = 1e9;
            for(int icol=0;icol<hit_image_v[iPlane].meta().cols();icol++){
                for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                    if(icol != endPix[iPlane].first && irow != endPix[iPlane].second)continue;
                    if(hit_image_v[iPlane].pixel(irow,icol) == 0)continue;
                    dist = sqrt(pow(pointCol-icol,2)+pow(pointRow-irow,2));
                    if(dist < minDist)minDist=dist;
                }
            }
        }
        return minDist;
    }
    //______________________________________________________
    std::vector<TVector3> AStarTracker::GetOpenSet(TVector3 newPoint, double dR){
        std::vector<TVector3> openSet;

        for(int ix = -1;ix<2;ix++){
            for(int iy = -1;iy<2;iy++){
                for(int iz = -1;iz<2;iz++){
                    if(ix == 0 && iy == 0 && iz == 0)continue;
                    TVector3 neighbour = newPoint;
                    neighbour.SetX(newPoint.X()+dR*ix);
                    neighbour.SetY(newPoint.Y()+dR*iy);
                    neighbour.SetZ(newPoint.Z()+dR*iz);
                    openSet.push_back(neighbour);
                }
            }
        }
        return openSet;
    }
    //______________________________________________________
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
    //______________________________________________________
    larlite::track AStarTracker::ComputedQdX(larlite::track newTrack, const std::vector<larcv::Image2D>& hit_image_v, const std::vector<larcv::Image2D>& chstatus_image_v){
        double X0,Y0,X1,Y1,X2,Y2;
        std::vector<double> planedQdx;
        double Npx = 5;
        double dist,xp,yp,alpha;
        double localdQdx;
        double nodeDirX,nodeDirY,nodeDirZ,norm;

        tellMe("computing dqdx",1);
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            //std::cout << "\t" << iPlane << std::endl;

            for(size_t iNode=0;iNode<newTrack.NumberTrajectoryPoints()-1;iNode++){
                localdQdx=0;
                int NpxLocdQdX = 0;
                X1 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(iNode  ),iPlane).w/0.3);
                Y1 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(iNode  ).X(),iPlane));

                X2 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(iNode+1),iPlane).w/0.3);
                Y2 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(iNode+1).X(),iPlane));

                if(X1 == X2){X1 -= 0.0001;}
                if(Y1 == Y2){Y1 -= 0.0001;}

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
                        X0 = icol;//+0.5;
                        Y0 = irow;//+0.5;
                        if(X0 == X1 || X0 == X2){X0 += 0.0001;}
                        if(Y0 == Y1 || Y0 == Y2){Y0 += 0.0001;}
                        alpha = ( (X2-X1)*(X0-X1)+(Y2-Y1)*(Y0-Y1) )/( pow(X2-X1,2)+pow(Y2-Y1,2) ); // normalized scalar product of the 2 vectors
                        //if(alpha <= -0.25 || alpha >= 1.25) continue;
                        xp = X1+alpha*(X2-X1);
                        yp = Y1+alpha*(Y2-Y1);
                        dist = sqrt(pow(xp-X0,2)+pow(yp-Y0,2));
                        double dist2point = std::min( sqrt(pow(X0-X1,2)+pow(Y0-Y1,2)), sqrt(pow(X0-X2,2)+pow(Y0-Y2,2)) );

                        if( alpha >= 0 && alpha <= 1 && dist > Npx ) continue;
                        if( (alpha < 0 || alpha > 1) && dist2point > Npx ) continue;

                        theta0 = std::abs(std::acos(((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0))/(sqrt(pow(X1-X0,2)+pow(Y1-Y0,2))*sqrt(pow(X2-X0,2)+pow(Y2-Y0,2)))));
                        // check all other nodes, if I find an angle bigger with a distance lower than Npx, reject the point for the current segment

                        bool bestCandidate = true;
                        for(size_t jNode=0;jNode<newTrack.NumberTrajectoryPoints()-1;jNode++){
                            if(jNode==iNode)continue;
                            double x1,y1,x2,y2, alpha2, xp2,yp2, dist2;
                            x1 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(jNode  ),iPlane).w/0.3);
                            x2 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(newTrack.LocationAtPoint(jNode+1),iPlane).w/0.3);

                            y1 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(jNode  ).X(),iPlane));
                            y2 = hit_image_v[iPlane].meta().row(X2Tick(newTrack.LocationAtPoint(jNode+1).X(),iPlane));

                            if(x1 == x2){x1 -= 0.0001;}
                            if(y1 == y2){y1 -= 0.0001;}

                            alpha2 = ( (x2-x1)*(X0-x1)+(y2-y1)*(Y0-y1) )/( pow(x2-x1,2)+pow(y2-y1,2) );
                            xp2 = x1+alpha2*(x2-x1);
                            yp2 = y1+alpha2*(y2-y1);
                            dist2 = sqrt(pow(xp2-X0,2)+pow(yp2-Y0,2));

                            double dist2point2 = std::min( sqrt(pow(X0-x1,2)+pow(Y0-y1,2)), sqrt(pow(X0-x2,2)+pow(Y0-y2,2)) );

                            if( alpha2 >= 0 && alpha2 <= 1 && dist2 > Npx ) continue;
                            if( (alpha2 < 0 || alpha2 > 1) && dist2point2 > Npx ) continue;

                            theta1 = std::abs(std::acos(((x1-X0)*(x2-X0)+(y1-Y0)*(y2-Y0))/(sqrt(pow(x1-X0,2)+pow(y1-Y0,2))*sqrt(pow(x2-X0,2)+pow(y2-Y0,2)))));
                            if(theta1 >= theta0){bestCandidate=false;}
                        }
                        if(!bestCandidate)continue;


                        nodeDirX =  newTrack.LocationAtPoint(iNode+1).X()-newTrack.LocationAtPoint(iNode).X();//Tick2X(RecoedPath.at(iNode-1).tyz[0],0)-Tick2X(RecoedPath.at(iNode).tyz[0],0);
                        nodeDirY =  newTrack.LocationAtPoint(iNode+1).Y()-newTrack.LocationAtPoint(iNode).Y();;//RecoedPath.at(iNode-1).tyz[1]-RecoedPath.at(iNode).tyz[1];
                        nodeDirZ =  newTrack.LocationAtPoint(iNode+1).Z()-newTrack.LocationAtPoint(iNode).Z();;//RecoedPath.at(iNode-1).tyz[2]-RecoedPath.at(iNode).tyz[2];
                        norm = sqrt(nodeDirX*nodeDirX+nodeDirY*nodeDirY+nodeDirZ*nodeDirZ);
                        localdQdx+=hit_image_v[iPlane].pixel(irow,icol)/norm;
                        NpxLocdQdX++;
                    }
                }
                if(NpxLocdQdX!=0)planedQdx.push_back(localdQdx/NpxLocdQdX);
                else planedQdx.push_back(localdQdx);
                //planedQdx.push_back(localdQdx);
            }
            newTrack.add_dqdx(planedQdx);
            planedQdx.clear();
        }
        return newTrack;
    }
    //______________________________________________________
    void AStarTracker::CreateDataImage(std::vector<larlite::wire> wire_v){
        tellMe("Entering CreateImages(wire)",1);
        hit_image_v.clear();
        //hit_image_v.reserve(3);
        const size_t num_planes = 3;

        // Using the range, construct Image2D

        for(size_t iPlane=0; iPlane<num_planes; ++iPlane) {
            // Retrieve boundaries
            auto const& time_bound = time_bounds[iPlane];
            auto const& wire_bound = wire_bounds[iPlane];
            // If no hit on this plane, make an empty image
            if(wire_bound.second <= wire_bound.first ||
               time_bound.second <= time_bound.first ) {
                tellMe("Adding Empty Image",0);
                hit_image_v[iPlane] = larcv::Image2D();
                continue;
            }
            // Construct meta
            //std::cout << "construct meta" << std::endl;
            size_t image_cols = (size_t)(wire_bound.second - wire_bound.first);
            size_t image_rows = (size_t)(time_bound.second - time_bound.first);
            tellMe(Form("%zu cols and %zu rows",image_cols,image_rows),1);
            larcv::ImageMeta hit_meta((double)image_cols, (double)image_rows,
                                      image_rows, image_cols,
                                      (size_t) wire_bound.first,  // origin x = min x
                                      (size_t) time_bound.second, // origin y = max time
                                      iPlane);


            // Prepare hit image data + fill

            std::vector<float> image_data(hit_meta.rows() * hit_meta.cols(), 0.);
            size_t row,col;
            for(auto wire : wire_v){
                unsigned int ch = wire.Channel();
                unsigned int detWire = larutil::Geometry::GetME()->ChannelToWire(ch);
                unsigned int detPlane = larutil::Geometry::GetME()->ChannelToPlane(ch);
                tellMe(Form("ch : %d, wire : %d, plane : %d", ch, detWire, detPlane),2);
                if(detPlane != iPlane) continue;
                if(detWire > wire_bound.first && detWire < wire_bound.second){
                    for (auto & iROI : wire.SignalROI().get_ranges()) {
                        int FirstTick = iROI.begin_index();
                        int time_tick = FirstTick;
                        for (float ADC : iROI) {
                            if ( time_tick > time_bound.first && time_tick < time_bound.second ) {
                                row = (size_t)(time_bound.second - time_tick + 0.5);
                                col = (size_t)(detWire - wire_bound.first + 0.5);
                                if(ADC >= _ADCthreshold){image_data[hit_meta.rows() * col + row] = ADC;}
                            }
                            time_tick++;
                        } // ADC
                    }// iROI
                } // detWire
            }
            tellMe("image_data OK",1);
            larcv::Image2D hit_image(std::move(hit_meta),std::move(image_data));
            tellMe("hit_image OK",1);
            tellMe(Form("hit_image_v.size() = %zu",hit_image_v.size()),1);
            // compress Images
            tellMe("compress Images",1);
            if(_compressionFactor_t > 1 || _compressionFactor_w > 1 ){
                tellMe(Form("plane %zu size : %zu rows, %zu cols for a %d x %d compression",iPlane,hit_image.meta().rows(),hit_image.meta().cols(),_compressionFactor_t,_compressionFactor_w),1);
                hit_image.compress(hit_image.meta().rows()/_compressionFactor_t, hit_image.meta().cols()/_compressionFactor_w);
                tellMe(Form("plane %zu %zu rows and %zu cols after compression", iPlane,hit_image.meta().rows(),hit_image.meta().cols()),1);
                tellMe("image compressed",1);
            }
            hit_image_v.push_back(hit_image);
        }
    }
    //______________________________________________________
    void AStarTracker::DrawTrack(const larlite::track recoTrack){
        double X0,Y0,X1,Y1,X2,Y2;
        double Npx = 5;
        double dist,xp,yp,alpha;
        TH2D *hImage[3];
        TH2D *hIntegral[3];
        TGraph *gTrack[3];
        TGraph *gStartNend[3];

        for(size_t iPlane=0;iPlane<3;iPlane++){
            tellMe(Form("computedQdx \t plane %zu", iPlane),1);
            gTrack[iPlane] = new TGraph();
            hImage[iPlane] = new TH2D(Form("hImage_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hImage_%d_%d_%d_%d_%zu;wire;time",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());
            hIntegral[iPlane] = new TH2D(Form("hIntegral_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hIntegral_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());
            for(int icol=0;icol<hit_image_v[iPlane].meta().cols();icol++){
                for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                    hImage[iPlane]->SetBinContent(icol+1,irow+1,hit_image_v[iPlane].pixel(irow,icol));
                }
            }

            gStartNend[iPlane] = new TGraph();
            gStartNend[iPlane]->SetPoint(0,hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w/0.3),hit_image_v[iPlane].meta().row(X2Tick(start_pt.X(),iPlane)));
            gStartNend[iPlane]->SetPoint(1,hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w/0.3),hit_image_v[iPlane].meta().row(X2Tick(end_pt.X(),iPlane)));

            for(size_t iNode=0;iNode<recoTrack.NumberTrajectoryPoints();iNode++){
                X1 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(recoTrack.LocationAtPoint(iNode  ),iPlane).w/0.3);
                Y1 = hit_image_v[iPlane].meta().row(X2Tick(recoTrack.LocationAtPoint(iNode  ).X(),iPlane));
                gTrack[iPlane]->SetPoint(iNode,X1+0.5,Y1+0.5);
            }
        }

        tellMe("computing dqdx",1);
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            //std::cout << "\t" << iPlane << std::endl;

            for(size_t iNode=0;iNode<recoTrack.NumberTrajectoryPoints()-1;iNode++){
                X1 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(recoTrack.LocationAtPoint(iNode  ),iPlane).w/0.3);
                Y1 = hit_image_v[iPlane].meta().row(X2Tick(recoTrack.LocationAtPoint(iNode  ).X(),iPlane));

                X2 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(recoTrack.LocationAtPoint(iNode+1),iPlane).w/0.3);
                Y2 = hit_image_v[iPlane].meta().row(X2Tick(recoTrack.LocationAtPoint(iNode+1).X(),iPlane));

                if(X1 == X2){X1 -= 0.0001;}
                if(Y1 == Y2){Y1 -= 0.0001;}

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
                        X0 = icol;//+0.5;
                        Y0 = irow;//+0.5;
                        if(X0 == X1 || X0 == X2){X0 += 0.0001;}
                        if(Y0 == Y1 || Y0 == Y2){Y0 += 0.0001;}
                        alpha = ( (X2-X1)*(X0-X1)+(Y2-Y1)*(Y0-Y1) )/( pow(X2-X1,2)+pow(Y2-Y1,2) ); // normalized scalar product of the 2 vectors
                        //if(alpha <= -0.25 || alpha >= 1.25) continue;
                        xp = X1+alpha*(X2-X1);
                        yp = Y1+alpha*(Y2-Y1);
                        dist = sqrt(pow(xp-X0,2)+pow(yp-Y0,2));
                        double dist2point = std::min( sqrt(pow(X0-X1,2)+pow(Y0-Y1,2)), sqrt(pow(X0-X2,2)+pow(Y0-Y2,2)) );

                        if( alpha >= 0 && alpha <= 1 && dist > Npx ) continue;
                        if( (alpha < 0 || alpha > 1) && dist2point > Npx ) continue;

                        theta0 = std::abs(std::acos(((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0))/(sqrt(pow(X1-X0,2)+pow(Y1-Y0,2))*sqrt(pow(X2-X0,2)+pow(Y2-Y0,2)))));
                        // check all other nodes, if I find an angle bigger with a distance lower than Npx, reject the point for the current segment

                        bool bestCandidate = true;
                        for(size_t jNode=0;jNode<recoTrack.NumberTrajectoryPoints()-1;jNode++){
                            if(jNode==iNode)continue;
                            double x1,y1,x2,y2, alpha2, xp2,yp2, dist2;
                            x1 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(recoTrack.LocationAtPoint(jNode  ),iPlane).w/0.3);
                            x2 = hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(recoTrack.LocationAtPoint(jNode+1),iPlane).w/0.3);

                            y1 = hit_image_v[iPlane].meta().row(X2Tick(recoTrack.LocationAtPoint(jNode  ).X(),iPlane));
                            y2 = hit_image_v[iPlane].meta().row(X2Tick(recoTrack.LocationAtPoint(jNode+1).X(),iPlane));

                            if(x1 == x2){x1 -= 0.0001;}
                            if(y1 == y2){y1 -= 0.0001;}

                            alpha2 = ( (x2-x1)*(X0-x1)+(y2-y1)*(Y0-y1) )/( pow(x2-x1,2)+pow(y2-y1,2) );
                            xp2 = x1+alpha2*(x2-x1);
                            yp2 = y1+alpha2*(y2-y1);
                            dist2 = sqrt(pow(xp2-X0,2)+pow(yp2-Y0,2));

                            double dist2point2 = std::min( sqrt(pow(X0-x1,2)+pow(Y0-y1,2)), sqrt(pow(X0-x2,2)+pow(Y0-y2,2)) );

                            if( alpha2 >= 0 && alpha2 <= 1 && dist2 > Npx ) continue;
                            if( (alpha2 < 0 || alpha2 > 1) && dist2point2 > Npx ) continue;

                            theta1 = std::abs(std::acos(((x1-X0)*(x2-X0)+(y1-Y0)*(y2-Y0))/(sqrt(pow(x1-X0,2)+pow(y1-Y0,2))*sqrt(pow(x2-X0,2)+pow(y2-Y0,2)))));
                            if(theta1 >= theta0){bestCandidate=false;}
                        }
                        if(!bestCandidate)continue;


                        hIntegral[iPlane]->Fill(X0,Y0);
                    }
                }
            }
        }

        TCanvas *c = new TCanvas(Form("c_%d_%d_%d_%d",_run,_subrun,_event,_track),Form("c_%d_%d_%d_%d",_run,_subrun,_event,_track),450,150);
        c->Divide(3,1);
        for(int iPlane=0;iPlane<3;iPlane++){
            c->cd(iPlane+1);
            //hImage[iPlane]->Draw("colz");
            hIntegral[iPlane]->Draw();
            //hIntegral[iPlane]->Draw("colz");
            hImage[iPlane]->Draw("same colz");
            gTrack[iPlane]->SetLineColor(2);
            gTrack[iPlane]->Draw("same LP");
            gStartNend[iPlane]->SetMarkerStyle(20);
            gStartNend[iPlane]->SetMarkerSize(0.25);
            gStartNend[iPlane]->Draw("same P");
        }
        c->SaveAs(Form("%s.pdf",c->GetName()));
        for(int iPlane = 0;iPlane<3;iPlane++){
            hImage[iPlane]->Delete();
            hIntegral[iPlane]->Delete();
            gTrack[iPlane]->Delete();
            gStartNend[iPlane]->Delete();
        }
    }
    //______________________________________________________
    void AStarTracker::DrawROI(){
        TH2D *hImage[3];
        TGraph *gTrack[3];
        TGraph *gStartNend[3];

        for(size_t iPlane=0;iPlane<3;iPlane++){
            gTrack[iPlane] = new TGraph();
            hImage[iPlane] = new TH2D(Form("hImage_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hImage_%d_%d_%d_%d_%zu;wire;time",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());

            for(int icol=0;icol<hit_image_v[iPlane].meta().cols();icol++){
                for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                    hImage[iPlane]->SetBinContent(icol+1,irow+1,hit_image_v[iPlane].pixel(irow,icol));
                }
            }

            gStartNend[iPlane] = new TGraph();
            gStartNend[iPlane]->SetPoint(0,hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w/0.3),hit_image_v[iPlane].meta().row(X2Tick(start_pt.X(),iPlane)));
            gStartNend[iPlane]->SetPoint(1,hit_image_v[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w/0.3),hit_image_v[iPlane].meta().row(X2Tick(end_pt.X(),iPlane)));
        }

        TCanvas *c = new TCanvas(Form("ROI_%d_%d_%d_%d",_run,_subrun,_event,_track),Form("ROI_%d_%d_%d_%d",_run,_subrun,_event,_track),450,150);
        c->Divide(3,1);
        for(int iPlane=0;iPlane<3;iPlane++){
            c->cd(iPlane+1);
            hImage[iPlane]->Draw("colz");
            gStartNend[iPlane]->SetMarkerStyle(20);
            gStartNend[iPlane]->SetMarkerSize(0.25);
            gStartNend[iPlane]->Draw("same P");
        }
        c->SaveAs(Form("%s.pdf",c->GetName()));
        for(int iPlane = 0;iPlane<3;iPlane++){
            hImage[iPlane]->Delete();
            gStartNend[iPlane]->Delete();
        }
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
    //______________________________________________________
    std::vector<TVector3> AStarTracker::CorrectSCE(larlite::mctrack thisTrack){
        std::vector<TVector3> newPath;
        for(size_t iNode = 0;iNode<thisTrack.size();iNode++){
            TVector3 node(thisTrack[iNode].X(),thisTrack[iNode].Y(),thisTrack[iNode].Z());
            newPath.push_back(node);
        }
        return CorrectSCE(newPath);
    }
    //______________________________________________________
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
        std::vector<int> trackinfo(16);
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
            trackinfo[4] = WireUMin;
            trackinfo[5] = TimeUMin;
            trackinfo[6] = WireUMax;
            trackinfo[7] = TimeUMax;
            trackinfo[8] = WireVMin;
            trackinfo[9] = TimeVMin;
            trackinfo[10] = WireVMax;
            trackinfo[11] = TimeVMax;
            trackinfo[12] = WireYMin;
            trackinfo[13] = TimeYMin;
            trackinfo[14] = WireYMax;
            trackinfo[15] = TimeYMax;

            _SelectableTracks.push_back(trackinfo);
            if(file.eof()){goOn=false;break;}
        }
    }
    //______________________________________________________
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
    //______________________________________________________

}
#endif
