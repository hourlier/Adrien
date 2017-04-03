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

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/SCE/SpaceChargeMicroBooNE.h"

namespace larlite {

    bool AStarTracker::initialize() {
        hdQdx = new TH1D("hdQdx","hdQdx"                     ,25,0,25);
        hdQdxEntries = new TH1D("hdQdxEntries","hdQdxEntries",25,0,25);
        hDistance2MC  = new TH1D("hDistance2MC" ,"hDistance2MC;#DeltaR_{reco=>MC} [cm]", 20,0,20);
        hDistance2MCX = new TH1D("hDistance2MCX","hDistance2MCX;#DeltaR_{reco=>MC} [cm]",20,0,20);
        hDistance2MCY = new TH1D("hDistance2MCY","hDistance2MCY;#DeltaR_{reco=>MC} [cm]",20,0,20);
        hDistance2MCZ = new TH1D("hDistance2MCZ","hDistance2MCZ;#DeltaR_{reco=>MC} [cm]",20,0,20);
        hDistance2Hit = new TH1D("hDistance2Hit","hDistance2Hit;distance RecoPoint=>closest Hit [wire,tick]",100,0,20);
        hDistanceMC2Hit = new TH1D("hDistanceMC2Hit","hDistanceMC2Hit;distace MC pt proj. => closest hit [wire,tick]",100,0,20);
        c2 = new TCanvas("c2","c2",800,600);
        return true;
    }

    bool AStarTracker::analyze(storage_manager* storage) {
        std::cout << "CAREFUL, next step is to limit ourselves to start and end points that are reasonably reconstructed so that we don't have obviously wrong tracks" << std::endl;
        gStyle->SetOptStat(0);
        //
        // Retrieve data products
        //
        auto ev_track    = storage->get_data<event_track>(_track_producer);       // tracks
        auto ev_chstatus = storage->get_data<event_chstatus>(_chstatus_producer); // chstatus
        auto ev_mct      = storage->get_data<event_mctrack>(_mctrack_producer);   // mctracks

        // Check validity
        if(!ev_track)    throw DataFormatException("Could not locate event_track data product!");
        if(!ev_chstatus) throw DataFormatException("Could not locate event_chstatus data product!");
        if(!ev_mct)      throw DataFormatException("Could not locate event_mctrack data product!");
        if(ev_track->empty()){std::cout << "ERROR, ev_track empty" << std::endl;return true;}

        _run = storage->run_id();
        _subrun = storage->subrun_id();
        _event = storage->event_id();

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

            std::cout << _run << "\t" << _subrun << "\t" << _event << "\t" << _track << "\t" << (*ev_track)[track_index].ID() <<  std::endl;
            if (track_to_hit[track_index].size() == 0){std::cout << "no hits associated to this track, stop here for now" << std::endl;continue;}

            //
            // 0) retrieve 3D start/end points
            //

            auto const& track = (*ev_track)[track_index];
            auto const& start_pt = track.Vertex();
            auto const& end_pt   = track.End();
            int thisTrackID = track.ID();

            //
            // 0-bis) find corresponding MC track
            //

            larlite::mctrack trueTrack;
            for(auto const& mct : *ev_mct) {
                double DRstart = sqrt(pow(start_pt.X()-mct.Start().X(),2)+pow(start_pt.Y()-mct.Start().Y(),2)+pow(start_pt.Z()-mct.Start().Z(),2));
                double DRend = sqrt(pow(end_pt.X()-mct.End().X(),2)+pow(end_pt.Y()-mct.End().Y(),2)+pow(end_pt.Z()-mct.End().Z(),2));
                if(DRstart > 20 || DRend > 20)continue;
                std::cout << "mc   track id " << mct.TrackID() << std::endl;
                std::cout << "reco track id " << track.ID() << std::endl;
                std::cout << "ΔR start   " << DRstart << std::endl;
                std::cout << "ΔR end     " << DRend << std::endl;
                std::cout << "PDG code   " << mct.PdgCode() << std::endl;
                trueTrack = mct;
            }
            std::cout << "Npoints in MC = " << trueTrack.size() << std::endl;
            if(trueTrack.size() == 0){std::cout << "Could not find MC track corresponding to provided start and end points => stop" << std::endl;continue;}

            //
            // 1) construct hit image2d ... start from figuring out the range
            // 2) construct chstatus image2d ... do this in the same loop
            //
            const size_t num_planes = 3;
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
            /*for(auto const& hit_index : track_to_hit[track_index]) {
                auto const& h = (*ev_hit)[hit_index];
                auto& time_bound = time_bounds[h.WireID().Plane];
                auto& wire_bound = wire_bounds[h.WireID().Plane];
                // make sure current hit is within the range
                if(h.PeakTime() < time_bound.first ) time_bound.first  = h.PeakTime();
                if(h.PeakTime() > time_bound.second) time_bound.second = h.PeakTime();
                if(h.WireID().Wire < wire_bound.first ) wire_bound.first  = h.WireID().Wire;
                if(h.WireID().Wire > wire_bound.second) wire_bound.second = h.WireID().Wire;
            }*/
            // make sure the start and end points are projected back into the range
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
                    if(time_bounds[jPlane].first  < time_bounds[iPlane].first) {time_bounds[iPlane].first  = time_bounds[jPlane].first;}
                    if(time_bounds[jPlane].second > time_bounds[iPlane].second){time_bounds[iPlane].second = time_bounds[jPlane].second;}
                }
            }
            // Update the range with margin
            double ImageMargin = 50.;
            for(size_t iPlane=0;iPlane<3;iPlane++){
                time_bounds[iPlane].first  -= ImageMargin;
                time_bounds[iPlane].second += ImageMargin;
                wire_bounds[iPlane].first  -= ImageMargin;
                wire_bounds[iPlane].second += ImageMargin;

                if(time_bounds[iPlane].first < 0) time_bounds[iPlane].first = 0;
                if(wire_bounds[iPlane].first < 0) wire_bounds[iPlane].first = 0;
                if(iPlane < 2 && wire_bounds[iPlane].second >= 2400) wire_bounds[iPlane].second = 2399;
                if(iPlane == 2 && wire_bounds[iPlane].second >= 3456) wire_bounds[iPlane].second = 3455;

                wire_bounds[iPlane].first  = (size_t)(wire_bounds[iPlane].first + 0.5);
                wire_bounds[iPlane].second = (size_t)(wire_bounds[iPlane].second + 0.5);
            }
            // make sure the number of row is divisible by 2
            if(!( (size_t)(time_bounds[0].second - time_bounds[0].first)%_rebinTime)){for(size_t iPlane=0;iPlane<3;iPlane++){time_bounds[iPlane].second+=1;}}
            // check if any associated hit is within the limits
            bool HitsPresent = false;
            for(auto const& hit_index : track_to_hit[track_index]) {
                auto const& h = (*ev_hit)[hit_index];
                auto& time_bound = time_bounds[h.WireID().Plane];
                auto& wire_bound = wire_bounds[h.WireID().Plane];
                // make sure current hit is within the range
                if(h.PeakTime() < time_bound.first )continue;
                if(h.PeakTime() > time_bound.second)continue;
                if(h.WireID().Wire < wire_bound.first ) continue;
                if(h.WireID().Wire > wire_bound.second) continue;
                HitsPresent=true;
            }
            if(HitsPresent == false){std::cout << "no hits in Image2D" << std::endl; continue;}

            // Using the range, construct Image2D
            std::vector<larcv::Image2D> hit_image_v;
            std::vector<larcv::Image2D> chstatus_image_v;
            hit_image_v.reserve(num_planes);
            chstatus_image_v.reserve(num_planes);
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
                size_t image_cols = (size_t)(wire_bound.second - wire_bound.first + 1);
                size_t image_rows = (size_t)(time_bound.second - time_bound.first + 1);
                larcv::ImageMeta hit_meta((double)image_cols, (double)image_rows,
                                          image_rows, image_cols,
                                          wire_bound.first,  // origin x = min x
                                          time_bound.second, // origin y = max time
                                          plane);
                auto chstatus_meta = hit_meta;


                // Prepare hit image data + fill
                std::vector<float> image_data(hit_meta.rows() * hit_meta.cols(), 0.);
                size_t row,col;
                for(auto const& hit_index : track_to_hit[track_index]) {
                    auto const& h = (*ev_hit)[hit_index];
                    if(h.WireID().Plane != plane) continue;
                    if(h.PeakTime()    < time_bound.first || h.PeakTime()    >= time_bound.second) continue;
                    if(h.WireID().Wire < wire_bound.first || h.WireID().Wire >= wire_bound.second) continue;
                    row = (size_t)(time_bound.second - h.PeakTime() + 0.5);
                    col = (size_t)(h.WireID().Wire - wire_bound.first + 0.5);
                    image_data[hit_meta.rows() * col + row] = h.Integral();
                }

                larcv::Image2D hit_image(std::move(hit_meta),std::move(image_data));
                hit_image_v.emplace_back(std::move(hit_image));

                // Prepare chstatus image data + fill
                std::vector<float> chstatus_data(chstatus_meta.rows() * chstatus_meta.cols(), 0.);
                for(auto const& chstatus_v : *ev_chstatus) {
                    if( chstatus_v.plane().Plane != plane) continue;
                    for(size_t ch=0; ch<chstatus_v.status().size(); ++ch) {
                        auto const& status = chstatus_v.status()[ch];
                        if(ch < wire_bound.first || ch >= wire_bound.second) continue;
                        for(size_t row=0; row<chstatus_meta.rows(); ++row){
                            // eveything that has a non-zero value is a bad ch for A* algo
                            if(status == larlite::larch::kGOOD){chstatus_data[(ch - wire_bound.first) * chstatus_meta.rows() + row] = 0;}
                            else {chstatus_data[(ch - wire_bound.first) * chstatus_meta.rows() + row] = 1;/*(float)(status);*/}
                        }
                    }
                    break;
                }
                larcv::Image2D chstatus_image(std::move(chstatus_meta),std::move(chstatus_data));
                chstatus_image_v.emplace_back(std::move(chstatus_image));
            }

            // compress Image2D
            for(size_t iPlane = 0;iPlane<3;iPlane++){
                hit_image_v.at(iPlane).compress(hit_image_v.at(iPlane).meta().rows()/_rebinTime,hit_image_v.at(iPlane).meta().cols());
            }

            //
            // 3) call reconstruct! and set the return track ovject to override the proto-track
            //
            (*ev_track)[track_index] = Reconstruct(start_pt, end_pt, hit_image_v, chstatus_image_v);
            (*ev_track)[track_index].set_track_id(thisTrackID);


            //
            // 4) Evaluate the differences between the true and reco'ed tracks
            //
            CompareReco2MC3D((*ev_track)[track_index],trueTrack);
            CompareReco2hits(hit_image_v);
            //CompareMC2hits(trueTrack,hit_image_v);

            std::cout << std::endl;
        }
        return true;
    }

    double AStarTracker::X2Tick(double x, size_t plane) const {

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

    double AStarTracker::Tick2X(double tick, size_t plane) const{

        auto ts = larutil::TimeService::GetME();
        auto larp = larutil::LArProperties::GetME();
        // remove tick offset due to plane
        if(plane == 0) {tick/=2.;}
        if(plane == 1) {tick += 0.3/larp->DriftVelocity(larp->Efield(1));tick/=2.;}
        if(plane == 3) {
            tick += 0.3/larp->DriftVelocity(larp->Efield(1));
            tick += 0.3/larp->DriftVelocity(larp->Efield(2));
            tick/=2.;
        }
        double x = (tick+ts->TriggerOffsetTPC()+_speedOffset)*larp->DriftVelocity();
        return x;
    }


    larlite::track AStarTracker::Reconstruct(const TVector3& start_pt, const TVector3& end_pt, const std::vector<larcv::Image2D>& hit_image_v, const std::vector<larcv::Image2D>& chstatus_image_v){
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
        else{std::cout << "error: 0 rows " << std::endl; return larlite::track();}

        for(size_t iPlane=0;iPlane<3;iPlane++){
            if(hit_image_v[iPlane].meta().width()==0)continue;
            double wireProjStartPt = larutil::GeometryHelper::GetME()->Point_3Dto2D(start_pt,iPlane).w / 0.3;
            double wireProjEndPt   = larutil::GeometryHelper::GetME()->Point_3Dto2D(end_pt,iPlane).w / 0.3;
            if(wireProjStartPt < 0) wireProjStartPt=0;
            if(wireProjEndPt < 0  ) wireProjEndPt = 0;
            start_cols[iPlane] = hit_image_v[iPlane].meta().col(wireProjStartPt);
            end_cols[iPlane]   = hit_image_v[iPlane].meta().col(wireProjEndPt);
        }

        //_______________________________________________
        // Get histograms and graphs ready
        //-----------------------------------------------
        TCanvas *cOutput = new TCanvas(Form("cOutput_%d_%d_%d_%d",_run,_subrun,_event,_track),Form("cOutput_%d_%d_%d_%d",_run,_subrun,_event,_track),600,200);
        cOutput->Divide(3,1);
        TH2D *hHitImages[3];
        TH2D *hDeadCh[3];
        TGraph *gStart[3];
        TGraph *gEnd[3];
        TGraph *gReco[3];
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            gStart[iPlane] = new TGraph();
            gEnd[iPlane] = new TGraph();
            gReco[iPlane] = new TGraph();
            gStart[iPlane]->SetPoint(0, start_cols[iPlane],start_row);
            gEnd[iPlane]->SetPoint(0, end_cols[iPlane],end_row);

            if(hit_image_v[iPlane].meta().rows() == 0 || hit_image_v[iPlane].meta().cols() == 0)continue;
            hHitImages[iPlane] = new TH2D(Form("hHitImages_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hHitImages_%d_%d_%d_%d_%zu;wires;ticks",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());
            hDeadCh[iPlane] = new TH2D(Form("hDeadCh_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),Form("hDeadCh_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane),hit_image_v[iPlane].meta().cols(),0,hit_image_v[iPlane].meta().cols(),hit_image_v[iPlane].meta().rows(),0,hit_image_v[iPlane].meta().rows());


            for(int row = 0;row<hit_image_v[iPlane].meta().rows();row++){
                for(int col = 0;col<hit_image_v[iPlane].meta().cols();col++){
                    hHitImages[iPlane]->SetBinContent(col+1,row+1,hit_image_v[iPlane].pixel(row,col));
                    hDeadCh[iPlane]->SetBinContent(col+1,row+1,chstatus_image_v[iPlane].pixel(row,col));
                }
            }
        }
        //if(hHitImages[0]->Integral()+hHitImages[1]->Integral()+hHitImages[2]->Integral() == 0){std::cout << "ERROR => EMPTY EVENT" << std::endl;return larlite::track();}

        //_______________________________________________
        // Configure A* algo
        //-----------------------------------------------
        larlitecv::AStar3DAlgoConfig config;
        config.accept_badch_nodes = true;
        config.astar_threshold.resize(3,5);
        //config.astar_threshold[2] = 10.0;
        config.astar_neighborhood.resize(3,5); //can jump over n empty pixels
        config.astar_start_padding = 15;// allowed region around the start point
        config.astar_end_padding = 15;  // allowed region around the end point
        config.lattice_padding = 5; // margin around the edges
        config.min_nplanes_w_hitpixel = 1;
        config.restrict_path = false; // do I want to restrict to a cylinder around the strainght line of radius defined bellow
        config.path_restriction_radius = 30.0;

        //_______________________________________________
        // Define A* algo
        //-----------------------------------------------
        larlitecv::AStar3DAlgo algo( config );
        algo.setVerbose(0);
        algo.setPixelValueEsitmation(true);

        std::cout << "starting A*" << std::endl;
        RecoedPath = algo.findpath( hit_image_v, chstatus_image_v, tag_image_v, start_row, end_row, start_cols, end_cols, goal_reached );

        //_______________________________________________
        // Make track out of 3D nodes
        //-----------------------------------------------
        larlite::track newTrack;
        int nNodes[3] = {0,0,0};
        double nodeX,nodeY,nodeZ,nodeDirX,nodeDirY,nodeDirZ,norm;
        //for(size_t iNode=0;iNode<RecoedPath.size();iNode++){
        for(size_t iNode=RecoedPath.size()-1;iNode>0;iNode--){
            nodeX = Tick2X(RecoedPath.at(iNode).tyz[0],0);
            nodeY = RecoedPath.at(iNode).tyz[1];
            nodeZ = RecoedPath.at(iNode).tyz[2];
            TVector3 node(nodeX,nodeY,nodeZ);
            newTrack.add_vertex(node);
            if(iNode<RecoedPath.size()-1){
                nodeDirX =  Tick2X(RecoedPath.at(iNode+1).tyz[0],0)-Tick2X(RecoedPath.at(iNode).tyz[0],0);
                nodeDirY =  RecoedPath.at(iNode+1).tyz[1]-RecoedPath.at(iNode).tyz[1];
                nodeDirZ =  RecoedPath.at(iNode+1).tyz[2]-RecoedPath.at(iNode).tyz[2];
                norm = sqrt(nodeDirX*nodeDirX+nodeDirY*nodeDirY+nodeDirZ*nodeDirZ);
                TVector3 nodeDir(nodeDirX/norm,nodeDirY/norm,nodeDirZ/norm);
                newTrack.add_direction(nodeDir);
            }
            else{
                nodeDirX = Tick2X(RecoedPath.at(iNode).tyz[0],0)-Tick2X(RecoedPath.at(0).tyz[0],0);
                nodeDirY = RecoedPath.at(iNode).tyz[1]-RecoedPath.at(0).tyz[1];
                nodeDirZ = RecoedPath.at(iNode).tyz[2]-RecoedPath.at(0).tyz[2];
                norm = sqrt(nodeDirX*nodeDirX+nodeDirY*nodeDirY+nodeDirZ*nodeDirZ);
                TVector3 nodeDir(nodeDirX/norm,nodeDirY/norm,nodeDirZ/norm);
                newTrack.add_direction(nodeDir);

            }
            //_______________________________________________
            // fill graph
            //-----------------------------------------------
            for(int iPlane=0;iPlane<3;iPlane++){
                gReco[iPlane]->SetPoint(nNodes[iPlane],RecoedPath.at(iNode).cols[iPlane],RecoedPath.at(iNode).row);
                nNodes[iPlane]++;
            }
        }

        //_______________________________________________
        // Compute dQdx for newTrack
        //-----------------------------------------------

        double a,b,c;
        double X0,Y0;
        double X1,Y1;
        std::vector<double> planedQdx;
        double Npx = 2.5;
        double dist,xp,yp;
        double localdQdx;

        std::cout << "computing dqdx" << std::endl;
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            for(size_t iNode = RecoedPath.size()-1; iNode>1 ;iNode--){
                localdQdx=0;
                //_______________________________________________
                // Get the line equation between node i and i+1
                //-----------------------------------------------
                X0 = RecoedPath.at(iNode).cols[iPlane];
                X1 = RecoedPath.at(iNode-1).cols[iPlane];
                Y0 = RecoedPath.at(iNode).row;
                Y1 = RecoedPath.at(iNode-1).row;
                //if(X1[iPlane] == X0[iPlane] && Y1[iPlane] == Y0[iPlane]){continue;}
                if(X1 == X0){X0 -= 0.001;}
                if(Y1 == Y0){Y0 -= 0.001;}
                a = 1.*(Y1-Y0)/(X1-X0);
                b = -1.;
                c = -1.*b*Y0-a*X0;

                //_______________________________________________
                // what we need to do now is to get all hits
                // that are within Npx pixels of the line
                //-----------------------------------------------
                for(int icol=0;icol<hit_image_v[iPlane].meta().cols();icol++){
                    for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                        if(hit_image_v[iPlane].pixel(irow,icol) == 0)continue;
                        double x = 1.*icol+0.5;
                        double y = 1.*irow+0.5;

                        // for point (x,y) and line ax+by+c=0, distance = |ax+by+c|/sqrt(a^2+b^2)
                        dist = sqrt(pow((a*x+b*y+c),2)/(a*a+b*b));
                        if(dist > Npx)continue;

                        // the projection of the point on the track line must also be within the track segment
                        xp = 1.*(b*( b*x-a*y)-a*c)/(a*a+b*b);
                        yp = 1.*(a*(-b*x+a*y)-b*c)/(a*a+b*b);
                        if(!((X0 <= xp && xp <= X1) || (X1 <= xp && xp <= X0)))continue;

                        nodeDirX =  Tick2X(RecoedPath.at(iNode-1).tyz[0],0)-Tick2X(RecoedPath.at(iNode).tyz[0],0);
                        nodeDirY =  RecoedPath.at(iNode-1).tyz[1]-RecoedPath.at(iNode).tyz[1];
                        nodeDirZ =  RecoedPath.at(iNode-1).tyz[2]-RecoedPath.at(iNode).tyz[2];
                        norm = sqrt(nodeDirX*nodeDirX+nodeDirY*nodeDirY+nodeDirZ*nodeDirZ);
                        localdQdx+=hit_image_v[iPlane].pixel(irow,icol)/norm;
                        //hHitImages[iPlane]->SetBinContent(icol+1,irow+1,1);
                    }
                }
                planedQdx.push_back(localdQdx);
            }
            newTrack.add_dqdx(planedQdx);
            planedQdx.clear();

        }



        //_______________________________________________
        // Draw outputs
        //-----------------------------------------------

        for(size_t iPlane = 0;iPlane<3;iPlane++){
            cOutput->cd(iPlane+1);
            cOutput->cd(iPlane+1)->SetGrid();
            hHitImages[iPlane]->Draw("colz");
            hDeadCh[iPlane]->SetMarkerColor(kGray);
            if(hDeadCh[iPlane]->Integral()>0)hDeadCh[iPlane]->Draw("same");
            gStart[iPlane]->SetMarkerStyle(7);
            gStart[iPlane]->SetMarkerColor(2);
            gEnd[iPlane]->SetMarkerStyle(7);
            gEnd[iPlane]->SetMarkerColor(6);
            gStart[iPlane]->Draw("same P");
            gEnd[iPlane]->Draw("same P");
            gReco[iPlane]->SetMarkerColor(4);
            gReco[iPlane]->SetLineColor(4);
            gReco[iPlane]->SetMarkerStyle(4);
            gReco[iPlane]->SetMarkerSize(0.75);
            gReco[iPlane]->SetName(Form("gReco_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane));
            gReco[iPlane]->SetTitle(Form("gReco_%d_%d_%d_%d_%zu",_run,_subrun,_event,_track,iPlane));
            if(gReco[iPlane]->GetN()>0)gReco[iPlane]->Draw("same LP");
            cOutput->Modified();
            cOutput->Update();
        }

        cOutput->SaveAs(Form("plot/cOutput_%d_%d_%d_%d.pdf",_run,_subrun,_event,_track));
        DrawdQdX(newTrack);
        return newTrack;
    }

    void AStarTracker::DrawdQdX(larlite::track thisTrack){
        for(size_t iNode = 0;iNode<thisTrack.NumberTrajectoryPoints()-1;iNode++){
            double dqdx = thisTrack.DQdxAtPoint(iNode,larlite::geo::kZ)+thisTrack.DQdxAtPoint(iNode,larlite::geo::kU)+thisTrack.DQdxAtPoint(iNode,larlite::geo::kV);
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
            hDistance2MCX->Fill(dPointSegmentX);
            hDistance2MCY->Fill(dPointSegmentY);
            hDistance2MCZ->Fill(dPointSegmentZ);
            if(dPointSegment > 20){std::cout << "STRANGE, TOO FAR FROM MC" << std::endl;}
        }

    }

    void AStarTracker::CompareReco2hits(const std::vector<larcv::Image2D> hit_image_v){
        for(size_t iNode = 0;iNode<RecoedPath.size();iNode++){
            int thisrow = RecoedPath.at(iNode).row;
            int thiscol[3];
            double dist=1e9;
            for(size_t iPlane=0;iPlane<3;iPlane++){
                thiscol[iPlane] = RecoedPath.at(iNode).cols[iPlane];
                //for a given Image2D (a given iPlane), we loop over the image and search for the closest non-zero pixel
                for(size_t icol = 0;icol<hit_image_v[iPlane].meta().cols();icol++){
                    for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                        if(hit_image_v[iPlane].pixel(irow,icol)==0)continue;
                        double thisdist = sqrt(pow(2*(thisrow-irow),2)+pow(thiscol[iPlane]-icol,2)); // factor of 2 in row to acount for the fact that my rows are 2 ticks long
                        if(dist > thisdist)dist = thisdist;
                    }
                }
            }
            hDistance2Hit->Fill(dist);
        }
    }

    void AStarTracker::CompareMC2hits(const larlite::mctrack mct, const std::vector<larcv::Image2D> hit_image_v){
        //std::vector<larcv::Image2D> MChits(3);
        /*for(size_t iNode = 0;iNode<mct.size();iNode++){
            int thisrow = hit_image_v[0].meta().row(X2Tick(mct[iNode].X(),0));
            int thiscol[3];
            double dist=1e9;
            for(size_t iPlane=0;iPlane<3;iPlane++){
                //MChits[iPlane].paint(0);
                double wireProj = larutil::GeometryHelper::GetMe()->Point_3Dto2D(mct[iNode],iPlane).w/0.3;
                thiscol[iPlane] = hit_image_v[iPlane].meta().col(wireProj);
                //if(iNode<mct.size()-1)MChits[iPlane].set_pixel(thisrow,thiscol[iPlane],mct.dQdx()[iNode][iPlane]);
                //else MChits[iPlane].set_pixel(thisrow,thiscol[iPlane],mct.dQdx()[iNode-1][iPlane]);
                //for a given Image2D (a given iPlane), we loop over the image and search for the closest non-zero pixel
                for(size_t icol = 0;icol<hit_image_v[iPlane].meta().cols();icol++){
                    for(int irow=0;irow<hit_image_v[iPlane].meta().rows();irow++){
                        if(hit_image_v[iPlane].pixel(irow,icol)<=0)continue;
                        double thisdist = sqrt(pow(2*(thisrow-irow),2)+pow(thiscol[iPlane]-icol,2)); // factor of 2 in row to acount for the fact that my rows are 2 ticks long
                        if(dist > thisdist)dist = thisdist;
                    }
                }
            }
            hDistanceMC2Hit->Fill(dist);
        }*/
    }
    
    bool AStarTracker::finalize() {
        gStyle->SetOptStat(1111);
        c2->cd();
        //TF1 *fGaus = new TF1("fGaus","gaus(0)",0,50);
        //fGaus->SetParameters(10,4,2);
        //hDistance2MC->Fit("fGaus","rl","",0,20);
        //hDistance2MC->SetLineWidth(3);
        //hDistance2MC->SetLineStyle(2);
        //hDistance2MCX->SetLineColor(1);
        //hDistance2MCX->SetFillColor(1);
        //hDistance2MCX->SetFillStyle(3003);
        //hDistance2MCY->SetLineColor(2);
        //hDistance2MCZ->SetLineColor(8);
        //if(hDistance2MC->GetBinContent(hDistance2MC->GetMaximumBin())>= hDistance2MCX->GetBinContent(hDistance2MCX->GetMaximumBin())
        //   && hDistance2MC->GetBinContent(hDistance2MC->GetMaximumBin())>= hDistance2MCY->GetBinContent(hDistance2MCY->GetMaximumBin())
        //   && hDistance2MC->GetBinContent(hDistance2MC->GetMaximumBin())>= hDistance2MCZ->GetBinContent(hDistance2MCZ->GetMaximumBin())){
            hDistance2MC->Draw();
        //    hDistance2MCX->Draw("same");
        //    hDistance2MCY->Draw("same");
        //    hDistance2MCZ->Draw("same");
        //}
        //else if(hDistance2MCX->GetBinContent(hDistance2MCX->GetMaximumBin())>= hDistance2MC->GetBinContent(hDistance2MC->GetMaximumBin())
        //        && hDistance2MCX->GetBinContent(hDistance2MCX->GetMaximumBin())>= hDistance2MCY->GetBinContent(hDistance2MCY->GetMaximumBin())
        //        && hDistance2MCX->GetBinContent(hDistance2MCX->GetMaximumBin())>= hDistance2MCZ->GetBinContent(hDistance2MCZ->GetMaximumBin())){
        //    hDistance2MCX->Draw();
        //    hDistance2MC->Draw("same");
        //    hDistance2MCY->Draw("same");
        //    hDistance2MCZ->Draw("same");
        //}
        /*else if(hDistance2MCY->GetBinContent(hDistance2MCY->GetMaximumBin())>= hDistance2MCX->GetBinContent(hDistance2MCX->GetMaximumBin())
                && hDistance2MCY->GetBinContent(hDistance2MCY->GetMaximumBin())>= hDistance2MC->GetBinContent(hDistance2MC->GetMaximumBin())
                && hDistance2MCY->GetBinContent(hDistance2MCY->GetMaximumBin())>= hDistance2MCZ->GetBinContent(hDistance2MCZ->GetMaximumBin())){
            hDistance2MCY->Draw();
            hDistance2MC->Draw("same");
            hDistance2MCX->Draw("same");
            hDistance2MCZ->Draw("same");
        }
        else{
            hDistance2MCZ->Draw();
            hDistance2MC->Draw("same");
            hDistance2MCX->Draw("same");
            hDistance2MCY->Draw("same");
        }*/
        c2->SaveAs("plot/distace2MC.pdf");

        hDistance2Hit->Draw();
        c2->SaveAs("plot/distance2Hit.pdf");

        hDistanceMC2Hit->Draw();
        c2->SaveAs("plot/distanceMC2Hit.pdf");

        hdQdx->Divide(hdQdxEntries);
        hdQdx->SetMarkerStyle(20);
        hdQdx->GetXaxis()->SetTitle("Distance to end of track");
        hdQdx->GetYaxis()->SetTitle("dQ/dx");
        hdQdx->Draw();
        c2->SaveAs("plot/dQdXprofile.pdf");
        return true;
    }

    //
    // Correct for Space-Charge Effects
    //
    larlite::track AStarTracker::CorrectSCE(larlite::track thisTrack){
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

}
#endif
