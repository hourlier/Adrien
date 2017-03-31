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
#include "TH2D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/SCE/SpaceChargeMicroBooNE.h"

namespace larlite {

    bool AStarTrackerDiagnostic::initialize() {
        return true;
    }

    bool AStarTrackerDiagnostic::analyze(storage_manager* storage) {
        auto ev_track    = storage->get_data<event_track>(_track_producer);       // tracks
        auto ev_chstatus = storage->get_data<event_chstatus>(_chstatus_producer); // chstatus
        auto ev_mct      = storage->get_data<event_mctrack>(_mctrack_producer);   // mctracks

        _run = storage->run_id();
        _subrun = storage->subrun_id();
        _event = storage->event_id();

        std::cout << _run << "\t" << _subrun << "\t" << _event << "\t";

        // Check validity
        if(!ev_track){      throw DataFormatException("Could not locate event_track data product!");}
        if(!ev_chstatus){   throw DataFormatException("Could not locate event_chstatus data product!");}
        if(ev_track->empty()){std::cout << "No tracks in ev_track" << std::endl; return true;}
        if(!ev_mct){        throw DataFormatException("Could not locate event_mctrack data product!");}

        std::cout << std::endl;

        for(auto track:*ev_track){
            if(track.NumberTrajectoryPoints() == 0)continue;
            std::cout << "\t\t" << track.ID() << "\t" << track.Length(0) << std::endl;
        }

        return true;
    }

    bool AStarTrackerDiagnostic::finalize() {
        return true;
    }

}
#endif
