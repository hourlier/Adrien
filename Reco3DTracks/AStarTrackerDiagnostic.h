/**
 * \file AStarTrackerDiagnostic.h
 *
 * \ingroup Adrien
 *
 * \brief Class def header for a class AStarTrackerDiagnostic
 *
 * @author hourlier
 */

/** \addtogroup Adrien

 @{*/

#ifndef LARLITE_AStarTrackerDiagnostic_H
#define LARLITE_AStarTrackerDiagnostic_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include "DataFormat/mctrack.h"
#include <TVector3.h>
#include "DataFormat/Image2D.h"
#include "TH1D.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"

namespace larlite {
    /**
     \class AStarTrackerDiagnostic
     User custom analysis class made by SHELL_USER_NAME
     */
    class AStarTrackerDiagnostic : public ana_base{

    public:

        /// Default constructor
        AStarTrackerDiagnostic(){
            _name="AStarTracker";
            _fout=0; _track_producer="dl";
            _chstatus_producer="chstatus";
            _mctrack_producer = "mcreco";
            _speedOffset=-2;
            _rebinTime = 2;
        }

        /// Default destructor
        virtual ~AStarTrackerDiagnostic(){}

        /** IMPLEMENT in AStarTrackerDiagnostic.cxx!
         Initialization method to be called before the analysis event loop.
         */
        virtual bool initialize();

        /** IMPLEMENT in AStarTrackerDiagnostic.cxx!
         Analyze a data event-by-event
         */
        virtual bool analyze(storage_manager* storage);

        /** IMPLEMENT in AStarTracker.cc!
         Finalize method to be called after all events processed.
         */
        virtual bool finalize();

        void set_producer(std::string track_producer,std::string chstatus_producer){ _track_producer = track_producer; _chstatus_producer = chstatus_producer; }

    protected:
        std::string _track_producer;
        std::string _chstatus_producer;
        std::string _mctrack_producer;
        int _run;
        int _subrun;
        int _event;
        int _track;
        int _rebinTime;
        double _speedOffset;
        
    };
}
#endif

//**************************************************************************
//
// For Analysis framework documentation, read Manual.pdf here:
//
// http://microboone-docdb.fnal.gov:8080/cgi-bin/ShowDocument?docid=3183
//
//**************************************************************************

/** @} */ // end of doxygen group 
