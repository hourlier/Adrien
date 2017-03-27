/**
 * \file AStarTracker.h
 *
 * \ingroup Adrien
 *
 * \brief Class def header for a class AStarTracker
 *
 * @author kazuhiro
 */

/** \addtogroup Adrien

 @{*/

#ifndef LARLITE_AStarTracker_H
#define LARLITE_AStarTracker_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"
#include <TVector3.h>
#include "DataFormat/Image2D.h"
#include "TH1D.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"

namespace larlite {
    /**
     \class AStarTracker
     User custom analysis class made by SHELL_USER_NAME
     */
    class AStarTracker : public ana_base{

    public:

        /// Default constructor
        AStarTracker(){
            _name="AStarTracker";
            _fout=0; _track_producer="dl";
            _chstatus_producer="chstatus";
            _mctrack_producer = "mcreco";
            _speedOffset=-2;
            _rebinTime = 2;
        }

        /// Default destructor
        virtual ~AStarTracker(){}

        /** IMPLEMENT in AStarTracker.cc!
         Initialization method to be called before the analysis event loop.
         */
        virtual bool initialize();

        /** IMPLEMENT in AStarTracker.cc!
         Analyze a data event-by-event
         */
        virtual bool analyze(storage_manager* storage);

        /** IMPLEMENT in AStarTracker.cc!
         Finalize method to be called after all events processed.
         */
        virtual bool finalize();

        void set_producer(std::string track_producer,
                          std::string chstatus_producer)
        { _track_producer = track_producer; _chstatus_producer = chstatus_producer; }

    protected:
        // X[cm] to TPC tick (waveform index) conversion
        double X2Tick(double x, size_t plane) const;
        //TPC tick (waveform index) to X[cm] conversion
        double Tick2X(double tick, size_t plane)const;

        larlite::track Reconstruct(const TVector3& start_pt, const TVector3& end_pt,
                                   const std::vector<larcv::Image2D>& hit_image_v,
                                   const std::vector<larcv::Image2D>& chstatus_image_v);

        void GetdQdX(const std::vector<larcv::Image2D>& hit_image_v, larlite::track thisTrack);

        std::string _track_producer;
        std::string _chstatus_producer;
        std::string _mctrack_producer;
        int _run;
        int _subrun;
        int _event;
        int _track;
        int _rebinTime;
        double _speedOffset;
        TH1D *hdQdx;
        TCanvas *c2;
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
