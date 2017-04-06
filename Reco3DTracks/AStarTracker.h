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
#include "DataFormat/mctrack.h"
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
            _fout=0;//TFile::Open("output.root","RECREATE");
            _track_producer="dl";
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

        void set_producer(std::string track_producer,std::string chstatus_producer){ _track_producer = track_producer; _chstatus_producer = chstatus_producer; }

    protected:
        // X[cm] to TPC tick (waveform index) conversion
        double X2Tick(double x, size_t plane) const;
        //TPC tick (waveform index) to X[cm] conversion
        double Tick2X(double tick, size_t plane)const;

        larlite::track Reconstruct(const TVector3& start_pt, const TVector3& end_pt,
                                   const std::vector<larcv::Image2D>& hit_image_v,
                                   const std::vector<larcv::Image2D>& chstatus_image_v);

        void CompareReco2MC3D(const larlite::track recoTrack, const larlite::mctrack trueTrack);
        void CompareReco2hits(const std::vector<larcv::Image2D> hitlist);
        void DrawdQdX(larlite::track thisTrack);

        larlite::track MakeTrack();
        larlite::track ComputedQdX(larlite::track newTrack, const std::vector<larcv::Image2D>& hit_image_v);
        larlite::track CorrectSCE(larlite::track thisTrack);
        std::vector<TVector3> CorrectSCE(larlite::mctrack thisTrack);
        std::vector<TVector3> CorrectSCE(std::vector<TVector3> thisTrack);

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
        TH1D *hdQdxEntries;
        TH1D *hDistance2MC;
        TH1D *hDistance2MCX;
        TH1D *hDistance2MCY;
        TH1D *hDistance2MCZ;
        TH1D *hDistance2Hit;
        TH1D *hDistanceMC2Hit;
        std::vector<larlitecv::AStar3DNode> RecoedPath;
        std::vector<TVector3> CorrectedPath;
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
