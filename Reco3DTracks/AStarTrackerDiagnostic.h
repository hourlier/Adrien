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
#include "DataFormat/hit.h"
#include "DataFormat/Image2D.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TVector3.h"
#include "TCanvas.h"

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
            _numPlanes = 3;
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

        void set_producer(std::string track_producer,std::string chstatus_producer){
            _track_producer = track_producer;
            _chstatus_producer = chstatus_producer;
        }

    protected:


        double X2Tick(double x, size_t plane) const;

        larlite::track CorrectSCE();

        std::vector<TVector3> CorrectSCE(std::vector<TVector3> originpath);

        void VisualizeTrack();
        void CompareRecoCorr2MC();
        void CompareReco2Hits();
        void CompareLengths();
        void DrawdQdX();
        void DrawdQdxSurface();
        void DrawPlots();

        std::string _track_producer;
        std::string _chstatus_producer;
        std::string _mctrack_producer;
        std::string _trackID;

        int _run;
        int _subrun;
        int _event;
        int _track;
        int _rebinTime;
        int _numPlanes;
        int chstat[3][3456];

        double _speedOffset;

        larlite::track thisTrack;
        larlite::track thisCorrectedTrack;
        larlite::mctrack thisTrueTrack;

        std::vector<larlite::hit> thisHit_v;
        std::vector<larlite::hit> thisGausHit_v;

        TH2D *hHitImage2D[3];
        TH2D *hLengthComparison;

        TH1D *hDistance2MC;
        TH1D *dist2MC_all;
        TH1D *hdQdx;
        TH1D *hdQdxEntries;

        TGraph *gRecoedTrack[3];
        TGraph *gCorrectedTrack[3];
        TGraph *gTrueTrack[3];

        TCanvas *Window;
        std::vector<TCanvas*> Window_v;
        std::vector<TCanvas*> dist2MC_v;
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
