#ifndef LARLITE_CHIMERAPATCHING_H
#define LARLITE_CHIMERAPATCHING_H
// std C++
//#include <string>
// larlite
#include "Analysis/ana_base.h"
#include "ChimeraTrackEvaluator.h"
//root
#include "TH1D.h"

namespace larlite {
    /**
     \class ChimeraPatching
     User custom analysis class made by SHELL_USER_NAME
     */
    class ChimeraPatching : public ana_base{

    public:

        /// Default constructor
        ChimeraPatching(){
            _name="ChimeraPatching";
            _fout=0;
            _track_producer="pandoraNu";
            _hit_producer="gaushit";

        }

        /// Default destructor
        virtual ~ChimeraPatching(){}

        /** Initialization method to be called before the analysis event loop.
         */
        virtual bool initialize();

        /** Analyze a data event-by-event
         */
        virtual bool analyze(storage_manager* storage);

        /** Finalize method to be called after all events processed.
         */
        virtual bool finalize();

    protected:
        std::string _track_producer;
        std::string _hit_producer;

        ChimeraTrackEvaluator _evaluator;

        TH1D *hX0;
        TH1D *hY0;
        TH1D *hZ0;
        TH1D *hL0;
        TH1D *hTheta0;
        TH1D *hPhi0;

        TH1D *hX;
        TH1D *hY;
        TH1D *hZ;
        TH1D *hL;
        TH1D *hTheta;
        TH1D *hPhi;

        TH1D *hScore;

    };
}
#endif
