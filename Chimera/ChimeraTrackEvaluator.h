#ifndef LARLITE_CHIMERATRACKEVALUATOR_H
#define LARLITE_CHIMERATRACKEVALUATOR_H

#include "Analysis/ana_base.h"
#include "DataFormat/track.h"

namespace larlite {
    class ChimeraTrackEvaluator{
    public:
        /// Default constructor
        ChimeraTrackEvaluator(){}

        /// Default destructor
        virtual ~ChimeraTrackEvaluator(){}

        /// Set Target parameters for a given track
        void SetTargets(double *par);
        /// Set allowed variation arround targets
        void SetSigmas(double *par);
        void SetEval(double *parTarget, double *parSigma);
        /// Evaluate Track
        double EvalTrack(larlite::track thisTrack);

        double GetScoreR()    {return _score_R;     }
        double GetScoreL()    {return _score_length;}
        double GetScoreTheta(){return _score_theta; }
        double GetScorePhi()  {return _score_phi;   }



    protected:
        larlite::track _protonCandidate;
        larlite::track _muonCandidate;
        std::vector<larlite::track> _coupleCandidate;

        /// track requirement
        double _X0;
        double _Y0;
        double _Z0;
        double _L0;
        double _theta0;
        double _phi0;
        /// how strict are we on the requirements?
        double _sigma_R;
        double _sigma_angles;
        double _sigma_length;
        // how well a given track do?
        double _score_R;
        double _score_theta;
        double _score_phi;
        double _score_length;
        double _score;
    };
}

#endif
