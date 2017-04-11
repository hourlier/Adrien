#ifndef LARLITE_CHIMERATRACKEVALUATOR_CXX
#define LARLITE_CHIMERATRACKEVALUATOR_CXX

#include "ChimeraTrackEvaluator.h"
#include "DataFormat/track.h"

namespace larlite {

    void ChimeraTrackEvaluator::SetTargets(double *par){
        _X0     = par[0];
        _Y0     = par[1];
        _Z0     = par[2];
        _L0     = par[3];
        _theta0 = par[4];
        _phi0   = par[5];
    }
    //------------------------------------------------------
    void ChimeraTrackEvaluator::SetSigmas(double *par){
        _sigma_R      = par[0];
        _sigma_length = par[1];
        _sigma_angles = par[2];
    }
    //______________________________________________________
    void ChimeraTrackEvaluator::SetEval(double *parTarget, double *parSigma){
        SetTargets(parTarget);
        SetSigmas(parSigma);
    }
    //______________________________________________________
    double ChimeraTrackEvaluator::EvalTrack(larlite::track thisTrack){
        double X0 = thisTrack.Vertex().X();
        double Y0 = thisTrack.Vertex().Y();
        double Z0 = thisTrack.Vertex().Z();
        double L0 = thisTrack.Length(0);
        double theta = thisTrack.Theta();
        double phi   = thisTrack.Phi();

        _score_R      = exp(-1.*(pow(X0-_X0,2)+pow(Y0-_Y0,2)+pow(Z0-_Z0,2))/(2*_sigma_R*_sigma_R));
        _score_length = exp(-1.*pow(L0-_L0,2)/(2*pow(_sigma_length*_L0,2)));
        _score_theta  = exp(-1.*pow(theta-_theta0,2)/(2*_sigma_angles*_sigma_angles));
        _score_phi    = exp(-1.*pow(phi-_phi0,2)/(2*_sigma_angles*_sigma_angles));
        _score = _score_R*_score_length*_score_theta*_score_phi;

        return _score;
    }
    //______________________________________________________
}

#endif
