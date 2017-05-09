#ifndef LARLITE_CHIMERATRACKEVALUATOR_CXX
#define LARLITE_CHIMERATRACKEVALUATOR_CXX

#include "ChimeraTrackEvaluator.h"
#include "DataFormat/track.h"

#include <TVector3.h>

namespace larlite {

    void ChimeraTrackEvaluator::SetTargets(double *par){
        _X0     = par[0];
        _Y0     = par[1];
        _Z0     = par[2];
        _L0     = par[3];
        _theta0 = par[4];
        _phi0   = par[5];
    }
    //______________________________________________________
    std::vector<double> ChimeraTrackEvaluator::GetTargets(){
        std::vector<double> targets(6);
        targets.push_back(_X0);
        targets.push_back(_Y0);
        targets.push_back(_Z0);
        targets.push_back(_L0);
        targets.push_back(_theta0);
        targets.push_back(_phi0);
        return targets;
    }
    //------------------------------------------------------
    void ChimeraTrackEvaluator::SetSigmas(double *par){
        _sigma_R      = par[0];
        _sigma_length = par[1];
        _sigma_angles = par[2];
    }
    //______________________________________________________
    std::vector<double> ChimeraTrackEvaluator::GetSigmas(){
        std::vector<double> sigmas(3);
        sigmas.push_back(_sigma_R);
        sigmas.push_back(_sigma_length);
        sigmas.push_back(_sigma_angles);
        return sigmas;
    }
    //______________________________________________________
    void ChimeraTrackEvaluator::SetEval(double *parTarget, double *parSigma){
        SetTargets(parTarget);
        SetSigmas(parSigma);
    }
    //______________________________________________________
    double ChimeraTrackEvaluator::EvalTrack(larlite::track thisTrack, std::string particleType){
        double score;
        if(particleType == "proton"){score = EvalTrackProton(thisTrack);}
        else if(particleType == "muon"){score = EvalTrackMuon(thisTrack);}
        else{std::cout << "[ChimeraTrackEvaluator::EvalTrack]...ERROR = unknown particle type" << std::endl;return 0;}
        return score;
    }
    double ChimeraTrackEvaluator::EvalTrackProton(larlite::track thisTrack){
        double X0 = thisTrack.Vertex().X();
        double Y0 = thisTrack.Vertex().Y();
        double Z0 = thisTrack.Vertex().Z();
        double L0 = thisTrack.Length(0);
        double theta = thisTrack.DirectionAtPoint(0).Theta();
        double phi   = thisTrack.DirectionAtPoint(0).Phi();

        if(_sigma_R!=0     ){_score_R      = exp(-1.*(pow(X0-_X0,2)+pow(Y0-_Y0,2)+pow(Z0-_Z0,2))/(2*_sigma_R*_sigma_R));}else{_score_R=1;}
        if(_sigma_length!=0){_score_length = exp(-1.*pow(L0-_L0,2)/(2*pow(_sigma_length*_L0,2)));                       }else{_score_length=1;}
        if(_sigma_angles!=0){_score_theta  = exp(-1.*pow(theta-_theta0,2)/(2*_sigma_angles*_sigma_angles));             }else{_score_theta=1;}
        if(_sigma_angles!=0){_score_phi    = exp(-1.*pow(phi-_phi0,2)/(2*_sigma_angles*_sigma_angles));                 }else{_score_phi=1;}
        _score = _score_R*_score_length*_score_theta*_score_phi;

        return _score;
    }
    double ChimeraTrackEvaluator::EvalTrackMuon(larlite::track thisTrack){
        //
        // This is assuming that the vertex in the muon track is actually the end-point
        double X0 = thisTrack.Vertex().X();
        double Y0 = thisTrack.Vertex().Y();
        double Z0 = thisTrack.Vertex().Z();
        double L0 = thisTrack.Length(0);
        double theta = thisTrack.DirectionAtPoint(0).Theta();
        double phi   = thisTrack.DirectionAtPoint(0).Phi();

        for(int itrackpoint=0;itrackpoint<thisTrack.NumberTrajectoryPoints();itrackpoint++){
            // 1) starting from vertex, follow the track until the length is what I want
            // 2) consider this the X0
            if(thisTrack.Length(0)-thisTrack.Length(itrackpoint)<_L0){
                X0 = thisTrack.LocationAtPoint(itrackpoint).X();
                Y0 = thisTrack.LocationAtPoint(itrackpoint).Y();
                Z0 = thisTrack.LocationAtPoint(itrackpoint).Z();
                L0 = thisTrack.Length(0)-thisTrack.Length(itrackpoint);
                TVector3 direction(-1.*thisTrack.DirectionAtPoint(itrackpoint).X(),-1.*thisTrack.DirectionAtPoint(itrackpoint).Y(),-1.*thisTrack.DirectionAtPoint(itrackpoint).Z());
                theta = direction.Theta();
                phi = direction.Phi();
            }
            else break;
        }

        if(_sigma_R!=0     ){_score_R      = exp(-1.*(pow(X0-_X0,2)+pow(Y0-_Y0,2)+pow(Z0-_Z0,2))/(2*_sigma_R*_sigma_R));}else{_score_R=1;}
        if(_sigma_length!=0){_score_length = exp(-1.*pow(L0-_L0,2)/(2*pow(_sigma_length*_L0,2)));                       }else{_score_length=1;}
        if(_sigma_angles!=0){_score_theta  = exp(-1.*pow(theta-_theta0,2)/(2*_sigma_angles*_sigma_angles));             }else{_score_theta=1;}
        if(_sigma_angles!=0){_score_phi    = exp(-1.*pow(phi-_phi0,2)/(2*_sigma_angles*_sigma_angles));                 }else{_score_phi=1;}
        _score = _score_R*_score_length*_score_theta*_score_phi;

        return _score;
    }
    //______________________________________________________
}

#endif
