#ifndef CHIMERAFINDING_CXX
#define CHIMERAFINDING_CXX

#include <string>
#include <fstream>

#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/storage_manager.h"

#include "ChimeraFinding.h"
#include "ChimeraTrackEvaluator.h"

bool ChimeraFinding::initialize(){

    ReadTrackFile();
    //ReadTargetFile();

    NrequestedEvents = FullTargetParameters.size();
    std::cout << NrequestedEvents << " events requested" << std::endl;
    _sigmaPar[0] =  10;                //σ_R
    _sigmaPar[1] = 0.1;                //σ_L
    _sigmaPar[2] = 10*TMath::Pi()/180.;//σ_angles

    _max_score.resize(NrequestedEvents,0);
    _best_run.resize(NrequestedEvents);
    _best_subrun.resize(NrequestedEvents);
    _best_event.resize(NrequestedEvents);
    _best_track_id.resize(NrequestedEvents);
    _best_Track.resize(NrequestedEvents);
    _best_HitCluster.resize(NrequestedEvents);

    return true;
}
//______________________________________________________________________________________________________
bool ChimeraFinding::analyze(larlite::storage_manager &storage){
    auto ev_track = storage.get_data<larlite::event_track>(_TrackGenerator); // tracks
    if(!ev_track)  throw larlite::DataFormatException("Could not locate event_track data product!");
    _run = storage.run_id();
    _subrun = storage.subrun_id();
    _event = storage.event_id();
    //std::cout << _run << "\t" << _subrun << "\t" << _event << std::endl;

    larlite::event_hit* ev_hit=nullptr;
    auto const& track_to_hit = storage.find_one_ass(ev_track->id(), ev_hit, ev_track->id().second);
    if(!ev_hit) throw larlite::DataFormatException("Could not find associated hit data product!");

    std::vector<double> score(NrequestedEvents,0);
    
    for(size_t track_index=0; track_index < ev_track->size(); ++track_index) {
        _track = track_index;
        if(_particleType == "proton"){
            if(!IsGoodTrack())continue;
        }

        auto const& track = (*ev_track)[track_index];
        if(track.NumberTrajectoryPoints() == 0)continue; // no poitns in track, not an interesting track
        //
        // Test that the track is a valid proton or muon track

        //
        // For a given track, loop over the number of events requested
        // and evaluate the track for each target event
        for(int itarget = 0;itarget<NrequestedEvents;itarget++){
            std::vector<larlite::hit> thistrackhits;
            //
            // Evaluate track score
            //std::cout << "Need to set the target parameters! cannot work without it!!!!" << std::endl;
            double parTargets[6];
            parTargets[0] = FullTargetParameters[itarget][0];
            parTargets[1] = FullTargetParameters[itarget][1];
            parTargets[2] = FullTargetParameters[itarget][2];
            if(_particleType == "proton"){
                parTargets[3] = FullTargetParameters[itarget][3];
                parTargets[4] = FullTargetParameters[itarget][4];
                parTargets[5] = FullTargetParameters[itarget][5];
            }
            else if(_particleType == "muon"){
                parTargets[3] = FullTargetParameters[itarget][6];
                parTargets[4] = FullTargetParameters[itarget][7];
                parTargets[5] = FullTargetParameters[itarget][8];
            }
            else{
                std::cout << "ERROR : Particle type not specified" << std::endl;
            }
            evaluator.SetTargets(parTargets);
            evaluator.SetSigmas(_sigmaPar);

            score[itarget] = evaluator.EvalTrack(track);
            //std::cout << itarget << "\t" << score[itarget] << "\t" << track.Vertex().X() << std::endl;
            if(score[itarget] > _max_score[itarget]){
                _max_score[itarget]=score[itarget];
                _best_run[itarget] = _run;
                _best_subrun[itarget] = _subrun;
                _best_event[itarget] = _event;
                _best_track_id[itarget] = _track;
                _best_Track[itarget] = track;
                for(auto const& hit_index : track_to_hit[track_index]) {
                    thistrackhits.push_back( (*ev_hit)[hit_index] );
                }
                //std::cout << thistrackhits.size() << std::endl;
                _best_HitCluster[itarget] = thistrackhits;
            }
        }
    }
    return true;
}
//______________________________________________________________________________________________________
bool ChimeraFinding::finalize(){
    /*for(int iTrack = 0; iTrack< _best_Track.size();iTrack++){
        std::cout << "Track #" << iTrack << "  : L = " << _best_Track[iTrack].Length(0) << " cm, X0 = (" << _best_Track[iTrack].Vertex().X() << "," << _best_Track[iTrack].Vertex().Y() << "," << _best_Track[iTrack].Vertex().Z() << "), theta = " << _best_Track[iTrack].Theta()*180/TMath::Pi() << ", Phi = " << _best_Track[iTrack].Phi()*180/TMath::Pi() << "... score = " << _max_score[iTrack] << std::endl;
    }*/
    return true;
}
//______________________________________________________________________________________________________
void ChimeraFinding::ReadTargetFile(){
    FullTargetParameters.clear();
    std::ifstream file(_TargetFile);
    if(!file){std::cout << "ERROR, could not open file" << std::endl;return;}
    int Run, SubRun, Event;
    double X0,Y0,Z0,L_p,theta_p,phi_p,L_u,theta_u,phi_u;
    char coma;
    std::vector<double> eventTarget(9);
    bool goOn = true;
    while(goOn){
        file >> Run >> coma >> SubRun >> coma >> Event >> coma >> X0 >> coma >> Y0 >> coma >> Z0 >> coma >> L_p >> coma >> theta_p >> coma >> phi_p >> coma >> L_u >> coma >> theta_u >> coma >> phi_u;
        eventTarget[0] = X0;
        eventTarget[1] = Y0;
        eventTarget[2] = Z0;
        eventTarget[3] = L_p;
        eventTarget[4] = theta_p;
        eventTarget[5] = phi_p;
        eventTarget[6] = L_u;
        eventTarget[7] = theta_u;
        eventTarget[8] = phi_u;
        FullTargetParameters.push_back(eventTarget);
        /*std::cout << Run << "\t" << SubRun << "\t" << Event << "\t";
        for(auto ipar:eventTarget){std::cout << ipar << "\t";}
        std::cout << std::endl;*/
        if(file.eof()){goOn=false;break;}
    }
    file.close();
}
//______________________________________________________________________________________________________
void ChimeraFinding::ReadTrackFile(){
    if(_particleType == "proton"){ReadProtonTrackFile();}
    else if(_particleType == "muon"){ReadMuonTrackFile();}
    else {std::cout << "ERROR, particle type is not recognized" << std::endl;}

}
//______________________________________________________________________________________________________
void ChimeraFinding::ReadProtonTrackFile(){
    _SelectableTracks.clear();
    std::vector<int> trackinfo(4);
    std::ifstream file(_TrackFile);
    if(!file){std::cout << "ERROR, could not open file of tracks to sort through" << std::endl;return;}
    std::string firstline;
    getline(file, firstline);
    bool goOn = true;
    int Run,SubRun,Event,TrackID,WireUMin,TimeUMin,WireUMax,TimeUMax,WireVMin,TimeVMin,WireVMax,TimeVMax,WireYMin,TimeYMin,WireYMax,TimeYMax;
    double BDTScore;
    while(goOn){
        file >> Run >> SubRun >> Event >> TrackID >> WireUMin >> TimeUMin >> WireUMax >> TimeUMax >> WireVMin >> TimeVMin >> WireVMax >> TimeVMax >> WireYMin >> TimeYMin >> WireYMax >> TimeYMax >> BDTScore;
        trackinfo[0] = Run;
        trackinfo[1] = SubRun;
        trackinfo[2] = Event;
        trackinfo[3] = TrackID;
        _SelectableTracks.push_back(trackinfo);
        if(file.eof()){goOn=false;break;}
    }
}
//______________________________________________________________________________________________________
void ChimeraFinding::ReadMuonTrackFile(){
    _SelectableTracks.clear();
    std::vector<int> trackinfo(4);
    std::ifstream file(_TrackFile);
    std::cout << _TrackFile << std::endl;
    if(!file){std::cout << "ERROR, could not open file of tracks to sort through" << std::endl;return;}

    bool goOn = true;

    int Run,SubRun,Event,decayIdx,Ntracks,trackIndex,Ymin,Tmin,Ymax,Tmax;
    char coma;

    while(goOn){
        file >> Run >> coma >> SubRun >> coma >> Event >> coma >> decayIdx >> coma >> Ntracks;
        if(Ntracks == 0){continue;}
        else{
            for(int iTrack = 0;iTrack < Ntracks;iTrack++){
                file >> coma >> trackIndex >> coma >> Ymin >> coma >> Tmin >> coma >> Ymax >> coma >> Tmax;
                if(file.eof()){goOn = false;return;}
                trackinfo[0] = Run;
                trackinfo[1] = SubRun;
                trackinfo[2] = Event;
                trackinfo[3] = trackIndex;
                _SelectableTracks.push_back(trackinfo);
            }
        }
        if(file.eof()){goOn=false;break;}
    }
}
//______________________________________________________________________________________________________
bool ChimeraFinding::IsGoodTrack(){
    for(auto trackinfo:_SelectableTracks){
        if(_run != trackinfo[0] && _subrun == trackinfo[1] && _event == trackinfo[2] && _track == trackinfo[3]){return true;}
    }
    return false;
}
//______________________________________________________________________________________________________
#endif
