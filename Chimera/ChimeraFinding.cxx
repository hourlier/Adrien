#ifndef CHIMERAFINDING_CXX
#define CHIMERAFINDING_CXX

#include <string>
#include <fstream>

#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/storage_manager.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"

#include "ChimeraFinding.h"
#include "ChimeraTrackEvaluator.h"
#include "ChimeraUtils.h"

#include "TVector3.h"



bool ChimeraFinding::initialize(){

    ReadTrackFile();

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

        auto /*const&*/ track = (*ev_track)[track_index];
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
            parTargets[3] = FullTargetParameters[itarget][3+3*_particleIndex];
            parTargets[4] = FullTargetParameters[itarget][4+3*_particleIndex];
            parTargets[5] = FullTargetParameters[itarget][5+3*_particleIndex];
            evaluator.SetTargets(parTargets);
            evaluator.SetSigmas(_sigmaPar);
            
            score[itarget] = evaluator.EvalTrack(track,_particleType);
            //std::cout << itarget << "\t" << score[itarget] << "\t" << track.Vertex().X() << std::endl;
            if(score[itarget] > _max_score[itarget]){

                _max_score[itarget]=score[itarget];
                _best_run[itarget] = _run;
                _best_subrun[itarget] = _subrun;
                _best_event[itarget] = _event;
                _best_track_id[itarget] = _track;

                if(_particleType == "muon"){
                    //std::cout << "new best muon track" << std::endl;
                    larlite::track muonTrack;
                    int newNpoints = 0;
                    while(track.Length(0)-track.Length(newNpoints)<parTargets[3] && newNpoints<track.NumberTrajectoryPoints()){newNpoints++;}
                    if(newNpoints >= track.NumberTrajectoryPoints()){newNpoints = track.NumberTrajectoryPoints()-1;}
                    //std::cout << "newNpoints = " << newNpoints  << " and original track has " << track.NumberTrajectoryPoints() << std::endl;

                    for(int ipoint = 0;ipoint<newNpoints;ipoint++){
                        //std::cout << newNpoints-ipoint << "/" << track.NumberTrajectoryPoints() << std::endl;
                        muonTrack.add_vertex(track.LocationAtPoint(newNpoints-ipoint));
                        TVector3 newDirection(-track.DirectionAtPoint(newNpoints-ipoint).X(),-track.DirectionAtPoint(newNpoints-ipoint).Y(),-track.DirectionAtPoint(newNpoints-ipoint).Z());
                        muonTrack.add_direction(newDirection);
                    }
                    _best_Track[itarget] = muonTrack;
                    //std::cout << "new muon tracks OK" << std::endl;
                }
                else{
                    _best_Track[itarget] = track;
                }
                if(_particleType == "muon"){
                    for (int ipoint = 1;ipoint<_best_Track[itarget].NumberTrajectoryPoints();ipoint++){
                        for(int iPlane = 0;iPlane<3;iPlane++){
                            double X0 = larutil::GeometryHelper::GetME()->Point_3Dto2D(_best_Track[itarget].LocationAtPoint(ipoint-1),iPlane).w / 0.3;
                            double Y0 = X2Tick(_best_Track[itarget].LocationAtPoint(ipoint-1).X(),iPlane);
                            double X1 = larutil::GeometryHelper::GetME()->Point_3Dto2D(_best_Track[itarget].LocationAtPoint(ipoint),iPlane).w / 0.3;
                            double Y1 = X2Tick(_best_Track[itarget].LocationAtPoint(ipoint).X(),iPlane);
                            for(auto const& hit_index : track_to_hit[track_index]) {
                                if((*ev_hit)[hit_index].WireID().Plane != iPlane)continue;
                                bool pixelProjectsOnNewTrack = false;
                                double X2 = (*ev_hit)[hit_index].WireID().Wire;
                                double Y2 = (*ev_hit)[hit_index].PeakTime();
                                if(((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0))/(pow((X1-X0),2)+pow((Y1-Y0),2)) >=0 && ((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0))/(pow((X1-X0),2)+pow((Y1-Y0),2)) <= 1){pixelProjectsOnNewTrack = true;}
                                if(pixelProjectsOnNewTrack == true)thistrackhits.push_back( (*ev_hit)[hit_index] );
                            }

                        }
                    }

                }
                else{
                    for(auto const& hit_index : track_to_hit[track_index]) {
                        thistrackhits.push_back( (*ev_hit)[hit_index] );
                    }
                }
                /*for(auto const& hit_index : track_to_hit[track_index]) {
                    if(_particleType == "muon"){
                        bool pixelProjectsOnNewTrack = false;
                        int iPlane =(*ev_hit)[hit_index].WireID().Plane;
                        for (int ipoint = 1;ipoint<_best_Track[itarget].NumberTrajectoryPoints();ipoint++){

                            double X1 = larutil::GeometryHelper::GetME()->Point_3Dto2D(_best_Track[itarget].LocationAtPoint(ipoint),iPlane).w / 0.3;
                            double Y1 = X2Tick(_best_Track[itarget].LocationAtPoint(ipoint).X(),iPlane);

                            double X2 = (*ev_hit)[hit_index].WireID().Wire;
                            double Y2 = (*ev_hit)[hit_index].PeakTime();

                            if(((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0))/(pow((X1-X0),2)+pow((Y1-Y0),2)) >=0 && ((X1-X0)*(X2-X0)+(Y1-Y0)*(Y2-Y0))/(pow((X1-X0),2)+pow((Y1-Y0),2)) <= 1){pixelProjectsOnNewTrack = true; break;}
                            }
                        if(pixelProjectsOnNewTrack == true)thistrackhits.push_back( (*ev_hit)[hit_index] );
                    }
                    else{
                        thistrackhits.push_back( (*ev_hit)[hit_index] );
                    }
                }*/
                //std::cout << thistrackhits.size() << std::endl;
                _best_HitCluster[itarget] = thistrackhits;
            }
        }
    }
    return true;
}
//______________________________________________________________________________________________________
bool ChimeraFinding::finalize(){
    return true;
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
