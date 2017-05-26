#ifndef CHIMERAPATCHING_H
#define CHIMERAPATCHING_H

#include <string>
#include <vector>

#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/storage_manager.h"

#include "TH2D.h"
#include "TCanvas.h"
#include "TVector3.h"

class ChimeraPatching{

public:
    ChimeraPatching(){}
    virtual ~ChimeraPatching(){};
    void Initialize();
    void AddTrack(larlite::track newTrack, std::vector<larlite::hit> newHitCluster, double score);
    void DrawEvent();
    void NewEvent(std::string evtID, std::vector<double> parVector);
    void TranslateClusters();
    void SetScoreLimit(double scoremin){_scoreLimit = scoremin;}
    void SetNpartperevent(int Npart){_Npartperevent = Npart;}
    void SetPartTypes(std::vector<std::string> partTypes){_particleTypes = partTypes;}

private:
    std::vector<double> _scores;
    std::vector<larlite::track> _Tracks;
    std::vector<larlite::track> _translatedTracks;
    std::vector<std::vector<larlite::hit> > _HitClusters;
    std::vector<std::vector<larlite::hit> > _translatedHitClusters;
    TCanvas *cEventImage;
    TH2D *hEventImage[3];
    std::string _evtID;
    std::vector<std::string> _particleTypes;
    TVector3 _X0;
    std::vector<double> _L0;
    std::vector<double> _theta0;
    std::vector<double> _phi0;
    double _scoreLimit;
    int _Npartperevent;
};
#endif
