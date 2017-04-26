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
    void AddTrack(larlite::track newTrack, std::vector<larlite::hit> newHitCluster);
    void DrawEvent();
    void NewEvent(std::string evtID, TVector3 X0);
    void TranslateClusters();

private:
    std::vector<larlite::track> _Tracks;
    std::vector<std::vector<larlite::hit> > _HitClusters;
    std::vector<std::vector<larlite::hit> > _translatedHitClusters;
    TCanvas *cEventImage;
    TH2D *hEventImage[3];
    std::string _evtID;
    TVector3 _X0;
};
#endif
