#include <iostream>
#include <string>
#include <vector>

#include "DataFormat/opflash.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/storage_manager.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"

#include "../ChimeraPatching.h"
#include "../ChimeraTrackEvaluator.h"
//#include "utils.cpp"

#include "TVector3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"

void DrawTrack(const std::vector< std::vector<larlite::hit> >& HitCluster, std::vector<TVector3> points, std::string evtID);
std::vector<larlite::hit> TranslateHit(std::vector<larlite::hit>& HitCluster, const TVector3 X0, const TVector3 thisvertex);
double X2Tick(double x, size_t plane);
std::vector<std::vector<double> >ReadTargetFile(std::string filename);
