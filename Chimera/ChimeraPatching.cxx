#ifndef CHIMERAPATCHING_CXX
#define CHIMERAPATCHING_CXX

#include <string>
#include <vector>

#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/storage_manager.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"

#include "TH2D.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "ChimeraPatching.h"

double X2Tick(double x, size_t plane);

void ChimeraPatching::Initialize(){
    cEventImage = new TCanvas("cEventImage","cEventImage",1200,300);
    cEventImage->Divide(3,1);
}
//_______________________________________________________________________
void ChimeraPatching::AddTrack(larlite::track newTrack, std::vector<larlite::hit> newHitCluster){
    _Tracks.push_back(newTrack);
    _HitClusters.push_back(newHitCluster);
}
//_______________________________________________________________________
void ChimeraPatching::NewEvent(std::string evtID, TVector3 X0){
    _evtID = evtID;
    _Tracks.clear();
    _HitClusters.clear();
    _X0.SetXYZ(X0.X(),X0.Y(),X0.Z());
}
//_______________________________________________________________________
void ChimeraPatching::DrawEvent(){
    gStyle->SetOptStat(0);
    int timebounds[2] = {1000000,0};
    int wirebounds[3][2];
    for(size_t iPlane=0;iPlane<3;iPlane++){
        wirebounds[iPlane][0] = 1000000 ;
        wirebounds[iPlane][1] = 0;
        for(size_t iCluster = 0;iCluster<_translatedHitClusters.size();iCluster++){
            for(auto h:_translatedHitClusters[iCluster]){
                if(h.WireID().Plane!=iPlane)continue;
                if(h.PeakTime()    < timebounds[0]){timebounds[0] = h.PeakTime();}
                if(h.PeakTime()    > timebounds[1]){timebounds[1] = h.PeakTime();}
                if(h.WireID().Wire < wirebounds[iPlane][0]){wirebounds[iPlane][0] = h.WireID().Wire;}
                if(h.WireID().Wire > wirebounds[iPlane][1]){wirebounds[iPlane][1] = h.WireID().Wire;}
            }
        }
        int marginwire = (int)(0.5*(wirebounds[iPlane][1]-wirebounds[iPlane][0]));
        wirebounds[iPlane][0] -= marginwire;
        wirebounds[iPlane][1] += marginwire;
    }
    int margintime = (int)(0.5*(timebounds[1]-timebounds[0]));
    timebounds[0]-=margintime;
    timebounds[1]+=margintime;

    for(size_t iPlane=0;iPlane<3;iPlane++){
        std::cout << wirebounds[iPlane][0] << " => " << wirebounds[iPlane][1] << "\t" << timebounds[0] << " => " << timebounds[1] << std::endl;
    }

    for(size_t iPlane = 0;iPlane<3;iPlane++){
        hEventImage[iPlane] = new TH2D(Form("hEventImage_%s_%zu",_evtID.c_str(),iPlane),Form("hEventImage_%s_%zu",_evtID.c_str(),iPlane),wirebounds[iPlane][1]-wirebounds[iPlane][0],wirebounds[iPlane][0],wirebounds[iPlane][1],timebounds[1]-timebounds[0],timebounds[0],timebounds[1]);
        /*for(int icol=0;icol<hEventImage[iPlane]->GetNbinsX();icol++){
            for(int irow=0;irow<hEventImage[iPlane]->GetNbinsY();irow++){
                hEventImage[iPlane]->SetBinContent(icol+1,irow+1,0.01);
            }
        }*/
        for(size_t iCluster = 0;iCluster<_translatedHitClusters.size();iCluster++){
            for(auto h:_translatedHitClusters[iCluster]){
                if(h.WireID().Plane!=iPlane)continue;
                hEventImage[iPlane]->SetBinContent(h.WireID().Wire+1-wirebounds[iPlane][0],h.PeakTime()+1-timebounds[0],h.Integral());
            }
        }
        cEventImage->cd(iPlane+1);
        //hEventImage[iPlane]->Rebin2D(5,5);
        hEventImage[iPlane]->Draw("colz");
    }
    cEventImage->SaveAs(Form("ChimeraEvent_%s.pdf",_evtID.c_str()));
}
//_______________________________________________________________________
void ChimeraPatching::TranslateClusters(){
    _translatedHitClusters.resize(_HitClusters.size());
    for(size_t iTrack = 0;iTrack<_HitClusters.size();iTrack++){
        _translatedHitClusters[iTrack].resize(_HitClusters[iTrack].size());
        std::vector<std::pair<double, double> > X0Proj(3);
        std::vector<std::pair<double, double> > VertexProj(3);
        TVector3 thisvertex = _Tracks[iTrack].Vertex();
        double dt[3],dw[3];
        for(size_t iPlane = 0;iPlane<3;iPlane++){
            X0Proj[iPlane].first  = X2Tick(_X0.X(),iPlane);
            X0Proj[iPlane].second = larutil::GeometryHelper::GetME()->Point_3Dto2D(_X0,iPlane).w / 0.3;
            VertexProj[iPlane].first  = X2Tick(thisvertex.X(),iPlane);
            VertexProj[iPlane].second = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisvertex,iPlane).w / 0.3;

            dt[iPlane] = VertexProj[iPlane].first-X0Proj[iPlane].first;
            dw[iPlane] = VertexProj[iPlane].second-X0Proj[iPlane].second;
        }
        // Translate hits
        double new_peaktime, new_sigmapeaktime;
        larlite::geo::WireID wid;
        unsigned int new_wire;
        for(size_t iHit = 0;iHit<_HitClusters[iTrack].size();iHit++){
            auto h = _HitClusters[iTrack][iHit];
            int iPlane = h.WireID().Plane;
            new_peaktime      = h.PeakTime()-dt[iPlane];
            new_sigmapeaktime = h.SigmaPeakTime();
            wid = h.WireID();
            new_wire = wid.Wire-dw[iPlane];
            wid.Wire = new_wire;
            h.set_time_peak(new_peaktime,new_sigmapeaktime);
            h.set_wire(wid);
            _translatedHitClusters[iTrack][iHit] = h;
        }
    }
}
//_______________________________________________________________________
double X2Tick(double x, size_t plane){
    auto ts = larutil::TimeService::GetME();
    auto larp = larutil::LArProperties::GetME();

    // (X [cm] / Drift Velocity [cm/us] - TPC waveform tick 0 offset) ... [us]
    double tick = (x / larp->DriftVelocity() - ts->TriggerOffsetTPC());
    // 1st plane's tick
    if(plane==0) return tick * 2;// 2 ticks/us
    // 2nd plane needs more drift time for plane0=>plane1 gap (0.3cm) difference
    tick -= 0.3 / larp->DriftVelocity(larp->Efield(1));
    // 2nd plane's tick
    if(plane==1) return tick * 2;
    // 3rd plane needs more drift time for plane1=>plane2 gap (0.3cm) difference
    tick -= 0.3 / larp->DriftVelocity(larp->Efield(2));
    return tick * 2;
}

#endif
