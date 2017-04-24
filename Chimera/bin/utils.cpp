#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "DataFormat/opflash.h"
#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/storage_manager.h"

#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"

#include "../ChimeraPatching.h"
#include "../ChimeraTrackEvaluator.h"
#include "utils.h"

#include "TVector3.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"

void DrawTrack(const std::vector< std::vector<larlite::hit> >& HitCluster, std::vector<TVector3> points, std::string evtID){

    TCanvas *cTrack = new TCanvas(Form("cTrack_%s",evtID.c_str()),Form("cTrack_%s",evtID.c_str()),1200,300);
    cTrack->Divide(3,1);
    TH2D *hTrack[3];
    int timebounds[2] = {(int)(1e9),0};
    int wirebounds[3][2];
    for(size_t iPlane=0;iPlane<3;iPlane++){
        wirebounds[iPlane][0] = (int)(1e9) ;
        wirebounds[iPlane][1] = 0;
        for(size_t iCluster = 0;iCluster < HitCluster.size();iCluster++){
            for(larlite::hit h:HitCluster[iCluster]){
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

    for(size_t iPlane = 0;iPlane<3;iPlane++){
        hTrack[iPlane] = new TH2D(Form("hTrack_%s_%zu",evtID.c_str(),iPlane),Form("hTrack_%s_%zu;wire;tick",evtID.c_str(),iPlane),wirebounds[iPlane][1]-wirebounds[iPlane][0],wirebounds[iPlane][0],wirebounds[iPlane][1],timebounds[1]-timebounds[0],timebounds[0],timebounds[1]);

        for(size_t iCluster = 0;iCluster < HitCluster.size();iCluster++){
            for(larlite::hit h:HitCluster[iCluster]){
                if(h.WireID().Plane!=iPlane)continue;
                hTrack[iPlane]->SetBinContent(h.WireID().Wire+1-wirebounds[iPlane][0],h.PeakTime()+1-timebounds[0],h.Integral());
            }
        }
        cTrack->cd(iPlane+1);
        hTrack[iPlane]->Draw("colz");
    }
    size_t nPoints = points.size();
    TGraph *gX0[nPoints][3];
    for(size_t iPoint = 0;iPoint < nPoints;iPoint++){
        auto X0=points[iPoint];
        for(size_t iPlane=0;iPlane<3;iPlane++){
            gX0[iPoint][iPlane] = new TGraph();
            gX0[iPoint][iPlane]->SetPoint(0,larutil::GeometryHelper::GetME()->Point_3Dto2D(X0,iPlane).w / 0.3,X2Tick(X0.X(),iPlane));
            //gX0[iPoint][iPlane]->SetMarkerStyle(7);
            cTrack->cd(iPlane+1);
            gX0[iPoint][iPlane]->Draw("P");
        }
    }

    cTrack->SaveAs(Form("BestTracks_%s.pdf",evtID.c_str()));
}

std::vector<larlite::hit> TranslateHit(std::vector<larlite::hit>& HitCluster, const TVector3 X0, const TVector3 thisvertex){
    std::vector<larlite::hit> newHitCluster(HitCluster.size());
    std::vector<std::pair<double, double> > X0Proj(3);
    std::vector<std::pair<double, double> > VertexProj(3);
    double dt[3],dw[3];

    for(size_t iPlane = 0;iPlane<3;iPlane++){
        X0Proj[iPlane].first  = X2Tick(X0.X(),iPlane);
        X0Proj[iPlane].second = larutil::GeometryHelper::GetME()->Point_3Dto2D(X0,iPlane).w / 0.3;
        VertexProj[iPlane].first  = X2Tick(thisvertex.X(),iPlane);
        VertexProj[iPlane].second = larutil::GeometryHelper::GetME()->Point_3Dto2D(thisvertex,iPlane).w / 0.3;

        dt[iPlane] = VertexProj[iPlane].first-X0Proj[iPlane].first;
        dw[iPlane] = VertexProj[iPlane].second-X0Proj[iPlane].second;
    }
    // Translate hits
    double new_peaktime, new_sigmapeaktime;
    larlite::geo::WireID wid;
    unsigned int new_wire;
    for(size_t iHit = 0;iHit<HitCluster.size();iHit++){
        auto h = HitCluster[iHit];
        int iPlane = h.WireID().Plane;
        new_peaktime      = h.PeakTime()-dt[iPlane];
        new_sigmapeaktime = h.SigmaPeakTime();
        wid = h.WireID();
        new_wire = wid.Wire-dw[iPlane];
        wid.Wire = new_wire;
        h.set_time_peak(new_peaktime,new_sigmapeaktime);
        h.set_wire(wid);
        newHitCluster[iHit] = h;
    }
    return newHitCluster;
}

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

std::vector<std::vector<double> >ReadTargetFile(std::string filename){
    std::vector<std::vector<double>> FullTargets;
    std::ifstream file(filename);
    if(!file){std::cout << "ERROR, could not open file" << std::endl; return FullTargets;}
    double X0,Y0,Z0,L_p,theta_p,phi_p,L_u,theta_u,phi_u;
    char coma;
    std::vector<double> eventTarget(9);
    bool goOn = true;
    while(goOn && FullTargets.size() < 10){
        file >> X0 >> coma >> Y0 >> coma >> Z0 >> coma >> L_p >> coma >> theta_p >> coma >> phi_p >> coma >> L_u >> coma >> theta_u >> coma >> phi_u;
        eventTarget[0] = X0;
        eventTarget[1] = Y0;
        eventTarget[2] = Z0;
        eventTarget[3] = L_p;
        eventTarget[4] = theta_p;
        eventTarget[5] = phi_p;
        eventTarget[6] = L_u;
        eventTarget[7] = theta_u;
        eventTarget[8] = phi_u;
        FullTargets.push_back(eventTarget);
        for(auto ipar:eventTarget){std::cout << ipar << "\t";}
        std::cout << std::endl;
        //std::cout << X0 <<" " << Y0 << std::endl;
        if(file.eof()){goOn=false;std::cout << "end reached" << std::endl;break;}
    }
    return FullTargets;
}
