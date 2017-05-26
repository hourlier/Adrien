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
#include "TLatex.h"
#include "TVector3.h"
#include "TGraph.h"

#include "ChimeraPatching.h"
#include "ChimeraUtils.h"

//double X2Tick(double x, size_t plane);

void ChimeraPatching::Initialize(){
    cEventImage = new TCanvas("cEventImage","cEventImage",2400,1200);
    cEventImage->Divide(1,2);
    cEventImage->cd(1)->Divide(3,1);
    //for(int i = 0;i<3;i++){cEventImage->cd(1)->cd(i+1)->SetLogz();}
    _scoreLimit = 1e-10;
}
//_______________________________________________________________________
void ChimeraPatching::AddTrack(larlite::track newTrack, std::vector<larlite::hit> newHitCluster, double score){
    _Tracks.push_back(newTrack);
    _HitClusters.push_back(newHitCluster);
    _scores.push_back(score);
}
//_______________________________________________________________________
void ChimeraPatching::NewEvent(std::string evtID, std::vector<double> parVector){
    TVector3 X0(parVector[0],parVector[1],parVector[2]);
    _evtID = evtID;
    _Tracks.clear();
    _HitClusters.clear();
    _translatedHitClusters.clear();
    _translatedTracks.clear();
    _X0.SetXYZ(X0.X(),X0.Y(),X0.Z());
    _L0.clear();
    _theta0.clear();
    _phi0.clear();
    for(int i=0;i<_Npartperevent;i++){
        _L0.push_back(parVector.at(3+i*3));
        _theta0.push_back(parVector.at(4+i*3));
        _phi0.push_back(parVector.at(5+i*3));
    }
    _scores.clear();
}
//_______________________________________________________________________
void ChimeraPatching::DrawEvent(){
    bool LowScore = false;
    for(int iscore=0;iscore<_scores.size();iscore++){
        if(_scores[iscore] < _scoreLimit){
            //std::cout << "score too low : " << _scores[iscore] << std::endl;
            LowScore = true;
            break;
        }
    }
    if(LowScore) return;
    else std::cout << "Drawing event " << _evtID << std::endl;
    for(int iscore=0;iscore<_scores.size();iscore++){
        std::cout << "Particle " << iscore+1 << " score : " << _scores[iscore] << std::endl;
    }

    if(_translatedTracks.size() == 0){_translatedTracks = _Tracks;}
    if(_translatedHitClusters.size() == 0){_translatedHitClusters = _HitClusters;}
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
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
        wirebounds[iPlane][0]-=1;
        wirebounds[iPlane][1]+=1;

        int marginwire = (int)(0.5*(wirebounds[iPlane][1]-wirebounds[iPlane][0]));
        wirebounds[iPlane][0] -= marginwire;
        wirebounds[iPlane][1] += marginwire;
    }
    timebounds[0]-=1;
    timebounds[1]+=1;
    int margintime = (int)(0.5*(timebounds[1]-timebounds[0]));
    timebounds[0]-=margintime;
    timebounds[1]+=margintime;

    for(size_t iPlane=0;iPlane<3;iPlane++){
        //std::cout << wirebounds[iPlane][0] << " => " << wirebounds[iPlane][1] << "\t" << timebounds[0] << " => " << timebounds[1] << std::endl;
        if(wirebounds[iPlane][0] >= wirebounds[iPlane][1] || timebounds[0] >= timebounds[1]){
            std::cout << "ERROR! problem with the boundaries" << std::endl;
            return;
        }
    }

    std::vector<std::pair<double, double> > X0Proj(3);
    TGraph gX0[3];
    std::vector<std::vector<TGraph> > gTrackPoint(3);

    for(size_t iPlane = 0;iPlane<3;iPlane++){

        gTrackPoint[iPlane].resize(_translatedTracks.size());

        X0Proj[iPlane].first  = X2Tick(_X0.X(),iPlane);
        X0Proj[iPlane].second = larutil::GeometryHelper::GetME()->Point_3Dto2D(_X0,iPlane).w / 0.3;
        gX0[iPlane] = TGraph(1);
        gX0[iPlane].SetPoint(0,X0Proj[iPlane].second,X0Proj[iPlane].first);
        gX0[iPlane].SetMarkerStyle(20);

        for(int iTrack = 0;iTrack<_translatedTracks.size();iTrack++){
            gTrackPoint[iPlane][iTrack] = TGraph();
            for(int ipoint = 0;ipoint<_translatedTracks[iTrack].NumberTrajectoryPoints();ipoint++){
                gTrackPoint[iPlane][iTrack].SetPoint(ipoint,larutil::GeometryHelper::GetME()->Point_3Dto2D(_translatedTracks[iTrack].LocationAtPoint(ipoint),iPlane).w / 0.3,X2Tick(_translatedTracks[iTrack].LocationAtPoint(ipoint).X(),iPlane));
            }
        }

        hEventImage[iPlane] = new TH2D(Form("hEventImage_%s_%zu",_evtID.c_str(),iPlane),Form("hEventImage_%s_%zu",_evtID.c_str(),iPlane),wirebounds[iPlane][1]-wirebounds[iPlane][0],wirebounds[iPlane][0],wirebounds[iPlane][1],timebounds[1]-timebounds[0],timebounds[0],timebounds[1]);
        /*for(int icol=0;icol<hEventImage[iPlane]->GetNbinsX();icol++){
            for(int irow=0;irow<hEventImage[iPlane]->GetNbinsY();irow++){
                hEventImage[iPlane]->SetBinContent(icol+1,irow+1,10);
            }
        }*/

        for(size_t iCluster = 0;iCluster<_translatedHitClusters.size();iCluster++){
            for(auto h:_translatedHitClusters[iCluster]){
                if(h.WireID().Plane!=iPlane)continue;
                //hEventImage[iPlane]->SetBinContent(h.WireID().Wire+1-wirebounds[iPlane][0],h.PeakTime()+1-timebounds[0],h.Integral());
                hEventImage[iPlane]->SetBinContent(h.WireID().Wire+1-wirebounds[iPlane][0],h.PeakTime()+1-timebounds[0],iCluster+0.1);
            }
        }
        cEventImage->cd(1)->cd(iPlane+1);
        hEventImage[iPlane]->Rebin2D(1,1);
        hEventImage[iPlane]->Draw("colz");
        gX0[iPlane].Draw("P");
        for(int iTrack = 0;iTrack<_translatedTracks.size();iTrack++){
            gTrackPoint[iPlane][iTrack].SetMarkerStyle(7);
            gTrackPoint[iPlane][iTrack].Draw("P");
        }
    }

    cEventImage->cd(2);
    cEventImage->cd(2)->Clear();
    TLatex tex;
    tex.SetTextSize(0.07);
    tex.SetTextAlign(11);
    std::vector<double> Xtext(_Npartperevent);
    double Y0 = 0.9;
    double dY = 0.09;
    for(int ipart = 0;ipart<_Npartperevent;ipart++){
        Xtext[ipart] = 0.05+ipart*(0.95-0.05)/_Npartperevent;

        tex.DrawLatex(Xtext[ipart],Y0,Form("%s track : ",_particleTypes[ipart].c_str()));

        tex.DrawLatex(Xtext[ipart],Y0-dY,Form("vertex: (%.1f, %.1f, %.1f)",_Tracks[ipart].Vertex().X(), _Tracks[ipart].Vertex().Y(),_Tracks[ipart].Vertex().Z()));
        tex.DrawLatex(Xtext[ipart],Y0-2.*0.9*dY,Form("#color[2]{vertex: (%.1f, %.1f, %.1f)}",_X0.X(), _X0.Y(),_X0.Z()));

        tex.DrawLatex(Xtext[ipart],Y0-3*dY,Form("(#theta, #phi) : (%.1f#circ, %.1f#circ)",_Tracks[ipart].Vertex().Theta()*180/3.14159, _Tracks[ipart].Vertex().Phi()*180/3.14159));
        tex.DrawLatex(Xtext[ipart],Y0-3.9*dY,Form("#color[2]{(#theta, #phi) : (%.1f#circ, %.1f#circ)}",_theta0[ipart]*180/3.14159,_phi0[ipart]*180/3.14159));

        tex.DrawLatex(Xtext[ipart],Y0-5*dY,Form("Length: %.1f cm",_Tracks[ipart].Length(0)));
        tex.DrawLatex(Xtext[ipart],Y0-5.9*dY,Form("#color[2]{Length: %.1f cm}",_L0[ipart]));

        tex.DrawLatex(Xtext[ipart],Y0-7*dY,Form("#DeltaR : %.f cm",sqrt(pow(_Tracks[ipart].Vertex().X()-_X0.X(),2)+pow(_Tracks[ipart].Vertex().Y()-_X0.Y(),2)+pow(_Tracks[ipart].Vertex().Z()-_X0.Z(),2))));
        tex.DrawLatex(Xtext[ipart],Y0-7.9*dY,Form("score : %.2e", _scores[ipart]));
    }


    cEventImage->SaveAs(Form("ChimeraEvent_%s.png",_evtID.c_str()));
}
//_______________________________________________________________________
void ChimeraPatching::TranslateClusters(){
    // I shoudl also translate 3D points in the track...
    _translatedHitClusters.resize(_HitClusters.size());
    _translatedTracks.clear();
    _translatedTracks.resize(_Tracks.size());

    for(size_t iTrack = 0;iTrack<_HitClusters.size();iTrack++){
        _translatedHitClusters[iTrack].resize(_HitClusters[iTrack].size());
        std::vector<std::pair<double, double> > X0Proj(3);
        std::vector<std::pair<double, double> > VertexProj(3);
        if(_Tracks.at(iTrack).NumberTrajectoryPoints() == 0)return;
        double dt[3],dw[3];
        TVector3 vertexPoint = _Tracks[iTrack].Vertex();

        for(size_t iPlane = 0;iPlane<3;iPlane++){
            X0Proj[iPlane].first  = X2Tick(_X0.X(),iPlane);
            X0Proj[iPlane].second = larutil::GeometryHelper::GetME()->Point_3Dto2D(_X0,iPlane).w / 0.3;
            VertexProj[iPlane].first  = X2Tick(vertexPoint.X(),iPlane);
            VertexProj[iPlane].second = larutil::GeometryHelper::GetME()->Point_3Dto2D(vertexPoint,iPlane).w / 0.3;

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
        // Translate tracks
        //TVector3 newLocationAtPoint;
        TVector3 dR = _X0-_Tracks[iTrack].Vertex();
        for(int ipoint = 0;ipoint<_Tracks[iTrack].NumberTrajectoryPoints();ipoint++){
            _translatedTracks[iTrack].add_vertex(_Tracks[iTrack].LocationAtPoint(ipoint)+dR);
            _translatedTracks[iTrack].add_direction(_Tracks[iTrack].DirectionAtPoint(ipoint));
        }
    }
}
//_______________________________________________________________________

#endif
