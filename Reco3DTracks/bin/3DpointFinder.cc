#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/wire.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/core/DataFormat/ChStatus.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgoProton.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/SCE/SpaceChargeMicroBooNE.h"

//______________________________________________________
// Prototypes
//______________________________________________________
double Tick2X(double tick, size_t plane);
double X2Tick(double x, size_t plane);
TVector3 GenerateOriginPoint();
std::vector<std::pair<int, int> > GetPointProjection(TVector3 point);
std::vector<std::pair<int, int> > GetEndPoint(std::vector<larcv::Image2D> viewOrigin);
std::vector<larcv::Image2D> GetViewOriginPoint(TVector3 point);
void DrawROI(std::vector<larcv::Image2D> viewOrigin, TVector3 PointOrigin, TVector3 newPoint);
TVector3 CheckEndPoints(TVector3 start_pt, std::vector<larcv::Image2D> viewOrigin, std::vector< std::pair<int,int> > endPix);
double EvalMinDist(TVector3 point, std::vector<larcv::Image2D> viewOrigin, std::vector< std::pair<int,int> > endPix);
std::vector<TVector3> GetOpenSet(TVector3 newPoint, double dR);
//______________________________________________________
// Implementations
//______________________________________________________
int main(){
    std::cout << "Hello World" << std::endl;
    TVector3 PointOrigin = GenerateOriginPoint();
    std::cout << "Original point defined : (" << PointOrigin.X() << "," << PointOrigin.Y() << "," << PointOrigin.Z() << ")" << std::endl;
    std::vector<std::pair<int, int> > pointProj = GetPointProjection(PointOrigin);
    for(int iPlane = 0;iPlane<pointProj.size();iPlane++){
        std::cout << "plane : " << iPlane << ", wire : " << pointProj[iPlane].first << ", tick : " << pointProj[iPlane].second << std::endl;
    }
    std::vector<larcv::Image2D> viewOrigin = GetViewOriginPoint(PointOrigin);
    std::vector<std::pair<int, int> > endPix = GetEndPoint(viewOrigin);
    TVector3 newPoint = CheckEndPoints(PointOrigin,viewOrigin,endPix);
    DrawROI(viewOrigin,PointOrigin, newPoint);
}
//______________________________________________________
TVector3 GenerateOriginPoint(){
    TVector3 point;
    TRandom3 *ran = new TRandom3();
    ran->SetSeed(0);

    point.SetX(ran->Uniform(0,260));
    point.SetY(ran->Uniform(-115,115));
    point.SetZ(ran->Uniform(1,1036));

    ran->Delete();
    return point;
}
//______________________________________________________
std::vector<std::pair<int, int> > GetPointProjection(TVector3 point){
    std::vector<std::pair<int, int> > pointProj(3);
    for(int iPlane = 0;iPlane<3;iPlane++){
        std::pair<int, int> planeProj;
        planeProj.first = larutil::GeometryHelper::GetME()->Point_3Dto2D(point,iPlane).w / 0.3;
        planeProj.second = X2Tick(point.X(),iPlane);
        pointProj[iPlane] = planeProj;
    }
    return pointProj;
}
//______________________________________________________
std::vector<larcv::Image2D> GetViewOriginPoint(TVector3 point){
    std::cout << "GetViewOriginPoint" << std::endl;
    std::vector<larcv::Image2D> views(3);
    //
    // Get wire and tick position from original point
    std::vector<std::pair<int, int> > pointProj = GetPointProjection(point);

    size_t ImageMargin = 10;
    for(int iPlane = 0;iPlane<3;iPlane++){
        //
        // Defint meta of the image
        larcv::ImageMeta meta(2.*ImageMargin,
                              2.*ImageMargin,
                              2*ImageMargin,
                              2*ImageMargin,
                              (size_t)(pointProj[iPlane].first-ImageMargin),
                              (size_t)(pointProj[iPlane].second+ImageMargin),
                              iPlane);
        std::cout <<  "plane : " << iPlane << ", meta.cols() : " << meta.cols()  << ", meta.rows() : " << meta.rows() << std::endl;
        std::cout <<  "\t   meta.min_x() : " << meta.min_x() << ", meta.max_x() : " << meta.max_x() << std::endl;
        std::cout <<  "\t   meta.min_y() : " << meta.min_y() << ", meta.max_y() : " << meta.max_y() << std::endl;

        //
        // Fill image data with the corresponding pixel
        std::vector<float> image_data(meta.rows() * meta.cols(), 0.);

        TRandom3 *ran = new TRandom3();
        ran->SetSeed(0);

        size_t col = meta.col(pointProj[iPlane].first)+ran->Gaus(0,2);
        size_t row = meta.row(pointProj[iPlane].second)+ran->Gaus(0,2);
        image_data[meta.rows()*col+row] = 1;
        //
        // Make Image2D
        larcv::Image2D Image(std::move(meta),std::move(image_data));
        views[iPlane] = Image;
    }
    return views;
}
//______________________________________________________
double X2Tick(double x, size_t plane){

    auto ts = larutil::TimeService::GetME();
    auto larp = larutil::LArProperties::GetME();
    double tick = (x / larp->DriftVelocity() - ts->TriggerOffsetTPC());
    return tick * 2;
}
//______________________________________________________
double Tick2X(double tick, size_t plane){
    auto ts = larutil::TimeService::GetME();
    auto larp = larutil::LArProperties::GetME();
    tick/=2.;
    double x = (tick+ts->TriggerOffsetTPC())*larp->DriftVelocity();

    return x;
}
//______________________________________________________
void DrawROI(std::vector<larcv::Image2D> viewOrigin, TVector3 PointOrigin, TVector3 newPoint){
    TH2D *hImage[3];
    TGraph *gPointOrigin[3];
    TGraph *gNewPoint[3];
    std::vector<std::pair<int, int> > pointProj    = GetPointProjection(PointOrigin);
    std::vector<std::pair<int, int> > newpointProj = GetPointProjection(newPoint);

    for(size_t iPlane=0;iPlane<3;iPlane++){
        gPointOrigin[iPlane] = new TGraph();
        gPointOrigin[iPlane]->SetPoint(0,viewOrigin[iPlane].meta().col(pointProj[iPlane].first)+0.5,viewOrigin[iPlane].meta().row(pointProj[iPlane].second)+0.5);
        gNewPoint[iPlane] = new TGraph();
        gNewPoint[iPlane]->SetPoint(0,viewOrigin[iPlane].meta().col(newpointProj[iPlane].first)+0.5,viewOrigin[iPlane].meta().row(newpointProj[iPlane].second)+0.5);

        hImage[iPlane] = new TH2D(Form("hImage_%zu",iPlane),Form("hImage_%zu;wire;time",iPlane),viewOrigin[iPlane].meta().cols(),0,viewOrigin[iPlane].meta().cols(),viewOrigin[iPlane].meta().rows(),0,viewOrigin[iPlane].meta().rows());

        for(int icol=0;icol<viewOrigin[iPlane].meta().cols();icol++){
            for(int irow=0;irow<viewOrigin[iPlane].meta().rows();irow++){
                hImage[iPlane]->SetBinContent(icol+1,irow+1,viewOrigin[iPlane].pixel(irow,icol));
            }
        }
    }

    TCanvas *c = new TCanvas("ROI","ROI",450,150);
    c->Divide(3,1);
    for(int iPlane=0;iPlane<3;iPlane++){
        c->cd(iPlane+1);
        hImage[iPlane]->Draw("colz");
        gPointOrigin[iPlane]->SetMarkerStyle(20);
        gPointOrigin[iPlane]->SetMarkerSize(0.3);
        gPointOrigin[iPlane]->Draw("same P");
        gNewPoint[iPlane]->SetMarkerColor(2);
        gNewPoint[iPlane]->SetMarkerStyle(20);
        gNewPoint[iPlane]->SetMarkerSize(0.3);
        gNewPoint[iPlane]->Draw("same P");
    }
    c->SaveAs(Form("%s.pdf",c->GetName()));
    for(int iPlane = 0;iPlane<3;iPlane++){
        hImage[iPlane]->Delete();
    }
}
//______________________________________________________
TVector3 CheckEndPoints(TVector3 start_pt, std::vector<larcv::Image2D> viewOrigin, std::vector< std::pair<int,int> > endPix){
    //tellMe(Form("start point : (%f,%f,%f)",start_pt.X(),start_pt.Y(),start_pt.Z()),2);
    double dR = 0.1;

    double closeEnough = 6;
    TVector3 newPoint = start_pt;
    double minDist = EvalMinDist(start_pt,viewOrigin,endPix);
    double PreviousMinDist = 2*minDist;
    if(minDist <= closeEnough){std::cout << "Point OK" << std::endl;return newPoint;}
    std::cout <<"Updating point"<<std::endl;
    int iter = 0;
    int iterMax = 100;
    while (minDist > closeEnough && iter<iterMax) {
        iter++;
        std::vector<TVector3> openSet = GetOpenSet(newPoint,dR);
        double minDist = 1e9;
        for(auto neighbour : openSet){
            double dist = EvalMinDist(neighbour,viewOrigin,endPix);
            std::cout << dist << std::endl;
            if(dist < minDist){minDist=dist;newPoint = neighbour;}
        }
        std::cout << "iteration #" << iter << "min dist : " << minDist << std::endl;
        if(minDist < closeEnough || iter >= iterMax || minDist >= PreviousMinDist) return newPoint;
        PreviousMinDist = minDist;
    }
    return newPoint;
}
//______________________________________________________
double EvalMinDist(TVector3 point, std::vector<larcv::Image2D> viewOrigin, std::vector< std::pair<int,int> > endPix){
    //tellMe(Form("EvalMinDist : (%f,%f,%f)",point.X(),point.Y(),point.Z()),2);
    double minDist = 1e9;
    double minDistplane[3] = {1e9,1e9,1e9};
    for(size_t iPlane = 0;iPlane<3;iPlane++){
        //tellMe(Form("plane %zu :",iPlane),2);
        int pointCol = viewOrigin[iPlane].meta().col(larutil::GeometryHelper::GetME()->Point_3Dto2D(point,iPlane).w/0.3);
        int pointRow = viewOrigin[iPlane].meta().row(X2Tick(point.X(),iPlane));
        //tellMe(Form("\t (%d,%d)",pointRow,pointCol),2);
        double dist = 1e9;
        for(int icol=0;icol<viewOrigin[iPlane].meta().cols();icol++){
            for(int irow=0;irow<viewOrigin[iPlane].meta().rows();irow++){
                if(icol != endPix[iPlane].first && irow != endPix[iPlane].second)continue;
                if(viewOrigin[iPlane].pixel(irow,icol) == 0)continue;
                dist = sqrt(pow(pointCol-icol,2)+pow(pointRow-irow,2));
                if(dist < minDistplane[iPlane])minDistplane[iPlane]=dist;
            }
        }
    }
    minDist = minDistplane[0]+minDistplane[1]+minDistplane[2];
    return minDist;
}
//______________________________________________________
std::vector<TVector3> GetOpenSet(TVector3 newPoint, double dR){
    std::vector<TVector3> openSet;

    for(int ix = -1;ix<2;ix++){
        for(int iy = -1;iy<2;iy++){
            for(int iz = -1;iz<2;iz++){
                //if(ix == 0 && iy == 0 && iz == 0)continue;
                TVector3 neighbour = newPoint;
                neighbour.SetX(newPoint.X()+dR*ix);
                neighbour.SetY(newPoint.Y()+dR*iy);
                neighbour.SetZ(newPoint.Z()+dR*iz);
                openSet.push_back(neighbour);
            }
        }
    }
    return openSet;
}
//______________________________________________________
std::vector<std::pair<int, int> > GetEndPoint(std::vector<larcv::Image2D> viewOrigin){
    std::vector<std::pair<int, int> > endPix(3);
    for(size_t iPlane = 0;iPlane<3;iPlane++){
        for(int icol=0;icol<viewOrigin[iPlane].meta().cols();icol++){
            for(int irow=0;irow<viewOrigin[iPlane].meta().rows();irow++){
                if(viewOrigin[iPlane].pixel(irow,icol) == 1){
                    endPix[iPlane].first = icol;
                    endPix[iPlane].second = irow;
                }
            }
        }
    }
    return endPix;
}

