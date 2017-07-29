#include "DataFormat/track.h"
#include "DataFormat/hit.h"
#include "DataFormat/Image2D.h"
#include "DataFormat/mctrack.h"
#include "DataFormat/wire.h"
#include "LArUtil/Geometry.h"
#include "LArUtil/GeometryHelper.h"
#include "LArUtil/TimeService.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TF1.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/core/DataFormat/ChStatus.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/DeepLearning/myLArCV/app/LArOpenCVHandle/LArbysUtils.h"

#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgo.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/ThruMu/AStar3DAlgoProton.h"


#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/myLArLiteCV/app/SCE/SpaceChargeMicroBooNE.h"
#include "/Users/hourlier/Documents/PostDocMIT/Research/MicroBooNE/mylarlite/UserDev/Adrien/Reco3DTracks/AStarTracker.h"

int main(){
    TFile *f = new TFile("../data/c_ana_0200_0299.root");
    TTree *T = (TTree*)f->Get("dQdSAnalysis");
    TTree *Tv =(TTree*)f->Get("VertexTree");

    UInt_t run, subrun, event, entry;
    Int_t roid, vtxid;
    Int_t run_v, subrun_v, event_v, entry_v;
    Int_t roid_v, vtxid_v;
    double x,y,z;

    double particle0_end_x;
    double particle0_end_y;
    double particle0_end_z;

    double particle1_end_x;
    double particle1_end_y;
    double particle1_end_z;

    std::vector< std::pair<int,int> > *particle0_end_point_v = 0;
    std::vector< std::pair<int,int> > *particle1_end_point_v = 0;
    std::vector<TVector3> *image_particle0_plane0_tmp = 0;
    std::vector<TVector3> *image_particle0_plane1_tmp = 0;
    std::vector<TVector3> *image_particle0_plane2_tmp = 0;

    std::vector<TVector3> *image_particle1_plane0_tmp = 0;
    std::vector<TVector3> *image_particle1_plane1_tmp = 0;
    std::vector<TVector3> *image_particle1_plane2_tmp = 0;

    std::vector<larcv::Image2D> *adc_img_v = 0;

    T->SetBranchAddress("run",&run);
    T->SetBranchAddress("subrun",&subrun);
    T->SetBranchAddress("event",&event);
    T->SetBranchAddress("entry",&entry);
    T->SetBranchAddress("roid",&roid);
    T->SetBranchAddress("vtxid",&vtxid);
    T->SetBranchAddress("x",&x);
    T->SetBranchAddress("y",&y);
    T->SetBranchAddress("z",&z);
    T->SetBranchAddress("particle0_end_point_v",&particle0_end_point_v);
    T->SetBranchAddress("particle1_end_point_v",&particle1_end_point_v);

    T->SetBranchAddress("particle0_end_x",&particle0_end_x);
    T->SetBranchAddress("particle0_end_y",&particle0_end_y);
    T->SetBranchAddress("particle0_end_z",&particle0_end_z);
    T->SetBranchAddress("particle1_end_x",&particle1_end_x);
    T->SetBranchAddress("particle1_end_y",&particle1_end_y);
    T->SetBranchAddress("particle1_end_z",&particle1_end_z);

    T->SetBranchAddress("image_particle0_plane0_tmp", &image_particle0_plane0_tmp);
    T->SetBranchAddress("image_particle0_plane1_tmp", &image_particle0_plane1_tmp);
    T->SetBranchAddress("image_particle0_plane2_tmp", &image_particle0_plane2_tmp);
    T->SetBranchAddress("image_particle1_plane0_tmp", &image_particle1_plane0_tmp);
    T->SetBranchAddress("image_particle1_plane1_tmp", &image_particle1_plane1_tmp);
    T->SetBranchAddress("image_particle1_plane2_tmp", &image_particle1_plane2_tmp);

    Tv->SetBranchAddress("run",&run_v);
    Tv->SetBranchAddress("subrun",&subrun_v);
    Tv->SetBranchAddress("event",&event_v);
    Tv->SetBranchAddress("entry",&entry_v);
    Tv->SetBranchAddress("roid",&roid_v);
    Tv->SetBranchAddress("vtxid",&vtxid_v);
    Tv->SetBranchAddress("adc_img_v", &adc_img_v);

    int nentries = T->GetEntries();

    TCanvas *cImage = new TCanvas("cImage","cImage",900,300);
    cImage->Divide(3,1);
    TH2D *hImage[3];
    TGraph *gVertex[3];
    TGraph *gEndPoint[3];
    TGraph *gROI[3];

    std::vector<std::vector<larcv::Image2D> > DataImages;
    std::vector<larcv::Image2D> SingleParticleImages(3);
    larlite::AStarTracker tracker;
    tracker.initialize();
    tracker.SetCompressionFactors(1,6);
    tracker.SetVerbose(0);

    TVector3 vertex, endPoint[2];
    for(int i = 0;i<nentries;i++){
        T->GetEntry(i);
        std::cout << "Entry " << i << std::endl;

        //
        // GetEstimated 3D Start and End Points
        vertex.SetXYZ(x,y,z);
        endPoint[0].SetXYZ(particle0_end_x,particle0_end_y,particle0_end_z);
        endPoint[1].SetXYZ(particle1_end_x,particle1_end_y,particle1_end_z);

        std::vector<TVector3> EndPoints;
        EndPoints.push_back(vertex);
        EndPoints.push_back(endPoint[0]);
        EndPoints.push_back(endPoint[1]);

        for(int iPoints = 0; iPoints < EndPoints.size(); iPoints++){
            bool thisPointOK = tracker.CheckEndPointsInVolume(EndPoints[iPoints]);
            if(!thisPointOK) std::cout << "point " << iPoints << " out of range" << std::endl;
        }

        tracker.SetTimeAndWireBounds(EndPoints);
        std::vector<std::pair<double,double> > time_bounds = tracker.GetTimeBounds();
        std::vector<std::pair<double,double> > wire_bounds = tracker.GetWireBounds();

        Tv->GetEntry(i);
        for(int iPlane = 0;iPlane<3;iPlane++){

            std::cout << "(*adc_img_v)[" << iPlane << "].meta().rows()         = " << (*adc_img_v)[iPlane].meta().rows() << std::endl;
            std::cout << "(*adc_img_v)[" << iPlane << "].meta().height()       = " << (*adc_img_v)[iPlane].meta().height() << std::endl;
            std::cout << "(*adc_img_v)[" << iPlane << "].meta().pixel_height() = " << (*adc_img_v)[iPlane].meta().pixel_height() << std::endl;
            std::cout << std::endl;
            std::cout << "(*adc_img_v)[" << iPlane << "].meta().cols()         = " << (*adc_img_v)[iPlane].meta().cols() << std::endl;
            std::cout << "(*adc_img_v)[" << iPlane << "].meta().width()        = " << (*adc_img_v)[iPlane].meta().width() << std::endl;
            std::cout << "(*adc_img_v)[" << iPlane << "].meta().pixel_width()  = " << (*adc_img_v)[iPlane].meta().pixel_width() << std::endl;
            std::cout << std::endl;
            std::cout << "time bounds : " << time_bounds[iPlane].first/((*adc_img_v)[iPlane].meta().pixel_height()) << " => " << time_bounds[iPlane].second/((*adc_img_v)[iPlane].meta().pixel_height()) << std::endl;
            std::cout << "wire bounds : " << wire_bounds[iPlane].first/((*adc_img_v)[iPlane].meta().pixel_width())  << " => " << wire_bounds[iPlane].second/((*adc_img_v)[iPlane].meta().pixel_width())  << std::endl;
            std::cout << std::endl << std::endl;

            double image_cols = wire_bounds[iPlane].second - wire_bounds[iPlane].first;
            double image_rows = time_bounds[iPlane].second - time_bounds[iPlane].first;

            double image_width  = wire_bounds[iPlane].second - wire_bounds[iPlane].first;
            double image_height = time_bounds[iPlane].second - time_bounds[iPlane].first;
            size_t col_count    = (size_t)(image_width / (*adc_img_v)[iPlane].meta().pixel_width()  );
            size_t row_count    = (size_t)(image_height/ (*adc_img_v)[iPlane].meta().pixel_height() );
            double origin_x     = wire_bounds[iPlane].first;
            double origin_y     = time_bounds[iPlane].second;

            std::cout << "newMeta definition:" << std::endl;
            std::cout << "width      = " << image_width << std::endl;
            std::cout << "height     = " << image_height << std::endl;

            std::cout << "row_count  = " << row_count << std::endl;
            std::cout << "col_count  = " << col_count << std::endl;

            std::cout << "origin_x   = " << origin_x << std::endl;
            std::cout << "origin_y   = " << origin_y << std::endl;
            std::cout << "plane      = " << iPlane << std::endl << std::endl;

            larcv::ImageMeta newMeta(image_width, image_height,
                                     row_count, col_count,
                                     origin_x,  // origin x = min wire
                                     origin_y, // origin y = max time
                                     iPlane);
            //larcv::Image2D newImage = (*adc_img_v)[iPlane].crop(newMeta);
            larcv::Image2D newImage = (*adc_img_v)[iPlane];

            hImage[iPlane] = new TH2D(Form("hImage_%d_%d",i,iPlane),Form("hImage_%d_%d",i,iPlane),newImage.meta().cols(),0,newImage.meta().width(),newImage.meta().rows(),0,newImage.meta().height());

            for(int icol = 0;icol<newImage.meta().cols();icol++){
                for(int irow = 0;irow<newImage.meta().rows();irow++){
                    if(newImage.pixel(irow, icol) < 10) continue;
                    hImage[iPlane]->SetBinContent(icol+1,irow+1,newImage.pixel(irow, icol));
                }
            }
            cImage->cd(iPlane+1);
            hImage[iPlane]->Draw("colz");
            /*gROI[iPlane] = new TGraph();
            gROI[iPlane]->SetPoint(0,newMeta.tl().x,newMeta.tl().y);
            gROI[iPlane]->SetPoint(1,newMeta.tr().x,newMeta.tr().y);
            gROI[iPlane]->SetPoint(2,newMeta.br().x,newMeta.br().y);
            gROI[iPlane]->SetPoint(3,newMeta.bl().x,newMeta.bl().y);
            gROI[iPlane]->SetPoint(4,newMeta.tl().x,newMeta.tl().y);
            gROI[iPlane]->SetLineWidth(1);
            gROI[iPlane]->Draw("same LP");*/

            gVertex[iPlane] = new TGraph();
            gVertex[iPlane]->SetMarkerStyle(20);
            gVertex[iPlane]->SetMarkerSize(0.25);
            double _parent_x = vertex.X();
            double _parent_y = vertex.Y();
            double _parent_z = vertex.Z();
            double _parent_t = 0;
            double x_pixel, y_pixel;
            Project3D(newImage.meta(),_parent_x,_parent_y,_parent_z,_parent_t,iPlane,x_pixel,y_pixel);

            gVertex[iPlane]->SetPoint(0,x_pixel,y_pixel);
            std::cout << iPlane << "\t vertex projection : " << tracker.GetWireTimeProjection(vertex)[iPlane].first << ", " << tracker.GetWireTimeProjection(vertex)[iPlane].second << std::endl;
            gVertex[iPlane]->Draw("same P");

            /*gEndPoint[iPlane] = new TGraph();
            gEndPoint[iPlane]->SetMarkerStyle(20);
            gEndPoint[iPlane]->SetMarkerSize(0.25);
            gEndPoint[iPlane]->SetPoint(0,particle0_end_point_v->at(iPlane).second,particle0_end_point_v->at(iPlane).first);
            gEndPoint[iPlane]->SetPoint(1,particle1_end_point_v->at(iPlane).second,particle1_end_point_v->at(iPlane).first);
             */

            std::cout << "new Meta bounds : " << std::endl;
            std::cout <<"("<< newMeta.tl().x << ","<< newMeta.tl().y << ")" << std::endl;
            std::cout <<"("<< newMeta.br().x << ","<< newMeta.br().y << ")" << std::endl;
            std::cout << std::endl << std::endl;

        }
        cImage->Modified();
        cImage->Update();
        cImage->SaveAs(Form("cImage_%d.png",i));

        continue;

        //tracker.SetImages((*adc_img_v));



        std::cout << run << "\t" << subrun << "\t" << event << "\t vertex (x,y,z) :" << x << " , " << y << " , " << z << "\t" << " end point particle 0 : " << std::endl;
        for(int iPlane = 0;iPlane<3;iPlane++){
            std::cout << "plane " << iPlane << " , col :" << particle0_end_point_v->at(iPlane).second << " , row : " << particle0_end_point_v->at(iPlane).first << std::endl;
        }
        std::cout << std::endl;
        std::cout << run << "\t" << subrun << "\t" << event << "\t vertex (x,y,z) : " << x << " , " << y << " , " << z << "\t" << " end point particle 1 : " << std::endl;
        for(int iPlane = 0;iPlane<3;iPlane++){
            std::cout << "plane " << iPlane << " , col :" << particle1_end_point_v->at(iPlane).second << " , row : " << particle1_end_point_v->at(iPlane).first << std::endl;
        }

        std::cout << std::endl;

        for(int iPlane = 0; iPlane<3;iPlane++){
            hImage[iPlane] = new TH2D(Form("hImage_%d_%d",i,iPlane),Form("hImage_%d_%d",i,iPlane),1000,0,1000,1000,0,1000);
            gEndPoint[iPlane] = new TGraph();
            gEndPoint[iPlane]->SetMarkerStyle(20);
            gEndPoint[iPlane]->SetMarkerSize(0.25);
            gEndPoint[iPlane]->SetPoint(0,particle0_end_point_v->at(iPlane).second,particle0_end_point_v->at(iPlane).first);
            gEndPoint[iPlane]->SetPoint(1,particle1_end_point_v->at(iPlane).second,particle1_end_point_v->at(iPlane).first);
            gVertex[iPlane] = new TGraph();
            gVertex[iPlane]->SetMarkerStyle(20);
            gVertex[iPlane]->SetMarkerSize(0.25);
            gVertex[iPlane]->SetPoint(0,tracker.GetWireTimeProjection(vertex)[iPlane].first/tracker.GetCompressionFactorWire(),tracker.GetWireTimeProjection(vertex)[iPlane].second/tracker.GetCompressionFactorTime());
            std::cout << iPlane << "\t" << tracker.GetWireTimeProjection(vertex)[iPlane].first/tracker.GetCompressionFactorWire() << ", " << tracker.GetWireTimeProjection(vertex)[iPlane].second/tracker.GetCompressionFactorTime() << std::endl;
            gROI[iPlane] = new TGraph();
        }

        //_______________________________________________________________
        //_______________________________________________________________
        // Process Particle 1
        std::cout << "particle 1" << std::endl;


        tracker.SetEndPoints(EndPoints[0],EndPoints[1]);


        std::cout << tracker.GetTimeBounds().size() << " planes?" << std::endl;

        for(int iPlane = 0; iPlane<3;iPlane++){
            std::cout << iPlane << std::endl;

            // Retrieve boundaries
            auto const& time_bound = time_bounds[iPlane];
            auto const& wire_bound = wire_bounds[iPlane];

            gROI[iPlane]->SetPoint(0, wire_bound.first,time_bound.first/tracker.GetCompressionFactorTime());
            gROI[iPlane]->SetPoint(1, wire_bound.first,time_bound.second/tracker.GetCompressionFactorTime());
            gROI[iPlane]->SetPoint(2, wire_bound.second,time_bound.second/tracker.GetCompressionFactorTime());
            gROI[iPlane]->SetPoint(3, wire_bound.second,time_bound.first/tracker.GetCompressionFactorTime());
            gROI[iPlane]->SetPoint(4, wire_bound.first,time_bound.first/tracker.GetCompressionFactorTime());

            tracker.tellMe("time_bound and wire_bound OK",0);

            // If no hit on this plane, make an empty image
            if(wire_bound.second <= wire_bound.first ||
               time_bound.second <= time_bound.first ) {
                tracker.tellMe("Adding Empty Image",0);
                SingleParticleImages[iPlane] = larcv::Image2D();
                continue;
            }
            // Construct meta
            tracker.tellMe("construct meta",0);
            size_t image_cols = (size_t)(wire_bound.second - wire_bound.first)/*/tracker.GetCompressionFactorWire()*/;
            size_t image_rows = (size_t)(time_bound.second - time_bound.first)/*/tracker.GetCompressionFactorTime()*/;
            std::cout << "wire bounds : " << wire_bound.first << " => " << wire_bound.second << " => " << image_cols << " columns" <<  std::endl;
            std::cout << "time bounds : " << time_bound.first << " => " << time_bound.second << " => " << image_rows << " rows" <<  std::endl;

            /*tracker.tellMe(Form("%zu cols and %zu rows",image_cols,image_rows),1);
            larcv::ImageMeta hit_meta((double)image_cols, (double)image_rows,
                                      image_rows, image_cols,
                                      (size_t) wire_bound.first,  // origin x = min wire
                                      (size_t) time_bound.second, // origin y = max time
                                      iPlane);*/


            // Prepare hit image data + fill
            
            //std::vector<float> image_data(hit_meta.rows() * hit_meta.cols(), 0.);

            std::vector<TVector3> *imageparticle = 0;

            if(iPlane == 0){imageparticle = image_particle0_plane0_tmp;}
            if(iPlane == 1){imageparticle = image_particle0_plane1_tmp;}
            if(iPlane == 2){imageparticle = image_particle0_plane2_tmp;}

            std::cout << "imageparticle->size() : " << imageparticle->size() << std::endl;
            for(int iPx = 0;iPx < imageparticle->size();iPx++){

                int col = imageparticle->at(iPx).Y();
                int row = imageparticle->at(iPx).X();
                //int tick = row*tracker.GetCompressionFactorTime();
                //int wire = col*tracker.GetCompressionFactorWire();
                int intensity = imageparticle->at(iPx).Z();

                //std::cout << col << "\t" << row << "\t" << intensity << std::endl;
                hImage[iPlane]->SetBinContent(col,row,intensity);
                //image_data[hit_meta.rows() * col + row] = intensity;

                //larcv::Image2D hit_image(std::move(hit_meta),std::move(image_data));
                //SingleParticleImages[iPlane] = hit_image;
            }
            cImage->cd(iPlane+1);
            hImage[iPlane]->Draw("colz");
        }

        //_______________________________________________________________
        //_______________________________________________________________
        // Process Particle2
        std::cout << "particle 2" << std::endl;

        tracker.SetTimeAndWireBounds(vertex, endPoint[1]);
        std::vector<std::pair<double,double> > time_bounds_particle1 = tracker.GetTimeBounds();
        std::vector<std::pair<double,double> > wire_bounds_particle1 = tracker.GetWireBounds();

        for(int iPlane = 0; iPlane<3;iPlane++){
            std::cout << iPlane << std::endl;

            std::vector<TVector3> *imageparticle = 0;

            if(iPlane == 0){imageparticle = image_particle1_plane0_tmp;}
            if(iPlane == 1){imageparticle = image_particle1_plane1_tmp;}
            if(iPlane == 2){imageparticle = image_particle1_plane2_tmp;}
            std::cout << "imageparticle->size() : " << imageparticle->size() << std::endl;
            for(int iPx = 0;iPx < imageparticle->size();iPx++){

                int col = imageparticle->at(iPx).Y();
                int row = imageparticle->at(iPx).X();
                //int tick = row*tracker.GetCompressionFactorTime();
                //int wire = col*tracker.GetCompressionFactorWire();
                int intensity = imageparticle->at(iPx).Z();

                //std::cout << iPx << "\t" << col << "\t" << row << "\t" << intensity << std::endl;
                hImage[iPlane]->SetBinContent(col,row,intensity);

            }
            cImage->cd(iPlane+1);
            hImage[iPlane]->Draw("colz");
            gEndPoint[iPlane]->Draw("same P");
            gROI[iPlane]->Draw("same PL");
            gVertex[iPlane]->Draw("same P");
        }

        cImage->Modified();
        cImage->Update();
        cImage->SaveAs(Form("cImage_%d.pdf",i));
    }
}
