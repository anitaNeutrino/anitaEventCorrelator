
const float TrueScaleLat=71;
const float CentralMeridian=0;
const float RadiusOfEarth=6378.1e3; //Metres
const float xOffest=375;
const float yOffset=312.5;
const float scale=271.5/2.19496e+06;
const float xSize=750;
const float ySize=625;


void getRelXYFromLatLong(float latitude, float longitude,
			 float &x, float &y)
{
    //Negative longitude is west
 //    //All latitudes assumed south
    float absLat=TMath::Abs(latitude);
    float r=RadiusOfEarth*TMath::Cos((90.-TrueScaleLat)*TMath::DegToRad())*TMath::Tan((90-absLat)*TMath::DegToRad());
    y=r*TMath::Cos(longitude*TMath::DegToRad());
    x=r*TMath::Sin(longitude*TMath::DegToRad());   

    y*=scale;
    y+=yOffset;
    y/=ySize;
    x*=scale;
    x+=xOffest;
    x/=xSize;
 
}


void quickDirtyInitialMapPlot() {
   
   TChain *sillyTree = new TChain("sillyTree");
   //   sillyTree->Add("noCutFirstDay/sillyDST*");
   sillyTree->Add("noCutInitial/sillyDST*");
   
   TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
   if(!canMap)
      canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
   canMap->Clear();
   canMap->SetLogz(1);
   canMap->SetTopMargin(0);
   canMap->SetBottomMargin(0);
   canMap->SetLeftMargin(0);
   canMap->SetRightMargin(0);
   TImage *img = TImage::Open("antarcticaIceMapBW.png");
   if (!img) {
      printf("Could not create an image... exit\n");
      return;
    }
   img->SetConstRatio(kFALSE);
   img->Draw("");
   Int_t numEvents=sillyTree->Draw("sourceY:sourceX","chiSq<20","same");
   cout << "Total:\t" << numEvents << endl;
   
   TCanvas *can = new TCanvas("can","can");
   sillyTree->Draw("chiSq","chiSq<100");

}
