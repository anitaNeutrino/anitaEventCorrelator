
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


void quickDirtyMapPlot() {
   
   TChain *sillyTree = new TChain("sillyTree");
   sillyTree->Add("newClockNoCutFirstDay/sillyDST*");
   //   sillyTree->Add("initial/sillyDST*");
   
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
   Int_t numEvents=sillyTree->Draw("sourceY:sourceX>>hMap(1000,,,1000,,)","chiSq<10","a col same");
   cout << "Total:\t" << numEvents << endl;
   
   TCanvas *canLat = new TCanvas("canLat","canLat");
   sillyTree->Draw("sourceLat:sourceLon>>hLat(1000,,,1000,,)","chiSq<10","col");

   TCanvas *canMcMurdo = new TCanvas("canMcMurdo","canMcMurdo",(int)xSize,(int)ySize);
   canMcMurdo->Clear();
   canMcMurdo->SetLogz(1);
   canMcMurdo->SetTopMargin(0);
   canMcMurdo->SetBottomMargin(0);
   canMcMurdo->SetLeftMargin(0);
   canMcMurdo->SetRightMargin(0);
   img->SetConstRatio(kFALSE);
   img->Draw("");
   Int_t numEventsMcMurdo=sillyTree->Draw("sourceY:sourceX>>hMap2(1000,,,1000,,)","chiSq<10   && sourceLon>164 && sourceLon<170","a col same");
   sillyTree->Draw("balloonY:balloonX","chiSq<10   && sourceLon>164 && sourceLon<170","same");
   cout << "McMurdo:\t" << numEventsMcMurdo << endl;


   TCanvas *canTaylor = new TCanvas("canTaylor","canTaylor",(int)xSize,(int)ySize);
   canTaylor->Clear();
   canTaylor->SetLogz(1);
   canTaylor->SetTopMargin(0);
   canTaylor->SetBottomMargin(0);
   canTaylor->SetLeftMargin(0);
   canTaylor->SetRightMargin(0);
   img->SetConstRatio(kFALSE);
   img->Draw("");
   Int_t numEventsTaylor=sillyTree->Draw("sourceY:sourceX>>hMap3(1000,,,1000,,)","chiSq<10   && sourceLon>156 && sourceLon<161","a col same");
   sillyTree->Draw("balloonY:balloonX","chiSq<10   && sourceLon>156 && sourceLon<161","same");
   cout << "Taylor:\t" << numEventsTaylor << endl;
   
   TCanvas *canTest= new TCanvas("canTest","canTest");
   sillyTree->Draw("triggerTime>>htrig","chiSq<10   && sourceLon>156 && sourceLon<161","");


   //   TCanvas *can = new TCanvas("can","can");
   //   sillyTree->Draw("log10(chiSq)","");
}
