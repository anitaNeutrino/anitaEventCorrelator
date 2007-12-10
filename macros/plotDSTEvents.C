
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



void plotDots() {
  
  TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
    if(!canMap)
	canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
    canMap->Clear();
    canMap->SetTopMargin(0);
    canMap->SetBottomMargin(0);
    canMap->SetLeftMargin(0);
    canMap->SetRightMargin(0);
    TImage *img = TImage::Open("antarcticaIceMap.png");
    if (!img) {
	printf("Could not create an image... exit\n");
	return;
    }
    img->SetConstRatio(kFALSE);
    img->Draw("");
    
    TMarker *mark = new TMarker() ;
    mark->SetMarkerColor(6);
    mark->SetMarkerStyle(29);
    mark->SetMarkerSize(2);

}



const float arrowLength=500;

void plotEventArrow(float latitude, float longitude, float direction) 
{
    
    TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
    if(!canMap)
	canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
    canMap->Clear();
    canMap->SetTopMargin(0);
    canMap->SetBottomMargin(0);
    canMap->SetLeftMargin(0);
    canMap->SetRightMargin(0);
    TImage *img = TImage::Open("antarcticaIceMap.png");
    if (!img) {
	printf("Could not create an image... exit\n");
	return;
    }
    img->SetConstRatio(kFALSE);
    img->Draw("");
    
    TMarker *mark = new TMarker() ;
    mark->SetMarkerColor(6);
    mark->SetMarkerStyle(29);
    mark->SetMarkerSize(2);
    float xBalloon,yBalloon;
    getRelXYFromLatLong(latitude,longitude,xBalloon,yBalloon);
    cout << xBalloon << " " << yBalloon << endl;
    mark->DrawMarker(xBalloon,yBalloon);

    float yLength=arrowLength/5000.;
    float xLength=arrowLength/6000.;
    float relAngle=(direction-(180-longitude));

    float xArrow=xBalloon-xLength*TMath::Sin(relAngle*TMath::DegToRad());
    float yArrow=yBalloon-yLength*TMath::Cos(relAngle*TMath::DegToRad());

//     cout << xBalloon << " " << yBalloon << endl;
//     cout << xArrow << " " << yArrow << endl;
    
    TArrow *arrow = new TArrow();
    arrow->SetLineWidth(2);
    arrow->SetLineColor(1);
    arrow->SetFillColor(1);
    arrow->DrawArrow(xBalloon,yBalloon,xArrow,yArrow,0.02);

}

void plotEventArrowAndDot(float latitude, float longitude, float direction, float sourceLat, float sourceLon) 
{
    
    TCanvas *canMap=(TCanvas*)gROOT->FindObject("canMap");
    if(!canMap)
	canMap = new TCanvas("canMap","canMap",(int)xSize,(int)ySize);
    canMap->Clear();
    canMap->SetTopMargin(0);
    canMap->SetBottomMargin(0);
    canMap->SetLeftMargin(0);
    canMap->SetRightMargin(0);
    TImage *img = TImage::Open("antarcticaIceMap.png");
    if (!img) {
	printf("Could not create an image... exit\n");
	return;
    }
    img->SetConstRatio(kFALSE);
    img->Draw("");
    
    TMarker *mark = new TMarker() ;
    mark->SetMarkerColor(6);
    mark->SetMarkerStyle(29);
    mark->SetMarkerSize(2);
    float xBalloon,yBalloon;
    getRelXYFromLatLong(latitude,longitude,xBalloon,yBalloon);
    cout << xBalloon << " " << yBalloon << endl;
    mark->DrawMarker(xBalloon,yBalloon);

    float yLength=arrowLength/5000.;
    float xLength=arrowLength/6000.;
    float relAngle=(direction-(180-longitude));

    float xArrow=xBalloon-xLength*TMath::Sin(relAngle*TMath::DegToRad());
    float yArrow=yBalloon-yLength*TMath::Cos(relAngle*TMath::DegToRad());

//     cout << xBalloon << " " << yBalloon << endl;
//     cout << xArrow << " " << yArrow << endl;
    
    TArrow *arrow = new TArrow();
    arrow->SetLineWidth(2);
    arrow->SetLineColor(1);
    arrow->SetFillColor(1);
    arrow->DrawArrow(xBalloon,yBalloon,xArrow,yArrow,0.02);

    float xSource,ySource;
    getRelXYFromLatLong(sourceLat,sourceLon,xSource,ySource);
    cout << xSource << "\t" << ySource << "\n";
    mark->SetMarkerColor(kRed);
    mark->DrawMarker(xSource,ySource);

}





void plotEventArrowWithPhi(float latitude, float longitude, float heading,float phi) 
{

    if(phi<0 || phi>16) return;
    float direction=(2*360./16)+heading-(phi*360./16);
//    float direction=heading+(phi*360./16);
    if(direction>360) direction-=360;
    cout << direction << endl;
    return plotEventArrow(latitude,longitude,direction);

}

TFile *fp;
TTree *adu5PatTree;
TTreeIndex *adu5TreeIndex;
Float_t fHeading,fLatitude,fLongitude;
UInt_t fUnixTime;


void plotEventArrowWithTimeAndPhi(char *filename,unsigned long unixTime, float phi)
{
    static int loadedFile=0;
    if(!loadedFile) {
	fp = new TFile(filename);
	if(!fp) {
	    cerr << "Couldn't open file: " << filename << endl;
	    return;
	}
	adu5PatTree = (TTree*) fp->Get("adu5PatTree");
	if(!adu5PatTree) {
	    cerr  << "Couldn't get adu5PatTree" << endl;
	    return;
	}

	adu5PatTree->SetBranchAddress("heading",&fHeading);
	adu5PatTree->SetBranchAddress("latitude",&fLatitude);
	adu5PatTree->SetBranchAddress("longitude",&fLongitude);
	adu5PatTree->SetBranchAddress("unixTime",&fUnixTime);

	adu5PatTree->SetBranchStatus("*",0);
	adu5PatTree->SetBranchStatus("unixTime",1);
	adu5PatTree->SetBranchStatus("heading",1);
	adu5PatTree->SetBranchStatus("latitude",1);
	adu5PatTree->SetBranchStatus("longitude",1);

	adu5PatTree->BuildIndex("unixTime");
	adu5TreeIndex=(TTreeIndex*)adu5PatTree->GetTreeIndex();
	loadedFile=1;
    }
    Long64_t entry=adu5TreeIndex->GetEntryNumberWithBestIndex(unixTime,0);
//    cout << entry << endl;
    adu5PatTree->GetEntry(entry);
    cout << fUnixTime << " " << fLatitude << " " << fLongitude << " " << fHeading << endl;
    if(TMath::Abs(float(unixTime)-float(fUnixTime))>100) {
	cout << "Warning your selected time and the match time are not close" << endl;
    }
    return plotEventArrowWithPhi(fLatitude, fLongitude, fHeading,phi);
}


















