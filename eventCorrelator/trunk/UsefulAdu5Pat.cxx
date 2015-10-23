//////////////////////////////////////////////////////////////////////////////
/////  UsefulAdu5Pat.cxx        ANITA event reading class                  /////
/////                                                                    /////
/////  Description:                                                      /////
/////    A simple class for utilising Adu5Pat objects                      /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "BedmapReader.h"
#include "RampdemReader.h"
#include "TProfile2D.h"
#include <iostream>
#include <fstream>

#include "TVector3.h"

//#define CMINCH 2.54

ClassImp(UsefulAdu5Pat);

AnitaGeomTool *fUPGeomTool=0;
//RampdemReader *fRampdemReader=0;
Int_t fRampdemOpen=0;


UsefulAdu5Pat::UsefulAdu5Pat() 
  : Adu5Pat()
{
  //Default Constructor
  fIncludeGroupDelay=0;
  fThetaWave=0;
  fPhiWave=0;
  fSourceLongitude=-1;
  fSourceLatitude=-1;
  fBalloonCoords[0]=0;
  fBalloonCoords[1]=0;
  fBalloonCoords[2]=0;
  fBalloonHeight=0;
  if(!fUPGeomTool){
    fUPGeomTool=AnitaGeomTool::Instance();
  }
  if(!fRampdemOpen){
    fRampdemReader=RampdemReader::Instance();
    fRampdemOpen=1;
  }

}

UsefulAdu5Pat::UsefulAdu5Pat(Adu5Pat *patPtr,double deltaR,double deltaRL,double deltaUD)
  : Adu5Pat(*patPtr)
{
  fIncludeGroupDelay=0;
  pitch=AnitaStaticAdu5Offsets::pitch;
  roll=AnitaStaticAdu5Offsets::roll;
  // pitch+=AnitaStaticAdu5Offsets::pitch;
  // roll+=AnitaStaticAdu5Offsets::roll;  
  heading+=AnitaStaticAdu5Offsets::heading;
  if(heading>=360) heading-=360;
  if(heading<0) heading+=360;
  fThetaWave=0;
  fPhiWave=0;
  fSourceLongitude=-1;
  fSourceLatitude=-1;
  if(!fUPGeomTool)
    fUPGeomTool=AnitaGeomTool::Instance();
  fUPGeomTool->updateAnt(deltaR,deltaRL,deltaUD);
  fUPGeomTool->getCartesianCoords(latitude,longitude,altitude,
				  fBalloonCoords);
  fBalloonPos.SetXYZ(fBalloonCoords[0],fBalloonCoords[1],fBalloonCoords[2]);
  fBalloonTheta=fBalloonPos.Theta();
  fBalloonPhi=fBalloonPos.Phi();
  if(fBalloonPhi<0) fBalloonPhi+=TMath::TwoPi();
  fBalloonHeight=fBalloonPos.Mag();
  if(!fRampdemOpen){
    fRampdemReader=RampdemReader::Instance();
    fRampdemOpen=1;
  }
}

UsefulAdu5Pat::UsefulAdu5Pat(Adu5Pat *patPtr)
  : Adu5Pat(*patPtr)
{
  fIncludeGroupDelay=0;

  pitch=AnitaStaticAdu5Offsets::pitch;
  roll=AnitaStaticAdu5Offsets::roll;
  // pitch+=AnitaStaticAdu5Offsets::pitch;
  // roll+=AnitaStaticAdu5Offsets::roll;
  heading+=AnitaStaticAdu5Offsets::heading;

  if(heading>=360) heading-=360;
  if(heading<0) heading+=360;
  fThetaWave=0;
  fPhiWave=0;
  fSourceLongitude=-1;
  fSourceLatitude=-1;
  if(!fUPGeomTool){
    fUPGeomTool=AnitaGeomTool::Instance();
    //     std::cout << "making the geom" << std::endl;
  }
  //std::cout << "LatLonAlt: " << latitude << "\t" << longitude << "\t" << altitude << "\n";
  fUPGeomTool->getCartesianCoords(latitude,longitude,altitude,
				  fBalloonCoords);
  // std::cout << "Balloon Coords: " << fBalloonCoords[0] << "\t" << fBalloonCoords[1] << "\t" << fBalloonCoords[2] << "\n";
  fBalloonPos.SetXYZ(fBalloonCoords[0],fBalloonCoords[1],fBalloonCoords[2]);
  fBalloonTheta=fBalloonPos.Theta();
  fBalloonPhi=fBalloonPos.Phi();
  if(fBalloonPhi<0) fBalloonPhi+=TMath::TwoPi();
  fBalloonHeight=fBalloonPos.Mag();
  if(!fRampdemOpen){
    fRampdemReader=RampdemReader::Instance();
    fRampdemOpen=1;
  }
}


UsefulAdu5Pat::~UsefulAdu5Pat() 
{
  //Default Destructor
}


void UsefulAdu5Pat::getThetaAndPhiWaveWillySeavey(Double_t &thetaWave, Double_t &phiWave)
{   
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_SURF_SEAVEY,AnitaLocations::LATITUDE_SURF_SEAVEY,AnitaLocations::ALTITUDE_SURF_SEAVEY,thetaWave,phiWave);
}

void UsefulAdu5Pat::getThetaAndPhiWaveWillyBorehole(Double_t &thetaWave,Double_t &phiWave)
{
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_BH,AnitaLocations::LATITUDE_BH,AnitaLocations::ALTITUDE_BH,thetaWave,phiWave);
}


void UsefulAdu5Pat::getThetaAndPhiWaveTaylorDome(Double_t &thetaWave, Double_t &phiWave)
{   
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_TD,AnitaLocations::LATITUDE_TD,AnitaLocations::ALTITUDE_TD,thetaWave,phiWave);
}


int UsefulAdu5Pat::getSourceLonAndLatAltZero(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat) {

  // std::cout << "getSourceLonAndLatAltZero " << phiWave << "\t" << thetaWave << "\n";
  return getSourceLonAndLatAtDesiredAlt(phiWave, thetaWave, sourceLon, sourceLat, 0.0);

}

int UsefulAdu5Pat::getSourceLonAndLatAtDesiredAlt(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat, Double_t desiredAlt = 0.0)
{

  //std::cout << "getSourceLonAndLatAtDesiredAlt " << phiWave << "\t" << thetaWave << "\n";

  if(fPhiWave!=phiWave) fPhiWave=phiWave;
  if(fThetaWave!=thetaWave) fThetaWave=thetaWave;

  //   Double_t thetaBalloon=fUPGeomTool->getThetaFromLat(TMath::Abs(latitude));
  //   Double_t phiBalloon=fUPGeomTool->getPhiFromLon(longitude);
  //   Double_t balloonHeight=fUPGeomTool->getGeoid(thetaBalloon)+altitude;
   
  Double_t tempPhiWave=phiWave;
  Double_t tempThetaWave=TMath::PiOver2()-thetaWave;
  //Now need to take account of balloon heading

  TVector3 rollAxis,pitchAxis;
  rollAxis=fUPGeomTool->fRollRotationAxis;
  pitchAxis=fUPGeomTool->fPitchRotationAxis;

  //std::cout << "Attitude: " << heading << "\t" << pitch << "\t" << roll 
  //	       << "\n";

  if(heading>=0 && heading<=360) {
    TVector3 arbDir;
    arbDir.SetMagThetaPhi(1,tempThetaWave,-1*tempPhiWave);

    //std::cout << "tempPhiWave3: " << arbDir.Phi() << "\t" << tempPhiWave << "\n";

    arbDir.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    rollAxis.Rotate((heading)*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    pitchAxis.Rotate((heading)*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    //     std::cout << arbDir.Phi() << "\t" << arbDir.Theta() << "\n";
    //     std::cout << "Pitch Axis\t" << pitchAxis.X() << " " 
    //	       << pitchAxis.Y() << " " << pitchAxis.Z() << "\t" << pitch << "\n";


    arbDir.Rotate(pitch*TMath::DegToRad(),pitchAxis);
    //     std::cout << arbDir.Phi() << "\t" << arbDir.Theta() << "\n";
    rollAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);
    arbDir.Rotate(roll*TMath::DegToRad(),rollAxis);//roll and pitch

    tempPhiWave=arbDir.Phi();
    //std::cout << "tempPhiWave2: " << tempPhiWave << "\n";
    if(tempPhiWave>TMath::TwoPi()) {
      tempPhiWave-=TMath::TwoPi();
    }
    if(tempPhiWave<0) {
      tempPhiWave+=TMath::TwoPi();
    }
    tempThetaWave=arbDir.Theta();
  }
  else std::cout << "heading bad" << std::endl;
   
  //   std::cout << "Get source angles: " <<  tempThetaWave << "\t" << tempPhiWave << "\n";

  //Double_t re=balloonHeight-altitude+2000;
  //   Double_t re=fBalloonHeight-altitude;
   
  Double_t reBalloon=fUPGeomTool->getDistanceToCentreOfEarth(latitude)+desiredAlt; // ACG mod
  //   std::cout << "Radius difference: " << re-re2 << "\t" << re << "\t" << re2 << "\n";
  Double_t re=reBalloon;
  Double_t nextRe=re;
  Double_t reh=reBalloon+altitude;
  do {
    //Okay so effectively what we do here is switch to cartesian coords with the balloon at RE_balloon + altitude and then try to fins where the line at thetaWave cuts the Earth's surface.
    //This is iterative because the Earth is cruel and isn't flat, and I couldn't be bothered to work out how to do it more elegantly.
    re=nextRe;

    Double_t sintw=TMath::Sin(tempThetaWave);
    Double_t costw=TMath::Cos(tempThetaWave);
    //     Double_t sqrtArg=(reh*reh*costw*costw - (reh*reh-re*re));
    Double_t sqrtArg(re*re-reh*reh*sintw*sintw);
    if(sqrtArg<0) {
      // No solution possible
      //     std::cout << "No solution possible\n";
      return 0;
    }
    Double_t L=reh*costw - TMath::Sqrt(sqrtArg);
    Double_t sinThetaL=L*sintw/re;
    Double_t sourceTheta=TMath::ASin(sinThetaL);
     
    //std::cout << "Source Theta: " << sourceTheta*TMath::RadToDeg() << "\t" << tempThetaWave << "\t" << reh << "\t" << re << "\t" << fUPGeomTool->getDistanceToCentreOfEarth(latitude) << "\n";
     
    //One possible improvement is to reiterate using the radius at the source location.
     
     
    //   Double_t sinSquiggle=(fBalloonHeight/tempRE)*TMath::Sin(tempThetaWave);
    //   Double_t squiggle=TMath::ASin(sinSquiggle);
    //   if(squiggle<0) squiggle+=TMath::Pi();
     
    //   Double_t sourceTheta=TMath::Pi()-tempThetaWave-squiggle;
    //   std::cout << thetaWave*TMath::RadToDeg() << "\t" << L << "\t" << sourceTheta*TMath::RadToDeg() << "\t" << phiWave*TMath::RadToDeg() << std::endl;
    //Start at ground below balloon
    fSourcePos.SetX(0);
    fSourcePos.SetY(0);
    fSourcePos.SetZ(re);
     
    //std::cout << "00ReLoc: " << fSourcePos.X() << "\t" << fSourcePos.Y() << "\t" << fSourcePos.Z() << "\n";
     
    //Rotate to latitude relative to balloon
    fSourcePos.RotateY(sourceTheta);   
    //Rotate to longitude relative to balloon
    fSourcePos.RotateZ(-1*tempPhiWave);
     
    //std::cout << "RelBalloonLoc: " << fSourcePos.X() << "\t" << fSourcePos.Y() << "\t" << fSourcePos.Z() << "\n";
     
    //Rotate to correct absolute values
    //std::cout << "Balloon angles: " << fBalloonTheta << "\t" << fBalloonPhi << "\n";
    fSourcePos.RotateY(fBalloonTheta);
    fSourcePos.RotateZ(fBalloonPhi);
    //Goofy sign thing
    //   fSourcePos.SetZ(-1*fSourcePos.Z());
     
    Double_t sourceVec[3];
    fSourcePos.GetXYZ(sourceVec);
    //std::cout << "BalloonPos: " << fBalloonPos.X() << "\t" << fBalloonPos.Y() << "\t" << fBalloonPos.Z() << "\n";
    //std::cout << "ZeroAltLoc: " << sourceVec[0] << "\t" << sourceVec[1] << "\t" << sourceVec[2] << "\n";
    Double_t sourceAlt;
    fUPGeomTool->getLatLonAltFromCartesian(sourceVec,sourceLat,sourceLon,sourceAlt);
    //std::cout << "SourceLatLonAlt: " << sourceLat << "\t" << sourceLon << "\t" 
    //	       << sourceAlt << "\n";
    nextRe=fUPGeomTool->getDistanceToCentreOfEarth(sourceLat);
    //std::cout << "Earth radius: " << nextRe << "\t" << re << "\n";
    //     break;
  } while(TMath::Abs(nextRe-re)>1);
  //   fUPGeomTool->getLonLat(fSourcePos,sourceLon,sourceLat);
  //   sourceLat*=-1;
  //   std::cout << "source lat " << sourceLat << " sourceLon " << sourceLon << std::endl;
  return 1;
}



int UsefulAdu5Pat::getSourceLonAndLatAtAlt(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat,Double_t &sourceAltitude)
{

  bool debug=false;

  if(fPhiWave!=phiWave) fPhiWave=phiWave;
  if(fThetaWave!=thetaWave) fThetaWave=thetaWave;

  //   Double_t thetaBalloon=fUPGeomTool->getThetaFromLat(TMath::Abs(latitude));
  //   Double_t phiBalloon=fUPGeomTool->getPhiFromLon(longitude);
  //   Double_t balloonHeight=fUPGeomTool->getGeoid(thetaBalloon)+altitude;

  bool endlessLoop = 0;
  bool endlessMess = 0;

  Double_t tempPhiWave=phiWave;
  Double_t tempThetaWave=TMath::PiOver2()-thetaWave;
  //Now need to take account of balloon heading

  TVector3 rollAxis,pitchAxis;
  rollAxis=fUPGeomTool->fRollRotationAxis;
  pitchAxis=fUPGeomTool->fPitchRotationAxis;

  //std::cout << "Attitude: " << heading << "\t" << pitch << "\t" << roll 
  //	       << "\n";

  if(heading>=0 && heading<=360) {
    TVector3 arbDir;
    arbDir.SetMagThetaPhi(1,tempThetaWave,-1*tempPhiWave);

    //std::cout << "tempPhiWave3: " << arbDir.Phi() << "\t" << tempPhiWave << "\n";

    arbDir.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    rollAxis.Rotate((heading)*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    pitchAxis.Rotate((heading)*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    //     std::cout << arbDir.Phi() << "\t" << arbDir.Theta() << "\n";
    //     std::cout << "Pitch Axis\t" << pitchAxis.X() << " " 
    //	       << pitchAxis.Y() << " " << pitchAxis.Z() << "\t" << pitch << "\n";


    arbDir.Rotate(pitch*TMath::DegToRad(),pitchAxis);
    //     std::cout << arbDir.Phi() << "\t" << arbDir.Theta() << "\n";
    rollAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);
    arbDir.Rotate(roll*TMath::DegToRad(),rollAxis);//roll and pitch

    tempPhiWave=arbDir.Phi();
    //std::cout << "tempPhiWave2: " << tempPhiWave << "\n";
    if(tempPhiWave>TMath::TwoPi()) {
      tempPhiWave-=TMath::TwoPi();
    }
    if(tempPhiWave<0) {
      tempPhiWave+=TMath::TwoPi();
    }
    tempThetaWave=arbDir.Theta();
  }
  else std::cout << "heading bad" << std::endl;
   
  //   std::cout << "Get source angles: " <<  tempThetaWave << "\t" << tempPhiWave << "\n";

  //Double_t re=balloonHeight-altitude+2000;
  //   Double_t re=fBalloonHeight-altitude;
   
  Double_t reBalloon=fUPGeomTool->getDistanceToCentreOfEarth(latitude);
  //   std::cout << "Radius difference: " << re-re2 << "\t" << re << "\t" << re2 << "\n";
  Double_t chosenAlt = fRampdemReader->SurfaceAboveGeoid(longitude,latitude);

  Double_t re=reBalloon+chosenAlt;
  Double_t nextRe=re;
  Double_t reh=reBalloon+altitude;

  //    if(debug){
  //      std::cout << "balloon lat " << latitude << " lon " << longitude << " re " << reBalloon+chosenAlt << " surface elevation " << chosenAlt << std::endl;
  //      std::cout << "phiW " << phiWave << " thetaW " << thetaWave << " fPhiW " << fPhiWave << " fPhiWave " << " fThetaWave " << fThetaWave << std::endl;
  //      std::cout << "head " << heading << " p " << pitch << " r " << roll << std::endl;
  //    }

  int inTheLoop=0;
  Double_t lastButOneRe=0; //in case we get stuck in the while loop - this is a bad method and needs to be improved!
  Double_t lastButTwoRe=0; //in case we get stuck in the while loop - this is a bad method and needs to be improved!
  Double_t lastButThreeRe=0; //in case we get stuck in the while loop - this is a bad method and needs to be improved!
  Double_t lastButFourRe=0; //in case we get stuck in the while loop - this is a bad method and needs to be improved!

  do {
    //Okay so effectively what we do here is switch to cartesian coords with the balloon at RE_balloon + altitude and then try to fins where the line at thetaWave cuts the Earth's surface.
    //This is iterative because the Earth is cruel and isn't flat, and I couldn't be bothered to work out how to do it more elegantly.
    if(inTheLoop>2) lastButFourRe = lastButThreeRe;
    if(inTheLoop>2) lastButThreeRe = lastButTwoRe;
    if(inTheLoop>1) lastButTwoRe = lastButOneRe;
    if(inTheLoop>0) lastButOneRe = re;
    re=nextRe;
    inTheLoop++;

    Double_t sintw=TMath::Sin(tempThetaWave);
    Double_t costw=TMath::Cos(tempThetaWave);
    //     Double_t sqrtArg=(reh*reh*costw*costw - (reh*reh-re*re));
    Double_t sqrtArg(re*re-reh*reh*sintw*sintw);
    if(sqrtArg<0) {
      // No solution possible
      //     std::cout << "No solution possible\n";
      return 0;
    }
    Double_t L=reh*costw - TMath::Sqrt(sqrtArg);
    Double_t sinThetaL=L*sintw/re;
    Double_t sourceTheta=TMath::ASin(sinThetaL);
     
    //std::cout << "Source Theta: " << sourceTheta*TMath::RadToDeg() << "\t" << tempThetaWave << "\t" << reh << "\t" << re << "\t" << fUPGeomTool->getDistanceToCentreOfEarth(latitude) << "\n";
     
    //One possible improvement is to reiterate using the radius at the source location.
     
     
    //   Double_t sinSquiggle=(fBalloonHeight/tempRE)*TMath::Sin(tempThetaWave);
    //   Double_t squiggle=TMath::ASin(sinSquiggle);
    //   if(squiggle<0) squiggle+=TMath::Pi();
     
    //   Double_t sourceTheta=TMath::Pi()-tempThetaWave-squiggle;
    //   std::cout << thetaWave*TMath::RadToDeg() << "\t" << L << "\t" << sourceTheta*TMath::RadToDeg() << "\t" << phiWave*TMath::RadToDeg() << std::endl;
    //Start at ground below balloon
    fSourcePos.SetX(0);
    fSourcePos.SetY(0);
    fSourcePos.SetZ(re);
     
    //std::cout << "00ReLoc: " << fSourcePos.X() << "\t" << fSourcePos.Y() << "\t" << fSourcePos.Z() << "\n";
     
    //Rotate to latitude relative to balloon
    fSourcePos.RotateY(sourceTheta);   
    //Rotate to longitude relative to balloon
    fSourcePos.RotateZ(-1*tempPhiWave);
     
    //std::cout << "RelBalloonLoc: " << fSourcePos.X() << "\t" << fSourcePos.Y() << "\t" << fSourcePos.Z() << "\n";
     
    //Rotate to correct absolute values
    //std::cout << "Balloon angles: " << fBalloonTheta << "\t" << fBalloonPhi << "\n";
    fSourcePos.RotateY(fBalloonTheta);
    fSourcePos.RotateZ(fBalloonPhi);
    //Goofy sign thing
    //   fSourcePos.SetZ(-1*fSourcePos.Z());
     
    Double_t sourceVec[3];
    fSourcePos.GetXYZ(sourceVec);
    //std::cout << "BalloonPos: " << fBalloonPos.X() << "\t" << fBalloonPos.Y() << "\t" << fBalloonPos.Z() << "\n";
    //std::cout << "ZeroAltLoc: " << sourceVec[0] << "\t" << sourceVec[1] << "\t" << sourceVec[2] << "\n";
    Double_t sourceAlt;
    fUPGeomTool->getLatLonAltFromCartesian(sourceVec,sourceLat,sourceLon,sourceAlt);
    //std::cout << "SourceLatLonAlt: " << sourceLat << "\t" << sourceLon << "\t" 
    //	       << sourceAlt << "\n";
    chosenAlt = fRampdemReader->SurfaceAboveGeoid(sourceLon,sourceLat);
    nextRe=fUPGeomTool->getDistanceToCentreOfEarth(sourceLat)+chosenAlt;
    //std::cout << "Earth radius: " << nextRe << "\t" << re << "\n";
    //     break;

    if(debug)
      if(inTheLoop>30 && inTheLoop<50) std::cout << "loop " << inTheLoop << " sourceLat " << sourceLat << " long " << sourceLon << " alt " << chosenAlt << " re " << re << " nextRe " << nextRe << " lastButOneRe " << lastButOneRe << " abs(last-next) " << TMath::Abs(lastButOneRe-nextRe) << std::endl;
   

    if(TMath::Abs(lastButOneRe-nextRe)<1.){
      if(debug)
	std::cout << "breaking out nextRe " << nextRe << " lastButTwoRe " << lastButOneRe << std::endl;
      endlessLoop=1;
      break;//WAY OF MAKING SURE WE DON'T GET STUCK IN THE LOOP WITH ALTERNATING RE VALUES - THIS IS NOT THE WAY TO DO THINGS!
    }
    if(TMath::Abs(lastButTwoRe-nextRe)<1.){
      if(debug)
	std::cout << "breaking out nextRe " << nextRe << " lastButOneRe " << lastButTwoRe << std::endl;
      endlessLoop=1;
      break;//WAY OF MAKING SURE WE DON'T GET STUCK IN THE LOOP WITH ALTERNATING RE VALUES - THIS IS NOT THE WAY TO DO THINGS!
    }
    if(TMath::Abs(lastButThreeRe-nextRe)<1.){
      if(debug)
	std::cout << "breaking out nextRe " << nextRe << " lastButThreeRe " << lastButThreeRe << std::endl;
      endlessLoop=1;
      break;//WAY OF MAKING SURE WE DON'T GET STUCK IN THE LOOP WITH ALTERNATING RE VALUES - THIS IS NOT THE WAY TO DO THINGS!
    }
    if(TMath::Abs(lastButFourRe-nextRe)<1.){
      if(debug)
	std::cout << "breaking out nextRe " << nextRe << " lastButFourRe " << lastButFourRe << std::endl;
      endlessLoop=1;
      break;//WAY OF MAKING SURE WE DON'T GET STUCK IN THE LOOP WITH ALTERNATING RE VALUES - THIS IS NOT THE WAY TO DO THINGS!
    }

    //fuck it, no sodding loops that go on too long!
    if(inTheLoop>50){
      endlessMess=1;
      break;
    }

  } while(TMath::Abs(nextRe-re)>1);

  //    if(debug)
  //      std::cout << "sourceLat " << sourceLat << " long " << sourceLon << " alt " << chosenAlt << " re " << re << " nextRe " << nextRe << std::endl;
  sourceAltitude = chosenAlt;
  //   fUPGeomTool->getLonLat(fSourcePos,sourceLon,sourceLat);
  //   sourceLat*=-1;

  if(endlessLoop==1) return 2;
  else if(endlessMess==1) return 3;
  else if(int(chosenAlt)==-9999) return 4;
  else
    return 1;
}



void UsefulAdu5Pat::getThetaAndPhiWaveWaisDivide(Double_t &thetaWave, Double_t &phiWave)
{   
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_WAIS,AnitaLocations::LATITUDE_WAIS,AnitaLocations::ALTITUDE_WAIS,thetaWave,phiWave);
}



void UsefulAdu5Pat::getThetaAndPhiWave(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t &thetaWave, Double_t &phiWave) {
  // Double_t thetaBalloon=fUPGeomTool->getThetaFromLat(TMath::Abs(latitude));
  // Double_t phiBalloon=fUPGeomTool->getPhiFromLon(longitude);
  // Double_t balloonHeight=fUPGeomTool->getGeoid(thetaBalloon)+altitude;
  // std::cout << "Theta " << thetaBalloon << "\t" << fBalloonTheta << "\n";
  // std::cout << "Phi " << phiBalloon << "\t" << fBalloonPhi << "\n";
   
  //    Double_t thetaSource=fUPGeomTool->getThetaFromLat(TMath::Abs(sourceLat));
  //    Double_t phiSource=fUPGeomTool->getPhiFromLon(sourceLon);
  //    Double_t radiusSource=fUPGeomTool->getGeoid(thetaSource)+sourceAlt;

  //    //Get vector from Earth's centre to source
  //    fSourcePos.SetX(radiusSource*TMath::Sin(thetaSource)*TMath::Cos(phiSource));
  //    fSourcePos.SetY(radiusSource*TMath::Sin(thetaSource)*TMath::Sin(phiSource));
  //    fSourcePos.SetZ(radiusSource*TMath::Cos(thetaSource));

  Double_t pSource[3]={0};
  fUPGeomTool->getCartesianCoords(TMath::Abs(sourceLat),sourceLon,sourceAlt,pSource);
  fSourcePos.SetXYZ(pSource[0],pSource[1],pSource[2]);

  //std::cout << "SourceLoc: " << pSource[0] << "\t" << pSource[1] << "\t" << pSource[2] << "\n";
   
  //Rotate such that balloon is at 0,0,fBalloonHeight
  fSourcePos.RotateZ(-1*fBalloonPhi);
  fSourcePos.RotateY(-1*fBalloonTheta);


  //    std::cout << "rotated source\tx " << fSourcePos.X() << "\ty " << fSourcePos.Y() << "\tz " << fSourcePos.Z() << std::endl;
  //    std::cout << "balloon\tx " << fBalloonPos.X() << "\ty " << fBalloonPos.Y() << "\tz " << fBalloonPos.Z() << std::endl;
  //    std::cout << "rotated source - balloon\tx " << fSourcePos.X()-fBalloonPos.X() << "\ty " << fSourcePos.Y()-fBalloonPos.Y() << "\tz " << fSourcePos.Z()-fBalloonPos.Z() << std::endl;

  //Now find thetaWave and phiWave
  thetaWave=TMath::ATan((fBalloonHeight-fSourcePos.Z())/TMath::Sqrt(fSourcePos.X()*fSourcePos.X() + fSourcePos.Y()*fSourcePos.Y()));
   
  //phiWave is just atan(yp/xp) only looks confusing to make sure I get the sign and 0-360 convention
  phiWave=0;
  if(fSourcePos.X()==0) {
    phiWave=TMath::PiOver2();
    if(fSourcePos.Y()<0)
      phiWave+=TMath::Pi();
  }
  else if(fSourcePos.X()<0) {
    phiWave=TMath::Pi()+TMath::ATan(fSourcePos.Y()/fSourcePos.X());
  }
  else {
    phiWave=TMath::ATan(fSourcePos.Y()/fSourcePos.X());
    if(fSourcePos.Y()<0) {
      phiWave+=TMath::TwoPi();
    }
  }   

  //Now need to take account of balloon heading
  //   //Will have to check heading at some point
  //    if(heading>=0 && heading<=360) {
  //       phiWave+=heading*TMath::DegToRad();
  //       if(phiWave>TMath::TwoPi())
  // 	 phiWave-=TMath::TwoPi();
  //    }

  //   std::cout << "Pre rotate: " << TMath::PiOver2()-thetaWave << "\t" << TMath::TwoPi()-phiWave << "\n";

  TVector3 rollAxis,pitchAxis;
  rollAxis=fUPGeomTool->fRollRotationAxis;
  pitchAxis=fUPGeomTool->fPitchRotationAxis;

  Double_t tempThetaWave=TMath::PiOver2()-thetaWave;

  if(heading>=0 && heading<=360) {
    TVector3 arbDir;
    arbDir.SetMagThetaPhi(1,tempThetaWave,-1*phiWave);
     
    rollAxis.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    pitchAxis.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    rollAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);
     
    arbDir.Rotate(-1.*roll*TMath::DegToRad(),rollAxis);
    arbDir.Rotate(-1*pitch*TMath::DegToRad(),pitchAxis);

    // std::cout << roll << "\t" << pitch << "\t in function " << __PRETTY_FUNCTION__ << std::endl;
     
    arbDir.Rotate(-1*heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);

    //need to rotate roll and pitch axes - in this function heading has been
    //chosen as +ve x axis, pitch and roll axes are defined in terms of x and
    //y axes from ANITA setup document
    //rollAxis.Rotate(-135.*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    //pitchAxis.Rotate(45.*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
    //pitchAxis.Rotate(-1*heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);

    phiWave=arbDir.Phi();
    phiWave*=-1;
    if(phiWave>TMath::TwoPi()) {
      phiWave-=TMath::TwoPi();
    }
    if(phiWave<0) {
      phiWave+=TMath::TwoPi();
    }
    thetaWave=TMath::PiOver2()-arbDir.Theta();
     
  }

  fPhiWave=phiWave;
  fThetaWave=thetaWave;
  fSourceLongitude=sourceLon;
  fSourceLatitude=sourceLat;
  fSourceAltitude=sourceAlt;
}

Double_t UsefulAdu5Pat::getGroupDelay(Double_t phiToAntBoresight)
{
  Double_t phiDeg=phiToAntBoresight*TMath::RadToDeg();
  Double_t delayTime=(phiDeg*phiDeg*phiDeg*phiDeg)*1.45676e-8;
  delayTime-=(phiDeg*phiDeg)*5.01452e-6;
  return delayTime;
}





void UsefulAdu5Pat::getThetaWaveAtBase(Double_t baseLon, Double_t baseLat, Double_t baseAlt, Double_t &thetaWave) {

  //gets the theta from a base to the balloon - can check if we are beyond horizon

  Double_t pSource[3]={0};
  fUPGeomTool->getCartesianCoords(TMath::Abs(baseLat),baseLon,baseAlt,pSource);
  fSourcePos.SetXYZ(pSource[0],pSource[1],pSource[2]);

  //make copy of balloon position and rotate it such that base is at 0,0,baseAlt
  TVector3 rotatedBalloonPos = fBalloonPos;
  rotatedBalloonPos.RotateZ(-1*fSourcePos.Phi());   
  rotatedBalloonPos.RotateY(-1*fSourcePos.Theta()); 

  TVector3 rotatedBase = fSourcePos;
  rotatedBase.RotateZ(-1*fSourcePos.Phi());   
  rotatedBase.RotateY(-1*fSourcePos.Theta()); 


  //    std::cout << "rotated balloon\tx " << rotatedBalloonPos.X() << "\ty " << rotatedBalloonPos.Y() << "\tz " << rotatedBalloonPos.Z() << std::endl;
  //    std::cout << "base\tx " << fSourcePos.X() << "\ty " << fSourcePos.Y() << "\tz " << fSourcePos.Z() << std::endl;
  //    std::cout << "rotated base\tx " << rotatedBase.X() << "\ty " << rotatedBase.Y() << "\tz " << rotatedBase.Z() << std::endl;
  //    std::cout << "rotated balloon - base\tx " << rotatedBalloonPos.X()-fSourcePos.X() << "\ty " << rotatedBalloonPos.Y()-fSourcePos.Y() << "\tz " << rotatedBalloonPos.Z()-fSourcePos.Z() << std::endl;

  //Now find thetaWave
  thetaWave=TMath::ATan((rotatedBase.Z()-rotatedBalloonPos.Z())/TMath::Sqrt(rotatedBalloonPos.X()*rotatedBalloonPos.X() + rotatedBalloonPos.Y()*rotatedBalloonPos.Y())); 




  //    Double_t pSource[3]={0};
  //    fUPGeomTool->getCartesianCoords(TMath::Abs(sourceLat),sourceLon,sourceAlt,pSource);
  //    fSourcePos.SetXYZ(pSource[0],pSource[1],pSource[2]);

  //    //std::cout << "SourceLoc: " << pSource[0] << "\t" << pSource[1] << "\t" << pSource[2] << "\n";
   
  //    //Rotate such that balloon is at 0,0,fBalloonHeight
  //    fSourcePos.RotateZ(-1*fBalloonPhi);
  //    fSourcePos.RotateY(-1*fBalloonTheta);

  //    //Now find thetaWave and phiWave
  //    thetaWave=TMath::ATan((fBalloonHeight-fSourcePos.Z())/TMath::Sqrt(fSourcePos.X()*fSourcePos.X() + fSourcePos.Y()*fSourcePos.Y()));

}


Double_t UsefulAdu5Pat::getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt)
{

  Double_t tempTheta,tempPhi;
  if(fSourceAltitude!=sourceAlt || fSourceLongitude!=sourceLon || fSourceLatitude!=sourceLat)
    getThetaAndPhiWave(sourceLon,sourceLat,sourceAlt,tempTheta,tempPhi);

  //Now fThetaWave and fPhiWave should be correctly set.
  Double_t phi1=fUPGeomTool->getAntPhiPositionRelToAftFore(ant1);
  Double_t r1=fUPGeomTool->getAntR(ant1);
  Double_t z1=fUPGeomTool->getAntZ(ant1);

  Double_t phi2=fUPGeomTool->getAntPhiPositionRelToAftFore(ant2);
  Double_t r2=fUPGeomTool->getAntR(ant2);
  Double_t z2=fUPGeomTool->getAntZ(ant2);

  Double_t tanThetaW=TMath::Tan(fThetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(fPhiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(fPhiWave-phi2);

  Double_t geomTime= 1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
  if(fIncludeGroupDelay) {
    Double_t phi1Diff=fUPGeomTool->getPhiDiff(fPhiWave,phi1);
    Double_t delay1=getGroupDelay(phi1Diff);
    Double_t phi2Diff=fUPGeomTool->getPhiDiff(fPhiWave,phi2);
    Double_t delay2=getGroupDelay(phi2Diff);
    geomTime+=(delay1-delay2);
  }
   

  return geomTime;
}



Double_t UsefulAdu5Pat::getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t phiWave, Double_t thetaWave)
{

  if(fPhiWave!=phiWave || fThetaWave!=thetaWave){
    fThetaWave=thetaWave;
    fPhiWave=phiWave;
  }

  Double_t phi1=fUPGeomTool->getAntPhiPositionRelToAftFore(ant1);
  Double_t r1=fUPGeomTool->getAntR(ant1);
  Double_t z1=fUPGeomTool->getAntZ(ant1);

  Double_t phi2=fUPGeomTool->getAntPhiPositionRelToAftFore(ant2);
  Double_t r2=fUPGeomTool->getAntR(ant2);
  Double_t z2=fUPGeomTool->getAntZ(ant2);

  Double_t tanThetaW=TMath::Tan(fThetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(fPhiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(fPhiWave-phi2);
   
  return  1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}


//does not take into account anita's position/heading...
Double_t UsefulAdu5Pat::getDeltaTExpected(Int_t ant1,Int_t ant2,Double_t cosPhi,Double_t sinPhi,Double_t cosTheta,Double_t sinTheta)
{
  Double_t x1,y1,z1;
  fUPGeomTool->getAntXYZ(ant1,x1,y1,z1);

  Double_t x2,y2,z2;
  fUPGeomTool->getAntXYZ(ant2,x2,y2,z2);

  Double_t part1 = x1*cosTheta*cosPhi + y1*cosTheta*sinPhi + z1*sinTheta;
  Double_t part2 = x2*cosTheta*cosPhi + y2*cosTheta*sinPhi + z2*sinTheta;
  return 1e9 * (part1 - part2) / C_LIGHT;

}


Double_t UsefulAdu5Pat::getDeltaTWillySeavey(Int_t ant1, Int_t ant2)
{
  return getDeltaTExpected(ant1, ant2,
			   AnitaLocations::LONGITUDE_SURF_SEAVEY,
			   AnitaLocations::LATITUDE_SURF_SEAVEY,
			   AnitaLocations::ALTITUDE_SURF_SEAVEY);
}

Double_t UsefulAdu5Pat::getDeltaTWillyBorehole(Int_t ant1, Int_t ant2)
{
  return getDeltaTExpected(ant1, ant2,
			   AnitaLocations::LONGITUDE_BH,
			   AnitaLocations::LATITUDE_BH,
			   AnitaLocations::ALTITUDE_BH);
}

Double_t UsefulAdu5Pat::getDeltaTTaylor(Int_t ant1, Int_t ant2)
{
  return getDeltaTExpected(ant1, ant2,
			   AnitaLocations::LONGITUDE_TD,
			   AnitaLocations::LATITUDE_TD,
			   AnitaLocations::ALTITUDE_TD);
}

Double_t UsefulAdu5Pat::getDeltaTTaylorOpt(Int_t ant1, Int_t ant2, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi)
{
  return getDeltaTExpectedOpt(ant1, ant2,
			      AnitaLocations::LONGITUDE_TD,
			      AnitaLocations::LATITUDE_TD,
			      AnitaLocations::ALTITUDE_TD,
			      deltaR, deltaZ,deltaPhi);
}

Double_t UsefulAdu5Pat::getDeltaTExpectedOpt(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi)
{
  Double_t tempTheta,tempPhi;
  if(fSourceAltitude!=sourceAlt || fSourceLongitude!=sourceLon || fSourceLatitude!=sourceLat)
    getThetaAndPhiWave(sourceLon,sourceLat,sourceAlt,tempTheta,tempPhi);

  //Now fThetaWave and fPhiWave should be correctly set.
  Double_t phi1=fUPGeomTool->getAntPhiPositionRelToAftFore(ant1)+deltaPhi[ant1];
  Double_t r1=fUPGeomTool->getAntR(ant1)+deltaR[ant1];
  Double_t z1=fUPGeomTool->getAntZ(ant1)+deltaZ[ant1];

  Double_t phi2=fUPGeomTool->getAntPhiPositionRelToAftFore(ant2)+deltaPhi[ant2];
  Double_t r2=fUPGeomTool->getAntR(ant2)+deltaR[ant2];
  Double_t z2=fUPGeomTool->getAntZ(ant2)+deltaZ[ant2];

  //   std::cout << ant1 << deltaPhi[ant1] << "  " << deltaPhi[ant2]<< std::endl;

  Double_t tanThetaW=TMath::Tan(fThetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(fPhiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(fPhiWave-phi2);
   
  return  1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}


Double_t UsefulAdu5Pat::getDeltaTSeaveyOpt(Int_t ant1, Int_t ant2, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi)
{
  return getDeltaTExpectedSeaveyOpt(ant1, ant2,
				    AnitaLocations::LONGITUDE_SURF_SEAVEY,
				    AnitaLocations::LATITUDE_SURF_SEAVEY,
				    AnitaLocations::ALTITUDE_SURF_SEAVEY,
				    deltaR, deltaZ,deltaPhi);
}

Double_t UsefulAdu5Pat::getDeltaTExpectedSeaveyOpt(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi)
{
  Double_t tempTheta,tempPhi;
  if(fSourceAltitude!=sourceAlt || fSourceLongitude!=sourceLon || fSourceLatitude!=sourceLat)
    getThetaAndPhiWave(sourceLon,sourceLat,sourceAlt,tempTheta,tempPhi);

  //Now fThetaWave and fPhiWave should be correctly set.
  Double_t phi1=fUPGeomTool->getAntPhiPositionRelToAftFore(ant1)+deltaPhi[ant1];
  Double_t r1=fUPGeomTool->getAntR(ant1)+deltaR[ant1];
  Double_t z1=fUPGeomTool->getAntZ(ant1)+deltaZ[ant1];

  Double_t phi2=fUPGeomTool->getAntPhiPositionRelToAftFore(ant2)+deltaPhi[ant2];
  Double_t r2=fUPGeomTool->getAntR(ant2)+deltaR[ant2];
  Double_t z2=fUPGeomTool->getAntZ(ant2)+deltaZ[ant2];

  //   std::cout << ant1 << deltaPhi[ant1] << "  " << deltaPhi[ant2]<< std::endl;

  Double_t tanThetaW=TMath::Tan(fThetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(fPhiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(fPhiWave-phi2);
   
  return  1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}


UInt_t UsefulAdu5Pat::getTaylorDomeTriggerTimeNs()
{
  return getTriggerTimeNsFromSource(AnitaLocations::LATITUDE_TD,
				    AnitaLocations::LONGITUDE_TD,
				    AnitaLocations::ALTITUDE_TD);
}

UInt_t UsefulAdu5Pat::getWaisDivideTriggerTimeNs()
{
  return getTriggerTimeNsFromSource(AnitaLocations::LATITUDE_WAIS,
				    AnitaLocations::LONGITUDE_WAIS,
				    AnitaLocations::ALTITUDE_WAIS);
}



UInt_t UsefulAdu5Pat::getTriggerTimeNsFromSource(Double_t sourceLat, Double_t sourceLong, Double_t sourceAlt)
{
   
  //    Double_t thetaBalloon=fUPGeomTool->getThetaFromLat(TMath::Abs(latitude));
  //    Double_t phiBalloon=fUPGeomTool->getPhiFromLon(longitude);
  //    Double_t balloonHeight=fUPGeomTool->getGeoid(thetaBalloon)+altitude;
   
  //    static Double_t thetaTaylor=fUPGeomTool->getThetaFromLat(TMath::Abs(AnitaLocations::LATITUDE_TD));
  //    static Double_t phiTaylor=fUPGeomTool->getPhiFromLon(AnitaLocations::LONGITUDE_TD);
  //    static Double_t radiusTaylor=fUPGeomTool->getGeoid(thetaTaylor)+AnitaLocations::ALTITUDE_TD;
   

  //    //Get vector from Earth's centre to taylor
  //    TVector3 fTaylorPos;
  //    fTaylorPos.SetX(radiusTaylor*TMath::Sin(thetaTaylor)*TMath::Cos(phiTaylor));
  //    fTaylorPos.SetY(radiusTaylor*TMath::Sin(thetaTaylor)*TMath::Sin(phiTaylor));
  //    fTaylorPos.SetZ(radiusTaylor*TMath::Cos(thetaTaylor));

  static Int_t firstTime=1;   
  static Double_t pTaylor[3]={0};
  if(firstTime) {
    fUPGeomTool->getCartesianCoords(TMath::Abs(sourceLat),
				    sourceLong,
				    sourceAlt,
				    pTaylor);
  }

  //   std::cout << "Old Geoid Model:\t" << fTaylorPos.x() << "\t" << fTaylorPos.y() << "\t" << fTaylorPos.z() << "\n";
  //    std::cout << "New Model:\t" << pTaylor[0] << "\t" << pTaylor[1] << "\t" << pTaylor << "\n";   

   
  //    //Rotate such that balloon is at 0,0,fBalloonHeight
  //    fTaylorPos.RotateZ(-1*phiBalloon);
  //    fTaylorPos.RotateY(-1*thetaBalloon);
  //    Double_t s2=(fTaylorPos.x()*fTaylorPos.x()) + (fTaylorPos.y()*fTaylorPos.y()) + TMath::Power(fBalloonHeight-fTaylorPos.z(),2);
  Double_t s2=(pTaylor[0]-fBalloonCoords[0])*(pTaylor[0]-fBalloonCoords[0])+
    (pTaylor[1]-fBalloonCoords[1])*(pTaylor[1]-fBalloonCoords[1])+
    (pTaylor[2]-fBalloonCoords[2])*(pTaylor[2]-fBalloonCoords[2]);

  // std::cout << fBalloonCoords[0] << "\t" << fBalloonCoords[1] << "\t" << fBalloonCoords[2] << std::endl;
   

  Double_t distanceToFly=TMath::Sqrt(s2);
  Double_t timeOfFlight=distanceToFly/C_LIGHT;
  timeOfFlight*=1e9;
  //   Double_t expTime=timeOfFlight-40e3;
  Double_t expTime=timeOfFlight;//-40e3;
  UInt_t expTrigTime=(UInt_t)expTime;
   
  //   std::cout << distanceToFly << "\t" << timeOfFlight << "\n";
  return expTrigTime;


}



Double_t UsefulAdu5Pat::getAngleBetweenPayloadAndSource(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt){ //ACG additional function

  Double_t thetaBalloon=fUPGeomTool->getThetaFromLat(TMath::Abs(latitude));
  Double_t phiBalloon=fUPGeomTool->getPhiFromLon(longitude);
  Double_t balloonHeight=fUPGeomTool->getGeoid(thetaBalloon)+altitude;

  Double_t thetaSource=fUPGeomTool->getThetaFromLat(TMath::Abs(sourceLat));
  Double_t phiSource=fUPGeomTool->getPhiFromLon(sourceLon);
  Double_t radiusSource=fUPGeomTool->getGeoid(thetaSource)+sourceAlt;

  //Get vector from Earth's centre to source

  fSourcePos.SetX(radiusSource*TMath::Sin(thetaSource)*TMath::Cos(phiSource));
  fSourcePos.SetY(radiusSource*TMath::Sin(thetaSource)*TMath::Sin(phiSource));
  fSourcePos.SetZ(radiusSource*TMath::Cos(thetaSource));

  //Rotate such that balloon is at 0,0,balloonHeight

  fSourcePos.RotateZ(-1*phiBalloon);
  fSourcePos.RotateY(-1*thetaBalloon);

  Double_t dotProduct=balloonHeight*fSourcePos.Z();//+0*X+0*Y
  Double_t magSource=TMath::Sqrt(fSourcePos.X()*fSourcePos.X() + fSourcePos.Y()*fSourcePos.Y()+fSourcePos.Z()*fSourcePos.Z());
  Double_t magBalloon=balloonHeight;

  Double_t phiEarthCenter=TMath::ACos(TMath::Abs(dotProduct/(magSource*magBalloon))); 

  return phiEarthCenter;

} 


// TProfile2D *UsefulAdu5Pat::rampMap(int coarseness,UInt_t &xBins,UInt_t &yBins){
//   TProfile2D *rampMap = fRampdemReader->rampMap(coarseness,0,xBins,yBins);
//   return rampMap;
// }
