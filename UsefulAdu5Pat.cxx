//////////////////////////////////////////////////////////////////////////////
/////  UsefulAdu5Pat.cxx        ANITA event reading class                  /////
/////                                                                    /////
/////  Description:                                                      /////
/////    A simple class for utilising Adu5Pat objects                      /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include <iostream>
#include <fstream>

#include "TVector3.h"

//#define CMINCH 2.54

ClassImp(UsefulAdu5Pat);

AnitaGeomTool *fUPGeomTool=0;


UsefulAdu5Pat::UsefulAdu5Pat() 
   : Adu5Pat()
{
   //Default Constructor
   fThetaWave=0;
   fPhiWave=0;
   fSourceLongitude=-1;
   fSourceLatitude=-1;
   if(!fUPGeomTool){
      fUPGeomTool=AnitaGeomTool::Instance();
   }
}

UsefulAdu5Pat::UsefulAdu5Pat(Adu5Pat *patPtr)
   : Adu5Pat(*patPtr)
{
   fThetaWave=0;
   fPhiWave=0;
   fSourceLongitude=-1;
   fSourceLatitude=-1;
   if(!fUPGeomTool){
      fUPGeomTool=AnitaGeomTool::Instance();
      std::cout << "making the geom" << std::endl;
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



int UsefulAdu5Pat::getSourceLonAndLatAltZero(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat)
{

  if(fPhiWave!=phiWave) fPhiWave=phiWave;
  if(fThetaWave!=thetaWave) fThetaWave=thetaWave;

   Double_t thetaBalloon=fUPGeomTool->getThetaFromLat(TMath::Abs(latitude));
   Double_t phiBalloon=fUPGeomTool->getPhiFromLon(longitude);
   Double_t balloonHeight=fUPGeomTool->getGeoid(thetaBalloon)+altitude;

   Double_t tempPhiWave=phiWave;
   Double_t tempThetaWave=TMath::PiOver2()-thetaWave;
   //Now need to take account of balloon heading
   //Will have to check heading at some point

   //double tempphi2 = tempPhiWave;
   //   if(heading>=0 && heading<=360) {
   // //     tempPhiWave-=heading*TMath::DegToRad();
   // tempphi2=heading*TMath::DegToRad()-tempphi2;
   // if(tempphi2<0)
   // tempphi2+=TMath::TwoPi();
   //}

   

   TVector3 rollAxis,pitchAxis;
   rollAxis=fUPGeomTool->fRollRotationAxis;
   pitchAxis=fUPGeomTool->fPitchRotationAxis;

   if(heading>=0 && heading<=360) {
     TVector3 arbDir;
     arbDir.SetMagThetaPhi(1,tempThetaWave,-1.*tempPhiWave);

     arbDir.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
     rollAxis.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
     pitchAxis.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);

     arbDir.Rotate(pitch*TMath::DegToRad(),pitchAxis);
     rollAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);

     arbDir.Rotate(roll*TMath::DegToRad(),rollAxis);//roll and pitch

     tempPhiWave=arbDir.Phi();
     if(tempPhiWave>TMath::TwoPi()) {
       tempPhiWave-=TMath::TwoPi();
     }
     if(tempPhiWave<0) {
       tempPhiWave+=TMath::TwoPi();
     }
     tempThetaWave=arbDir.Theta();
   }
   else std::cout << "heading bad" << std::endl;
   
   //Double_t re=balloonHeight-altitude+2000;
   Double_t re=balloonHeight-altitude;
   Double_t reh=balloonHeight;
   Double_t costw=TMath::Cos(tempThetaWave);
   Double_t sqrtArg=(reh*reh*costw*costw - (reh*reh-re*re));
   if(sqrtArg<0) {
     // No solution possible
     return 0;
   }
   Double_t L=reh*costw - TMath::Sqrt(sqrtArg);
   Double_t sinThetaL=L*TMath::Sin(tempThetaWave)/re;
   Double_t sourceTheta=TMath::ASin(sinThetaL);

   

   //   Double_t sinSquiggle=(balloonHeight/tempRE)*TMath::Sin(tempThetaWave);
   //   Double_t squiggle=TMath::ASin(sinSquiggle);
   //   if(squiggle<0) squiggle+=TMath::Pi();
   
   //   Double_t sourceTheta=TMath::Pi()-tempThetaWave-squiggle;
   //   std::cout << thetaWave*TMath::RadToDeg() << "\t" << L << "\t" << sourceTheta*TMath::RadToDeg() << "\t" << phiWave*TMath::RadToDeg() << std::endl;
   //Start at ground below balloon
   fSourcePos.SetX(0);
   fSourcePos.SetY(0);
   fSourcePos.SetZ(re);
   
   //Rotate to latitude relative to balloon
   fSourcePos.RotateY(sourceTheta);   
   //Rotate to longitude relative to balloon
   fSourcePos.RotateZ(-1*tempPhiWave);

   //Rotate to correct absolute values
   fSourcePos.RotateY(thetaBalloon);
   fSourcePos.RotateZ(phiBalloon);
   //Goofy sign thing
   //   fSourcePos.SetZ(-1*fSourcePos.Z());
   
   fUPGeomTool->getLonLat(fSourcePos,sourceLon,sourceLat);
   sourceLat*=-1;
   //std::cout << "source lat " << sourceLat << " sourceLon " << sourceLon << std::endl;
   return 1;
}



void UsefulAdu5Pat::getThetaAndPhiWave(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t &thetaWave, Double_t &phiWave) {
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

   //Now find thetaWave and phiWave
   thetaWave=TMath::ATan((balloonHeight-fSourcePos.Z())/TMath::Sqrt(fSourcePos.X()*fSourcePos.X() + fSourcePos.Y()*fSourcePos.Y()));
   
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

   TVector3 rollAxis,pitchAxis;
   rollAxis=fUPGeomTool->fRollRotationAxis;
   pitchAxis=fUPGeomTool->fPitchRotationAxis;

   if(heading>=0 && heading<=360) {
     TVector3 arbDir;
     arbDir.SetMagThetaPhi(1,thetaWave,phiWave);

     arbDir.Rotate(heading*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);

     //need to rotate roll and pitch axes - in this function heading has been
     //chosen as +ve x axis, pitch and roll axes are defined in terms of x and
     //y axes from ANITA setup document
     rollAxis.Rotate(-135.*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);
     pitchAxis.Rotate(45.*TMath::DegToRad(),fUPGeomTool->fHeadingRotationAxis);

     arbDir.Rotate(pitch*TMath::DegToRad(),pitchAxis);
     arbDir.Rotate(roll*TMath::DegToRad(),rollAxis);


     phiWave=arbDir.Phi();
     if(phiWave>TMath::TwoPi()) {
       phiWave-=TMath::TwoPi();
     }
     if(phiWave<0) {
       phiWave+=TMath::TwoPi();
     }
     thetaWave=arbDir.Theta();
   }

   fPhiWave=phiWave;
   fThetaWave=thetaWave;
   fSourceLongitude=sourceLon;
   fSourceLatitude=sourceLat;
   fSourceAltitude=sourceAlt;
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
   
   return  1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
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
   return getDeltaTExpected(ant1,ant2,AnitaLocations::LONGITUDE_SURF_SEAVEY,AnitaLocations::LATITUDE_SURF_SEAVEY,AnitaLocations::ALTITUDE_SURF_SEAVEY);
}

Double_t UsefulAdu5Pat::getDeltaTWillyBorehole(Int_t ant1, Int_t ant2)
{
   return getDeltaTExpected(ant1,ant2,AnitaLocations::LONGITUDE_BH,AnitaLocations::LATITUDE_BH,AnitaLocations::ALTITUDE_BH);
}

Double_t UsefulAdu5Pat::getDeltaTTaylor(Int_t ant1, Int_t ant2)
{
   return getDeltaTExpected(ant1,ant2,AnitaLocations::LONGITUDE_TD,AnitaLocations::LATITUDE_TD,AnitaLocations::ALTITUDE_TD);
}


Double_t UsefulAdu5Pat::getDeltaTTaylorOpt(Int_t ant1, Int_t ant2, Double_t *deltaR, Double_t *deltaZ, Double_t *deltaPhi)
{
  return getDeltaTExpectedOpt(ant1,ant2,AnitaLocations::LONGITUDE_TD,AnitaLocations::LATITUDE_TD,AnitaLocations::ALTITUDE_TD, deltaR, deltaZ,deltaPhi);
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

   Double_t tanThetaW=TMath::Tan(fThetaWave);
   Double_t part1=z1*tanThetaW - r1 * TMath::Cos(fPhiWave-phi1);
   Double_t part2=z2*tanThetaW - r2 * TMath::Cos(fPhiWave-phi2);
   
   return  1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}


