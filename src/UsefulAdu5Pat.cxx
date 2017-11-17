//////////////////////////////////////////////////////////////////////////////
/////  UsefulAdu5Pat.cxx        ANITA event reading class                  /////
/////                                                                    /////
/////  Description:                                                      /////
/////    A simple class for utilising Adu5Pat objects                      /////
/////  Author: Ryan Nichol (rjn@hep.ucl.ac.uk)                           /////
//////////////////////////////////////////////////////////////////////////////

#include "AnitaVersion.h"
#include "UsefulAdu5Pat.h"
#include "AnitaGeomTool.h"
#include "AnitaConventions.h"
#include "BedmapReader.h"
#include "RampdemReader.h"
#include "TProfile2D.h"
#include <iostream>
#include <fstream>

#include "TVector3.h"

ClassImp(UsefulAdu5Pat);

UsefulAdu5Pat::UsefulAdu5Pat()
  : Adu5Pat()
{
  //Default Constructor
  fIncludeGroupDelay=0;
  pitch = AnitaStaticAdu5Offsets::pitch;
  roll = AnitaStaticAdu5Offsets::roll;
  heading += AnitaStaticAdu5Offsets::heading;
  longitude = 0;
  latitude = 0;
  altitude = 0;
  // fThetaWave=0;
  // fPhiWave=0;
  // fSourceLongitude=-1;
  // fSourceLatitude=-1;
  fBalloonCoords[0]=0;
  fBalloonCoords[1]=0;
  fBalloonCoords[2]=0;
  fBalloonHeight=0;
  fBalloonLonCache=0;
  fBalloonLatCache=0;
  fBalloonAltCache=0;
  fDebug = false;

}

UsefulAdu5Pat::UsefulAdu5Pat(const Adu5Pat *patPtr)
  : Adu5Pat(*patPtr)
{
  fIncludeGroupDelay=0;

  pitch = AnitaStaticAdu5Offsets::pitch;
  roll = AnitaStaticAdu5Offsets::roll;
  heading += AnitaStaticAdu5Offsets::heading;

  if(heading>=360){
    heading-=360;
  }
  if(heading<0){
    heading+=360;
  }
  // fThetaWave=0;
  // fPhiWave=0;
  // fSourceLongitude=-1;
  // fSourceLatitude=-1;
  fBalloonLonCache=0;
  fBalloonLatCache=0;
  fBalloonAltCache=0;
  fDebug = false;
  updateCartesianBalloonInfo();
}


UsefulAdu5Pat::~UsefulAdu5Pat()
{

}


void UsefulAdu5Pat::updateCartesianBalloonInfo(){
  
  if(longitude!=fBalloonLonCache ||
     latitude!=fBalloonLatCache ||
     altitude!=fBalloonAltCache){

    AnitaGeomTool::Instance()->getCartesianCoords(latitude,
						  longitude,
						  altitude,
						  fBalloonCoords);
    fBalloonPos.SetXYZ(fBalloonCoords[0],fBalloonCoords[1],fBalloonCoords[2]);
    fBalloonTheta=fBalloonPos.Theta();
    fBalloonPhi=fBalloonPos.Phi();
    if(fBalloonPhi<0){
      fBalloonPhi+=TMath::TwoPi();
    }
    fBalloonHeight=fBalloonPos.Mag();

    fBalloonLonCache = longitude;
    fBalloonLatCache = latitude;
    fBalloonAltCache = altitude;
  }
}


void UsefulAdu5Pat::getThetaAndPhiWaveWillySeavey(Double_t &thetaWave, Double_t &phiWave)
{
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_SURF_SEAVEY,
			    AnitaLocations::LATITUDE_SURF_SEAVEY,
			    AnitaLocations::ALTITUDE_SURF_SEAVEY,
			    thetaWave, phiWave);
}

void UsefulAdu5Pat::getThetaAndPhiWaveWillyBorehole(Double_t &thetaWave,Double_t &phiWave)
{
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_BH,
			    AnitaLocations::LATITUDE_BH,
			    AnitaLocations::ALTITUDE_BH,
			    thetaWave, phiWave);
}


void UsefulAdu5Pat::getThetaAndPhiWaveTaylorDome(Double_t &thetaWave, Double_t &phiWave)
{
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_TD,
			    AnitaLocations::LATITUDE_TD,
			    AnitaLocations::ALTITUDE_TD,
			    thetaWave, phiWave);
}


int UsefulAdu5Pat::getSourceLonAndLatAltZero(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat) {
  return getSourceLonAndLatAtDesiredAlt(phiWave, thetaWave, sourceLon, sourceLat, 0.0);
}





void UsefulAdu5Pat::accountForPitchAndRollInPhiWaveThetaWave(Double_t& phiWave, Double_t& thetaWave) const {


  AnitaGeomTool * geom = AnitaGeomTool::Instance();
  Double_t tempPhiWave=phiWave;
  Double_t tempThetaWave=TMath::PiOver2()-thetaWave;

  TVector3 rollAxis = geom->fRollRotationAxis;
  TVector3 pitchAxis = geom->fPitchRotationAxis;

  TVector3 arbDir;
  arbDir.SetMagThetaPhi(1,tempThetaWave,-1*tempPhiWave);

  arbDir.Rotate(heading*TMath::DegToRad(),geom->fHeadingRotationAxis);
  rollAxis.Rotate((heading)*TMath::DegToRad(),geom->fHeadingRotationAxis);
  pitchAxis.Rotate((heading)*TMath::DegToRad(),geom->fHeadingRotationAxis);

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


  thetaWave = tempThetaWave;
  phiWave = tempPhiWave;  
}




int UsefulAdu5Pat::getSourceLonAndLatAtDesiredAlt(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat, Double_t desiredAlt = 0.0)
{
  updateCartesianBalloonInfo();

  // if(fPhiWave!=phiWave) fPhiWave=phiWave;
  // if(fThetaWave!=thetaWave) fThetaWave=thetaWave;

  Double_t tempPhiWave=phiWave;
  Double_t tempThetaWave=thetaWave;

  if(heading>=0 && heading<=360) {
    accountForPitchAndRollInPhiWaveThetaWave(tempPhiWave, tempThetaWave);
  }
  else{
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ": Bad heading = " << heading << std::endl;
    sourceLon = -9999;
    sourceLat = -9999;
    desiredAlt = -9999;
    return -1;
  }

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  Double_t reBalloon=geom->getDistanceToCentreOfEarth(latitude)+desiredAlt; // ACG mod
  Double_t re=reBalloon;
  Double_t nextRe=re;
  Double_t reh=reBalloon+altitude;

  const Double_t sintw = TMath::Sin(tempThetaWave);
  const Double_t costw = TMath::Cos(tempThetaWave);

  do {
    //Okay so effectively what we do here is switch to cartesian coords with the balloon at RE_balloon + altitude and then try to find where the line at thetaWave cuts the Earth's surface.
    //This is iterative because the Earth is cruel and isn't flat, and I couldn't be bothered to work out how to do it more elegantly.
    re=nextRe;

    Double_t sqrtArg(re*re-reh*reh*sintw*sintw);
    if(sqrtArg<0) {
      // No solution possible
      sourceLon = -9999;
      sourceLat = -9999;
      desiredAlt = -9999;
      return 0;
    }
    Double_t L=reh*costw - TMath::Sqrt(sqrtArg);
    Double_t sinThetaL=L*sintw/re;
    Double_t sourceTheta=TMath::ASin(sinThetaL);

    TVector3 sourcePos(0, 0, re);
    sourcePos.RotateY(sourceTheta);
    sourcePos.RotateZ(-1*tempPhiWave);

    sourcePos.RotateY(fBalloonTheta);
    sourcePos.RotateZ(fBalloonPhi);

    Double_t sourceVec[3];
    sourcePos.GetXYZ(sourceVec);

    Double_t sourceAlt;
    geom->getLatLonAltFromCartesian(sourceVec,sourceLat,sourceLon,sourceAlt);
    nextRe=geom->getDistanceToCentreOfEarth(sourceLat);

  } while(TMath::Abs(nextRe-re)>1);

  return 1;
}



int UsefulAdu5Pat::getSourceLonAndLatAtAlt(Double_t phiWave, Double_t thetaWave, Double_t &sourceLon, Double_t &sourceLat,Double_t &sourceAltitude) const
{
  // updateCartesianBalloonInfo();

  // if(fPhiWave!=phiWave) fPhiWave=phiWave;
  // if(fThetaWave!=thetaWave) fThetaWave=thetaWave;

  Double_t tempPhiWave=phiWave;
  Double_t tempThetaWave=thetaWave;

  if(heading>=0 && heading<=360) {
    accountForPitchAndRollInPhiWaveThetaWave(tempPhiWave, tempThetaWave);
  }
  else{
    std::cerr << "Error in " << __PRETTY_FUNCTION__ << ": Bad heading = " << heading << std::endl;
    sourceLon = -9999;
    sourceLat = -9999;
    sourceAltitude = -9999;
    return -1;
  }

  AnitaGeomTool* geom = AnitaGeomTool::Instance();
  
  Double_t reBalloon = geom->getDistanceToCentreOfEarth(latitude);		// reBalloon, radius of earth where balloon is
  Double_t chosenAlt = RampdemReader::SurfaceAboveGeoid(longitude,latitude);	// ChosenAlt, the ice surface below balloon
  Double_t re = reBalloon+chosenAlt;						// re, radius of earth at balloon plus ice height above geoid
  Double_t nextRe = re;								// nextRe, radius of geoid+ice at next iteration
  Double_t reh = reBalloon+altitude;						// reh, radius of earth at balloon plus balloon altitude

  if(fDebug){
    std::cout << "balloon: lat = " << latitude << ", lon = " << longitude << ", re = " << reBalloon+chosenAlt << ", surface elevation = " << chosenAlt << std::endl;
    // std::cout << "phiWave = " << phiWave << ", thetaWave = " << thetaWave << ", fPhiWave = " << fPhiWave << ", fThetaWave = " << fThetaWave << std::endl;
    std::cout << "heading = " << heading << ", pitch = " << pitch << ", roll = " << roll << std::endl;
  }

  int inTheLoop = 0;
  Double_t lastButOneRe   = 0;	// in case we get stuck in the while loop - this is a bad method and needs to be improved!
  Double_t lastButTwoRe   = 0;	// in case we get stuck in the while loop - this is a bad method and needs to be improved!
  Double_t lastButThreeRe = 0;	// in case we get stuck in the while loop - this is a bad method and needs to be improved!
  Double_t lastButFourRe  = 0;	// in case we get stuck in the while loop - this is a bad method and needs to be improved!

  const int maxLoopIterations = 50; // maybe make this configurable?
  const double deltaReCloseEnough = 1.0;  
  Bool_t ocillatingAroundSolution = false;  
  Bool_t tooManyLoops = false;

  const Double_t sintw = TMath::Sin(tempThetaWave);
  const Double_t costw = TMath::Cos(tempThetaWave);
  
  do {
    /**
     * Okay so effectively what we do here is switch to cartesian coords
     * with the balloon at RE_balloon + altitude and then try to find where the line at thetaWave cuts the Earth's surface.
     * This is iterative because the Earth is cruel and isn't flat, and I couldn't be bothered to work out how to do it more elegantly.
     */
    if(inTheLoop>3){lastButFourRe = lastButThreeRe;}
    if(inTheLoop>2){lastButThreeRe = lastButTwoRe;}
    if(inTheLoop>1){lastButTwoRe = lastButOneRe;}
    if(inTheLoop>0){lastButOneRe = re;}
    re = nextRe;
    inTheLoop++;

    Double_t sqrtArg = re*re - reh*reh*sintw*sintw;
    // clearly we're going for some kind of trigonometric relation here
    // re = radius of earth under balloon + chosenAlt
    // reh = radius of earth under balloon + balloon altitude
    // sintw = sin(thetaWave)

    if(sqrtArg<0) {
      if(fDebug){
	std::cerr << "Error in " << __PRETTY_FUNCTION__ << ", no solution possible! Returning 0\n";
	std::cerr << "Balloon alt = " << altitude << "\n";	
	std::cerr << "re = " << re << ",  reh = " << reh << ", sintw = " << sintw << ", tempThetaWave = " << tempThetaWave << "\n";
        std::cerr << "re*re = " << re*re << ",  reh*reh*sintw*sintw = " << reh*reh*sintw*sintw << "\n";
        std::cerr << "sqrtArg = " << sqrtArg << std::endl;
      }
      return 0;
    }

    Double_t L = reh*costw - TMath::Sqrt(sqrtArg);
    Double_t sinThetaL = L*sintw/re;
    Double_t sourceTheta = TMath::ASin(sinThetaL);

    TVector3 sourcePos(0, 0, re); // Put the source at the south pole at height, re

    /* Move the source to near the balloon? */
    sourcePos.RotateY(sourceTheta);
    sourcePos.RotateZ(-1*tempPhiWave);
    sourcePos.RotateY(fBalloonTheta);
    sourcePos.RotateZ(fBalloonPhi);

    Double_t sourceVec[3];
    sourcePos.GetXYZ(sourceVec);
    Double_t sourceAlt;
    geom->getLatLonAltFromCartesian(sourceVec,sourceLat,sourceLon,sourceAlt);
    chosenAlt = RampdemReader::SurfaceAboveGeoid(sourceLon,sourceLat);
    nextRe = geom->getDistanceToCentreOfEarth(sourceLat) + chosenAlt;

    if(fDebug){
      std::cerr << "Debug in " << __PRETTY_FUNCTION__ << ": inTheLoop = " << inTheLoop
		<< ", sourceLon = " << sourceLon << ", sourceLat = " << sourceLat << ", sourceAlt = " << sourceAlt
		<< ", re = " << re << ", lastButOneRe = " << lastButOneRe << ", lastButTwoRe = " << lastButTwoRe
		<< ", lastButThreeRe = " << lastButThreeRe << ", lastButFourRe = " << lastButFourRe
		<< ", TMath::Abs(lastButOneRe-nextRe) = " << TMath::Abs(lastButOneRe-nextRe) << "\n";
    }

    const int numLoopChecks = 4;
    const double deltaRes[numLoopChecks] = {TMath::Abs(lastButOneRe-nextRe),
					    TMath::Abs(lastButTwoRe-nextRe),
					    TMath::Abs(lastButThreeRe-nextRe),
					    TMath::Abs(lastButFourRe-nextRe)};
    int smallestRe = TMath::LocMin(numLoopChecks, deltaRes);
    if(deltaRes[smallestRe] <= deltaReCloseEnough){
      if(fDebug){
	const double Res[numLoopChecks]    = {lastButOneRe,    lastButTwoRe,   lastButThreeRe,   lastButFourRe};
	const char* whichRe[numLoopChecks] = {"lastButOneRe", "lastButTwoRe", "lastButThreeRe", "lastButFourRe"};
	std::cerr << "Breaking out due to small enough difference over several loop iterations, nextRe = "
		  << nextRe << ", " << whichRe[smallestRe] << " = " << Res[smallestRe]
		  << ", abs(" << whichRe[smallestRe] << " - nextRe) = " << deltaRes[smallestRe] << std::endl;
      }
      ocillatingAroundSolution = true;
      break;
    }

    
    if(inTheLoop>maxLoopIterations){
      if(fDebug){
	std::cerr << "Debug in " << __PRETTY_FUNCTION__ << ": Breaking out due to too many loop iterations! inTheLoop > " << maxLoopIterations << std::endl;	
      }
      tooManyLoops=1;
      break;
    }

  } while(TMath::Abs(nextRe-re)>deltaReCloseEnough);

  sourceAltitude = chosenAlt;

  if(ocillatingAroundSolution){
    // It converged although the loop oscillated around the true solution a bit,
    // so let's trust this one and not overwrite the source lon/lat/alt
    return 2;
  }
  else if(tooManyLoops){
    // We couldn't get to the true answer in quite a lot of iterations, so we probably don't trust the answer in this case.
    sourceLon = -9999;
    sourceLat = -9999;
    sourceAltitude = -9999;
    return 3;
  }
  else if(TMath::Nint(chosenAlt)==-9999){
    // I think this means the source altitude is off the edge of the Rampdem map, so it's probably nonsense
    sourceLon = -9999;
    sourceLat = -9999;
    sourceAltitude = -9999;
    return 4;
  }
  else{
    // It converged on the first attempt
    return 1;
  }
}



void UsefulAdu5Pat::getThetaAndPhiWaveWaisDivide(Double_t &thetaWave, Double_t &phiWave)
{  
  return getThetaAndPhiWave(AnitaLocations::getWaisLongitude(), AnitaLocations::getWaisLatitude(),AnitaLocations::getWaisAltitude(),thetaWave,phiWave);
}

void UsefulAdu5Pat::getThetaAndPhiWaveLDB(Double_t &thetaWave, Double_t &phiWave)
{
  return getThetaAndPhiWave(AnitaLocations::LONGITUDE_LDB,AnitaLocations::LATITUDE_LDB,AnitaLocations::ALTITUDE_LDB,thetaWave,phiWave);
}



void UsefulAdu5Pat::getThetaAndPhiWave(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt, Double_t &thetaWave, Double_t &phiWave) {
  // updateCartesianBalloonInfo();
  
  Double_t pSource[3]={0};
  AnitaGeomTool * geom = AnitaGeomTool::Instance();
  geom->getCartesianCoords(TMath::Abs(sourceLat),sourceLon,sourceAlt,pSource);
  TVector3 sourcePos(pSource);//[0],pSource[1],pSource[2]);

  //Rotate such that balloon is at 0,0,fBalloonHeight
  sourcePos.RotateZ(-1*fBalloonPhi);
  sourcePos.RotateY(-1*fBalloonTheta);

  //Now find thetaWave and phiWave
  thetaWave = TMath::ATan((fBalloonHeight-sourcePos.Z())/TMath::Sqrt(sourcePos.X()*sourcePos.X() + sourcePos.Y()*sourcePos.Y()));

  // phiWave is just atan(yp/xp)
  // It only looks confusing to make sure I get the sign and 0-360 convention
  phiWave = 0;
  if(sourcePos.X()==0) {
    phiWave = TMath::PiOver2();
    if(sourcePos.Y()<0)
      phiWave += TMath::Pi();
  }
  else if(sourcePos.X()<0) {
    phiWave = TMath::Pi()+TMath::ATan(sourcePos.Y()/sourcePos.X());
  }
  else {
    phiWave = TMath::ATan(sourcePos.Y()/sourcePos.X());
    if(sourcePos.Y()<0) {
      phiWave += TMath::TwoPi();
    }
  }

  TVector3 rollAxis = geom->fRollRotationAxis;
  TVector3 pitchAxis = geom->fPitchRotationAxis;

  Double_t tempThetaWave=TMath::PiOver2()-thetaWave;

  if(heading>=0 && heading<=360) {
    TVector3 arbDir; //
    arbDir.SetMagThetaPhi(1,tempThetaWave,-1*phiWave);

    rollAxis.Rotate(heading*TMath::DegToRad(),geom->fHeadingRotationAxis);
    pitchAxis.Rotate(heading*TMath::DegToRad(),geom->fHeadingRotationAxis);
    rollAxis.Rotate(pitch*TMath::DegToRad(),pitchAxis);

    arbDir.Rotate(-1.*roll*TMath::DegToRad(),rollAxis);
    arbDir.Rotate(-1*pitch*TMath::DegToRad(),pitchAxis);
    arbDir.Rotate(-1*heading*TMath::DegToRad(),geom->fHeadingRotationAxis);

    phiWave = arbDir.Phi();
    phiWave *= -1;
    if(phiWave>TMath::TwoPi()) {
      phiWave-=TMath::TwoPi();
    }
    if(phiWave<0) {
      phiWave+=TMath::TwoPi();
    }
    thetaWave=TMath::PiOver2()-arbDir.Theta();
  }

  // fPhiWave = phiWave;
  // fThetaWave = thetaWave;
  // fSourceLongitude = sourceLon;
  // fSourceLatitude = sourceLat;
  // fSourceAltitude = sourceAlt;
}

Double_t UsefulAdu5Pat::getGroupDelay(Double_t phiToAntBoresight)
{
  Double_t phiDeg=phiToAntBoresight*TMath::RadToDeg();
  Double_t delayTime=(phiDeg*phiDeg*phiDeg*phiDeg)*1.45676e-8;
  delayTime-=(phiDeg*phiDeg)*5.01452e-6;
  return delayTime;
}





void UsefulAdu5Pat::getThetaWaveAtBase(Double_t baseLon, Double_t baseLat, Double_t baseAlt, Double_t &thetaWave) {
  updateCartesianBalloonInfo();
  
  //gets the theta from a base to the balloon - can check if we are beyond horizon

  Double_t pSource[3]={0};
  AnitaGeomTool::Instance()->getCartesianCoords(TMath::Abs(baseLat),baseLon,baseAlt,pSource);
  TVector3 sourcePos(pSource);

  //make copy of balloon position and rotate it such that base is at 0,0,baseAlt
  TVector3 rotatedBalloonPos = fBalloonPos;
  rotatedBalloonPos.RotateZ(-1*sourcePos.Phi());
  rotatedBalloonPos.RotateY(-1*sourcePos.Theta());

  TVector3 rotatedBase = sourcePos;
  rotatedBase.RotateZ(-1*sourcePos.Phi());
  rotatedBase.RotateY(-1*sourcePos.Theta());

  //Now find thetaWave
  thetaWave=TMath::ATan((rotatedBase.Z()-rotatedBalloonPos.Z())/TMath::Sqrt(rotatedBalloonPos.X()*rotatedBalloonPos.X() + rotatedBalloonPos.Y()*rotatedBalloonPos.Y()));

}


Double_t UsefulAdu5Pat::getDeltaTExpected(Int_t ant1, Int_t ant2,Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt)
{
  updateCartesianBalloonInfo();

  Double_t tempTheta,tempPhi;
  if(fSourceAltitude!=sourceAlt || fSourceLongitude!=sourceLon || fSourceLatitude!=sourceLat){
    getThetaAndPhiWave(sourceLon,sourceLat,sourceAlt,tempTheta,tempPhi);    
  }
  

  AnitaGeomTool *geom = AnitaGeomTool::Instance();
  //Now fThetaWave and fPhiWave should be correctly set.
  Double_t phi1=geom->getAntPhiPositionRelToAftFore(ant1);
  Double_t r1=geom->getAntR(ant1);
  Double_t z1=geom->getAntZ(ant1);

  Double_t phi2=geom->getAntPhiPositionRelToAftFore(ant2);
  Double_t r2=geom->getAntR(ant2);
  Double_t z2=geom->getAntZ(ant2);

  Double_t tanThetaW=TMath::Tan(fThetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(fPhiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(fPhiWave-phi2);

  Double_t geomTime= 1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
  if(fIncludeGroupDelay) {
    Double_t phi1Diff=geom->getPhiDiff(fPhiWave,phi1);
    Double_t delay1=getGroupDelay(phi1Diff);
    Double_t phi2Diff=geom->getPhiDiff(fPhiWave,phi2);
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

  AnitaGeomTool *geom = AnitaGeomTool::Instance();
  Double_t phi1=geom->getAntPhiPositionRelToAftFore(ant1);
  Double_t r1=geom->getAntR(ant1);
  Double_t z1=geom->getAntZ(ant1);

  Double_t phi2=geom->getAntPhiPositionRelToAftFore(ant2);
  Double_t r2=geom->getAntR(ant2);
  Double_t z2=geom->getAntZ(ant2);

  Double_t tanThetaW=TMath::Tan(fThetaWave);
  Double_t part1=z1*tanThetaW - r1 * TMath::Cos(fPhiWave-phi1);
  Double_t part2=z2*tanThetaW - r2 * TMath::Cos(fPhiWave-phi2);

  return  1e9*((TMath::Cos(fThetaWave) * (part1 - part2))/C_LIGHT);    //returns time in ns
}


//does not take into account anita's position/heading...
Double_t UsefulAdu5Pat::getDeltaTExpected(Int_t ant1,Int_t ant2,Double_t cosPhi,Double_t sinPhi,Double_t cosTheta,Double_t sinTheta)
{  
  Double_t x1,y1,z1;
  AnitaGeomTool *geom = AnitaGeomTool::Instance();
  geom->getAntXYZ(ant1,x1,y1,z1);

  Double_t x2,y2,z2;
  geom->getAntXYZ(ant2,x2,y2,z2);

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
  updateCartesianBalloonInfo();

  Double_t tempTheta,tempPhi;
  if(fSourceAltitude!=sourceAlt || fSourceLongitude!=sourceLon || fSourceLatitude!=sourceLat)
    getThetaAndPhiWave(sourceLon,sourceLat,sourceAlt,tempTheta,tempPhi);

  //Now fThetaWave and fPhiWave should be correctly set.
  AnitaGeomTool *geom = AnitaGeomTool::Instance();
  Double_t phi1=geom->getAntPhiPositionRelToAftFore(ant1)+deltaPhi[ant1];
  Double_t r1=geom->getAntR(ant1)+deltaR[ant1];
  Double_t z1=geom->getAntZ(ant1)+deltaZ[ant1];

  Double_t phi2=geom->getAntPhiPositionRelToAftFore(ant2)+deltaPhi[ant2];
  Double_t r2=geom->getAntR(ant2)+deltaR[ant2];
  Double_t z2=geom->getAntZ(ant2)+deltaZ[ant2];

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
  updateCartesianBalloonInfo();
  
  Double_t tempTheta,tempPhi;
  if(fSourceAltitude!=sourceAlt || fSourceLongitude!=sourceLon || fSourceLatitude!=sourceLat){
    getThetaAndPhiWave(sourceLon,sourceLat,sourceAlt,tempTheta,tempPhi);
  }

  //Now fThetaWave and fPhiWave should be correctly set.
  AnitaGeomTool *geom = AnitaGeomTool::Instance();
  Double_t phi1=geom->getAntPhiPositionRelToAftFore(ant1)+deltaPhi[ant1];
  Double_t r1=geom->getAntR(ant1)+deltaR[ant1];
  Double_t z1=geom->getAntZ(ant1)+deltaZ[ant1];

  Double_t phi2=geom->getAntPhiPositionRelToAftFore(ant2)+deltaPhi[ant2];
  Double_t r2=geom->getAntR(ant2)+deltaR[ant2];
  Double_t z2=geom->getAntZ(ant2)+deltaZ[ant2];

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
  return getTriggerTimeNsFromSource(AnitaLocations::getWaisLatitude(),
				    AnitaLocations::getWaisLongitude(),
				    AnitaLocations::getWaisAltitude());
}

UInt_t UsefulAdu5Pat::getSipleTriggerTimeNs()
{
  return getTriggerTimeNsFromSource(AnitaLocations::LATITUDE_SIPLE,
				    AnitaLocations::LONGITUDE_SIPLE,
				    AnitaLocations::ALTITUDE_SIPLE);
}

UInt_t UsefulAdu5Pat::getLDBTriggerTimeNs()
{
  return getTriggerTimeNsFromSource(AnitaLocations::LATITUDE_LDB,
				    AnitaLocations::LONGITUDE_LDB,
				    AnitaLocations::ALTITUDE_LDB);
}



Double_t UsefulAdu5Pat::getDistanceFromSource(Double_t sourceLat, Double_t sourceLong, Double_t sourceAlt) const 
{
  // updateCartesianBalloonInfo();
  // static Double_t pTaylor[3]={0};
  Double_t pTaylor[3]={0};  
  AnitaGeomTool::Instance()->getCartesianCoords(TMath::Abs(sourceLat),
						sourceLong,
						sourceAlt,
						pTaylor);
  Double_t s2 = ((pTaylor[0]-fBalloonCoords[0])*(pTaylor[0]-fBalloonCoords[0]) +
		 (pTaylor[1]-fBalloonCoords[1])*(pTaylor[1]-fBalloonCoords[1]) +
		 (pTaylor[2]-fBalloonCoords[2])*(pTaylor[2]-fBalloonCoords[2]));
  Double_t distanceToFly=TMath::Sqrt(s2);
  return distanceToFly;
}





UInt_t UsefulAdu5Pat::getTriggerTimeNsFromSource(Double_t sourceLat, Double_t sourceLong, Double_t sourceAlt)
{
  updateCartesianBalloonInfo();  
  Double_t distanceToFly = getDistanceFromSource(sourceLat, sourceLong, sourceAlt);
  Double_t timeOfFlight=distanceToFly/C_LIGHT;
  timeOfFlight*=1e9;
  //   Double_t expTime=timeOfFlight-40e3;
  Double_t expTime=timeOfFlight;//-40e3;
  UInt_t expTrigTime=(UInt_t)expTime;

  //   std::cout << distanceToFly << "\t" << timeOfFlight << "\n";
  return expTrigTime;


}



Double_t UsefulAdu5Pat::getAngleBetweenPayloadAndSource(Double_t sourceLon, Double_t sourceLat, Double_t sourceAlt){ //ACG additional function
  updateCartesianBalloonInfo();
  AnitaGeomTool * geom = AnitaGeomTool::Instance();
  Double_t thetaBalloon=geom->getThetaFromLat(TMath::Abs(latitude));
  Double_t phiBalloon=geom->getPhiFromLon(longitude);
  Double_t balloonHeight=geom->getGeoid(thetaBalloon)+altitude;

  Double_t thetaSource=geom->getThetaFromLat(TMath::Abs(sourceLat));
  Double_t phiSource=geom->getPhiFromLon(sourceLon);
  Double_t radiusSource=geom->getGeoid(thetaSource)+sourceAlt;

  //Get vector from Earth's centre to source
  TVector3 sourcePos;
  sourcePos.SetX(radiusSource*TMath::Sin(thetaSource)*TMath::Cos(phiSource));
  sourcePos.SetY(radiusSource*TMath::Sin(thetaSource)*TMath::Sin(phiSource));
  sourcePos.SetZ(radiusSource*TMath::Cos(thetaSource));

  //Rotate such that balloon is at 0,0,balloonHeight

  sourcePos.RotateZ(-1*phiBalloon);
  sourcePos.RotateY(-1*thetaBalloon);

  Double_t dotProduct=balloonHeight*sourcePos.Z();//+0*X+0*Y
  Double_t magSource=TMath::Sqrt(sourcePos.X()*sourcePos.X() + sourcePos.Y()*sourcePos.Y()+sourcePos.Z()*sourcePos.Z());
  Double_t magBalloon=balloonHeight;

  Double_t phiEarthCenter=TMath::ACos(TMath::Abs(dotProduct/(magSource*magBalloon)));

  return phiEarthCenter;

}


//---------------------------------------------------------------------------------------------------------
/**
 * @brief Calculates the azimuth position of the sun if ANITA was facing directly north.
 *
 * @return azimuth (Degrees) relative to heading=0.
 *
 * Uses the timeOfDay variable from inside the Adu5Pat to get the sun's azimuth position relative to north.
 * i.e. returns -1*heading of the sun, as payload phi increases anti-clockwise and heading increases clockwise.
 * The more useful quantity is probably returned by UsefulAdu5Pat::getAzimuthOfSun().
 */
Double_t UsefulAdu5Pat::getAzimuthOfSunRelativeToNorth(){
  updateCartesianBalloonInfo();
  //work out the angle to the Sun for heading 0
  TTimeStamp theTime = this->realTime;

  UInt_t hour = 0;
  UInt_t min = 0;
  UInt_t sec = 0;
  theTime.GetTime(true,0,&hour,&min,&sec);
  Int_t timeOfDay = hour*3600 + min*60 + sec;


  // this would make 12pm at 180 degrees to sun - so we have to subtract 180 from the result
  Int_t secondsPerDay = (24*60*60);
  Double_t zeroLonAngleToSun = 360*Double_t(timeOfDay)/secondsPerDay;


  zeroLonAngleToSun-=180.;
  Double_t thisLonAngleToSun = zeroLonAngleToSun+this->longitude;

  while(thisLonAngleToSun<-180) thisLonAngleToSun+=360.;
  while(thisLonAngleToSun>180) thisLonAngleToSun-=360.;

  return thisLonAngleToSun;
}





//---------------------------------------------------------------------------------------------------------
/**
 * Get's the sun's azimuth position in payload phi relative to adu5 aft-fore
 *
 * @return the azimuth position of the sun taking into account ANITA's heading.
 */
Double_t UsefulAdu5Pat::getAzimuthOfSun(){
  updateCartesianBalloonInfo();
  Double_t sunAngle = UsefulAdu5Pat::getAzimuthOfSunRelativeToNorth();
  Double_t sunAngleAtHeading = sunAngle + this->heading;

  while(sunAngleAtHeading<-180) sunAngleAtHeading+=360;
  while(sunAngleAtHeading>180)  sunAngleAtHeading-=360;

  return sunAngleAtHeading;

}




//---------------------------------------------------------------------------------------------------------
/**
 * @brief Quick way of getting difference between any azimuth phi and the sun's azimuth phi.
 *
 * @param phiAngle is the input azimuth in either degrees or radians
 * @param inputInDegrees is a boolian, which should be true is the phiAngle is in degrees or false if the phiAngle is in radians.
 *
 * @return the difference between the two angles in Degrees (regardless of whether the input is in radians or degrees).
 */
Double_t UsefulAdu5Pat::getDifferencePointingToSun(Double_t phiAngle, Bool_t inputInDegrees){
  updateCartesianBalloonInfo();

  if(phiAngle>9999 || phiAngle<-9999){
    std::cerr << "UsefulAdu5Pat::getDifferencePointingToSun() " << phiAngle << " ignore\n";
    return -9999.;
  }

  Double_t sunAngle = UsefulAdu5Pat::getAzimuthOfSun();

  Double_t angleFromSun=0;
  if(inputInDegrees==false){
    angleFromSun = sunAngle - phiAngle*TMath::RadToDeg();
  }
  else{
    angleFromSun = sunAngle - phiAngle;
  }

  while(angleFromSun<-180){
    angleFromSun+=360;
  }
  while(angleFromSun>=180){
    angleFromSun-=360;
  }

  //always returns degrees
  return angleFromSun;

}



/**
 * @brief Uses realTime, latitude, longitude to calculate the sun's position.
 *
 * @param phiDeg will contain the phi position in ANITA's coordinates after the function call
 * @param thetaDeg will contain the theta position in ANITA's coordinates after the function call
 *
 * The meat of this function was taken from here http://www.psa.es/sdg/sunpos.htm
 */
void UsefulAdu5Pat::getSunPosition(Double_t& phiDeg, Double_t& thetaDeg){
  updateCartesianBalloonInfo();

  // Declaration of some constants
  const Double_t dEarthMeanRadius = 6371.01;	// In km
  const Double_t dAstronomicalUnit =149597890;	// In km


  TTimeStamp timeStamp = realTime;
  UInt_t year, month, day;
  timeStamp.GetDate(true, 0, &year, &month, &day);
  UInt_t hour, min, sec;
  timeStamp.GetTime(true, 0, &hour, &min, &sec);

  // std::cerr << year << "-" << month << "-" << day << "\t" << hour << ":" << min << ":" << sec << std::endl;

  // need to convert from unsigned to signed ints otherwise this bastard breaks
  Int_t iDay = day;
  Int_t iMonth = month;
  Int_t iYear = year;

  // Main variables
  Double_t dAzimuth;
  Double_t dZenithAngle;
  Double_t dElapsedJulianDays;
  Double_t dDecimalHours;
  Double_t dEclipticLongitude;
  Double_t dEclipticObliquity;
  Double_t dRightAscension;
  Double_t dDeclination;

  // Auxiliary variables
  Double_t dY;
  Double_t dX;

  // Calculate difference in days between the current Julian Day
  // and JD 2451545.0, which is noon 1 January 2000 Universal Time
  {
    Double_t dJulianDate;
    long int liAux1;
    long int liAux2;
    // Calculate time of the day in UT decimal hours
    dDecimalHours = (double)hour + ((double)min + (double)sec / 60.0 ) / 60.0;
    // Calculate current Julian Day
    liAux1 =(iMonth-14)/12;
    liAux2=(1461*(iYear + 4800 + liAux1))/4 + (367*(iMonth
						    - 2-12*liAux1))/12- (3*((iYear + 4900
									     + liAux1)/100))/4+iDay-32075;
    dJulianDate=(double)(liAux2)-0.5+dDecimalHours/24.0;
    // Calculate difference between current Julian Day and JD 2451545.0
    dElapsedJulianDays = dJulianDate-2451545.0;
  }

  // Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
  // ecliptic in radians but without limiting the angle to be less than 2*Pi
  // (i.e., the result may be greater than 2*Pi)
  {
    Double_t dMeanLongitude;
    Double_t dMeanAnomaly;
    Double_t dOmega;
    dOmega=2.1429-0.0010394594*dElapsedJulianDays;
    dMeanLongitude = 4.8950630+ 0.017202791698*dElapsedJulianDays; // Radians
    dMeanAnomaly = 6.2400600+ 0.0172019699*dElapsedJulianDays;
    dEclipticLongitude = dMeanLongitude + 0.03341607*sin( dMeanAnomaly )
      + 0.00034894*sin( 2*dMeanAnomaly )-0.0001134
      -0.0000203*sin(dOmega);
    dEclipticObliquity = 0.4090928 - 6.2140e-9*dElapsedJulianDays
      +0.0000396*cos(dOmega);
  }

  // Calculate celestial coordinates ( right ascension and declination ) in radians
  // but without limiting the angle to be less than 2*Pi (i.e., the result may be
  // greater than 2*Pi)
  {
    Double_t dSin_EclipticLongitude;
    dSin_EclipticLongitude= sin( dEclipticLongitude );
    dY = cos( dEclipticObliquity ) * dSin_EclipticLongitude;
    dX = cos( dEclipticLongitude );
    dRightAscension = atan2( dY,dX );
    if( dRightAscension < 0.0 ) dRightAscension = dRightAscension + TMath::TwoPi();
    dDeclination = asin( sin( dEclipticObliquity )*dSin_EclipticLongitude );
  }

  // Calculate local coordinates ( azimuth and zenith angle ) in degrees
  {
    Double_t dGreenwichMeanSiderealTime;
    Double_t dLocalMeanSiderealTime;
    Double_t dLatitudeInRadians;
    Double_t dHourAngle;
    Double_t dCos_Latitude;
    Double_t dSin_Latitude;
    Double_t dCos_HourAngle;
    Double_t dParallax;
    dGreenwichMeanSiderealTime = 6.6974243242 +
      0.0657098283*dElapsedJulianDays
      + dDecimalHours;
    dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime*15
			      + longitude)*TMath::DegToRad();
    dHourAngle = dLocalMeanSiderealTime - dRightAscension;
    dLatitudeInRadians = latitude*TMath::DegToRad();
    dCos_Latitude = cos( dLatitudeInRadians );
    dSin_Latitude = sin( dLatitudeInRadians );
    dCos_HourAngle= cos( dHourAngle );
    dZenithAngle = (acos( dCos_Latitude*dCos_HourAngle
			  *cos(dDeclination) + sin( dDeclination )*dSin_Latitude));
    dY = -sin( dHourAngle );
    dX = tan( dDeclination )*dCos_Latitude - dSin_Latitude*dCos_HourAngle;
    dAzimuth = atan2( dY, dX );
    if ( dAzimuth < 0.0 )
      dAzimuth = dAzimuth + TMath::TwoPi();
    dAzimuth = dAzimuth/TMath::DegToRad();
    // Parallax Correction
    dParallax=(dEarthMeanRadius/dAstronomicalUnit)
      *sin(dZenithAngle);
    dZenithAngle=(dZenithAngle
		  + dParallax)/TMath::DegToRad();
  }

  // std::cerr << dAzimuth << "\t" << dZenithAngle << std::endl;

  // finally we need to convert to ANITA's coordinates
  // from my testing...
  // dAzimuth = pat->heading - recoPhiDeg;
  // therefore

  phiDeg = heading - dAzimuth;
  while(phiDeg < -180) phiDeg += 360;
  while(phiDeg >= 180) phiDeg -= 360;

  // and theta for these guys runs from
  thetaDeg = dZenithAngle - 90;

  // thetaDeg = dZenithAngle;
  // phiDeg = dAzimuth;

 }

int UsefulAdu5Pat::traceBackToContinent(Double_t phiWave, Double_t thetaWave,
                                        Double_t * lon_ptr, Double_t * lat_ptr, Double_t * alt_ptr,
                                        Double_t * adj_ptr, Double_t max_adjust, Int_t max_iter) const {
  // updateCartesianBalloonInfo();
  
  if(fDebug){
    std::cerr << "Debug info in " << __PRETTY_FUNCTION__ << ":\n";
  }
  
  
  if (thetaWave + max_adjust < 0){
    if(fDebug){
      std::cerr << "thetaWave + max_adjust < 0... no chance! returning 0!" << std::endl;
    }
    return 0; // no chance
  }

  Double_t lon,lat,alt;
  Double_t last_theta_tried =0;
  Double_t last_successful_theta = 0;
  Double_t last_failed_theta = 0;

  double last_good_lat = -9999; 
  double last_good_lon = -9990; 
  double last_good_alt = -9999; 
  Int_t iter = 0;

  while(iter++ < max_iter)
  {
    Double_t theta_try = last_theta_tried == 0 ? thetaWave
                       : last_theta_tried == thetaWave
                          ? thetaWave+max_adjust
                          : (last_successful_theta + last_failed_theta)/2;

    if(fDebug){
      std::cerr << "theta_try  = " << theta_try << "\n";
    }
    
    last_theta_tried = theta_try;

    int success = getSourceLonAndLatAtAlt(phiWave, theta_try, lon, lat, alt); 
    if(fDebug){
      std::cerr << "getSourceLonAndLatAtAlt: success  = " << success << "\n";
    }
    
    if (success == 1 || success == 2)
    {
      if (iter == 1) //this was the first try and it was a success
      {
        if (lon_ptr) *lon_ptr = lon; 
        if (lat_ptr) *lat_ptr = lat; 
        if (alt_ptr) *alt_ptr = alt; 
        if (adj_ptr) *adj_ptr = 0;
	if(fDebug){
	  std::cerr << "success on first try! Returning 1" << std::endl;
	}
        return 1; 
      }


      last_good_lat = lat; 
      last_good_lon = lon; 
      last_good_alt = alt; 
      last_successful_theta = theta_try; 

    }
    else 
    {
      if (theta_try == thetaWave + max_adjust)
      {
	if(fDebug){	
	  std::cerr << "theta_try = thetaWave + max_adjust = " << theta_try << "... it's hopeless! returning 0" << std::endl;
	}
        //it's hopeless. 
        return 0; 
      }

      last_failed_theta = theta_try; 
    }
  }

  if (lat_ptr) *lat_ptr = last_good_lat;
  if (lon_ptr) *lon_ptr = last_good_lon;
  if (alt_ptr) *alt_ptr = last_good_alt;
  if (adj_ptr) *adj_ptr = last_successful_theta - thetaWave;

  if(fDebug){
    std::cerr << "Got here! Returning 2 " << std::endl;
  }

  return 2;
}


Double_t UsefulAdu5Pat::getReflectionAngle(Double_t plAlt, Double_t el, Double_t imAlt) {
//float reflectionAngle(float h, float el) {
  Double_t result = -9999;
  const Double_t EARTH_POLAR_RADIUS = 6378000;
  if (el < -5.0) {
    Double_t theta = M_PI / 2.0 - (el * M_PI / 180.0);
    Double_t er = EARTH_POLAR_RADIUS;
    Double_t plH = er + plAlt;
    Double_t costr = cos(theta);
    Double_t d = -plH*costr - sqrt(plH*plH*costr*costr - 2*er*(plAlt-imAlt) - plAlt*plAlt + imAlt*imAlt);
    // alpha = angle over earth curvature from payload to image; < 3deg so don't bother w/asin
    Double_t alpha = d * sin(theta) / (er+imAlt);  
    //Double_t alpha = asin(d * sin(theta) / (er+imAlt));
    Double_t td = M_PI + 2*alpha - theta;
    result = (90.0 - 180.0*td/M_PI);
  }
  return result;
}

Double_t UsefulAdu5Pat::getReflectionAngle(Double_t el, Double_t imAlt) {
  return getReflectionAngle(this->altitude, el, imAlt);
}


int UsefulAdu5Pat::astronomicalCoordinates(Double_t phiWave, Double_t thetaWave, Double_t * RA_ptr, Double_t * dec_ptr, Double_t * l_ptr, Double_t * b_ptr) const
{
#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
  std::cerr << "Your  version of ROOT does not support all the features of eventReaderRoot\n";
  return 0;
#else


  // I used coordinate conversions based on Az going north
  double Az = phiWave - heading;
  double el = -thetaWave;


  //to radians
  Az *= M_PI/180;
  el *= M_PI/180;
  double lat = latitude * M_PI/180;


  //declination
  double dec = asin( sin(lat) * sin(el) + cos(lat) *cos(el)  *cos(Az));

  //hour angle
  double h  = atan2 ( sin(Az) , -cos(Az)* sin(lat)  + tan(el) * cos(lat));

//  printf("%f\n", h * 12 / M_PI);
  TTimeStamp ts ( (time_t) realTime, (timeOfDay % 1000)  * 1e6);

  //should this be AsLAST or AsLMST? Need to find a real astronomer
  // Also, the UT1 offset seems like a royal pain
  double lst = ts.AsLMST(longitude);

 // printf("lst: %f\n", lst);

  lst*= ( M_PI / 12);  //hour -> radians
  double RA = lst -h;

  if (RA_ptr) *RA_ptr = RA * 12 / M_PI;
  if (dec_ptr) *dec_ptr = dec* 180 / M_PI;

  if (!l_ptr && ! b_ptr) return 0;

  //we gotta precess to B1950

  double b1950 = 2433282.4235;  // julian day of B1950
  double jd = ts.AsJulianDate();
//  printf("%f\n",jd);
  double T = ( jd - 2451545) / (36525);  // offset from J2000 in centidays... ?
  double t = (b1950 - jd) / 36525; //offset from B1950... in centidays?
  double sectodeg = 1./3600 ;

  //tons of magic numbers! stolen from DMTPC code. Sorry Asher.
  double zeta = (2306.2181*sectodeg+1.39656*sectodeg*T-0.000139*sectodeg*T*T)*t+
               (0.30188*sectodeg-0.000344*sectodeg*T)*t*t+
               0.017998*sectodeg*t*t*t;
  double z = (2306.2181*sectodeg+1.39656*sectodeg*T-0.000139*sectodeg*T*T)*t+
            (1.09468*sectodeg+0.000066*sectodeg*T)*t*t+
             0.018203*sectodeg*t*t*t;
  double theta = (2004.3109*sectodeg-0.85330*sectodeg*T-0.000217*sectodeg*T*T)*t-
                  (0.42665*sectodeg+0.000217*sectodeg*T)*t*t+0.041833*sectodeg*t*t*t;
  zeta *=M_PI/180;
  z *=M_PI/180;
  theta *=M_PI/180;


  double A = cos(dec)*sin(RA+zeta);
  double B = cos(theta)*cos(dec)*cos(RA+zeta)-sin(theta)*sin(dec);
  double C = sin(theta)*cos(dec)*cos(RA+zeta)+cos(theta)*sin(dec);

  double ra1950 = (atan2(A,B)+z);
  double dec1950 = asin(C);


  if (l_ptr)
  {
    *l_ptr = 303 - 180 / M_PI *
             atan2(   sin(192.25 * M_PI/180 - ra1950),
                      cos (192.25 * M_PI/180  - ra1950) * sin(27.5 * M_PI/180) - tan(dec1950) * cos (27.4 * M_PI/180)
                   );

    while (*l_ptr > 180) *l_ptr -= 360 ;
    while (*l_ptr < -180) *l_ptr += 360 ;


  }


  if (b_ptr)
  {
    *b_ptr =  180 / M_PI * asin(
                                 sin(dec1950) * sin(27.4 * M_PI/180) +
                                 cos(dec1950) * cos(27.4 * M_PI/180) * cos (192.25 * M_PI/180 - ra1950)
                                ) ;

  }

  return 0;
#endif

}

int UsefulAdu5Pat::fromRADec(Double_t RA, Double_t dec, Double_t * phi, Double_t *theta) const
{
#if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
  std::cerr << "Your  version of ROOT does not support all the features of eventReaderRoot\n";
  return 0;
#else

  TTimeStamp ts ( (time_t) realTime, (timeOfDay % 1000)  * 1e6);

  //should this be AsLAST or AsLMST? Need to find a real astronomer
  // Also, the UT1 offset seems like a royal pain
  double lst = ts.AsLMST(longitude);
  double h = lst - RA; //hour angle 
  h*= ( M_PI / 12);  //hour -> radians
  dec *= (M_PI/180);  //deg to radians

  double lat = latitude * (M_PI /180); 
  double Az = atan2( cos(dec) * sin(h), -cos(h) * sin(lat)* cos(dec) + sin(dec) * cos(lat)); 
  double el = asin(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(h)); 

  if (phi) *phi = heading + Az * (180/M_PI); 
  if (theta) *theta = - el * (180/M_PI); 
  return 0; 
#endif

}
