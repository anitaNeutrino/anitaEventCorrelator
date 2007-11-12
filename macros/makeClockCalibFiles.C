
void makeClockCalibFiles() {
  int firstSurf,secondSurf,chip;
  double value,error;
  
  double diffs[8][4];

  ifstream SlowClock("diffSurfFastClock.txt");
  int count=0;
  double sumDiffs[4]={0};
  while(SlowClock >> firstSurf >> secondSurf >> chip >> value >> error) {
    diffs[firstSurf][chip]=value;
    sumDiffs[chip]+=value;
    count++;
    //    cout << chip << "\t" << sumDiffs[chip] << endl;
  }
  cout << "SumDiff:\t";
  for(int chip=0;chip<4;chip++) {
    cout << sumDiffs[chip] << "\t";
  }
  cout << "\n";
  
  double diffToSurf0[8][4];
  double diffToSurf0_2[8][4];
  double diffToSurf0Ave[8][4];
  
  for(int chip=0;chip<4;chip++) {
    for(int surf=0;surf<8;surf++) {
      diffToSurf0[surf][chip]=0;
      for(int i=0;i<surf;i++) {
	diffToSurf0[surf][chip]+=diffs[i][chip];
      }
      diffToSurf0_2[surf][chip]=0;
      if(surf)
	diffToSurf0_2[surf][chip]=diffToSurf0[surf][chip]-sumDiffs[chip];    
      
      diffToSurf0Ave[surf][chip]=(diffToSurf0[surf][chip]+diffToSurf0_2[surf][chip])/2.;
    }
  }
  
  ofstream ClockCalib("newFastClockCalibNums.dat");
  ClockCalib << "#SURF\tChip\tCalib\n";
  for(int surf=0;surf<8;surf++) {
    for(int chip=0;chip<4;chip++) {
      ClockCalib << surf << "\t" << chip << "\t";
      if(surf<4)
	ClockCalib << diffToSurf0[surf][chip] << "\n";
      else
      	ClockCalib << diffToSurf0_2[surf][chip] << "\n";
      //      else if(surf<6) 
      //      ClockCalib << diffToSurf0Ave[surf][chip] << "\n";
    
      
      cout << surf << "\t" << chip << "\t" << diffToSurf0[surf][chip] << "\t" << diffToSurf0_2[surf][chip] 
	   << "\t" << diffToSurf0Ave[surf][chip] << "\n";
    }
  }



}
