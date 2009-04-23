#!/bin/sh
cd /Users/simonbevan/ANITA/eventCorrelator/macros
echo "starting optimisation"
root -l -b -q 'runAntOpt.C(0,1)' 
cat /Users/simonbevan/ANITA/eventCorrelator/macros/optFile.txt >> /Users/simonbevan/ANITA/eventCorrelator/macros/optFileAll.txt 
#root -l -b -q 'runAntOpt.C(1,2)' >> dumpOut.txt
#cat /Users/simonbevan/ANITA/eventCorrelator/macros/optFile.txt >> /Users/simonbevan/ANITA/eventCorrelator/macros/optFileAll.txt 
#root -l -b -q 'runAntOpt.C(2,3)' >> dumpOut.txt
#cat /Users/simonbevan/ANITA/eventCorrelator/macros/optFile.txt >> /Users/simonbevan/ANITA/eventCorrelator/macros/optFileAll.txt 