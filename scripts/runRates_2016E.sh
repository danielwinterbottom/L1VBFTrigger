#/bin/bash
EventsPassedL1 --isData --input filelists/ZeroBias_temp.dat --outputFilename output/DoubleJet_90_30_Mj30j30_580/Run2016E_output.root --maxEvents -1 --runs 277069,277070,277071,277072,277073,277075,277076,277081,277086,277087,277093,277094,277096

EventsPassedL1 --input filelists/HToInvisible.dat --outputFilename output/DoubleJet_94_30_Mj30j30_620/HToInvisible.root --maxEvents -1 --cutsInputLine 4