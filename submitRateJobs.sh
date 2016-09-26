for i in 1 2 3 4
do
  #python scripts/makejobs.py --outputFolder jobs$i --fileList filelists/ZeroBias_RunBAndE.dat --splitting 30 --outputFilename output/PureRate$i --options "--runs 277069,277070,277071,277072,277073,277076,277087,277094,277096,277112,277126,277127,277148 --isData --cutsInputLine $i"
  #./jobs$i/submitJobs.sh
  
  python scripts/makejobs.py --outputFolder jobsMuoHad$i --fileList filelists/HToTauTau.dat --splitting 30 --outputFilename output/MuoTau$i --options "--cutsInputLine $i --decayType MuoHad --makeHistograms"
  ./jobsMuoHad$i/submitJobs.sh

   python scripts/makejobs.py --outputFolder jobsHadHad$i --fileList filelists/HToTauTau.dat --splitting 30 --outputFilename output/TauTau$i --options "--cutsInputLine $i --decayType HadHad --makeHistograms"
  ./jobsHadHad$i/submitJobs.sh

  python scripts/makejobs.py --outputFolder jobsEleHad$i --fileList filelists/HToTauTau.dat --splitting 30 --outputFilename output/EleTau$i --options "--cutsInputLine $i --decayType EleHad --makeHistograms"
  ./jobsEleHad$i/submitJobs.sh
 
  python scripts/makejobs.py --outputFolder jobsInvis$i --fileList filelists/HToInvisible.dat --splitting 30 --outputFilename output/Invis$i --options "--cutsInputLine $i --makeHistograms"
  ./jobsInvis$i/submitJobs.sh
  
  python scripts/makejobs.py --outputFolder jobsBB$i --fileList filelists/HToBB.dat --splitting 30 --outputFilename output/BB$i --options "--cutsInputLine $i --makeHistograms"
  ./jobsBB$i/submitJobs.sh
done
