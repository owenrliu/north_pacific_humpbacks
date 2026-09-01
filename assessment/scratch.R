# Run all this to seed the DoRun function
Code="B2F1";SensCase="BC";StrayBase=0;
IAmat=8;SA=0.96;SC=0.8;TimeLag=0;DensDepOpt=0;
SF=c(0,1,0,1,1,1); WithMirror=1;
UseKPrior=1; Kmax=60000;
YrSDevs=2000;
rvars="epsEnv";
SetNew=0;
SigmaDevS=6;SigmaDevF=0.01;Yr1=1970;Yr2=2023;StochSopt=1;
AddCV=T;MixWeights=c(1,1);CatchSer="B";AllPlots=F;DoBoot=F;
ByCatchFile="BycatchActual_2024_04_24.csv";
DataFileName= here('data',"Hump.dat");
FullDiag=T;
WghtTotal=1;Idirichlet=1;MaxN=100;seed=19101;
Init=NULL;Nage=11;DoBayes=F