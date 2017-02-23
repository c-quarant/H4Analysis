#include "TFile.h" 
#include "TTree.h" 
#include "TStyle.h" 
#include "TH1.h" 
#include "TH2.h" 
#include "TF1.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h" 
#include <iostream>

using namespace std;

TString AddSelection(TTree*, TString, TString, TString);

/*
Double_t fitfunc(Double_t *v, Double_t *par) {

  Double_t arg = 0;
  if (par[3] != 0) arg = v[0]/par[3];

  Double_t fitval;
  if (arg>1)
    fitval = par[0] + (par[1]*(1.- (1./(par[2]*log(arg)+1))));
  else
    fitval = 0;
  
  return fitval;
}
*/

void ComputeEfficiency() {

  gStyle->SetOptStat(0);

  // inputs - here each mcp is in 1st position (for 2corr it is not)
  //TFile* inFileBinp1 = new TFile("/cmsrm/pc28_2/crovelli/data/imcp/binps/newRecoWindowV4/btf2016_HVscanBinp1.root");
  //TFile* inFileBinp2 = new TFile("/cmsrm/pc28_2/crovelli/data/imcp/binps/newRecoWindowV4/btf2016_HVscanBinp4.root");   // in btf2016_HVscanBinp2.root binp2 was probably partly off
  //TFile* inFileBinp3 = new TFile("/cmsrm/pc28_2/crovelli/data/imcp/binps/newRecoWindowV4/btf2016_HVscanBinp3.root");
  //TFile* inFileBinp4 = new TFile("/cmsrm/pc28_2/crovelli/data/imcp/binps/newRecoWindowV4/btf2016_HVscanBinp4.root");
  TFile* inFileBinp1 = new TFile("/cmsrm/pc28_2/bmarzocc/btf2016_HVscanBinp1.root");
  TFile* inFileBinp2 = new TFile("/cmsrm/pc28_2/bmarzocc/btf2016_HVscanBinp4.root");   // in btf2016_HVscanBinp2.root binp2 was probably partly off
  TFile* inFileBinp3 = new TFile("/cmsrm/pc28_2/bmarzocc/btf2016_HVscanBinp3.root");
  TFile* inFileBinp4 = new TFile("/cmsrm/pc28_2/bmarzocc/btf2016_HVscanBinp4.root");
  //
  TFile* inFileTiming14 = new TFile("/cmsrm/pc28_2/crovelli/data/imcp/binps/newRecoWindowV4/btf2016_RU5_2378.root");
  TFile* inFileTiming2  = new TFile("/cmsrm/pc28_2/crovelli/data/imcp/binps/newRecoWindowV4/btf2016_RU6_2547.root");  
  TFile* inFileTiming3  = new TFile("/cmsrm/pc28_2/crovelli/data/imcp/binps/newRecoWindowV4/btf2016_RU4_merged_2363_2366.root");  
  
  TTree* h4binp1 = (TTree*)inFileBinp1->Get("h4");
  TTree* h4binp2 = (TTree*)inFileBinp2->Get("h4");
  TTree* h4binp3 = (TTree*)inFileBinp3->Get("h4");
  TTree* h4binp4 = (TTree*)inFileBinp4->Get("h4");
  //
  TTree* h4timing14 = (TTree*)inFileTiming14->Get("h4");
  TTree* h4timing2  = (TTree*)inFileTiming2->Get("h4");
  TTree* h4timing3  = (TTree*)inFileTiming3->Get("h4");


  // -------------------------------------
  // Efficiency runs
  //
  //
  // Nominal histos
  TH1F* num1 = new TH1F("num1","",2500,0.,2500.);
  TH1F* num2 = new TH1F("num2","",2500,0.,2500.);
  TH1F* num3 = new TH1F("num3","",2500,0.,2500.);
  TH1F* num4 = new TH1F("num4","",2500,0.,2500.);
  TH1F* den1 = new TH1F("den1","",2500,0.,2500.);
  TH1F* den2 = new TH1F("den2","",2500,0.,2500.);
  TH1F* den3 = new TH1F("den3","",2500,0.,2500.);
  TH1F* den4 = new TH1F("den4","",2500,0.,2500.);

  // To be quoted only - with timing infos as well
  TH1F* num1b = new TH1F("num1b","",2500,0.,2500.);
  TH1F* num2b = new TH1F("num2b","",2500,0.,2500.);
  TH1F* num3b = new TH1F("num3b","",2500,0.,2500.);
  TH1F* num4b = new TH1F("num4b","",2500,0.,2500.);
  TH1F* den1b = new TH1F("den1b","",2500,0.,2500.);
  TH1F* den2b = new TH1F("den2b","",2500,0.,2500.);
  TH1F* den3b = new TH1F("den3b","",2500,0.,2500.);
  TH1F* den4b = new TH1F("den4b","",2500,0.,2500.);



  // Nominal denominator: selection based on ref-mcp amplitudes and scintillators
  TString commonDenAmp = "adc_data>200 && adc_data<700 && n_hitsX>0 && n_hitsY>0 && n_hitsX<3 && n_hitsY<3 && amp_max[MiB2]>200 && fabs(time_max[MiB2])<150";
  //TString commonDenAmp = "adc_data<150 && n_hitsX==0 && n_hitsY==0";
  //TString commonDenAndRm2Amp = commonDenAmp + " && amp_max[Rm2]>200 && amp_max[Rm2]<1200 && fabs(time_max[Rm2])<150";
  TString commonDenAndRm2Amp = commonDenAmp + "";
  TString den1s = commonDenAmp;            // +"&& HVBINP1==2400";
  TString den2s = commonDenAndRm2Amp;      // +"&& HVBINP2==1400";
  TString den3s = commonDenAndRm2Amp;      // +"&& HVBINP3==1400";
  TString den4s = commonDenAndRm2Amp;      // +"&& HVBINP4==2400";
  h4binp1->Project("den1","HVBINP1",den1s);
  h4binp2->Project("den2","HVBINP2",den2s);
  h4binp3->Project("den3","HVBINP3",den3s);
  h4binp4->Project("den4","HVBINP4",den4s);

  // Compute mean timing distance between Mib2 and Rm2 (not for Binp1 runs, RM2 OFF) in the run at plateau to be added to denom - to be quoted only
  TH1F* tim2_rm2mib2 = new TH1F("tim2_rm2mib2","tim2_rm2mib2",4000,-20.,20.);
  TH1F* tim3_rm2mib2 = new TH1F("tim3_rm2mib2","tim3_rm2mib2",4000,-20.,20.);
  TH1F* tim4_rm2mib2 = new TH1F("tim4_rm2mib2","tim4_rm2mib2",4000,-20.,20.);
  h4binp2->Project("tim2_rm2mib2","time_max[Rm2]-time_max[MiB2]",commonDenAndRm2Amp+" && HVBINP2==1400");
  h4binp3->Project("tim3_rm2mib2","time_max[Rm2]-time_max[MiB2]",commonDenAndRm2Amp+" && HVBINP3==1400");
  h4binp4->Project("tim4_rm2mib2","time_max[Rm2]-time_max[MiB2]",commonDenAndRm2Amp+" && HVBINP4==2400");

  char mean2_rm2mib2[100];
  char mean3_rm2mib2[100];
  char mean4_rm2mib2[100];
  sprintf(mean2_rm2mib2,"%f",tim2_rm2mib2->GetMean());
  sprintf(mean3_rm2mib2,"%f",tim3_rm2mib2->GetMean());
  sprintf(mean4_rm2mib2,"%f",tim4_rm2mib2->GetMean());
  std::string sMean2_rm2mib2 = std::string(mean2_rm2mib2);
  std::string sMean3_rm2mib2 = std::string(mean3_rm2mib2);
  std::string sMean4_rm2mib2 = std::string(mean4_rm2mib2);


  TString denAndRm2AmpTim1 = commonDenAmp + " && HVBINP1==2400";

  TString denAndRm2AmpTim2;
  if(tim2_rm2mib2->GetMean() < 0.){
    sMean2_rm2mib2.erase(sMean2_rm2mib2.begin(),sMean2_rm2mib2.begin()+1);
    denAndRm2AmpTim2 = commonDenAndRm2Amp+std::string(" && HVBINP2==1400 && fabs(time_max[Rm2]-time_max[MiB2]+")+std::string(sMean2_rm2mib2)+std::string(")<1");
    } else { 
    denAndRm2AmpTim2 = commonDenAndRm2Amp+std::string(" && HVBINP2==1400 && fabs(time_max[Rm2]-time_max[MiB2]-")+std::string(sMean2_rm2mib2)+std::string(")<1");
  }

  TString denAndRm2AmpTim3;
  if(tim3_rm2mib2->GetMean() < 0.){
    sMean3_rm2mib2.erase(sMean3_rm2mib2.begin(),sMean3_rm2mib2.begin()+1);
    denAndRm2AmpTim3 = commonDenAndRm2Amp+std::string(" && HVBINP3==1400 && fabs(time_max[Rm2]-time_max[MiB2]+")+std::string(sMean3_rm2mib2)+std::string(")<1");
    } else { 
    denAndRm2AmpTim3 = commonDenAndRm2Amp+std::string(" && HVBINP3==1400 && fabs(time_max[Rm2]-time_max[MiB2]-")+std::string(sMean3_rm2mib2)+std::string(")<1");
  }

  TString denAndRm2AmpTim4;
  if(tim4_rm2mib2->GetMean() < 0.){
    sMean4_rm2mib2.erase(sMean4_rm2mib2.begin(),sMean4_rm2mib2.begin()+1);
    denAndRm2AmpTim4 = commonDenAndRm2Amp+std::string(" && HVBINP4==2400 && fabs(time_max[Rm2]-time_max[MiB2]+")+std::string(sMean4_rm2mib2)+std::string(")<1");
    } else { 
    denAndRm2AmpTim4 = commonDenAndRm2Amp+std::string(" && HVBINP4==2400 && fabs(time_max[Rm2]-time_max[MiB2]-")+std::string(sMean4_rm2mib2)+std::string(")<1");
  }

  h4binp1->Project("den1b","HVBINP1",denAndRm2AmpTim1);
  h4binp2->Project("den2b","HVBINP2",denAndRm2AmpTim2);
  h4binp3->Project("den3b","HVBINP3",denAndRm2AmpTim3);
  h4binp4->Project("den4b","HVBINP4",denAndRm2AmpTim4);




  // Compute mean timing distance between Binps and Mib2 - to be added to numerator (to be quoted only)
  TH1F* tim1mib2binp = new TH1F("tim1mib2binp", "tim1mib2binp",  4000,-20.,20.);
  TH1F* tim2mib2binp = new TH1F("tim2mib2binp", "tim2mib2binp",  4000,-20.,20.);
  TH1F* tim3mib2binp = new TH1F("tim3mib2binp", "tim3mib2binp",  4000,-20.,20.);
  TH1F* tim4mib2binp = new TH1F("tim4mib2binp", "tim4mib2binp",  4000,-20.,20.);
  h4binp1->Project("tim1mib2binp","time_max[BINP1]-time_max[MiB2]",denAndRm2AmpTim1+" && amp_max[BINP1]<3000 && amp_max[BINP1]>20");
  h4binp2->Project("tim2mib2binp","time_max[BINP2]-time_max[MiB2]",denAndRm2AmpTim2+" && amp_max[BINP2]<3000 && amp_max[BINP2]>20");
  h4binp3->Project("tim3mib2binp","time_max[BINP3]-time_max[MiB2]",denAndRm2AmpTim3+" && amp_max[BINP3]<3000 && amp_max[BINP3]>20");
  h4binp4->Project("tim4mib2binp","time_max[BINP4]-time_max[MiB2]",denAndRm2AmpTim4+" && amp_max[BINP4]<3000 && amp_max[BINP4]>20");

  char mean1_mib2binp[100];
  char mean2_mib2binp[100];
  char mean3_mib2binp[100];
  char mean4_mib2binp[100];
  sprintf(mean1_mib2binp,"%f",tim1mib2binp->GetMean());
  sprintf(mean2_mib2binp,"%f",tim2mib2binp->GetMean());
  sprintf(mean3_mib2binp,"%f",tim3mib2binp->GetMean());
  sprintf(mean4_mib2binp,"%f",tim4mib2binp->GetMean());
  std::string sMean1_mib2binp = std::string(mean1_mib2binp);
  std::string sMean2_mib2binp = std::string(mean2_mib2binp);
  std::string sMean3_mib2binp = std::string(mean3_mib2binp);
  std::string sMean4_mib2binp = std::string(mean4_mib2binp);

  TString numAmpTim1;
  if(tim1mib2binp->GetMean() < 0.){
    sMean1_mib2binp.erase(sMean1_mib2binp.begin(),sMean1_mib2binp.begin()+1);
    numAmpTim1 = denAndRm2AmpTim1+std::string(" && (amp_max[BINP1]>3000 || (amp_max[BINP1]<3000 && fabs(time_max[BINP1]-time_max[MiB2]+")+std::string(sMean1_mib2binp)+std::string(")<1))");
    } else {
    numAmpTim1 = denAndRm2AmpTim1+std::string(" && (amp_max[BINP1]>3000 || (amp_max[BINP1]<3000 && fabs(time_max[BINP1]-time_max[MiB2]-")+std::string(sMean1_mib2binp)+std::string(")<1))");
  }

  TString numAndRm2AmpTim2;
  if(tim2mib2binp->GetMean() < 0.){
    sMean2_mib2binp.erase(sMean2_mib2binp.begin(),sMean2_mib2binp.begin()+1);
    numAndRm2AmpTim2 = denAndRm2AmpTim2+std::string(" && (amp_max[BINP2]>3000 || (amp_max[BINP2]<3000 && fabs(time_max[BINP2]-time_max[MiB2]+")+std::string(sMean2_mib2binp)+std::string(")<1))");
    } else { 
    numAndRm2AmpTim2 = denAndRm2AmpTim2+std::string(" && (amp_max[BINP2]>3000 || (amp_max[BINP2]<3000 && fabs(time_max[BINP2]-time_max[MiB2]-")+std::string(sMean2_mib2binp)+std::string(")<1))");
  }

  TString numAndRm2AmpTim3;
  if(tim3mib2binp->GetMean() < 0.){
    sMean3_mib2binp.erase(sMean3_mib2binp.begin(),sMean3_mib2binp.begin()+1);
    numAndRm2AmpTim3 = denAndRm2AmpTim3+std::string(" && (amp_max[BINP3]>3000 || (amp_max[BINP3]<3000 && fabs(time_max[BINP3]-time_max[MiB2]+")+std::string(sMean3_mib2binp)+std::string(")<1))");
    } else { 
    numAndRm2AmpTim3 = denAndRm2AmpTim3+std::string(" && (amp_max[BINP3]>3000 || (amp_max[BINP3]<3000 && fabs(time_max[BINP3]-time_max[MiB2]-")+std::string(sMean3_mib2binp)+std::string(")<1))");
  }

  TString numAndRm2AmpTim4;
  if(tim4mib2binp->GetMean() < 0.){
    sMean4_mib2binp.erase(sMean4_mib2binp.begin(),sMean4_mib2binp.begin()+1);
    numAndRm2AmpTim4 = denAndRm2AmpTim4+std::string(" && (amp_max[BINP4]>3000 || (amp_max[BINP4]<3000 && fabs(time_max[BINP4]-time_max[MiB2]+")+std::string(sMean4_mib2binp)+std::string(")<1))");
    } else { 
    numAndRm2AmpTim4 = denAndRm2AmpTim4+std::string(" && (amp_max[BINP4]>3000 || (amp_max[BINP4]<3000 && fabs(time_max[BINP4]-time_max[MiB2]-")+std::string(sMean4_mib2binp)+std::string(")<1))");
  }


  // Nominal numerators
  TString snum1 = den1s + " && amp_max[BINP1]>17";
  TString snum2 = den2s + " && amp_max[BINP2]>14";
  TString snum3 = den3s + " && amp_max[BINP3]>14";
  TString snum4 = den4s + " && amp_max[BINP4]>15";

   //snum4 = snum4 + " && ((run == 2352 && fabs(time[MiB2]-time[BINP4]-2.316)<2.) || (run == 2353 && fabs(time[MiB2]-time[BINP4]-3.366)<2.) || (run == 2354 && fabs(time[MiB2]-time[BINP4]-2.520)<2.) || (run == 2355 && fabs(time[MiB2]-time[BINP4]-2.743)<2.) || (run == 2356 && fabs(time[MiB2]-time[BINP4]-2.884)<2.) || (run == 2357 && fabs(time[MiB2]-time[BINP4]-2.797)<2.) || (run == 2358 && fabs(time[MiB2]-time[BINP4]-2.886)<2.) || (run == 2359 && fabs(time[MiB2]-time[BINP4]-3.084)<2.) || (run == 2360 && fabs(time[MiB2]-time[BINP4]-2.885)<2.) || (run == 2361 && fabs(time[MiB2]-time[BINP4]-3.428)<2.) || (run == 2362 && fabs(time[MiB2]-time[BINP4]-3.648)<2.))";
   //snum4 = snum4 + " && ((run == 2352 && fabs(time[MiB2]-time[BINP4]-2.316)<5.) || (run == 2353 && fabs(time[MiB2]-time[BINP4]-3.366)<5.) || (run == 2354 && fabs(time[MiB2]-time[BINP4]-2.520)<5.) || (run == 2355 && fabs(time[MiB2]-time[BINP4]-2.743)<5.) || (run == 2356 && fabs(time[MiB2]-time[BINP4]-2.884)<5.) || (run == 2357 && fabs(time[MiB2]-time[BINP4]-2.797)<5.) || (run == 2358 && fabs(time[MiB2]-time[BINP4]-2.886)<5.) || (run == 2359 && fabs(time[MiB2]-time[BINP4]-3.084)<5.) || (run == 2360 && fabs(time[MiB2]-time[BINP4]-2.885)<5.) || (run == 2361 && fabs(time[MiB2]-time[BINP4]-3.428)<5.) || (run == 2362 && fabs(time[MiB2]-time[BINP4]-3.648)<5.))";

  //trigeref sample 
  //snum4 = snum4 + " && ((run == 2352 && fabs(time[MiB2]-time[BINP4]-2.365)<2.) || (run == 2353 && fabs(time[MiB2]-time[BINP4]-3.417)<2.) || (run == 2354 && fabs(time[MiB2]-time[BINP4]-2.960)<2.) || (run == 2355 && fabs(time[MiB2]-time[BINP4]-3.433)<2.) || (run == 2356 && fabs(time[MiB2]-time[BINP4]-3.018)<2.) || (run == 2357 && fabs(time[MiB2]-time[BINP4]-2.941)<2.) || (run == 2358 && fabs(time[MiB2]-time[BINP4]-2.025)<2.) || (run == 2359 && fabs(time[MiB2]-time[BINP4]-3.242)<2.) || (run == 2360 && fabs(time[MiB2]-time[BINP4]-3.520)<2.) || (run == 2361 && fabs(time[MiB2]-time[BINP4]-3.695)<2.) || (run == 2362 && fabs(time[MiB2]-time[BINP4]-3.900)<2.))";
  //snum4 = snum4 + " && ((run == 2352 && fabs(time[MiB2]-time[BINP4]-2.365)<5.) || (run == 2353 && fabs(time[MiB2]-time[BINP4]-3.417)<5.) || (run == 2354 && fabs(time[MiB2]-time[BINP4]-2.960)<5.) || (run == 2355 && fabs(time[MiB2]-time[BINP4]-3.433)<5.) || (run == 2356 && fabs(time[MiB2]-time[BINP4]-3.018)<5.) || (run == 2357 && fabs(time[MiB2]-time[BINP4]-2.941)<5.) || (run == 2358 && fabs(time[MiB2]-time[BINP4]-2.025)<5.) || (run == 2359 && fabs(time[MiB2]-time[BINP4]-3.242)<5.) || (run == 2360 && fabs(time[MiB2]-time[BINP4]-3.520)<5.) || (run == 2361 && fabs(time[MiB2]-time[BINP4]-3.695)<5.) || (run == 2362 && fabs(time[MiB2]-time[BINP4]-3.900)<5.))";
  snum4 = snum4 + " && ((run == 2352 && fabs(time[MiB2]-time[BINP4]-2.316)<10.) || (run == 2353 && fabs(time[MiB2]-time[BINP4]-3.366)<10.) || (run == 2354 && fabs(time[MiB2]-time[BINP4]-2.520)<10.) || (run == 2355 && fabs(time[MiB2]-time[BINP4]-2.743)<10.) || (run == 2356 && fabs(time[MiB2]-time[BINP4]-2.884)<10.) || (run == 2357 && fabs(time[MiB2]-time[BINP4]-2.797)<10.) || (run == 2358 && fabs(time[MiB2]-time[BINP4]-2.886)<10.) || (run == 2359 && fabs(time[MiB2]-time[BINP4]-3.084)<10.) || (run == 2360 && fabs(time[MiB2]-time[BINP4]-2.885)<10.) || (run == 2361 && fabs(time[MiB2]-time[BINP4]-3.428)<10.) || (run == 2362 && fabs(time[MiB2]-time[BINP4]-3.648)<10.))";

  h4binp1->Project("num1","HVBINP1",snum1); 
  h4binp2->Project("num2","HVBINP2",snum2); 
  h4binp3->Project("num3","HVBINP3",snum3); 
  h4binp4->Project("num4","HVBINP4",snum4); 

  // Numerators with timing infos - to be quoted only
  TString snum1b = numAmpTim1       + " && amp_max[BINP1]>17";
  TString snum2b = numAndRm2AmpTim2 + " && amp_max[BINP2]>14";
  TString snum3b = numAndRm2AmpTim3 + " && amp_max[BINP3]>14";
  TString snum4b = numAndRm2AmpTim4 + " && amp_max[BINP4]>15";
  h4binp1->Project("num1b","HVBINP1",snum1b); 
  h4binp2->Project("num2b","HVBINP2",snum2b); 
  h4binp3->Project("num3b","HVBINP3",snum3b); 
  h4binp4->Project("num4b","HVBINP4",snum4b); 

  cout << endl;
  cout << "-------------------------------------" << endl;
  cout << "check chiara: nominal " << endl;
  cout << "den1 = " << den1s << endl;
  cout << "den2 = " << den2s << endl;
  cout << "den3 = " << den3s << endl;
  cout << "den4 = " << den4s << endl;
  cout << endl;
  cout << "num1 = " << snum1 << endl;
  cout << "num2 = " << snum2 << endl;
  cout << "num3 = " << snum3 << endl;
  cout << "num4 = " << snum4 << endl;
  cout << "-------------------------------------" << endl;
  cout << endl;
  cout << endl;
  cout << "-------------------------------------" << endl;
  cout << "check chiara: with timing onfos " << endl;
  cout << "den1b = " << denAndRm2AmpTim1 << endl;
  cout << "den2b = " << denAndRm2AmpTim2 << endl;
  cout << "den3b = " << denAndRm2AmpTim3 << endl;
  cout << "den4b = " << denAndRm2AmpTim4 << endl;
  cout << endl;
  cout << "num1b = " << snum1b << endl;
  cout << "num2b = " << snum2b << endl;
  cout << "num3b = " << snum3b << endl;
  cout << "num4b = " << snum4b << endl;
  cout << "-------------------------------------" << endl;


  // efficiencies
  TGraphAsymmErrors* eff1 = new TGraphAsymmErrors(num1,den1);
  TGraphAsymmErrors* eff2 = new TGraphAsymmErrors(num2,den2);
  TGraphAsymmErrors* eff3 = new TGraphAsymmErrors(num3,den3);
  TGraphAsymmErrors* eff4 = new TGraphAsymmErrors(num4,den4);
  eff1->SetMarkerStyle(20);
  eff1->SetMarkerColor(1);
  eff1->SetLineColor(1);
  eff2->SetMarkerStyle(21);
  eff2->SetMarkerColor(2);
  eff2->SetLineColor(2);
  eff3->SetMarkerStyle(22);
  eff3->SetMarkerColor(3);
  eff3->SetLineColor(3);
  eff4->SetMarkerStyle(23);
  eff4->SetMarkerColor(4);
  eff4->SetLineColor(4);

  /*
  // Fitting function
  TF1 *func = new TF1("fit",fitfunc,0.,2500.,4);
  func->SetParLimits(1,0,1);

  cout << endl;
  cout << "Binp1 " << endl;
  func->SetParameters(0.2,1.,1,1200);
  eff1->Fit("fit","","",1200,2400);
  eff1->GetFunction("fit")->SetLineColor(1);
  float p3 = eff1->GetFunction("fit")->GetParameter(3);
  int VthrI = (int)p3 + 1;
  func->SetParameters(0.2,1.,1,1200);
  eff1->Fit("fit","","",VthrI,2400);
  eff1->GetFunction("fit")->SetLineColor(1);

  cout << endl;
  cout << "Binp2 " << endl;
  func->SetParameters(0.2,1.,1,500);
  eff2->Fit("fit","","",600,1300);
  eff2->GetFunction("fit")->SetLineColor(1);
  p3 = eff2->GetFunction("fit")->GetParameter(3);
  VthrI = (int)p3 + 1;
  func->SetParameters(0.2,1,1,500);
  eff2->Fit("fit","","",VthrI,1300);
  eff2->GetFunction("fit")->SetLineColor(1);

  cout << endl;
  cout << "Binp3 " << endl;
  func->SetParameters(0.2,1.,1,500);
  eff3->Fit("fit","","",600,1400);
  eff3->GetFunction("fit")->SetLineColor(1);
  p3 = eff3->GetFunction("fit")->GetParameter(3);
  VthrI = (int)p3 + 1;
  cout << "Binp3, VthrI = " << VthrI << endl;
  func->SetParameters(0.2,1.,1,500);
  eff3->Fit("fit","","",VthrI,1400);
  eff3->GetFunction("fit")->SetLineColor(1);

  cout << endl;
  cout << "Binp4 " << endl;
  func->SetParameters(0.2,1.,1,1300);
  eff4->Fit("fit","","",1300,2400);
  eff4->GetFunction("fit")->SetLineColor(1);
  p3 = eff4->GetFunction("fit")->GetParameter(3);
  VthrI = (int)p3 + 1;
  func->SetParameters(0.2,1.,1,1300);
  eff4->Fit("fit","","",VthrI,2400);
  eff4->GetFunction("fit")->SetLineColor(1);
  */

  TLegend *leg;
  leg = new TLegend(0.1,0.5,0.4,0.8);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.05);
  leg->SetFillColor(0);
  leg->AddEntry(eff1, "40x2+40x2 (Binp1)", "pl");
  leg->AddEntry(eff2, "40x1+40x2 (Binp2)", "pl");
  leg->AddEntry(eff3, "90x1+40x2 (Binp3)", "pl");
  leg->AddEntry(eff4, "90x2+40x2 (Binp4)", "pl");

  TCanvas* c1c = new TCanvas("c1c","c",1);
  c1c->SetGrid();
  TH2F* H2 = new TH2F("H2","",2500,0.,2500.,100,0.,1.);
  H2->GetXaxis()->SetTitle("Bias voltage [V]");
  H2->GetYaxis()->SetTitle("Efficiency");
  H2->Draw();
  eff1->Draw("Plsame"); 
  eff2->Draw("Plsame"); 
  eff3->Draw("Plsame"); 
  eff4->Draw("Plsame");
  leg->Draw();
  c1c->SaveAs("eff_binp2corr.png");
  c1c->SaveAs("eff_binp2corr.pdf");
  
  TFile myOutputF("myOutputFileBinps.root","RECREATE");
  myOutputF.cd();
  eff1->Write("binp1");
  eff2->Write("binp2");
  eff3->Write("binp3");
  eff4->Write("binp4");
  myOutputF.Close();
  
  // Efficiencies with additional timing infos
  cout << endl;
  cout << endl;
  cout << "----------------------------------------------------------------" << endl;
  cout << "Efficiencies in HV scan runs with timing infos as well: " << endl;
  cout << "Binp1 = " << num1b->Integral()/den1b->Integral()  << endl;
  cout << "Binp2 = " << num2b->Integral()/den2b->Integral()  << endl;
  cout << "Binp3 = " << num3b->Integral()/den3b->Integral()  << endl;
  cout << "Binp4 = " << num4b->Integral()/den4b->Integral()  << endl;
  cout << "----------------------------------------------------------------" << endl;
  /*
  cout << "Efficiencies in HV scan runs without timing infos: " << endl;
  cout << "Binp1 = " << num1->Integral()/den1->Integral()  << endl;
  cout << "Binp2 = " << num2->Integral()/den2->Integral()  << endl;
  cout << "Binp3 = " << num3->Integral()/den3->Integral()  << endl;
  cout << "Binp4 = " << num4->Integral()/den4->Integral()  << endl;
  */
  cout << endl;
  cout << endl;

  // Spectra comparison
  TH1F *amp1 = new TH1F("amp1","amp1",100,0.,4000);
  TH1F *amp2 = new TH1F("amp2","amp2",100,0.,4000);
  TH1F *amp3 = new TH1F("amp3","amp3",100,0.,4000);
  TH1F *amp4 = new TH1F("amp4","amp4",100,0.,4000);
  amp1->SetLineColor(1);
  amp2->SetLineColor(2);
  amp3->SetLineColor(3);
  amp4->SetLineColor(4);
  amp1->SetLineWidth(2);
  amp2->SetLineWidth(2);
  amp3->SetLineWidth(2);
  amp4->SetLineWidth(2);
  h4binp1 -> Project("amp1","amp_max[BINP1]",den1s+" && HVBINP1==2400");
  h4binp2 -> Project("amp2","amp_max[BINP2]",den2s+" && HVBINP2==1400"); 
  h4binp3 -> Project("amp3","amp_max[BINP3]",den3s+" && HVBINP3==1400");
  h4binp4 -> Project("amp4","amp_max[BINP4]",den4s+" && HVBINP4==2400");

  TCanvas* c1b = new TCanvas("cb","c",1);
  c1b->SetGrid();
  amp1->DrawNormalized(); 
  amp2->DrawNormalized("same"); 
  amp3->DrawNormalized("same"); 
  amp4->DrawNormalized("same"); 
  leg->Draw();
  c1b->SetLogy();
  c1b->SaveAs("amplitude.png");
  c1b->SaveAs("amplitude.pdf");

  TH1F *amp1z = new TH1F("amp1z","amp1z",20,0.,100);
  TH1F *amp2z = new TH1F("amp2z","amp2z",20,0.,100);
  TH1F *amp3z = new TH1F("amp3z","amp3z",20,0.,100);
  TH1F *amp4z = new TH1F("amp4z","amp4z",20,0.,100);
  amp1z->SetLineColor(1);
  amp2z->SetLineColor(2);
  amp3z->SetLineColor(3);
  amp4z->SetLineColor(4);
  amp1z->SetLineWidth(2);
  amp2z->SetLineWidth(2);
  amp3z->SetLineWidth(2);
  amp4z->SetLineWidth(2);
  h4binp1 -> Project("amp1z","amp_max[BINP1]",den1s+" && HVBINP1==2400");
  h4binp2 -> Project("amp2z","amp_max[BINP2]",den2s+" && HVBINP2==1400"); 
  h4binp3 -> Project("amp3z","amp_max[BINP3]",den3s+" && HVBINP3==1400");
  h4binp4 -> Project("amp4z","amp_max[BINP4]",den4s+" && HVBINP4==2400");

  TCanvas* c1bz = new TCanvas("cbz","c",1);
  c1bz->SetGrid();
  amp1z->DrawNormalized(); 
  amp2z->DrawNormalized("same"); 
  amp3z->DrawNormalized("same"); 
  amp4z->DrawNormalized("same"); 
  c1bz->SetLogy();
  c1bz->SaveAs("amplitudeZoom.png");
  c1bz->SaveAs("amplitudeZoom.pdf");



  // -------------------------------------------------------------------------------
  // Timing runs
  //
  //
  // Histos 
  TH1F* num14tim1 = new TH1F("num14tim1","",2500,-2500.,2500.);
  TH1F* num2tim   = new TH1F("num2tim",  "",2500,-2500.,2500.);
  TH1F* num3tim   = new TH1F("num3tim",  "",2500,-2500.,2500.);
  TH1F* num14tim4 = new TH1F("num14tim4","",2500,-2500.,2500.);

  TH1F* den14tim = new TH1F("den14tim","", 2500,-2500.,2500.);
  TH1F* den2tim  = new TH1F("den2tim", "", 2500,-2500.,2500.);
  TH1F* den3tim  = new TH1F("den3tim", "", 2500,-2500.,2500.);


  // denominators
  TString den3_timingRuns  = commonDenAmp;
  TString den2_timingRuns  = commonDenAndRm2Amp;
  TString den14_timingRuns = commonDenAndRm2Amp;
  h4timing14->Project("den14tim","X",den14_timingRuns);
  h4timing2 ->Project("den2tim", "X",den2_timingRuns);
  h4timing3 ->Project("den3tim", "X",den3_timingRuns);

  // numerators
  TString snum1_timing = den14_timingRuns + " && amp_max[BINP1]>20";
  TString snum2_timing = den2_timingRuns  + " && amp_max[BINP2]>20";
  TString snum3_timing = den3_timingRuns  + " && amp_max[BINP3]>20";
  TString snum4_timing = den14_timingRuns + " && amp_max[BINP4]>20";
  h4timing14->Project("num14tim1","X",snum1_timing);
  h4timing2 ->Project("num2tim",  "X",snum2_timing);
  h4timing3 ->Project("num3tim",  "X",snum3_timing);
  h4timing14->Project("num14tim4","X",snum4_timing);

  cout << "-----------------------------------------" << endl;
  cout << "Timing runs" << endl;
  cout << "den14: "  << den14_timingRuns << endl;
  cout << "den2: "   << den2_timingRuns  << endl;
  cout << "den3: "   << den3_timingRuns  << endl;
  cout << "num1: "   << snum1_timing << endl;
  cout << "num2: "   << snum2_timing << endl;
  cout << "num3: "   << snum3_timing << endl;
  cout << "num4: "   << snum4_timing << endl;

  float numT1 = num14tim1->Integral();
  float numT2 = num2tim->Integral();
  float numT3 = num3tim->Integral();
  float numT4 = num14tim4->Integral();
  //
  float denT14 = den14tim->Integral();
  float denT2  = den2tim->Integral();
  float denT3  = den3tim->Integral();

  cout << "Efficiencies in timing runs: " << endl;
  cout << "Binp1 = " << numT1/denT14  << endl;
  cout << "Binp2 = " << numT2/denT2   << endl;
  cout << "Binp3 = " << numT3/denT3   << endl;
  cout << "Binp4 = " << numT4/denT14  << endl;


  // Spectra comparison on timing runs
  TH1F *amp1T = new TH1F("amp1T","amp1T",100,0.,4000);
  TH1F *amp2T = new TH1F("amp2T","amp2T",100,0.,4000);
  TH1F *amp3T = new TH1F("amp3T","amp3T",100,0.,4000);
  TH1F *amp4T = new TH1F("amp4T","amp4T",100,0.,4000);
  amp1T->SetLineColor(1);
  amp2T->SetLineColor(2);
  amp3T->SetLineColor(3);
  amp4T->SetLineColor(4);
  amp1T->SetLineWidth(2);
  amp2T->SetLineWidth(2);
  amp3T->SetLineWidth(2);
  amp4T->SetLineWidth(2);
  h4timing14 ->Project("amp1T","amp_max[BINP1]",den14_timingRuns);
  h4timing2  ->Project("amp2T","amp_max[BINP2]",den2_timingRuns);
  h4timing3  ->Project("amp3T","amp_max[BINP3]",den3_timingRuns);
  h4timing14 ->Project("amp4T","amp_max[BINP4]",den14_timingRuns);

  TCanvas* c2a = new TCanvas("c2a","c",1);
  c2a->SetGrid();
  amp1T->DrawNormalized(); 
  amp2T->DrawNormalized("same"); 
  amp3T->DrawNormalized("same"); 
  amp4T->DrawNormalized("same"); 
  leg->Draw();
  c2a->SetLogy();
  c2a->SaveAs("amplitude_timingRuns.png");
  c2a->SaveAs("amplitude_timingRuns.pdf");

}

TString AddSelection(TTree* h4, TString Var, TString Selection, TString Cut = "0.")
{
    TH1F* h = new TH1F("h","h",40000,-20.,20.);
    h4->Draw(Var+" >> h",Selection);
    string sMean;

    char Mean [100];
    sprintf(Mean,"%f",h->GetMean());
    sMean = std::string(Mean);

    if(h->GetMean() < 0.){
       sMean.erase(sMean.begin(),sMean.begin()+1);
       Selection = Selection+" && fabs("+Var+"+"+sMean.c_str()+")<"+Cut;
    }else{
       Selection = Selection+" && fabs("+Var+"-"+sMean.c_str()+")<"+Cut;
    }
    
    delete h;
    return Selection;
}

	
