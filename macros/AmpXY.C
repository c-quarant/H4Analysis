#include "TLeaf.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile2D.h"
#include <string>
#include <fstream>

void AmpXYP(std::string FileIN, std::string detector)
{
  int i, Nentries, CH; 
  bool Icond=false, IIcond=false;
  Float_t Ampl[6], X[2], Y[2], AmpTemp;
  gStyle->SetOptStat(0);
  
  TFile *f = TFile::Open(FileIN.c_str());
  TTree *h4 = (TTree*)f->Get("h4");
  h4->GetEntry(0);
  CH=h4->GetLeaf(detector.c_str())->GetValue(0);
 
  //2DHist definition 
  auto *AmpXY0 = new TProfile2D("AmpXY0","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  auto *AmpXY1 = new TProfile2D("AmpXY1","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  auto *AmpXYM = new TProfile2D("AmpXYM","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  
  Nentries = h4->GetEntries();
  
  h4->SetBranchAddress("amp_max", Ampl);
  h4->SetBranchAddress("X", X);
  h4->SetBranchAddress("Y", Y);	
  
  //Filling TProfile2D histograms for p0, p1, pAVG
  for(i=0; i<Nentries; i++)
    {
      h4->GetEntry(i);
      
      if(X[0]>-800 && Y[0]>-800)
	{ 
	  AmpXY0->Fill(X[0], Y[0], Ampl[CH]);
	  Icond = true;
	}	  
      
      if(X[1]>-800 && Y[1]>-800)
	{
	  AmpXY1->Fill(X[1], Y[1], Ampl[CH]);
	  IIcond = true;
	}

      if(Icond && IIcond)
	{
	  AmpXYM->Fill( 0.5*(X[0]+X[1]), 0.5*(Y[0]+Y[1]), Ampl[CH]);
	}
    }
  
  //Drawing p0 histogram
  TCanvas* c1 = new TCanvas("c1","c1");
  TH2F* H1 = new TH2F("H1","", 32, -16, 16, 32, -16, 16);     
  H1->GetXaxis()->SetTitle("X[0]");    
  H1->GetYaxis()->SetTitle("Y[0]");
  
  H1->Draw();	
  AmpXY0->Draw("COLZ SAME");
  
  std::string sel;
  if(CH==0)sel =  "C3Gain";
  else if(CH==3 || CH==4) sel = "C0Gain";

  std::string Gain = std::to_string((int)h4->GetLeaf(sel.c_str())->GetValue(0));

  std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
 
  std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmpXY/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p0" + ".pdf";
  std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmpXY/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p0" + ".png";
  
  
  c1->SaveAs(fileOutpdf.c_str());
  c1->SaveAs(fileOutpng.c_str());
  
  //Drawing p1 histogram
  TCanvas* c2 = new TCanvas("c2","c2");
  H1->GetXaxis()->SetTitle("X[1]");    
  H1->GetYaxis()->SetTitle("Y[1]");
  
  H1->Draw();	
  AmpXY1->Draw("COLZ SAME");
  
  fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmpXY/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p1" + ".pdf";
  fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmpXY/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p1" + ".png";
  

  c2->SaveAs(fileOutpdf.c_str());
  c2->SaveAs(fileOutpng.c_str());
  
  //Drawving pAVG histogram
  TCanvas* c3 = new TCanvas("c3","c3");
  H1->GetXaxis()->SetTitle("X_AVG");    
  H1->GetYaxis()->SetTitle("Y_AVG");
  
  H1->Draw();	
  AmpXY1->Draw("COLZ SAME");
  
  fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmpXY/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_pAVG" + ".pdf";
  fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmpXY/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_pAVG" + ".png";
  
  
  c3->SaveAs(fileOutpdf.c_str());
  c3->SaveAs(fileOutpng.c_str());
}
