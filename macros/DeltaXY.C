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
#include "TProfile.h"
#include <string>
#include <fstream>

void DeltaXY(std::string FileIn, std::string detector)
{
	int Nentries, i, CH;
	Float_t X[2], Y[2], Xavg, Yavg;
	gStyle->SetOptStat(0);

	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	h4->GetEntry(0);
	
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	
	TH1F *h1 = new TH1F("h1", "", 80, -5, 2);
	auto *h2 = new TH1F("h2", "", 128, -16, 16, -10, 10);

	//DeltaX histogram
	h4->Draw("(X[0]-X[1])>>h1", "X[0]>-800 && X[1]>-800");
		
	TCanvas* c1 = new TCanvas("c1","c1");
 	TH2F* H1 = new TH2F("H1","", 80, -5, 2, 50, 0, (h1->GetMaximum())*1.05);     
	H1->GetXaxis()->SetTitle("X[0]-X[1]");    
 	H1->GetYaxis()->SetTitle("events");
  
	H1->Draw();	
	h1->Draw("SAME");
  
  	std::string sel;
  	if(CH==0)sel =  "C3Gain";
  	else if(CH==3 || CH==4) sel = "C0Gain";

  	std::string Gain = std::to_string((int)h4->GetLeaf(sel.c_str())->GetValue(0));
  	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
 
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/DeltaX/DX_" + detector + "_G" + Gain + "_" + Energy + "Gev.pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/DeltaX/DX_" + detector + "_G" + Gain + "_" + Energy + "Gev.png";
    
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());

	//Drawing Deltax vs X0 histogram
	h4->Draw("(X[0]-X[1]):X[0]>>h2", "X[0]>-800 && X[1]>-800");
	
	TCanvas *c2 = new TCanvas("c2", "c2");
	TH2F* H2 = new TH2F("H2","", 128, -16, 16, 50, 0, (h2->GetMaximum())*1.05);     
	H2->GetXaxis()->SetTitle("X[0]");    
 	H2->GetYaxis()->SetTitle("X[0]-X[1]");
  
	H2->Draw();	
	h2->Draw("SAME");
  
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/DeltaX/DXvsX_" + detector + "_G" + Gain + "_" + Energy + "Gev.pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/DeltaX/DXvsX_" + detector + "_G" + Gain + "_" + Energy + "Gev.png";
    
  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());

	//Drawing DeltaY histogram
	h4->Draw("(Y[0]-Y[1])>>h1", "Y[0]>-800 && Y[1]>-800");
	
	TCanvas* c3 = new TCanvas("c3","c3");
 	H1->GetXaxis()->SetTitle("Y[0]-Y[1]");    
  
	H1->Draw();	
	h1->Draw("SAME");
  
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/DeltaY/DY_" + detector + "_G" + Gain + "_" + Energy + "Gev.pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/DeltaY/DY_" + detector + "_G" + Gain + "_" + Energy + "Gev.png";
  
  
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());
}	
