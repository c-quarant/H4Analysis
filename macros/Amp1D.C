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

void Amp1D(std::string FileIn, std::string detector, float bound)
{
	int Nentries, i, CH;
	Float_t X[2], Y[2], Ampl[6], Xavg, Yavg;
	gStyle->SetOptStat(0);

	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	h4->GetEntry(0);
	
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	
	auto *AmpX0 = new TProfile("AmpX0", "", 32, -16, 16, 0, 10000);
	auto *AmpY0 = new TProfile("AmpY0", "", 32, -16, 16, 0, 10000);
	auto *AmpX1 = new TProfile("AmpX1", "", 32, -16, 16, 0, 10000);
	auto *AmpY1 = new TProfile("AmpY1", "", 32, -16, 16, 0, 10000);
	auto *AmpXavg = new TProfile("AmpXavg", "", 32, -16, 16, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 32, -16, 16, 0, 10000);
	
	h4->SetBranchAddress("amp_max", Ampl);	
	h4->SetBranchAddress("X", X);	
	h4->SetBranchAddress("Y", Y);

	for(i=0; i<Nentries; i++)
	   {
		h4->GetEntry(i);
		
		if(Y[0]>-bound && Y[0]<bound) AmpX0->Fill(X[0], Ampl[CH]);
		if(X[0]>-bound && X[0]<bound) AmpY0->Fill(Y[0], Ampl[CH]);

		if(Y[1]>-bound && Y[1]<bound) AmpX1->Fill(X[1], Ampl[CH]);
		if(X[1]>-bound && X[1]<bound) AmpY1->Fill(Y[1], Ampl[CH]);

		if(X[0]>-800 && X[1]>-800 && Y[0]>-800 && Y[1]>-800)
		   {
			Xavg = 0.5*(X[0]+X[1]);
			Yavg = 0.5*(Y[0]+Y[1]);
			if(Yavg>-bound && Yavg<bound) AmpXavg->Fill(X[0], Ampl[CH]);
			if(Xavg>-bound && Xavg<bound) AmpYavg->Fill(Y[0], Ampl[CH]); 
		   }
	   }
	
	//Drawing X0 histogram
	TCanvas* c1 = new TCanvas("c1","c1");
 	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 50, 0, (AmpX0->GetMaximum())*1.05);     
	H1->GetXaxis()->SetTitle("X[0]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();	
	AmpX0->Draw("SAME");
  
  	std::string sel;
  	if(CH==0)sel =  "C3Gain";
  	else if(CH==3 || CH==4) sel = "C0Gain";

  	std::string Gain = std::to_string((int)h4->GetLeaf(sel.c_str())->GetValue(0));
	cout << Gain << endl;
  	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
 
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p0/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_X0[+-" + std::to_string((int)bound) + "].pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p0/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_X0[+-" + std::to_string((int)bound) + "].png";
  
  
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
	
	//Drawing Y0 histogram
	TCanvas* c2 = new TCanvas("c2","c2");
 	H1->GetXaxis()->SetTitle("Y[0]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();	
	AmpY0->Draw("SAME");
  
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p0/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Y0[+-" + std::to_string((int)bound) + "].pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p0/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Y0[+-" + std::to_string((int)bound) + "].png";
  
  
  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());

	//Drawing X1 histogram
	TCanvas* c3 = new TCanvas("c3","c3");
 	H1->GetXaxis()->SetTitle("X[1]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();	
	AmpX1->Draw("SAME");
  
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p1/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_X1[+-" + std::to_string((int)bound) + "].pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p1/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_X1[+-" + std::to_string((int)bound) + "].png";
  
  
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());

	//Drawing Y1 histogram
	TCanvas* c4 = new TCanvas("c4","c4");
 	H1->GetXaxis()->SetTitle("Y[1]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();	
	AmpY1->Draw("SAME");
  
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p1/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Y1[+-" + std::to_string((int)bound) + "].pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amp1D/p1/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Y1[+-" + std::to_string((int)bound) + "].png";
  
  
  	c4->SaveAs(fileOutpdf.c_str());
  	c4->SaveAs(fileOutpng.c_str());

	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
 	H1->GetXaxis()->SetTitle("Xavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();	
	AmpXavg->Draw("SAME");
  
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amp1D/pAVG/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Xavg[+-" + std::to_string((int)bound) + "].pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amp1D/pAVG/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Xavg[+-" + std::to_string((int)bound) + "].png";
  
  
  	c5->SaveAs(fileOutpdf.c_str());
  	c5->SaveAs(fileOutpng.c_str());
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();	
	AmpYavg->Draw("SAME");
  
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amp1D/pAVG/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Yavg[+-" + std::to_string((int)bound) + "].pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amp1D/pAVG/Amp_" + detector + "_G" + Gain + "_" + Energy + "Gev_Yavg[+-" + std::to_string((int)bound) + "].png";
  
  
  	c6->SaveAs(fileOutpdf.c_str());
  	c6->SaveAs(fileOutpng.c_str());
}

















