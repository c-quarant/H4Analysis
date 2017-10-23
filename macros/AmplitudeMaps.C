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
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "FPCanvasStyle.C"
#include "TProfile2D.h"
#include "MyLib.h"
#include <string>
#include <fstream>

void AmplitudeMaps(std::string FileIN, std::string detector)
{	
  	int i, Nentries, CH, runNum; 
  	float AmpTemp;
	std::string Xshift_str, Yshift_str;
  	gStyle->SetOptStat(0);
  	
  	TFile *f = TFile::Open(FileIN.c_str());
  	TTree *h4 = (TTree*)f->Get("h4");
  	h4->GetEntry(0);
 	CH=h4->GetLeaf(detector.c_str())->GetValue(0);

	std::string pathToLinFitOutput = "/afs/cern.ch/user/c/cquarant/www/";

	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	runNum = h4->GetLeaf("run")->GetValue(0);

	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = Energy+"Gev_G"+Gain+"_"+std::to_string(runNum);
	
	//Calculating shift of hodo plane 1 respect to hodo plane 0 faro X and Y axis
	Xshift_str = std::to_string(HodoPlaneShift(h4, detector, pathToLinFitOutput, RunStats, "X"));
	Yshift_str = std::to_string(HodoPlaneShift(h4, detector, pathToLinFitOutput, RunStats, "Y"));
 
  	//2DHist definition 
  	auto *AmpXY0 = new TProfile2D("AmpXY0","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	auto *AmpXY1 = new TProfile2D("AmpXY1","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	auto *AmpXYM = new TProfile2D("AmpXYM","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	
  	Nentries = h4->GetEntries();
  	
	
   	h4->Draw(("amp_max["+detector+"]:Y[0]:X[0]>>AmpXY0").c_str(), "amp_max[MCP1]>200 && amp_max[MCP1]<2000 && X[0]>-800 && Y[0]>-800");
  	h4->Draw(("amp_max["+detector+"]:Y[1]-("+Yshift_str+"):X[1]-("+Xshift_str+")>>AmpXY1").c_str(), "amp_max[MCP1]>200 && amp_max[MCP1]<2000 && X[1]>-800 && Y[1]>-800");
  	h4->Draw(("amp_max["+detector+"]:(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>AmpXYM").c_str(), "amp_max[MCP1]>200 && amp_max[MCP1]<2000 && X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800");

  	//Drawing p0 histogram
  	TCanvas* c1 = new TCanvas("c1","c1");
	FPCanvasStyle(c1, "", "", 0, "", 0, 1);
  	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");
  	AmpXY0->GetZaxis()->SetTitle(("amp_max["+detector+"] (ADC counts)").c_str());
  	
  	H1->Draw();	
  	AmpXY0->Draw("COLZ SAME");
  	
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p0" + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p0" + ".png";
  	
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  	
  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
  	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p1" + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_p1" + ".png";
  	
	
  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	
  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");
  	
  	H1->Draw();	
  	AmpXYM->Draw("COLZ SAME");
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_NoMCPCut_" + detector + "_G" + Gain + "_" + Energy + "Gev_pAVG" + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_NoMCPCUt_" + detector + "_G" + Gain + "_" + Energy + "Gev_pAVG" + ".png";
  	
  	
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());

	auto *EventsXYM = new TH2D("EventsXYM","", 32, -16, 16, 32, -16, 16); 
  	
  	h4->Draw(("(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>EventsXYM").c_str(), "amp_max[MCP1]>200 && amp_max[MCP1]<2000 && X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800");

	//Drawving pAVG histogram
  	TCanvas* c4 = new TCanvas("c4","c4");
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");
  	
  	H1->Draw();	
  	EventsXYM->Draw("COLZ SAME");
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/EventsXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_pAVG" + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/EventsXY_" + detector + "_G" + Gain + "_" + Energy + "Gev_pAVG" + ".png";
  	
  	
  	c4->SaveAs(fileOutpdf.c_str());
  	c4->SaveAs(fileOutpng.c_str());
}				
