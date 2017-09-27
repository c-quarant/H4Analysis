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

#include <string>
#include <fstream>

void timeFit(std::string FileIN, std::string detector)
{
	//gStyle->SetOptStat(0);
	

	TFile *f = TFile::Open(FileIN.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	TH1D *tD1 = new TH1D("tD1", "", 100, -2, -0.5);

	std::string sel = "time[" + detector + "]-time[MCP1]>>tD1";
	h4->Draw(sel.c_str());
	gStyle->SetOptFit();	
	tD1->Fit("gaus");
	
	//setStyle();
	TCanvas* c2 = new TCanvas("c2","c2");
	TH2F* H1 = new TH2F("H3","", 100, -2, -0.5, 20, 0, (tD1->GetMaximum()*1.05));     
        H1->GetXaxis()->SetTitle("time ns");    
	H1->GetYaxis()->SetTitle("events");
	
	tD1->Draw();
	
	c2->SaveAs("/afs/cern.ch/user/c/cquarant/www/time_dist/FitProva.pdf");	
	c2->SaveAs("/afs/cern.ch/user/c/cquarant/www/time_dist/FitProva.png");
}	
