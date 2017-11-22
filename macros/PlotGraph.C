#include "TFile.h"
#include "TGraphErrors.h"
#include "TH2.h"
#include "TCanvas.h"
#include <string>
#include <iostream>
#include "stdio.h"
#include "TStyle.h"
#include "TGaxis.h"

void PlotGraph(std::string RootFile, std::string detector, std::string TimeRef)
{
	TFile *f = TFile::Open(RootFile.c_str());
	
	TGraphErrors *g = (TGraphErrors*)f->Get("Graph;1");
	TH2F *Hset = (TH2F*)f->Get("Hset");
	
	TGaxis::SetMaxDigits(2);
	gStyle->SetOptStat(0);

	Hset->GetXaxis()->SetTitle("A/#sigma_{n}");
	Hset->GetXaxis()->SetTitleSize(0.055);
	Hset->GetXaxis()->SetTitleOffset(0.75);

	Hset->GetYaxis()->SetTitle(("#bar{t_{"+detector+"}-t_{"+TimeRef+"}} (ps)").c_str());
	Hset->GetYaxis()->SetTitleSize(0.055);
	Hset->GetYaxis()->SetTitleOffset(0.9);

	TCanvas *c0 = new TCanvas("c0", "c0");
	c0->cd();
	Hset->Draw();	
	g->Draw("PSAME");

	RootFile.replace(RootFile.find(".root"), std::string(".root").length(), ".png");
	c0->SaveAs(RootFile.c_str());

	RootFile.replace(RootFile.find(".png"), std::string(".png").length(), ".pdf");
	c0->SaveAs(RootFile.c_str());
}
