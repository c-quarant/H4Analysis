#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

void TwoGraphPlotter(std::string RootFile1, std::string RootFile2, std::string PathAndName)
{

	TFile *f1 = TFile::Open(RootFile1.c_str());
	TFile *f2 = TFile::Open(RootFile2.c_str());

	TGraphErrors* g1 = (TGraphErrors*)f1->Get("Graph");
	TGraphErrors* g2 = (TGraphErrors*)f2->Get("Graph");

	g1->SetMarkerStyle(kFullCircle);
	g1->SetMarkerColor(kRed);
	g2->SetMarkerStyle(kFullSquare);
	g2->SetMarkerColor(kBlue);

	TCanvas* cOut = new TCanvas("cOut", "cOut");
	g1->Draw("AP");
	g2->Draw("SAMEP");
	cOut->SaveAs(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".png").c_str());
	cOut->SaveAs(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".pdf").c_str());
	cOut->~TCanvas();

	TFile* f3 = new TFile(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".root").c_str(), "RECREATE");
	f3->cd();
	g1->Write();
	g2->Write();
	f3->Write();
	f3->Close();
	f3->~TFile();
}
