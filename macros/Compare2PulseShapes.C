#include "TLeaf.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "TGaxis.h"
#include "TAxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <string>
#include <fstream>
#include <math.h>

void Compare2PulseShapes(std::string Waveform1, std::string Waveform2)
{
	gStyle->SetOptStat(0);

	TFile *f1 = TFile::Open(Waveform1.c_str());
	TFile *f2 = TFile::Open(Waveform2.c_str());

	TF1 *w1 = (TF1*)f1->Get("Waveform_");
	w1->SetTitle("");
	TH1D *w2 = (TH1D*)f2->Get("Waveform_");
	w2->SetTitle("");

	TCanvas *c0 = new TCanvas("c0","c0");
	c0->cd();

	//w1->GetXaxis()->SetTitle("time (ns)");
	//lw1->GetYaxis()->SetTitle("ADC counts");
	w1->Draw("L");
	w2->SetLineColor(2);
	w2->SetMarkerColor(2);
	w2->Draw("SAME");

	TLegend* legend = new TLegend(0.60, 0.7, 0.88, 0.88);
	legend->AddEntry(w1, "C0APD1 100 Gev Gain 50");
	legend->AddEntry(w2, "C3 100 Gev Gain 50");
	legend->Draw();
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Comp_C0APD1_C3_100Gev_G50.pdf"));
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Comp_C0APD1_C3_100Gev_G50.png"));


}
