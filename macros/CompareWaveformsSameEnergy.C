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
//#include "FPCanvasStyle.C"
//#include "setStyle.C"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <string>
#include <vector>
#include <fstream>
#include <math.h>

void CompareWaveformsSameEnergy(std::string detector, std::string Energy, float xmin, float xmax)
{
	unsigned int i, j;
	std::vector <TFile*> WaveFiles;
	std::vector <std::string> Gains;
	std::vector <TH1D*> Ratios;
	gStyle->SetOptStat(0);

	std::string WF_G50 = "WaveForms/"+detector+"_"+Energy+"Gev_G50_Waveform.root";
	std::string WF_G100 = "WaveForms/"+detector+"_"+Energy+"Gev_G100_Waveform.root";
	std::string WF_G200 = "WaveForms/"+detector+"_"+Energy+"Gev_G200_Waveform.root";
	
	TFile *f0 = TFile::Open(WF_G50.c_str());
	TFile *f2 = TFile::Open(WF_G100.c_str());
	TFile *f4 = TFile::Open(WF_G200.c_str());
	
	if(f2)
	{
		WaveFiles.push_back(f2);
		Gains.push_back("Gain 100");
		Ratios.push_back(new TH1D("R2", "", 300, -10, 40));
	}
			
	if(f4)
	{
		WaveFiles.push_back(f4);
		Gains.push_back("Gain 200");
		Ratios.push_back(new TH1D("R4", "", 300, -10, 40));
	}		
	if(f0)
	{
		WaveFiles.push_back(f0);	
		Gains.push_back("Gain 50");
		Ratios.push_back(new TH1D("R-REF", "", 300, -10, 40));
	}

	TCanvas *c0 = new TCanvas("c0","c0");
	c0->cd();

	TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.5,1.00,1.00);
	TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.5);

	cUp->SetBottomMargin(0.01); 
	cDown->SetTopMargin(0.01); 
    
    	cUp->Draw();
	cDown->Draw();

	TLegend* legend = new TLegend(0.65, 0.4, 0.89, 0.88);
	legend->SetHeader(("Detector "+detector+" Energy "+Energy).c_str());
	
	Ratios[0]->SetTitle("");
	Ratios[0]->GetXaxis()->SetTitle("Time_MCP1 (ns)");
	Ratios[0]->GetXaxis()->SetRangeUser(xmin, xmax);	
	Ratios[0]->GetYaxis()->SetRangeUser(0.8, 1.1);

	TH1D* w0 = new TH1D("w0", "", 300, -10, 40);
	w0->GetYaxis()->SetTitleSize(0.05);	
	w0->GetXaxis()->SetRangeUser(xmin, xmax);	
	w0->GetYaxis()->SetRangeUser(-0.05, 1.05);

	TH1D *wREF = new TH1D("wREF", "", 300, -10, 40); 
	TLegend *Rlegend = new TLegend(0.5, 0.88, 0.88, 0.96);
	Rlegend->SetHeader(("Ratio plot whit respect to " + Gains.back()).c_str());

	cUp->cd();
	wREF = (TH1D*)WaveFiles[WaveFiles.size()-1]->Get("Waveform_");
	wREF->SetLineColor(1);
	wREF->SetTitle("");
	wREF->GetXaxis()->SetRangeUser(xmin, xmax);	
	
	legend->AddEntry(wREF, (Gains.back() + " Ref value").c_str());

	wREF->Draw("R");
	
	for(i=0; i<WaveFiles.size()-1; i++)
	{
		cUp->cd();

		w0 = (TH1D*)WaveFiles[i]->Get("Waveform_");
		w0->SetLineColor(i+2);
		w0->SetTitle("");
			
		legend->AddEntry(w0, Gains[i].c_str());
	
		w0->Draw("SAME");
		legend->Draw("SAME");
	
		cDown->cd();
	
		for(j=0; j<300; j++)
		{
			Ratios[i]->SetBinContent(j, w0->GetBinContent(j));
			Ratios[i]->SetBinError(j, w0->GetBinError(j));
		}
		
		Ratios[i]->Divide(wREF);
		
		Ratios[i]->SetLineColor(i+2);
	
		//Rlegend->AddEntry(Ratios[i], (Energies[i]+"/100Gev").c_str());
	
		if(i==0)Ratios[i]->Draw();
		else Ratios[i]->Draw("SAME");			
	}

	TLine *line = new TLine(xmin, 1, xmax, 1);
	line->SetLineColor(1);
  	line->Draw("SAME");

	Rlegend->Draw("SAME");

	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Energy"+Energy+".pdf").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/PulseShapes/ComparePulseShape/Ratio_"+detector+"_Energy"+Energy+".png").c_str());

}
