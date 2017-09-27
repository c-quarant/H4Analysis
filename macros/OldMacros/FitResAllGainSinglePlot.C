#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "MyLib.h"
#include "TVirtualFitter.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void FitRes(std::string FileName, std::string detector, int NumberGains){
	gStyle->SetOptStat(0);	

	TCanvas* c0 = new TCanvas("c0", "c0");
	c0->cd();
	
	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + [1]*[1] )", 15, 1000);

	fitFunc->SetParLimits(0, 0, 100000);
	fitFunc->SetParLimits(1, 0, 1000);
	//fitFunc->SetParLimits(2, 0, 2000);  

	fitFunc->SetParameter(0, 10000);
	fitFunc->SetParameter(1, 28);
	//fitFunc->SetParameter(2, 0);

	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");
	//fitFunc->SetParName(2, "~stochastic");

	TLegend* l = new TLegend(0.55, 0.6, 0.88, 0.88);
	
	TGraphErrors* g50 = new TGraphErrors(FileName.c_str(), "%lg %lg %lg %lg");
	g50->SetTitle("");

	g50->GetXaxis()->SetTitle("A/#sigma(Noise)");
	g50->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
	g50->SetLineColor(kBlue);
	g50->SetMarkerStyle(kFullCircle);
	g50->SetMarkerSize(1);
	g50->SetMarkerColor(kBlue);	
	g50->GetXaxis()->SetRangeUser(0, 1000);
	//g50->Fit("fitFunc");
	l->AddEntry(g50, "C0APD1-C0APD2 Gain 50", "lep");

	TH2D* Hset = new TH2D("Hset", "", 100, 0, 135, 100, 0, 700);	
	Hset->GetXaxis()->SetTitle("A/#sigma(Noise)");
	Hset->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
	Hset->Draw();
	g50->Draw("LPSAME");


	TGraphErrors* g100 = new TGraphErrors(FileName.c_str(), "%*lg %*lg %*lg %*lg %lg %lg %lg %lg");
	g100->SetTitle("");
	
	g100->GetXaxis()->SetTitle("A/#sigma(Noise)");
	g100->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
	g100->SetLineColor(kGreen);
	g100->SetMarkerStyle(kFullSquare);
	g100->SetMarkerSize(1);
	g100->SetMarkerColor(kGreen);
	//g100->Fit("fitFunc");
	l->AddEntry(g100, "C0APD1-C0APD2 Gain 50", "lep");
       	g100->Draw("LPSAME");

	if(NumberGains>=3)
	{
		TGraphErrors* g200 = new TGraphErrors(FileName.c_str(), "%*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg %lg %lg");
		g200->SetTitle("");	

		g200->GetXaxis()->SetTitle("A/#sigma(Noise)");
		g200->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
		g200->SetLineColor(kViolet);
		g200->SetMarkerStyle(kFullTriangleDown);
		g200->SetMarkerSize(1);
		g200->SetMarkerColor(kViolet);
		//g200->Fit("fitFunc");
		l->AddEntry(g200, "C0APD1-C0APD2 Gain 50", "lep");
	        g200->Draw("LPSAME");
	}

	l->Draw("SAME");

	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/TimeRes_vs_ANoise_"+detector+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/TimeRes_vs_ANoise_"+detector+".pdf").c_str());

}
	 
