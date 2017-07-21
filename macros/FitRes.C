#include "TStyle.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TMath.h"
#include "MyLib.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void FitRes(std::string FileName, std::string detector, int NumberGains){
	gStyle->SetOptFit();	

	TCanvas* c0 = new TCanvas("c0", "c0");
	c0->cd();
	
	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + [1]*[1])", 15, 1000);

	fitFunc->SetParLimits(0, 0, 10000);
	fitFunc->SetParLimits(1, 0, 1000); 

	fitFunc->SetParameter(0, 100);
	fitFunc->SetParameter(1, 40);

	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");

	if(NumberGains==1)
	{
		TGraphErrors* g50 = new TGraphErrors(FileName.c_str(), "%lg %lg %lg %lg");
		g50->SetTitle("");

		g50->GetXaxis()->SetTitle("A/#sigma(Noise)");
		g50->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
		g50->SetMarkerStyle(kFullCircle);
		g50->SetMarkerSize(1);
		g50->Fit("fitFunc");
		g50->Draw("AP");
	
		detector += "_G50";
	}
	if(NumberGains==2)
	{
		TGraphErrors* g100 = new TGraphErrors(FileName.c_str(), "%*lg %*lg %*lg %*lg %lg %lg %lg %lg");
		g100->SetTitle("");

		g100->GetXaxis()->SetTitle("A/#sigma(Noise)");
		g100->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
		g100->SetLineColor(kBlue);
		g100->SetMarkerStyle(kFullSquare);
		g100->SetMarkerSize(1);
		g100->SetMarkerColor(kBlue);
		g100->Fit("fitFunc");
        	g100->Draw("AP");

		detector += "_G100";
	}
	if(NumberGains>=3)
	{
		TGraphErrors* g200 = new TGraphErrors(FileName.c_str(), "%*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg %lg %lg");
		g200->SetTitle("");	

		g200->GetXaxis()->SetTitle("A/#sigma(Noise)");
		g200->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
		g200->SetLineColor(kViolet);
		g200->SetMarkerStyle(kFullTriangleDown);
		g200->SetMarkerSize(1);
		g200->SetMarkerColor(6);
		g200->Fit("fitFunc");
	        g200->Draw("AP");

		detector += "_G200";
	}

	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/TimeRes_vs_ANoise_"+detector+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/TimeRes_vs_ANoise_"+detector+".pdf").c_str());


}
	 
