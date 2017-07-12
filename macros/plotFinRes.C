#include "TStyle.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

void plotRes(std::string FileName, std::string detector, int NumberGains){
	gStyle->SetOptStat(0);	

	if(NumberGains>3 || NumberGains<1)
	{
		cout << "Choose a number of gains between 1 and 3" << endl;
	}

	TCanvas* c0 = new TCanvas("c0", "c0");
	c0->cd();
	TLegend *legend = new TLegend(0.65, 0.6, 0.88, 0.88);


	if(NumberGains>=1)
	{
		TGraphErrors* g50 = new TGraphErrors(FileName.c_str(), "%lg %lg %lg");
		g50->SetTitle("");
		legend->AddEntry(g50, "Gain 50", "lep");
		g50->GetYaxis()->SetRangeUser(0, 0.32);
		g50->GetXaxis()->SetTitle("Energy (GeV)");
		g50->GetYaxis()->SetTitle("#sigma(APD-MCP) (ns)");
		g50->SetMarkerStyle(kFullCircle);
		g50->SetMarkerSize(1);
		g50->Draw();
	}
	if(NumberGains>=2)
	{
		TGraphErrors* g100 = new TGraphErrors(FileName.c_str(), "%lg %*lg %*lg %lg %lg");
		g100->SetTitle("");
		legend->AddEntry(g100, "Gain 100", "lep");
		g100->SetLineColor(kBlue);
		g100->SetMarkerStyle(kFullSquare);
		g100->SetMarkerSize(1);
		g100->SetMarkerColor(kBlue);
        	g100->Draw("PLSAME");
	}
	if(NumberGains>=3)
	{
		TGraphErrors* g200 = new TGraphErrors(FileName.c_str(), "%lg %*lg %*lg %*lg %*lg %lg %lg");
		g200->SetTitle("");	
		legend->AddEntry(g200, "Gain 200", "lep");
		g200->SetLineColor(kViolet);
		g200->SetMarkerStyle(kFullTriangleDown);
		g200->SetMarkerSize(1);
		g200->SetMarkerColor(6);
	        g200->Draw("PLSAME");
	}

	
	legend->Draw("SAME");
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/TimeRes_vs_Energy_"+detector+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/TimeRes_vs_Energy_"+detector+".pdf").c_str());


}
	 
