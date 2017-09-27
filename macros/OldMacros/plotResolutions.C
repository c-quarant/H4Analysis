#include "TStyle.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"

void plotRes(std::string FileName){
	gStyle->SetOptStat(0);	

	TGraphErrors* g20 = new TGraphErrors(FileName.c_str(), "%lg %lg %lg");
	g20->SetTitle("");
	TGraphErrors* g50 = new TGraphErrors(FileName.c_str(), "%lg %*lg %*lg %lg %lg");
	g50->SetTitle("");
	TGraphErrors* g100 = new TGraphErrors(FileName.c_str(), "%lg %*lg %*lg %*lg %*lg %lg %lg");
	g100->SetTitle("");	
	TGraphErrors* g150 = new TGraphErrors(FileName.c_str(), "%lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
	g150->SetTitle("");
	TGraphErrors* g200 = new TGraphErrors(FileName.c_str(), "%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %lg %lg");
	g200->SetTitle("");

	TCanvas* c0 = new TCanvas("c0", "c0");
	c0->cd();
	g20->GetYaxis()->SetRangeUser(0, 0.25);
	g20->GetXaxis()->SetTitle("fit window width (# samples)");
	g20->GetYaxis()->SetTitle("time resolution (ns)");
	g20->SetMarkerStyle(20);
	g20->SetMarkerSize(1);
	g20->Draw();

	g50->SetLineColor(2);
	g50->SetMarkerStyle(21);
	g50->SetMarkerSize(1);
	g50->SetMarkerColor(2);
        g50->Draw("PL SAME");

	g100->SetLineColor(3);
	g100->SetMarkerStyle(22);
	g100->SetMarkerSize(1);
	g100->SetMarkerColor(3);
        g100->Draw("PLSAME");

	g150->SetLineColor(4);
	g150->SetMarkerStyle(kFullDiamond);
	g150->SetMarkerSize(1);
	g150->SetMarkerColor(4);
        g150->Draw("PLSAME");

	g200->SetLineColor(6);
	g200->SetMarkerStyle(kFullTriangleDown);
	g200->SetMarkerSize(1);
	g200->SetMarkerColor(6);
        g200->Draw("PLSAME");

	TLegend *legend = new TLegend(0.65, 0.6, 0.88, 0.88);
	legend->SetHeader("Energy");
	legend->AddEntry(g20, "20 Gev", "lep");
	legend->AddEntry(g50, "50 Gev", "lep");
	legend->AddEntry(g100, "100 Gev", "lep");
	legend->AddEntry(g150, "150 Gev", "lep");
	legend->AddEntry(g200, "200 Gev", "lep");
	legend->Draw("SAME");
	
	c0->SaveAs("/afs/cern.ch/user/c/cquarant/www/fitWindow/TimeRes_vs_fitWin_C3.png");
	c0->SaveAs("/afs/cern.ch/user/c/cquarant/www/fitWindow/TimeRes_vs_fitWin_C3.pdf");

}
	 
