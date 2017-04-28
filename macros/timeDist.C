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

void timeDist(std::string FileIN, std::string detector, std::string cutValue, std::string run)
{
	gStyle->SetOptStat(0);
	

	TFile *f = TFile::Open(FileIN.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	
	//plot eff
	TH2F *h_num = new TH2F("h_num", "", 64, -32, 32, 64, -32, 32);
	TH2F *h_denum = new TH2F("h_denum", "", 64, -32, 32, 64, -32, 32);
	TH2F *h_eff = new TH2F("h_eff", "", 64, -32, 32, 64, -32, 32);
	
	std::string cut="amp_max[" + detector + "]>" + cutValue + " && X[0]>-800 && Y[0]>-800";
	
	h4->Draw("Y[0]:X[0] >> h_num", cut.c_str());
        h4->Draw("Y[0]:X[0] >> h_denum","X[0]>-800 && Y[0]>-800");  
	h_eff->Divide(h_num,h_denum,1,1,"B");

	TCanvas* c3 = new TCanvas("c3","c3");
	TH2F* H3 = new TH2F("H3","", 64, -32, 32, 64, -32, 32);     
        H3->GetXaxis()->SetTitle("X");    
	H3->GetYaxis()->SetTitle("Y");
	
	H3->Draw();	
	h_eff->Draw("COLZ SAME");
	
	std::string plotEffpdf = "Eff_" + detector + "_" + cutValue + "_" + run  + ".pdf";
	c3->SaveAs(plotEffpdf.c_str());
	
	std::string plotEffpng = "Eff_" + detector + "_" + cutValue + "_" + run  + ".png";
	c3->SaveAs(plotEffpng.c_str());

	//plot time distribution ref MCP1 
	TH1F *h1 = new TH1F("h1", "", 100, -5, 0.5);
	std::string sel1="(time[" + detector + "]-time[MCP1])>>h1";
	cut = "amp_max[" + detector + "]>" + cutValue;
	cout << cut << endl;	
	h4->Draw(sel1.c_str(), cut.c_str());
	
	setStyle();      
	TH2F* H1 = new TH2F("H1","", 100, -5, 0.5, 20, 0, (h1->GetMaximum()*(1.05)));     
        H1->GetXaxis()->SetTitle("time ns");    
	H1->GetYaxis()->SetTitle("events");

	TCanvas* c1 = new TCanvas("c1", "c1");
	FPCanvasStyle(c1);
	
	H1->Draw();
	h1->Draw("h1 SAME");

	std::string plotnamepdf = "TimeDist_" + run + "_" + detector + "-MPC1.pdf"; 
	std::string plotnamepng = "TimeDist_" + run + "_" + detector + "-MPC1.png"; 
	
	c1->SaveAs(plotnamepdf.c_str());	
	c1->SaveAs(plotnamepng.c_str());	
	
	//plot time distribution ref MCP2
	TH1F *h2 = new TH1F("h2", "", 100, -0.5, 5);
	std::string sel2="(time[" + detector + "]-time[MCP2])>>h2";
	h4->Draw(sel2.c_str());

	setStyle();
	TH2F* H2 = new TH2F("H2","", 100, -0.5, 5, 20, 0, (h1->GetMaximum()*1.05));     
        H2->GetXaxis()->SetTitle("time ns");    
	H2->GetYaxis()->SetTitle("events");

	TCanvas* c2 = new TCanvas("c2", "c2");
	FPCanvasStyle(c2);

	H2->Draw();
	h2->Draw("h2 SAME");

	plotnamepdf = "TimeDist_" + run + "_" + detector + "-MPC2.pdf"; 
	plotnamepng = "TimeDist_" + run + "_" + detector + "-MPC2.png"; 
	
	c2->SaveAs(plotnamepdf.c_str());	
	c2->SaveAs(plotnamepng.c_str());
	
	
	//plot amplitude
	std::string ampl = "amp_max[" + detector + "]>>hAmpl";
	TH1F *hAmpl = new TH1F("hAmpl","", 800, -50, 10000);
	h4->Draw(ampl.c_str());
	
	setStyle();
	TCanvas* c4 = new TCanvas("c4","c4");
	TH2F* H4 = new TH2F("H4","", 100, -50, (hAmpl->GetXaxis()->GetBinCenter(hAmpl->GetMaximumBin())*1.3), 64, 0, (hAmpl->GetMaximum()*1.1));     
	std::string xtit = "amp_max[" + detector + "]";        
	H4->GetXaxis()->SetTitle(xtit.c_str());    
	H4->GetYaxis()->SetTitle("events");
	FPCanvasStyle(c4);	

	H4->Draw();
	hAmpl->Draw("H4 SAME");

	std::string plotAmplpdf = "Amp_" + run + "_" + detector  + ".pdf";
	c4->SaveAs(plotAmplpdf.c_str());
	
	std::string plotAmplpng = "Amp_" + run + "_" + detector  + ".png";
	c4->SaveAs(plotAmplpng.c_str());


}
	
