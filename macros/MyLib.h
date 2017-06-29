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
#include <iostream>
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <string>
#include <fstream>
#include <math.h>

//returns the mean time shift between detector and MCP selected
float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOutput, std::string RunStats, std::string MCP)
{
	TH1F* raw_time_dist = new TH1F("raw_time_dist", "", 200, -3, 3);

	h4->Draw((std::string("time["+detector+"]-time["+MCP+"]>>raw_time_dist")).c_str(), Selection.c_str());

	raw_time_dist->GetXaxis()->SetTitle((std::string("time["+detector+"]-time["+MCP+"] (ns)")).c_str());
	raw_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* c0 = new TCanvas();
    	c0->cd();
	raw_time_dist->Fit("gaus", "", "", -3, 3);
    	raw_time_dist->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"RawTimeDistribution/RawTimeDist"+MCP+"_"+detector+"_"+RunStats+".png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"RawTimeDistribution/RawTimeDist"+MCP+"_"+detector+"_"+RunStats+".pdf").c_str());

	cout << raw_time_dist->GetMean() << endl;
	cout << raw_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return raw_time_dist->GetFunction("gaus")->GetMaximumX();
}

//returns the mean time measured by MCP selected
float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOutput, std::string RunStats, std::string MCP)
{
	TH1F* MCP_time_dist = new TH1F("MCP_time_dist", "", 200, 0, 50);

	h4->Draw((std::string("time["+MCP+"]>>MCP_time_dist")).c_str(), Selection.c_str());
	
	MCP_time_dist->GetXaxis()->SetTitle(("time["+MCP+"] (ns)").c_str());
	MCP_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* ca = new TCanvas();
    	ca->cd();
	MCP_time_dist->Fit("gaus", "", "", 0, 50);
    	MCP_time_dist->Draw();
    	ca -> SaveAs(std::string(pathToOutput+"MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".png").c_str());
    	ca -> SaveAs(std::string(pathToOutput+"MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".pdf").c_str());

	cout << MCP_time_dist->GetMean() << endl;
	cout << MCP_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return MCP_time_dist->GetFunction("gaus")->GetMaximumX();
}

//Draw Amplitude histogram and fit it with gaus
void AmplitudeHist(TTree *h4, std::string detector, std::string Selection, std::string pathToOutput,  std::string RunStats, float* AmpMean, float* AmpSigma)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -50, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.80), (int)(HAmp->GetMaximumBin()*1.10));

	HAmp->Fit("gaus");
	TCanvas* c1 = new TCanvas("c1", "c1");
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Draw();

	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOutput+"Amp_plot/Amp_"+detector+"_"+RunStats+".pdf").c_str());
	
	*AmpMean = HAmp->GetFunction("gaus")->GetParameter(1);
	*AmpSigma = HAmp->GetFunction("gaus")->GetParameter(2);
}

//returns the shift of plane 1 from plane 0 of hodoscope along selected axis
float HodoPlaneShift(TTree* h4, std::string detector, std::string pathToOutput, std::string RunStats, std::string axis)
{
	auto *DXvsX = new TProfile(("D"+axis+"vs"+axis+"").c_str(), "", 128, -16, 16, -10, 10);

	//Filling DeltaX histogram
	h4->Draw(("("+axis+"[0]-"+axis+"[1]):"+axis+"[0]>>D"+axis+"vs"+axis+"").c_str(), (axis+"[0]>-800 && "+axis+"[1]>-800").c_str());

	//Fitting and Drawing DeltaX
	TCanvas *cDeltaX = new TCanvas(("cDelta"+axis+"").c_str(), ("cDelta"+axis+"").c_str());
	TH2F* Hset = new TH2F(("Hset"+axis).c_str(),"", 128, -16, 16, 100, -5, 5);     
	Hset->GetXaxis()->SetTitle((axis+"[0]").c_str());    
 	Hset->GetYaxis()->SetTitle((axis+"[0]-"+axis+"[1]").c_str());
  
	Hset->Draw();	
	DXvsX->Fit("pol1", "Q", "", -4, 4);
	DXvsX->Draw("SAME");
	
	std::string fileOutpdf = pathToOutput + "Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".pdf";
  	std::string fileOutpng = pathToOutput + "Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".png";
    
  	cDeltaX->SaveAs(fileOutpdf.c_str());
  	cDeltaX->SaveAs(fileOutpng.c_str());

	return DXvsX->GetFunction("pol1")->GetParameter(0);
}

//Draw amplitude profiles in X and Y calculate maximum X
void AmplitudeProfilesFit(TTree *h4, std::string detector, std::string pathToOutput, std::string RunStats, Float_t bound, float* XMax, float* YMax)
{
	Int_t Nentries, i, CH, C3=0, C0APD1=0, C0APD2=0, hodoCfg;
	Float_t X1true, Y1true, Xshift, Yshift, Xavg, Yavg;
	Float_t Gain_val, Energy_val;
	float fitRange=5;
	
	h4->GetEntry(0);
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string runNum = std::to_string((int)h4->GetLeaf("run")->GetValue(0));

	gStyle->SetOptStat();
	gStyle->SetOptFit();

	auto *AmpXavg = new TProfile("AmpXavg", "", 128, -16, 16, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 128, -16, 16, 0, 10000);
	
	//Calculating shift of hodo plane 1 respect to hodo plane 0 faro X and Y axis
	Xshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "X");
	Yshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "Y");

	//Filling Amplitude profile histograms Draw method
	std::string varexp, selection;	
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);

	varexp = "amp_max[" + detector + "]:0.5*(X[0]+X[1]" + Xshift_str + ")>>AmpXavg";
	selection = "0.5*(Y[0]+Y[1]" + Yshift_str + ")>-" + bound_str + " && 0.5*(Y[0]+Y[1]" + Yshift_str + ")<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:0.5*(Y[0]+Y[1]" + Yshift_str + ")>>AmpYavg";
	selection = "0.5*(X[0]+X[1]" + Xshift_str + ")>-" + bound_str + " && 0.5*(X[0]+X[1]" + Xshift_str + ")<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());
	

	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
 	TH2F* H1 = new TH2F("H1","", 128, -16, 16, 50, 0, (AmpXavg->GetMaximum())*1.25);
 	H1->GetXaxis()->SetTitle("Xavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpXavg->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpXavg->Draw("SAME");

	TF1 *fitResX = AmpXavg->GetFunction("pol2");
  
	std::string fileOutpdf = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpXAVG_" + runNum + "_" + detector + "_G" + Gain + ".pdf";
  	std::string fileOutpng = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpXAVG_" + runNum + "_" + detector + "_G" + Gain + ".png";
  
  
  	c5->SaveAs(fileOutpdf.c_str());
  	c5->SaveAs(fileOutpng.c_str());
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw();  
	AmpYavg->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpYavg->Draw("SAME");
	
	TF1* fitResY = AmpYavg->GetFunction("pol2");
	
  	fileOutpdf = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpYAVG_" + runNum + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "Amplitude_profiles/pAVG/" + Energy + "Gev/" + "AmpYAVG_" + runNum + "_" + detector + "_G" + Gain + ".png";

    
  	c6->SaveAs(fileOutpdf.c_str());
  	c6->SaveAs(fileOutpng.c_str());

	*XMax = fitResX->GetMaximumX();
	*YMax = fitResY->GetMaximumX();
}

void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats, std::string runNum, std::string MCP)
{
	int i;
	float AmpMean, AmpSigma;

	std::string TimeShift;
	std::string TimeMCP;
	std::string AmpMean_str;
	std::string AmpSigma_str;

	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5);

	std::string Selection = "fabs(X[" + std::to_string(plane) + "]-(" + std::to_string(XMax) + "))<" + std::to_string(range) + " && fabs(Y[" + std::to_string(plane) + "]-(" + std::to_string(YMax) + "))<" + std::to_string(range) + " && amp_max["+MCP+"]>100";	

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<7";

	AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats, &AmpMean, &AmpSigma);
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	Selection = Selection + " && fabs(amp_max["+detector+"]-("+AmpMean_str+"))<5*"+AmpSigma_str;	
	
	TimeShift = std::to_string(MeanTimeShift(h4, detector, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = "WF_ch == " + detector + " && " + Selection;	
	
	cout << Selection << endl;
        h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+")>>p2D_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 1" << endl;
	
	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+") >> h2_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

	TObjArray aSlices;
	h2_amp_vs_time->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	
	TProfile *waveForm = new TProfile("waveForm", "", 300, -10, 40);
	waveForm = (TProfile*)aSlices[1];

	p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
    	h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
    	waveForm->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	h2_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	waveForm->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
   	h2_amp_vs_time->GetZaxis()->SetTitle("amp_max");	
    	    	
	gStyle->SetOptStat(0);
	
    	TCanvas* c1 = new TCanvas();
    	c1->cd();
    	h2_amp_vs_time->Draw("COLZ");
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_h2.png").c_str());
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+detector+"_h2.pdf").c_str());
	
    	TCanvas* c2 = new TCanvas();
    	c2->cd();
    	p2D_amp_vs_time->Draw("COLZ");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.png").c_str(),"png");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.pdf").c_str(),"pdf");
	
	TCanvas* c0 = new TCanvas();
    	c0->cd();
	waveForm->GetYaxis()->SetRangeUser(0, 1.05);
	waveForm->SetName("Waveform_");
    	waveForm->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.pdf").c_str());    	

    	TFile* output_Waveform = new TFile(std::string("WaveForms/"+detector+"_"+RunStats+"_Waveform.root").c_str(),"RECREATE");
    	output_Waveform->cd();
	
	p2D_amp_vs_time->SetName("Profile2DWaveforms_");
	h2_amp_vs_time->SetName("H2Waveforms_");

	p2D_amp_vs_time->Write();
    	h2_amp_vs_time->Write();
    	waveForm->Write();
    	output_Waveform->Close(); 
}


