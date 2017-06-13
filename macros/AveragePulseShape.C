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
#include "FPCanvasStyle.C"
#include "setStyle.C"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "AmplitudeHist.C"
#include <string>
#include <fstream>
#include <math.h>

void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats);
float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOut, std::string RunStats);
float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOut, std::string RunStats);

void AveragePulseShape(std::string FileIn, std::string detector, Float_t bound)
{
	Int_t Nentries, i, CH, C3=0, C0APD1=0, C0APD2=0, runNum, hodoCfg;
	Float_t X[2], Y[2], Ampl[6];
	Float_t X1true, Y1true, Xshift, Yshift, Xavg, Yavg;
	Float_t Gain_val, Energy_val;

	float fitRange=5;
	//setStyle();
	gStyle->SetOptStat();
	gStyle->SetOptFit();

	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/";
	std::string pathToLinFitOutput = "/afs/cern.ch/user/c/cquarant/www/Delta";

	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	runNum = h4->GetLeaf("run")->GetValue(0);

	auto *DXvsX = new TProfile("DXvsX", "", 128, -16, 16, -10, 10);
	auto *DYvsY = new TProfile("DYvsY", "", 128, -16, 16, -10, 10);
	
	auto *AmpX0 = new TProfile("AmpX0", "", 32, -16, 16, 0, 10000);
	auto *AmpY0 = new TProfile("AmpY0", "", 32, -16, 16, 0, 10000);
	auto *AmpX1 = new TProfile("AmpX1", "", 32, -16, 16, 0, 10000);
	auto *AmpY1 = new TProfile("AmpY1", "", 32, -16, 16, 0, 10000);
	auto *AmpXavg = new TProfile("AmpXavg", "", 32, -16, 16, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 32, -16, 16, 0, 10000);
	
	std::string sel;
  	if(CH==0)
	{
		sel =  "C3Gain";
		C3 = 1;
	}
  	else if(CH==3)
	{
		sel = "C0Gain";
		C0APD1 = 1;
	}
	else if(CH==4)
	{
		sel = "C0Gain";
		C0APD2 = 1;
	}

	std::string Gain = std::to_string((int)h4->GetLeaf(sel.c_str())->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string fileOutpdf;
	std::string fileOutpng;
	
	//Filling DeltaX & DeltaY histograms
	h4->Draw("(X[0]-X[1]):X[0]>>DXvsX", "X[0]>-800 && X[1]>-800");
	h4->Draw("(Y[0]-Y[1]):Y[0]>>DYvsY", "Y[0]>-800 && Y[1]>-800");

	//Fitting and Drawing DeltaX e DeltaY
	TCanvas *cDeltaX = new TCanvas("cDeltaX", "cDeltaX");
	TH2F* Hset = new TH2F("Hset","", 128, -16, 16, 100, -5, 5);     
	Hset->GetXaxis()->SetTitle("X[0]");    
 	Hset->GetYaxis()->SetTitle("X[0]-X[1]");
  
	Hset->Draw();	
	DXvsX->Fit("pol1", "", "", -4, 4);
	DXvsX->Draw("SAME");
	
	TF1 *fitResX = DXvsX->GetFunction("pol1");
	Xshift = fitResX->GetParameter(0);
	
	fileOutpdf = pathToLinFitOutput + "X/DXvsX_" + std::to_string(runNum) + "_" + detector + "_" + Energy + "Gev.pdf";
  	fileOutpng = pathToLinFitOutput + "X/DXvsX_" + std::to_string(runNum) + "_" + detector + "_" + Energy + "Gev.png";
    
  	cDeltaX->SaveAs(fileOutpdf.c_str());
  	cDeltaX->SaveAs(fileOutpng.c_str());

	TCanvas *cDeltaY = new TCanvas("cDeltaY", "cDeltaY");
	Hset->GetXaxis()->SetTitle("Y[0]");    
 	Hset->GetYaxis()->SetTitle("Y[0]-Y[1]");
  
	Hset->Draw();	
  	DYvsY->Fit("pol1", "", "", -4, 4);
	DYvsY->Draw("SAME");

	TF1 *fitResY = DYvsY->GetFunction("pol1");
	Yshift = fitResY->GetParameter(0);

	fileOutpdf = pathToLinFitOutput + "Y/DYvsY_" + std::to_string(runNum) + "_" + detector + "_" + Energy + "Gev.pdf";
  	fileOutpng = pathToLinFitOutput + "Y/DYvsY_" + std::to_string(runNum) + "_" + detector + "_" + Energy + "Gev.png";
    
  	cDeltaY->SaveAs(fileOutpdf.c_str());
  	cDeltaY->SaveAs(fileOutpng.c_str());

	//Filling Amplitude profile histograms Draw method
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);

	std::string varexp = "amp_max[" + detector + "]:X[0]>>AmpX0";
	std::string selection = "Y[0]>-" + bound_str + " && Y[0]<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());
	
	varexp = "amp_max[" + detector + "]:Y[0]>>AmpY0"; 
	selection = "X[0]>-" + bound_str + " && X[0]<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:(X[1]" + Xshift_str + ")>>AmpX1";
	selection = "Y[1]" + Yshift_str + ">-" + bound_str + " && Y[1]" + Yshift_str + "<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:(Y[1]" + Yshift_str + ")>>AmpY1";
	selection = "X[1]" + Xshift_str + ">-" + bound_str + " && X[1]" + Xshift_str + "<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:0.5*(X[0]+X[1]" + Xshift_str + ")>>AmpXavg";
	selection = "0.5*(Y[0]+Y[1]" + Yshift_str + ")>-" + bound_str + " && 0.5*(Y[0]+Y[1]" + Yshift_str + ")<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:0.5*(Y[0]+Y[1]" + Yshift_str + ")>>AmpYavg";
	selection = "0.5*(X[0]+X[1]" + Xshift_str + ")>-" + bound_str + " && 0.5*(X[0]+X[1]" + Xshift_str + ")<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());
	

	//Drawing and fitting X0 histogram
	TCanvas* c1 = new TCanvas("c1","c1");
 	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 50, 0, (AmpX0->GetMaximum())*1.25);
	H1->GetXaxis()->SetTitle("X[0]");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw("AXIS");
	AmpX0->Fit("pol2", "", "", -fitRange, fitRange);
	AmpX0->Draw("SAME");

	fileOutpdf = pathToOutput + "p0/" + Energy + "Gev/" + "AmpX0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p0/" + Energy + "Gev/" + "AmpX0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
	
	
	//Drawing and fitting Y0 histogram
	TCanvas* c2 = new TCanvas("c2","c2");
 	H1->GetXaxis()->SetTitle("Y[0]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw("AXIS");	
	AmpY0->Fit("pol2", "", "", -fitRange, fitRange);	
	AmpY0->Draw("SAME");

	fileOutpdf = pathToOutput + "p0/" + Energy + "Gev/" + "AmpY0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p0/" + Energy + "Gev/" + "AmpY0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());


	//Drawing X1 histogram
	TCanvas* c3 = new TCanvas("c3","c3");
 	H1->GetXaxis()->SetTitle("X[1]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpX1->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpX1->Draw("SAME");
  
	fileOutpdf = pathToOutput + "p1/" + Energy + "Gev/" + "AmpX1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p1/" + Energy + "Gev/" + "AmpX1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());

	//Drawing Y1 histogram
	TCanvas* c4 = new TCanvas("c4","c4");
 	H1->GetXaxis()->SetTitle("Y[1]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpY1->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpY1->Draw("SAME");

	fileOutpdf = pathToOutput + "p1/" + Energy + "Gev/" + "AmpY1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p1/" + Energy + "Gev/" + "AmpY1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c4->SaveAs(fileOutpdf.c_str());
  	c4->SaveAs(fileOutpng.c_str());

	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
 	H1->GetXaxis()->SetTitle("Xavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpXavg->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpXavg->Draw("SAME");

	fitResX = AmpXavg->GetFunction("pol2");
  
	fileOutpdf = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpXAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpXAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c5->SaveAs(fileOutpdf.c_str());
  	c5->SaveAs(fileOutpng.c_str());
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw();  
	AmpYavg->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpYavg->Draw("SAME");
	
	fitResY = AmpYavg->GetFunction("pol2");
	
  	fileOutpdf = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpYAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpYAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
    
  	c6->SaveAs(fileOutpdf.c_str());
  	c6->SaveAs(fileOutpng.c_str());

	cout << fitResX->GetMaximumX() << "        " << fitResY->GetMaximumX() << endl;

	std::string RunStats = Energy+"Gev_G"+Gain+"_"+std::to_string(runNum);
	float XMax = fitResX->GetMaximumX();
	float YMax = fitResY->GetMaximumX();

	if(fabs(XMax)<3 && fabs(YMax)<3) PulseShapes(h4, detector, 0, XMax, YMax, 2, "/afs/cern.ch/user/c/cquarant/www/", RunStats);
	else cout<< "SBALLATO!" << endl;
	
}


void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats)
{
	std::string TimeShift;
	std::string TimeMCP;
	std::string AmpThresh;

	//TProfile2D* p2D_amp_vs_time_no_cut = new TProfile2D("p2D_amp_vs_time_no_cut", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	//TH2F* h2_amp_vs_time_no_cut = new TH2F("h2_amp_vs_time_no_cut", "", 300, -10, 40, 300, -0.5, 1.5);
	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5);

	std::string Selection = "fabs(X[" + std::to_string(plane) + "]-(" + std::to_string(XMax) + "))<" + std::to_string(range) + " && fabs(Y[" + std::to_string(plane) + "]-(" + std::to_string(YMax) + "))<" + std::to_string(range) + " && amp_max[MCP1]>100";	

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"PulseShapes/", RunStats));
	Selection = Selection + " && fabs(time[MCP1]-("+TimeMCP+"))<7";

	AmpThresh = std::to_string(0.5*AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats));	
	Selection = Selection + " && amp_max["+detector+"]>"+AmpThresh;	
	
	TimeShift = std::to_string(MeanTimeShift(h4, detector, Selection, pathToOutput+"PulseShapes/", RunStats));
	Selection = "WF_ch == " + detector + " && " + Selection;	

	
	cout << Selection << endl;
        h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+")>>p2D_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 1" << endl;
	
	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+") >> h2_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

    	TProfile* waveForm = h2_amp_vs_time->ProfileX(); 
    	waveForm->SetName("Waveform");

	/*
	//Drawing WITHOUT CUTS on HODOSCOPE
	Selection = "amp_max[MCP1]>200 && fabs(time[MCP1]-("+TimeMCP+"))<8 && X[0]>-800 && Y[0]>-800";

	h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+")>>p2D_amp_vs_time_no_cut").c_str(),Selection.c_str());
        //std::cout << (std::string("amp_max[")+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+") >> p2D_amp_vs_time") << std::endl;	
	//cout << Selection << endl;
	cout << "Draw 1" << endl;

	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time[MCP1]"+TimeShift+" >> h2_amp_vs_time_no_cut").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

    	TProfile* waveForm_no_cut = h2_amp_vs_time_no_cut->ProfileX(); 
    	waveForm_no_cut->SetName(std::string(detector+std::string("_waveform_prof_no_cut")).c_str());

	p2D_amp_vs_time_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	h2_amp_vs_time_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	waveForm_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	p2D_amp_vs_time_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	h2_amp_vs_time_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	waveForm_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	p2D_amp_vs_time_no_cut->GetZaxis()->SetTitle("amp_max");
   	h2_amp_vs_time_no_cut->GetZaxis()->SetTitle("amp_max");	
	*/
	p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
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
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_h2.png").c_str());
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_h2.pdf").c_str());
	
    	TCanvas* c2 = new TCanvas();
    	c2->cd();
    	p2D_amp_vs_time->Draw("COLZ");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_Amp.png").c_str(),"png");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_Amp.pdf").c_str(),"pdf");
	
    	TCanvas* c3 = new TCanvas();
    	c3->cd();
    	waveForm->Draw("P");
    	c3 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_profile.png").c_str(),"png");
    	c3 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_profile.pdf").c_str(),"pdf");
	/*
	TCanvas* c4 = new TCanvas();
    	c5->cd();
    	h2_amp_vs_time_no_cut->Draw("COLZ");
    	c5 -> SaveAs(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_h2_no_cut.png").c_str());
    	c5 -> SaveAs(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_h2_no_cut.pdf").c_str());

    	TCanvas* c5 = new TCanvas();
    	c6->cd();
    	p2D_amp_vs_time_no_cut->Draw("COLZ");
    	c6 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_Amp_no_cut.png").c_str(),"png");
    	c6 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_Amp_no_cut.pdf").c_str(),"pdf");

    	TCanvas* c6 = new TCanvas();
    	c7->cd();
    	waveForm_no_cut->Draw("P");
    	c7 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_profile_no_cut.png").c_str(),"png");
    	c7 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_profile_no_cut.pdf").c_str(),"pdf");
	*/
    	TFile* output_Waveform = new TFile(std::string("WaveForms/"+detector+"_"+RunStats+"_Waveform.root").c_str(),"RECREATE");
    	output_Waveform->cd();
	
	p2D_amp_vs_time->SetName("Profile2DWaveforms");
	h2_amp_vs_time->SetName("H2Waveforms");

	p2D_amp_vs_time->Write();
    	h2_amp_vs_time->Write();
    	waveForm->Write();
    	output_Waveform->Close(); 
}

float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOut, std::string RunStats)
{
	TH1F* raw_time_dist = new TH1F("raw_time_dist", "", 200, -3, 3);

	h4->Draw((std::string("time["+detector+"]-time[MCP1]>>raw_time_dist")).c_str(), Selection.c_str());

	raw_time_dist->GetXaxis()->SetTitle((std::string("time["+detector+"]-time[MCP1] (ns)")).c_str());
	raw_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* c0 = new TCanvas();
    	c0->cd();
	raw_time_dist->Fit("gaus", "", "", -3, 3);
    	raw_time_dist->Draw();
    	c0 -> SaveAs(std::string(pathToOut+"RawTimeDistribution/RawTimeDist_"+detector+"_"+RunStats+".png").c_str());
    	c0 -> SaveAs(std::string(pathToOut+"RawTimeDistribution/RawTimeDist_"+detector+"_"+RunStats+".pdf").c_str());

	cout << raw_time_dist->GetMean() << endl;
	cout << raw_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return raw_time_dist->GetFunction("gaus")->GetMaximumX();
}

float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOut, std::string RunStats)
{
	TH1F* MCP_time_dist = new TH1F("MCP_time_dist", "", 200, 0, 50);

	h4->Draw((std::string("time[MCP1]>>MCP_time_dist")).c_str(), Selection.c_str());
	
	MCP_time_dist->GetXaxis()->SetTitle("time[MCP1] (ns)");
	MCP_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* ca = new TCanvas();
    	ca->cd();
	MCP_time_dist->Fit("gaus", "", "", 0, 50);
    	MCP_time_dist->Draw();
    	ca -> SaveAs(std::string(pathToOut+"MCPTimeDistribution/MCP_TimeDist_"+RunStats+".png").c_str());
    	ca -> SaveAs(std::string(pathToOut+"MCPTimeDistribution/MCP_TimeDist_"+RunStats+".pdf").c_str());

	cout << MCP_time_dist->GetMean() << endl;
	cout << MCP_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return MCP_time_dist->GetFunction("gaus")->GetMaximumX();
}










