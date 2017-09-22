#include "TFile.h"
#include "TStyle.h"
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
#include "TPaletteAxis.h"
#include "MyLib.h"
#include <string>
#include <fstream>
#include "Riostream.h"

void TimeResolutionAnalyzer(std::string NtupleList)
{
	ifstream in;
	in.open(NtupleList);
	int nline=0;
	TFile *f;
	TTree *h4;
	float XCentre, YCentre, Xshift, Yshift, AmpMean, AmpSigma, Xfirst, Xlast, TMax, TMaxValue;
	std::string detector, Gain, Energy, RunStats;
	std::string TimeMCP1, TimeMCP2, AmpMean_str, AmpSigma_str;
	std::string XCentre_str, YCentre_str, bound_str, Xshift_str, Yshift_str;
	std::string PosSel, MCP1Sel, MCP2Sel, AmpSel, tD_APD_MCP1_Sel, tD_APD_MCP2_Sel, tD_APD_MCP_Mean_Sel;
	GaussPar TimePar;

	TH1F* tD_APD_MCP1, tD_APD_MCP2, tD_APD_MCP_Mean, tD_APD_MCP1_Chi2, tD_APD_MCP2_Chi2;
	TCanvas* c0, c1, c2, c3, c_Chi2_1, c_Chi2_2;
	
	std::string path = (std::string)gSystem->UnixPathName(__FILE__);
	size_t found = path.find("macros/./TimeResScript.C");
	path.replace(found, std::string("macros/./TimeResScript.C").length(), "ntuples/");
	std::string ntuple="";


	while(ntuple!="END"){
		in >> ntuple;
		if(ntuple.find("C3")!=std::string::npos) detector = "C3";
		f = TFile::Open((path+ntuple).c_str());
		h4 = (TTree*)f->Get("h4");
	
		std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";
	
		h4->GetEntry(0);
		Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
		Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
		RunStats = Energy+"Gev_G"+Gain;

		//Shift between Hodo planes
		Xshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "X");
		Yshift = HodoPlaneShift(h4, detector, pathToOutput, RunStats, "Y");

		//Selection on distance of hitting position from the center of the detector
		AmplitudeProfilesFit(h4, detector, pathToOutput, RunStats, bound, &XCentre, &YCentre);

		XCentre_str = std::to_string(XCentre);
		YCentre_str = std::to_string(YCentre);
		bound_str = std::to_string(bound);
		Xshift_str = std::to_string(Xshift);
		Yshift_str = std::to_string(Yshift);

		PosSel = "(fabs(X[0]-("+XCentre_str+"))<"+bound_str+" || fabs(X[1]-("+XCentre_str+")-("+Xshift_str+"))<"+bound_str+") && (fabs(Y[0]-("+YCentre_str+"))<"+bound_str+" || fabs(Y[1]-("+Centre_str+")-("+Yshift_str+"))<"+bound_str+") && (fabs(X[0])<5 || fabs(X[1]-"+Xshift_str+")<5) && (fabs(Y[0])<5 || fabs(Y[1]-"+Yshift_str+")<5)";

		//Selection on MCP1 time & amplitude
		TimeMCP1 = std::to_string(MeanTimeMCP(h4, PosSel, pathToOutput+"fitTimeDist/", RunStats, "MCP1"));
		MCP1Sel = "amp_max[MCP1]>200 && amp_max[MCP1]<900 && fabs(time[MCP1]-("+TimeMCP1+"))<7";

		//Drawing Amplitude Maps
		AmplitudeMaps(h4, detector, MCP1Sel, Xshift_str, Yshift_str, XCentre, YCentre, bound, RunStats);	

		//Selection on Amplitude
		AmplitudeHist(h4, detector, PosSel+" && "+MCP1Sel, pathToOutput, RunStats, &AmpMean, &AmpSigma);
		AmpMean_str = std::to_string(AmpMean);
		AmpSigma_str = std::to_string(AmpSigma);	
		AmpSel = "fabs(amp_max["+detector+"]-("+AmpMean_str+"))<1*"+AmpSigma_str;	
	
		tD_APD_MCP1_Sel = PosSel + " && " + MCP1Sel + " && " + AmpSel;
		//cout << tD_APD_MCP1_Sel << endl;

		//Template fit Chi2 distribution
		tD_APD_MCP1_Chi2 = new TH1F("tD_APD_MCP1_Chi2", "", 50, 0, 8);
		h4->Draw(("fit_chi2["+detector+"]>>tD_APD_MCP1_Chi2").c_str(), tD_APD_MCP1_Sel.c_str());
	
		c_Chi2_1 = new TCanvas("c_Chi2", "c_Chi2");
		tD_APD_MCP1_Chi2->GetXaxis()->SetTitle("#chi^2");
		tD_APD_MCP1_Chi2->GetYaxis()->SetTitle("events");
		tD_APD_MCP1_Chi2->Draw();
		c_Chi2_1->SaveAs((pathToOutput+"fitTimeDist/TemplateChi2/Chi2_"+detector+"-MCP1_"+RunStats+".png").c_str());

		//define APD_MCP1 time distribution
		tD_APD_MCP1 = new TH1F("tD_APD_MCP1", "", 2000, -20, 20);
		if(Energy == "20" && Gain!="200") tD_APD_MCP1->SetBins(1000, -40, 40);
		if(Energy == "20" && Gain=="200") tD_APD_MCP1->SetBins(1500, -40, 40);
		h4->Draw(("fit_time["+detector+"]-time[MCP1]>>tD_APD_MCP1").c_str(), tD_APD_MCP1_Sel.c_str());
	
		Xfirst = tD_APD_MCP1->GetXaxis()->GetBinCenter(tD_APD_MCP1->GetMaximumBin())-1;
		Xlast = tD_APD_MCP1->GetXaxis()->GetBinCenter(tD_APD_MCP1->GetMaximumBin())+1;
		TMax = tD_APD_MCP1->GetBinCenter(tD_APD_MCP1->GetMaximumBin());
		TMaxValue = tD_APD_MCP1->GetMaximum();
	
		//plot and fit APD_MCP1 time distribution
		c0 = new TCanvas("c0", "c0");
		tD_APD_MCP1->GetXaxis()->SetRangeUser(Xfirst, Xlast);
		tD_APD_MCP1->GetXaxis()->SetTitle(("time["+detector+"]-time[MCP1] (ns)").c_str());
		tD_APD_MCP1->GetYaxis()->SetTitle("events");
		tD_APD_MCP1->Fit("gaus", "", "", Xfirst, Xlast);
		tD_APD_MCP1->Draw();
	
		c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP1_"+RunStats+".png").c_str());
		c0->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP1_"+RunStats+".pdf").c_str());

		//Drawing Time Maps MCP1
		TimeMaps(h4, detector, "MCP1", MCP1Sel+" && "+AmpSel, std::to_string(Xshift), std::to_string(Yshift), XMax, YMax, bound, RunStats, tD_APD_MCP1->GetFunction("gaus")->GetParameter(1), tD_APD_MCP1->GetFunction("gaus")->GetParameter(2));
	
	
	



		//Selection on MCP2 time && amplitude
		TimeMCP2 = std::to_string(MeanTimeMCP(h4, PosSel, pathToOutput+"fitTimeDist/", RunStats, "MCP2"));
		MCP2Sel = "amp_max[MCP2]>200 && amp_max[MCP2]<3500 && fabs(time[MCP2]-("+TimeMCP2+"))<7";

		tD_APD_MCP2_Sel = PosSel + " && " + MCP2Sel + " && " + AmpSel;
		//cout << tD_APD_MCP2_Sel << endl;
	
		//Drawing Time Maps MCP1
		TimeMaps(h4, detector, "MCP2", MCP2Sel+" && "+AmpSel, std::to_string(Xshift), std::to_string(Yshift), XMax, YMax, bound, RunStats,0,0);

		//Template fit Chi2 distribution
		tD_APD_MCP2_Chi2 = new TH1F("tD_APD_MCP2_Chi2", "", 50, 0, 8);
		h4->Draw(("fit_chi2["+detector+"]>>tD_APD_MCP2_Chi2").c_str(), tD_APD_MCP2_Sel.c_str());
	
		c_Chi2_2 = new TCanvas("c_Chi2_2", "c_Chi2_2");
		tD_APD_MCP2_Chi2->GetXaxis()->SetTitle("#chi^2");
		tD_APD_MCP2_Chi2->GetYaxis()->SetTitle("events");
		tD_APD_MCP2_Chi2->Draw();
		c_Chi2_2->SaveAs((pathToOutput+"fitTimeDist/TemplateChi2/Chi2_"+detector+"-MCP2_"+RunStats+".png").c_str());
	
		//plot and fit APD_MCP2 time distribution
		tD_APD_MCP2 = new TH1F("tD_APD_MCP2", "", 2000, -20, 20);
		if(Energy == "20" && Gain!="200") tD_APD_MCP2->SetBins(1500, -40, 40);
		if(Energy == "20" && Gain=="200") tD_APD_MCP2->SetBins(1500, -40, 40);
		h4->Draw(("fit_time["+detector+"]-time[MCP2]>>tD_APD_MCP2").c_str(), tD_APD_MCP2_Sel.c_str());

		//setting fit function	
		Xfirst = tD_APD_MCP2->GetXaxis()->GetBinCenter(tD_APD_MCP2->GetMaximumBin())-1;
		Xlast = tD_APD_MCP2->GetXaxis()->GetBinCenter(tD_APD_MCP2->GetMaximumBin())+1;
		TMax = tD_APD_MCP2->GetBinCenter(tD_APD_MCP2->GetMaximumBin());
		TMaxValue = tD_APD_MCP2->GetMaximum();

		c2 = new TCanvas("c2", "c2");
		tD_APD_MCP2->GetXaxis()->SetRangeUser(Xfirst, Xlast);
		tD_APD_MCP2->GetXaxis()->SetTitle(("time["+detector+"]-time[MCP2] (ns)").c_str());
		tD_APD_MCP2->GetYaxis()->SetTitle("events");
		tD_APD_MCP2->Fit("gaus", "", "", Xfirst, Xlast);
		tD_APD_MCP2->Draw();

		c2->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP2_"+RunStats+".png").c_str());
		c2->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP2_"+RunStats+".pdf").c_str());

		//Drawing Time Maps MCP2
		TimeMaps(h4, detector, "MCP2", MCP1Sel+" && "+AmpSel, std::to_string(Xshift), std::to_string(Yshift), XMax, YMax, bound, RunStats, tD_APD_MCP2->GetFunction("gaus")->GetParameter(1), tD_APD_MCP2->GetFunction("gaus")->GetParameter(2));




		//Selection on both MCP time & amplitude
		std::string tD_APD_MCP_Mean_Sel = PosSel + " && " + MCP1Sel + " && " + MCP2Sel + " && " + AmpSel;
		//cout << tD_APD_MCP_Mean_Sel
	
		//plot and fit APD_MCP_Mean time distribution
	
		tD_APD_MCP_Mean = new TH1F("tD_APD_MCP_Mean", "", 4000, -40, 40);
		if(Energy == "20" && Gain!="200") tD_APD_MCP_Mean->SetBins(750, -40, 40);
		if(Energy == "20" && Gain=="200") tD_APD_MCP_Mean->SetBins(1500, -40, 40);
		h4->Draw(("fit_time["+detector+"]-0.5*(time[MCP1]+time[MCP2])>>tD_APD_MCP_Mean").c_str(), tD_APD_MCP_Mean_Sel.c_str());
	
		Xfirst = tD_APD_MCP_Mean->GetXaxis()->GetBinCenter(tD_APD_MCP_Mean->GetMaximumBin())-1;
		Xlast = tD_APD_MCP_Mean->GetXaxis()->GetBinCenter(tD_APD_MCP_Mean->GetMaximumBin())+1;

		c3 = new TCanvas("c3", "c3");
		tD_APD_MCP_Mean->GetXaxis()->SetRangeUser(Xfirst, Xlast);
		tD_APD_MCP_Mean->GetXaxis()->SetTitle(("time["+detector+"]-timeMean (ns)").c_str());
		tD_APD_MCP_Mean->GetYaxis()->SetTitle("events");
		tD_APD_MCP_Mean->Fit("gaus", "", "", Xfirst, Xlast);
		tD_APD_MCP_Mean->Draw();

		c3->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP_Mean_"+RunStats+".png").c_str());
		c3->SaveAs((pathToOutput+"fitTimeDist/FinalTimeDistribution/Time_"+detector+"-MCP_Mean_"+RunStats+".pdf").c_str());
		}

	in.close();
}
