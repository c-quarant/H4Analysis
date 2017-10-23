#include "XtalXtalTimeTools.h"

XtalXtalTimeTools::XtalXtalTimeTools(std::string RunNtuple)
{
	TFile* f = TFile::Open(RunNtuple.c_str());
	h4 = (TTree*)f->Get("h4");
	Nentries = h4->GetEntries();

	h4->GetEntry(0);
	Run	= (int)h4->GetLeaf("run")->GetValue(0);
	Gain	= (int)h4->GetLeaf("CHGain")->GetValue(0);
	Energy	= (int)h4->GetLeaf("Energy")->GetValue(0);
	RunStats = std::to_string((int)Energy)+"Gev_G"+std::to_string((int)Gain)+"_"+std::to_string(Run);

	HodoShift.X	= -999; 
	HodoShift.Y	= -999;
	HodoShift.X_str	= std::to_string(-999); 
	HodoShift.Y_str	= std::to_string(-999);

	
	TTree* digi = (TTree*)f->Get("digi");
	TObjArray* BrList = digi->GetListOfBranches();
	std::string BrName;			
	for(int i=0; i<BrList->GetEntriesFast(); i++)
	{
		BrName = BrList->At(i)->GetName();
		if(BrName=="C0APD1" || BrName=="C0APD2" || BrName=="C2" || BrName=="C3" || BrName=="C4" || BrName=="D3")
			APDList.push_back(BrName);
	}

	for(auto const& APD : APDList)
	{
		Center[APD].X= -999; 
		Center[APD].Y= -999;
		GaussParInit(&NoiseAmplitude[APD]);
	}

	GaussParInit(&Default);
	GaussParInit(&SummedAmplitude);
}

void XtalXtalTimeTools::GaussParInit(GaussPar* gP)
{
	gP->Mean=-999; gP->MeanErr=-999; gP->Sigma=-999; gP->SigmaErr=-999;
}

//returns the shift of plane 1 from plane 0 of hodoscope along selected axis
float XtalXtalTimeTools::HodoPlaneShift(std::string axis)
{
	if(axis == "X" && HodoShift.X!=-999) return HodoShift.X;
	if(axis == "Y" && HodoShift.Y!=-999) return HodoShift.Y;
		
	//Filling DeltaX histogram
	auto *DXvsX = new TProfile(("D"+axis+"vs"+axis+"").c_str(), "", 32, -16, 16, -10, 10);
	h4->Draw(("("+axis+"[1]-"+axis+"[0]):"+axis+"[0]>>D"+axis+"vs"+axis+"").c_str(), (axis+"[0]>-800 && "+axis+"[1]>-800").c_str());
	DXvsX->Fit("pol1", "Q", "", -4, 4);

	TH2F* Hset = new TH2F(("Hset"+axis).c_str(),"", 32, -16, 16, 100, -5, 5);     
	Hset->GetXaxis()->SetTitle((axis+"[0]").c_str());    
 	Hset->GetYaxis()->SetTitle((axis+"[1]-"+axis+"[0]").c_str());
	
	//Fitting and Drawing DeltaX
	TCanvas *cDeltaX = new TCanvas(("cDelta"+axis+"").c_str(), ("cDelta"+axis+"").c_str());
  	Hset->Draw();	
	DXvsX->Draw("SAME");

	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Delta"+axis+"/D"+axis+"vs"+axis+"_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Delta"+axis+"/D"+axis+"vs"+axis+"_" + RunStats + ".png";
    	
  	cDeltaX->SaveAs(fileOutpdf.c_str(), "Q");
  	cDeltaX->SaveAs(fileOutpng.c_str(), "Q");

	if(axis == "X")
	{
		HodoShift.X	= DXvsX->GetFunction("pol1")->GetParameter(0);
		HodoShift.X_str	= std::to_string(HodoShift.X);
	}
	if(axis == "Y")
	{
		HodoShift.Y	= DXvsX->GetFunction("pol1")->GetParameter(0);
		HodoShift.Y_str	= std::to_string(HodoShift.Y);
	}

	cDeltaX->~TCanvas();
	Hset->~TH2F();
	
	return DXvsX->GetFunction("pol1")->GetParameter(0);
}

PlaneCoord XtalXtalTimeTools::GetXtalCenterEdge(std::string APD, std::string edge)
{
	HodoPlaneShift("X");
	HodoPlaneShift("Y");

	PlaneCoord Output;
	float fitRange=7;
	
	gStyle->SetOptStat();
	gStyle->SetOptFit();

	auto *AmpXavg = new TProfile("AmpXavg", "", 36, -18, 18, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 36, -18, 18, 0, 10000);

	//Filling Amplitude profile histograms Draw method
	std::string varexp = "fit_ampl[" + APD + "]:0.5*(X[0]+X[1]-(" + HodoShift.X_str + "))>>AmpXavg";	
	std::string Selection ="0.5*(Y[0]+Y[1]-("+HodoShift.Y_str+"))>-18 && 0.5*(Y[0]+Y[1]-("+HodoShift.Y_str+"))<18";
	h4->Draw(varexp.c_str(), Selection.c_str());

	varexp = "fit_ampl[" + APD + "]:0.5*(Y[0]+Y[1]-(" + HodoShift.Y_str + "))>>AmpYavg";
	Selection = "0.5*(X[0]+X[1]-(" + HodoShift.X_str + "))>-18 && 0.5*(X[0]+X[1]-(" + HodoShift.X_str + "))<18";
	h4->Draw(varexp.c_str(), Selection.c_str());

	TH2F* H1 = new TH2F("H1","", 36, -18, 18, 50, 0, (AmpXavg->GetMaximum())*1.25);
	H1->GetXaxis()->SetTitle("Xavg");    
	H1->GetYaxis()->SetTitle("amp_max");

	//Fitting Amplitude profile histograms
	if(edge == "E") AmpXavg->Fit("pol2", "Q", "", -16, -5);
	else AmpXavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	if(edge=="N") AmpYavg->Fit("pol2", "Q", "", 0, 15);
	else if(edge=="S") AmpYavg->Fit("pol2", "Q", "", -15, 0);
	else AmpYavg->Fit("pol2", "Q", "", -fitRange, fitRange);
	
	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
  
	H1->Draw();
	AmpXavg->Draw("SAME");

	TF1 *fitResX = AmpXavg->GetFunction("pol2");
  
	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+std::to_string(Energy)+"Gev/"+"AmpXAVG_"+std::to_string(Run)+"_"+APD+"_G"+std::to_string(Gain)+".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+std::to_string(Energy)+"Gev/"+"AmpXAVG_"+std::to_string(Run)+"_"+APD+"_G"+std::to_string(Gain)+".png";
  
	cout << "\n\nX Center Position = " << fitResX->GetMaximumX() << "\t";  
	
  	c5->SaveAs(fileOutpdf.c_str(), "Q");
  	c5->SaveAs(fileOutpng.c_str(), "Q");
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw();  
	AmpYavg->Draw("SAME");
	
	TF1* fitResY = AmpYavg->GetFunction("pol2");
	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+std::to_string(Energy)+"Gev/"+"AmpYAVG_"+std::to_string(Run)+"_"+APD+"_G"+std::to_string(Gain)+"_NoMCPSel.pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+std::to_string(Energy)+"Gev/"+"AmpYAVG_"+std::to_string(Run)+"_"+APD+"_G"+std::to_string(Gain)+"_NoMCPSel.png";

	cout << "Y Center Position = " << fitResY->GetMaximumX() << "\n" << endl;      

  	c6->SaveAs(fileOutpdf.c_str());
  	c6->SaveAs(fileOutpng.c_str());

	Output.X = fitResX->GetMaximumX();
	Output.Y = fitResY->GetMaximumX();
	Output.X_str = std::to_string(fitResX->GetMaximumX());
	Output.Y_str = std::to_string(fitResY->GetMaximumX());
	

	H1->~TH2F();
	c5->~TCanvas();
	c6->~TCanvas();

	return Output;
}

void XtalXtalTimeTools::AmplitudeMapsEdge(std::string APD, std::string edge)
{	
	HodoPlaneShift("X");
	HodoPlaneShift("Y");

   	
  	//2DHist definition 
  	auto *AmpXY0 = new TProfile2D("AmpXY0","", 32, -17, 15, 32, -16, 16, 0, 10000); 
	h4->Draw(("fit_ampl["+APD+"]:Y[0]:X[0]>>AmpXY0").c_str(), "X[0]>-800 && Y[0]>-800");
  	AmpXY0->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());

  	auto *AmpXY1 = new TProfile2D("AmpXY1","", 32, -17, 15, 32, -16, 16, 0, 10000); 
  	h4->Draw(("fit_ampl["+APD+"]:Y[1]-("+HodoShift.Y_str+"):X[1]-("+HodoShift.X_str+")>>AmpXY1").c_str(), "X[1]>-800 && Y[1]>-800");
  	AmpXY1->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());

  	auto *AmpXYM = new TProfile2D("AmpXYM","", 32, -17, 15, 32, -16, 16, 0, 10000); 
  	h4->Draw(("fit_ampl["+APD+"]:(0.5*(Y[0]+Y[1]-("+HodoShift.Y_str+"))):(0.5*(X[0]+X[1]-("+HodoShift.X_str+")))>>AmpXYM").c_str(), "X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800");	
  	AmpXYM->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());

	TH2F* H1 = new TH2F("H1","", 32, -17, 15, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");


	//Drawing p0 histogram
	gStyle->SetOptStat(0);
  	TCanvas* c1 = new TCanvas("c1","c1");
  	H1->Draw();	
  	AmpXY0->Draw("COLZ SAME");

  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + APD + "_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + APD + "_" + RunStats + ".png";
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  

	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
	gStyle->SetOptStat(0);
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + APD + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + APD + "_" + RunStats + ".png";

  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	

  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");
  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
	gStyle->SetOptStat(0);
  	H1->Draw();	
  	AmpXYM->Draw("COLZ SAME");

  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + APD + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + APD + "_" + RunStats + ".png";  
  	
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());
	
	H1->~TH2F();
	c1->~TCanvas();
	c2->~TCanvas();
	c3->~TCanvas();
}

//Draw Amplitude histogram and fit it with single gaus
GaussPar XtalXtalTimeTools::NoiseAmplitudeDistributionFit(std::string APD) 
{
	if(NoiseAmplitude[APD].Mean!=-999) return NoiseAmplitude[APD];

	gStyle->SetOptStat();	
	
	std::string NoiseNtuple;
	if(Gain == 50) NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5896.root";
	else if(Gain == 100) NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5894.root";
	else if(Gain == 200) NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5892.root";

	TFile* fNoise	= TFile::Open(NoiseNtuple.c_str());
	TTree* h4Noise 	= (TTree*)fNoise->Get("h4");
	TH1F* HAmp 	= new TH1F("HAmp", "", 45, -30, 60);
	
	if(APD=="C2" || APD=="D3") h4Noise->Draw("fit_ampl[C3]>>HAmp");
	else h4Noise->Draw(("fit_ampl["+APD+"]>>HAmp").c_str());
	HAmp->GetXaxis()->SetRangeUser(HAmp->GetBinCenter(HAmp->GetMaximumBin())-30, HAmp->GetBinCenter(HAmp->GetMaximumBin())+30);
	
	TCanvas* c0 = new TCanvas("c0", "c0");
	HAmp->Fit("gaus", "Q");	
	HAmp->Draw();
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/NoiseAmp_"+APD+"_"+std::to_string((int)Gain)+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/NoiseAmp_"+APD+"_"+std::to_string((int)Gain)+".pdf").c_str());
	
	NoiseAmplitude[APD].Mean = HAmp->GetFunction("gaus")->GetParameter(1);
	NoiseAmplitude[APD].MeanErr = HAmp->GetFunction("gaus")->GetParError(1);
	NoiseAmplitude[APD].Sigma = HAmp->GetFunction("gaus")->GetParameter(2);
	NoiseAmplitude[APD].SigmaErr = HAmp->GetFunction("gaus")->GetParError(2);
	
	c0->~TCanvas();
	HAmp->~TH1F();
	h4Noise->~TTree();
	fNoise->Close();
	fNoise->~TFile();

	return NoiseAmplitude[APD];
}

GaussPar XtalXtalTimeTools::SummedAmplitudeFit()
{
	if(SummedAmplitude.Mean!=Default.Mean) 
		return SummedAmplitude;

	//Define and fill histogram of the summed amplitude of two crystals 
	TH1F *SumAmp = new TH1F("SumAmp", "", 2000, 0, 5000);
	if(find(APDList.begin(), APDList.end(), "D3")!=APDList.end())
		h4->Draw(("fit_ampl["+APDList[0]+"]+fit_ampl["+APDList[1]+"]>>SumAmp").c_str(), "fabs(Y[0])<10 || fabs(Y[1])<10");
	else h4->Draw(("fit_ampl["+APDList[0]+"]+fit_ampl["+APDList[1]+"]>>SumAmp").c_str(), "fabs(X[0])<10 || fabs(X[1])<10");
	SumAmp->GetXaxis()->SetTitle(("fit_ampl["+APDList[0]+"]+fit_ampl["+APDList[1]+"]").c_str());
	SumAmp->GetYaxis()->SetTitle("events");

	//Fit histogram with gaus around maximum bin
	SumAmp->GetXaxis()->SetRangeUser(350, 5000);
	float XMax = SumAmp->GetBinCenter(SumAmp->GetMaximumBin());
	float Xfirst = XMax*0.6;
	float Xlast = XMax*1.2;
	SumAmp->Fit("gaus", "RQ", "", Xfirst, Xlast);
	//SumAmp->GetXaxis()->SetRangeUser(0, 5000);

	//Drawing and saving histogram as pdf, png and root file
	TCanvas* cAMP = new TCanvas("cAMP", "cAMP");
	cAMP->cd();
	SumAmp->Draw();
	cAMP->SaveAs(("/afs/cern.ch/user/c/cquarant/www/AmplitudeSum/SumAmp_"+APDList[0]+"+"+APDList[1]+"_"+RunStats+".png").c_str());
	cAMP->SaveAs(("/afs/cern.ch/user/c/cquarant/www/AmplitudeSum/SumAmp_"+APDList[0]+"+"+APDList[1]+"_"+RunStats+".pdf").c_str());

	TFile* f = new TFile(("/afs/cern.ch/user/c/cquarant/www/AmplitudeSum/SumAmp_"+APDList[0]+"+"+APDList[1]+"_"+RunStats+".root").c_str(), "RECREATE");
	f->cd();	
	SumAmp->Write();
	f->Write();
	f->Close();

	SummedAmplitude.Mean = SumAmp->GetFunction("gaus")->GetParameter(1);
	SummedAmplitude.MeanErr = SumAmp->GetFunction("gaus")->GetParError(1);
	SummedAmplitude.Sigma = SumAmp->GetFunction("gaus")->GetParameter(2);
	SummedAmplitude.SigmaErr = SumAmp->GetFunction("gaus")->GetParError(2);

	SumAmp->~TH1F();
	f->~TFile();
	cAMP->~TCanvas();

	return SummedAmplitude;
}
	
void XtalXtalTimeTools::TimeXTALvsXTALAeff()
{
	int NBins=30;
	if(Energy==50)NBins=15; 

	//Preliminary calculus for setting selection
	for(auto& APD : APDList)
		NoiseAmplitudeDistributionFit(APD);

	cout << ">>>>> Drawing spectrum of the summed amplitude of two neighbour crystals..." << endl; 
	SummedAmplitudeFit();
	cout << ">>>>> Drawing spectrum of the effective amplitude of two neighbour crystals..." << endl;
	TH1F* Aeff = AeffDistribution();	
	AeffMeanDistribution(NBins);

	//setting selections
	std::string PosSel;
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && (fabs(Y[0])<800 || fabs(Y[1])<800)";		

	std::string SumAmpl	= "fit_ampl["+APDList[0]+"]+fit_ampl["+APDList[1]+"]";

	std::string AmplSel;
	if(APDList[0]=="D3" || APDList[1]=="D3") 	
		AmplSel 	= "("+SumAmpl+"-"+std::to_string(SummedAmplitude.Mean)+")<"+std::to_string(SummedAmplitude.Sigma)+" && fit_ampl["+APDList[0]+"]<fit_ampl["+APDList[1]+"]";
	else
		AmplSel 	= "("+SumAmpl+"-"+std::to_string(SummedAmplitude.Mean)+")<"+std::to_string(SummedAmplitude.Sigma)+" && fit_ampl["+APDList[0]+"]>fit_ampl["+APDList[1]+"]";
	
	std::string tD_XTAL_XTAL_Sel = PosSel+" && "+AmplSel;
	
	//define XTAL_XTAL time distribution wrt Aeff
	TH2F* tD_XTAL_XTAL_Aeff = new TH2F("tD_XTAL_XTAL_Aeff", "", NBins, 0, Aeff->GetBinCenter(Aeff->FindLastBinAbove(10))+0.5*Aeff->GetBinWidth(1), 100, -3, 3);
	std::string Aeff_SigmaNoise = "1/TMath::Sqrt( pow("+std::to_string(NoiseAmplitude[APDList[0]].Sigma)+"/fit_ampl["+APDList[0]+"], 2) + pow("+std::to_string(NoiseAmplitude[APDList[1]].Sigma)+"/fit_ampl["+APDList[1]+"], 2) )";
	h4->Draw(("fit_time["+APDList[0]+"]-fit_time["+APDList[1]+"]:"+Aeff_SigmaNoise+">>tD_XTAL_XTAL_Aeff").c_str(), tD_XTAL_XTAL_Sel.c_str());
	tD_XTAL_XTAL_Aeff->GetXaxis()->SetTitle("Aeff/#sigma");
	tD_XTAL_XTAL_Aeff->GetYaxis()->SetTitle(SumAmpl.c_str());

	//plot Xtal-Xtal time distribution
	gStyle->SetOptStat(0);
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_XTAL_XTAL_Aeff->Draw("COLZ");
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/Time_"+APDList[0]+"-"+APDList[1]+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/Time_"+APDList[0]+"-"+APDList[1]+"_"+RunStats+".pdf").c_str());
	c0->~TCanvas();

	TFile* f = new TFile(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/Time_"+APDList[0]+"-"+APDList[1]+"_"+RunStats+".root").c_str(), "RECREATE");
	f->cd();
	tD_XTAL_XTAL_Aeff->Write();
	f->Write();
	f->Close();
	f->~TFile();

	//Divide TH2 in slices along Y and fit each slice with gaus, results stored in SliceResolutions
	MyFitSlicesY(tD_XTAL_XTAL_Aeff);
	
	//plot Slice Resolution vs mean effective amplitude of the slices and fit it with FitFunction
	std::vector<float> SliceResolutions, SliceResolutionsErr;
	for(auto& SliceRes : SliceGaussFit)
	{
		SliceResolutions.push_back(1000*SliceRes.Sigma);
		SliceResolutionsErr.push_back(1500*SliceRes.SigmaErr);
	}

	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + TMath::Sqrt(2)*[1]*[1] )", 15, 1000);

	fitFunc->SetParLimits(0, 10, 10000);
	fitFunc->SetParLimits(1, 0, 50);

	fitFunc->SetParameter(0, 5000);
	fitFunc->SetParameter(1, 33);

	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");

	
	TGraphErrors *ResvsAeff = new TGraphErrors(SliceResolutions.size(), &SliceAeffMean[0], &SliceResolutions[0], &SliceAeffMeanErr[0], &SliceResolutionsErr[0]);
	ResvsAeff->GetXaxis()->SetRangeUser(SliceAeffMean[5]-1, SliceAeffMean[NBins-1]+10);
	ResvsAeff->GetYaxis()->SetRangeUser(0, SliceResolutions[5]*1.05);	
	ResvsAeff->Fit("fitFunc", "", "", SliceAeffMean[10]-1, SliceAeffMean[NBins-2]+1);
	
	gStyle->SetOptFit();
	TCanvas* cResvsAeff = new TCanvas("ResvsAeff", "ResvsAeff");
	ResvsAeff->Draw("AP");
	cResvsAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimeResvsAeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".png").c_str());
	cResvsAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimeResvsAeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".pdf").c_str());
	cResvsAeff->~TCanvas();

	TFile* f1 = new TFile(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimeResvsAeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".root").c_str(), "RECREATE");
	f1->cd();
	ResvsAeff->Write();
	f1->Write();
	f1->Close();
	f1->~TFile();

	std::vector<float> SliceTimeMean, SliceTimeMeanErr;
	for(auto& SliceTime : SliceGaussFit)
	{
		SliceTimeMean.push_back(1000*SliceTime.Mean);
		SliceTimeMeanErr.push_back(1500*SliceTime.MeanErr);
	}
	
	TGraphErrors *TimevsAeff = new TGraphErrors(SliceTimeMean.size(), &SliceAeffMean[0], &SliceTimeMean[0], &SliceAeffMeanErr[0], &SliceTimeMeanErr[0]);
	TimevsAeff->GetXaxis()->SetRangeUser(SliceAeffMean[5]-1, SliceAeffMean[NBins-1]+10);
	TimevsAeff->GetYaxis()->SetRangeUser(*std::min_element(SliceTimeMean.begin()+5, SliceTimeMean.end())-15, *std::max_element(SliceTimeMean.begin()+5, SliceTimeMean.end())+15);	
	
	TCanvas* cTimevsAeff = new TCanvas("TimevsAeff", "TimevsAeff");
	TimevsAeff->Draw("AP");
	cTimevsAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimevsAeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".png").c_str());
	cTimevsAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimevsAeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".pdf").c_str());
	cTimevsAeff->~TCanvas();

	TFile* f2 = new TFile(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimevsAeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".root").c_str(), "RECREATE");
	f2->cd();
	TimevsAeff->Write();
	f2->Write();
	f2->Close();
	f2->~TFile();

}

TH1F* XtalXtalTimeTools::AeffDistribution()
{
	for(auto& APD : APDList)
		NoiseAmplitudeDistributionFit(APD);

	//setting selections
	std::string PosSel;
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && (fabs(Y[0])<800 || fabs(Y[1])<800)";		

	std::string SumAmpl	= "fit_ampl["+APDList[0]+"]+fit_ampl["+APDList[1]+"]";
	std::string AmplSel 	= "("+SumAmpl+"-"+std::to_string(SummedAmplitude.Mean)+")<"+std::to_string(SummedAmplitude.Sigma);

	std::string Aeff_Sel	= PosSel+" && "+AmplSel;

	//Define and fill Aeff distribution histogram
	TH1F* Aeff = new TH1F("Aeff", "", 400, 0, 140);
	std::string Aeff_Sigma = "1/TMath::Sqrt( pow("+std::to_string(NoiseAmplitude[APDList[0]].Sigma)+"/fit_ampl["+APDList[0]+"], 2) + pow("+std::to_string(NoiseAmplitude[APDList[1]].Sigma)+"/fit_ampl["+APDList[1]+"], 2) )";
	h4->Draw((Aeff_Sigma+">>Aeff").c_str(), (Aeff_Sel).c_str());

	Aeff->GetXaxis()->SetRangeUser(0, Aeff->GetBinCenter(Aeff->FindLastBinAbove(10)));
	Aeff->GetXaxis()->SetTitle("Aeff/#sigma(Noise)");
	Aeff->GetYaxis()->SetTitle("events");
 
	//Drawing and saving Aeff distribution histogram
	TCanvas *cH = new TCanvas("cH", "cH");
	Aeff->Draw();
	cH->SaveAs(("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/Aeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".png").c_str());
	cH->SaveAs(("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/Aeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".pdf").c_str());

	TFile* f = new TFile(("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/Aeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".root").c_str(), "RECREATE");
	f->cd();
	Aeff->Write();
	f->Write();
	f->Close();

	cH->~TCanvas();
	f->~TFile();

	return Aeff;
}

std::vector<float>* XtalXtalTimeTools::AeffMeanDistribution(int NSlices)
{
	if (SliceAeffMean.size()!=0) return &SliceAeffMean;

	int i;
	TH1F* Aeff = AeffDistribution();
	int AeffNBins = Aeff->FindLastBinAbove(10);
	//Calculate mean value of Aeff Slice by Slice and storing it in a TH1
	TH1F *hAeffMean = new TH1F("hAeffMean", "", NSlices, 0, Aeff->GetBinCenter(AeffNBins)+0.5*Aeff->GetBinWidth(AeffNBins));

	float BinInASlice = AeffNBins/NSlices;
	//cout << endl << "AeffBins "<< AeffNBins << "  AeffMax " << Aeff->GetBinCenter(AeffNBins) << "   Bin in a Slice " << (float)AeffNBins/NSlices << endl; 
	for(i=1; i<=NSlices; i++)
	{
		Aeff->GetXaxis()->SetRange((int)i*BinInASlice, (int)(i+1)*BinInASlice);
		SliceAeffMean.push_back(Aeff->GetMean());
		SliceAeffMeanErr.push_back(Aeff->GetMeanError());
		hAeffMean->SetBinContent(i, Aeff->GetMean());
		hAeffMean->SetBinError(i, Aeff->GetMeanError());
		//cout << "Aeff Low " << Aeff->GetBinCenter(i*BinInASlice) << "  Aeff high " << Aeff->GetBinCenter((i+1)*BinInASlice) << " AeffMean " << Aeff->GetMean() << endl;
	}
	
	TCanvas *cAeffMean = new TCanvas("cAeffMean", "cAeffMean");
	cAeffMean->cd();
	hAeffMean->Draw();
	cAeffMean->SaveAs(("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/AeffMean_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".png").c_str());
	cAeffMean->SaveAs(("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/AeffMean_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".pdf").c_str());
	cAeffMean->~TCanvas();

	TFile* fAeffMean = new TFile(("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/AeffMean_"+APDList[0]+"_"+APDList[1]+"_"+RunStats+".root").c_str(), "RECREATE");
	fAeffMean->cd();
	hAeffMean->Write();
	fAeffMean->Write();
	fAeffMean->Close();

	fAeffMean->~TFile();

	return &SliceAeffMean;
}

//Divide a TH2 in slices along Y (one slice for each bin), make a TH1 for each slcice and fit each one of them with gaus, storing each TH1 and fit	
std::vector<GaussPar>* XtalXtalTimeTools::MyFitSlicesY(TH2F* Xtal_Xtal_Time)
{      
	if (SliceGaussFit.size()!=0) return &SliceGaussFit;

	GaussPar ResTemp;
	int i,j;
	int NBinsX = Xtal_Xtal_Time->GetNbinsX();
	int NBinsy = Xtal_Xtal_Time->GetNbinsY();

	gStyle->SetOptStat();
	gStyle->SetOptFit();

	gSystem->mkdir(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+RunStats).c_str());
	gSystem->CopyFile("/afs/cern.ch/user/c/cquarant/www/Amp_plot/index.php", ("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+RunStats+"/index.php").c_str());

	TH1D* hTemp;
	TH1D* hSlices[1000];
	TCanvas* cOUT = new TCanvas("cOUT", "cOUT");
	TObjArray* f = new TObjArray();

	std::string hSliceName;
	for(i=1; i<NBinsX+1; i++)
	{
		if(i<10)hSliceName	= "TimeResSlice_0"+std::to_string(i);
		else	hSliceName	= "TimeResSlice_"+std::to_string(i);
		hSlices[i]		= new TH1D(hSliceName.c_str(), "", 200, -2, 2);
		hSlices[i]		= Xtal_Xtal_Time->ProjectionY(hSliceName.c_str(), i, i, "e");
		hSlices[i]->Fit("gaus", "Q");
		
		cOUT->cd();
		hSlices[i]->Draw();
		cOUT->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+RunStats+"/"+hSliceName+"_"+APDList[0]+"_"+APDList[1]+".pdf").c_str());
		cOUT->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+RunStats+"/"+hSliceName+"_"+APDList[0]+"_"+APDList[1]+".png").c_str());
		
		f->Add(TFile::Open(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+RunStats+"/"+hSliceName+"_"+APDList[0]+"_"+APDList[1]+".root").c_str(), "RECREATE"));	
		((TFile*)f->FindObject(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+RunStats+"/"+hSliceName+"_"+APDList[0]+"_"+APDList[1]+".root").c_str()))->cd();
		hSlices[i]->Draw();
		hSlices[i]->Write();
		((TFile*)f->FindObject(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+RunStats+"/"+hSliceName+"_"+APDList[0]+"_"+APDList[1]+".root").c_str()))->Close();

		ResTemp.Mean		= hSlices[i]->GetFunction("gaus")->GetParameter(1);
		ResTemp.MeanErr		= hSlices[i]->GetFunction("gaus")->GetParError(1); 
		ResTemp.Sigma		= hSlices[i]->GetFunction("gaus")->GetParameter(2);
		ResTemp.SigmaErr	= hSlices[i]->GetFunction("gaus")->GetParError(2); 

		SliceGaussFit.push_back(ResTemp);
		hSlices[i]->~TH1D();
	}

	//hSlice->~TH1D();
	//f->~TFile();
	cOUT->~TCanvas();

	return &SliceGaussFit;	
}
