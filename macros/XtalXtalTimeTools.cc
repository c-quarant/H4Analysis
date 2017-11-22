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
	DXvsX->GetXaxis()->SetTitle((axis+"[0]").c_str());    
 	DXvsX->GetYaxis()->SetTitle((axis+"[1]-"+axis+"[0]").c_str());

	//Drawing DeltaX
	SaveAsPngPdfRoot(DXvsX, "Delta"+axis+"/D"+axis+"vs"+axis+"_" + RunStats, "");

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
	
	return DXvsX->GetFunction("pol1")->GetParameter(0);
}

PlaneCoord XtalXtalTimeTools::GetXtalCenterEdge(std::string APD, std::string edge)
{
	HodoPlaneShift("X");
	HodoPlaneShift("Y");

	PlaneCoord Output;
	float fitRange=7;
	
	gStyle->SetOptStat();

	auto *AmpXavg = new TProfile("AmpXavg", "", 36, -18, 18, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 36, -18, 18, 0, 10000);

	//Filling Amplitude profile histograms Draw method
	std::string varexp = "fit_ampl[" + APD + "]:0.5*(X[0]+X[1]-(" + HodoShift.X_str + "))>>AmpXavg";	
	std::string Selection ="0.5*(Y[0]+Y[1]-("+HodoShift.Y_str+"))>-18 && 0.5*(Y[0]+Y[1]-("+HodoShift.Y_str+"))<18";
	h4->Draw(varexp.c_str(), Selection.c_str());
	AmpXavg->GetXaxis()->SetTitle("Xavg");    
	AmpXavg->GetYaxis()->SetTitle("amp_max");

	varexp = "fit_ampl[" + APD + "]:0.5*(Y[0]+Y[1]-(" + HodoShift.Y_str + "))>>AmpYavg";
	Selection = "0.5*(X[0]+X[1]-(" + HodoShift.X_str + "))>-18 && 0.5*(X[0]+X[1]-(" + HodoShift.X_str + "))<18";
	h4->Draw(varexp.c_str(), Selection.c_str());
	AmpYavg->GetXaxis()->SetTitle("Yavg");    
	AmpYavg->GetYaxis()->SetTitle("amp_max");

	//Fitting Amplitude profile histograms
	if(edge == "E") AmpXavg->Fit("pol2", "Q", "", -16, -5);
	else AmpXavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	if(edge=="N") AmpYavg->Fit("pol2", "Q", "", 0, 15);
	else if(edge=="S") AmpYavg->Fit("pol2", "Q", "", -15, 0);
	else AmpYavg->Fit("pol2", "Q", "", -fitRange, fitRange);
	
	//Drawing Xavg histogram
	SaveAsPngPdfRoot(AmpXavg, "Amplitude_profiles/pAVG/"+std::to_string(Energy)+"Gev/"+"AmpXAVG_"+std::to_string(Run)+"_"+APD+"_G"+std::to_string(Gain), "");
	SaveAsPngPdfRoot(AmpYavg, "Amplitude_profiles/pAVG/"+std::to_string(Energy)+"Gev/"+"AmpYAVG_"+std::to_string(Run)+"_"+APD+"_G"+std::to_string(Gain), "");

	TF1 *fitResX = AmpXavg->GetFunction("pol2");
	TF1* fitResY = AmpYavg->GetFunction("pol2");

	Output.X = fitResX->GetMaximumX();
	Output.Y = fitResY->GetMaximumX();
	Output.X_str = std::to_string(fitResX->GetMaximumX());
	Output.Y_str = std::to_string(fitResY->GetMaximumX());
	
	return Output;
}

void XtalXtalTimeTools::HodoPlaneMaps()
{	
	HodoPlaneShift("X");
	HodoPlaneShift("Y");
   	
  	//2DHist definition 
  	TH2F* HodoXY0 = new TH2F("HodoXY0","", 32, -17, 15, 32, -16, 16); 
	h4->Draw("Y[0]:X[0]>>HodoXY0", "X[0]>-800 && Y[0]>-800");
  	HodoXY0->GetXaxis()->SetTitle("X[0]");    
  	HodoXY0->GetYaxis()->SetTitle("Y[0]");
  	HodoXY0->GetZaxis()->SetTitle("events");

  	TH2F* HodoXY1 = new TH2F("HodoXY1","", 32, -17, 15, 32, -16, 16); 
  	h4->Draw(("Y[1]-("+HodoShift.Y_str+"):X[1]-("+HodoShift.X_str+")>>HodoXY1").c_str(), "X[1]>-800 && Y[1]>-800");
  	HodoXY1->GetXaxis()->SetTitle("X[1]");    
  	HodoXY1->GetYaxis()->SetTitle("Y[1]");
  	HodoXY1->GetZaxis()->SetTitle("events");

  	TH2F* HodoXYM = new TH2F("HodoXYM","", 32, -17, 15, 32, -16, 16); 
  	h4->Draw(("(0.5*(Y[0]+Y[1]-("+HodoShift.Y_str+"))):(0.5*(X[0]+X[1]-("+HodoShift.X_str+")))>>HodoXYM").c_str(), "X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800");	
  	HodoXYM->GetXaxis()->SetTitle("Xavg");    
  	HodoXYM->GetYaxis()->SetTitle("Yavg");
  	HodoXYM->GetZaxis()->SetTitle("events");

	//Drawing histograms
	gStyle->SetOptStat(0);
	SaveAsPngPdfRoot(HodoXY0, "HodoMaps/Plane0/AmpXY_" + RunStats, "COLZ");
	SaveAsPngPdfRoot(HodoXY1, "HodoMaps/Plane1/AmpXY_" + RunStats, "COLZ");
	SaveAsPngPdfRoot(HodoXYM, "HodoMaps/PlaneAVG/AmpXY_" + RunStats, "COLZ");

	HodoXY0->~TH2F();
	HodoXY1->~TH2F();
	HodoXYM->~TH2F();
}

void XtalXtalTimeTools::AmplitudeMapsEdge(std::string APD)
{	
	HodoPlaneShift("X");
	HodoPlaneShift("Y");
   	
  	//2DHist definition 
  	auto *AmpXY0 = new TProfile2D("AmpXY0","", 32, -17, 15, 32, -16, 16, 0, 10000); 
	h4->Draw(("fit_ampl["+APD+"]:Y[0]:X[0]>>AmpXY0").c_str(), "X[0]>-800 && Y[0]>-800");
  	AmpXY0->GetXaxis()->SetTitle("X[0]");    
  	AmpXY0->GetYaxis()->SetTitle("Y[0]");
  	AmpXY0->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());

  	auto *AmpXY1 = new TProfile2D("AmpXY1","", 32, -17, 15, 32, -16, 16, 0, 10000); 
  	h4->Draw(("fit_ampl["+APD+"]:Y[1]-("+HodoShift.Y_str+"):X[1]-("+HodoShift.X_str+")>>AmpXY1").c_str(), "X[1]>-800 && Y[1]>-800");
  	AmpXY1->GetXaxis()->SetTitle("X[1]");    
  	AmpXY1->GetYaxis()->SetTitle("Y[1]");
  	AmpXY1->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());

  	auto *AmpXYM = new TProfile2D("AmpXYM","", 32, -17, 15, 32, -16, 16, 0, 10000); 
  	h4->Draw(("fit_ampl["+APD+"]:(0.5*(Y[0]+Y[1]-("+HodoShift.Y_str+"))):(0.5*(X[0]+X[1]-("+HodoShift.X_str+")))>>AmpXYM").c_str(), "X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800");	
  	AmpXYM->GetXaxis()->SetTitle("Xavg");    
  	AmpXYM->GetYaxis()->SetTitle("Yavg");
  	AmpXYM->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());


	//Drawing histograms
	gStyle->SetOptStat(0);
	SaveAsPngPdfRoot(AmpXY0, "AmplitudeMaps/Plane0/AmpXY_" + APD + "_" + RunStats, "COLZ");
	SaveAsPngPdfRoot(AmpXY1, "AmplitudeMaps/Plane1/AmpXY_" + APD + "_" + RunStats, "COLZ");
	SaveAsPngPdfRoot(AmpXYM, "AmplitudeMaps/PlaneAVG/AmpXY_" + APD + "_" + RunStats, "COLZ");
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
	HAmp->Fit("gaus", "Q");	

	NoiseAmplitude[APD].Mean = HAmp->GetFunction("gaus")->GetParameter(1);
	NoiseAmplitude[APD].MeanErr = HAmp->GetFunction("gaus")->GetParError(1);
	NoiseAmplitude[APD].Sigma = HAmp->GetFunction("gaus")->GetParameter(2);
	NoiseAmplitude[APD].SigmaErr = HAmp->GetFunction("gaus")->GetParError(2);

	SaveAsPngPdfRoot(HAmp, "Amp_plot/NoiseAmp_"+APD+"_"+std::to_string((int)Gain), "");
	
	HAmp->~TH1F();
	h4Noise->~TTree();

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

	SummedAmplitude.Mean = SumAmp->GetFunction("gaus")->GetParameter(1);
	SummedAmplitude.MeanErr = SumAmp->GetFunction("gaus")->GetParError(1);
	SummedAmplitude.Sigma = SumAmp->GetFunction("gaus")->GetParameter(2);
	SummedAmplitude.SigmaErr = SumAmp->GetFunction("gaus")->GetParError(2);

	//Drawing and saving histogram as pdf, png and root file
	SaveAsPngPdfRoot(SumAmp, "AmplitudeSum/SumAmp_"+APDList[0]+"+"+APDList[1]+"_"+RunStats, "");
	SumAmp->~TH1F();

	return SummedAmplitude;
}

TH1F* XtalXtalTimeTools::AeffDistribution()
{
	for(auto& APD : APDList)
		NoiseAmplitudeDistributionFit(APD);

	//setting selections
	std::string PosSel;
	
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "((fabs(X[0])<800 && (X[0]<-8 || X[0]>-2)) || (fabs(X[1])<800 && (X[1]<-8 || X[1]>-2))) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else if(APDList[0]=="C2" || APDList[1]=="C2") 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>3) || (fabs(Y[1])<800 && fabs(Y[1])>3))";	
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>2) || (fabs(Y[1])<800 && fabs(Y[1])>2))";	
	/*
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && (fabs(Y[0])<800 || fabs(Y[1])<800)";			
	*/
	std::string SumAmpl	= "fit_ampl["+APDList[0]+"]+fit_ampl["+APDList[1]+"]";
	std::string MCPAmplSel	= "amp_max[MCP1]<1000";	

	std::string AmplSel 	= "("+SumAmpl+"-"+std::to_string(SummedAmplitude.Mean)+")<"+std::to_string(SummedAmplitude.Sigma)+" && fit_ampl["+APDList[0]+"]>fit_ampl["+APDList[1]+"] && fit_ampl["+APDList[0]+"]>0.1*"+std::to_string(SummedAmplitude.Mean)+"&& fit_ampl["+APDList[1]+"]>0.1*"+std::to_string(SummedAmplitude.Mean);
	std::string Aeff_Sel = PosSel+" && "+AmplSel;

	//Define and fill Aeff distribution histogram
	TH1F* Aeff = new TH1F("Aeff", "", 525, 0, 140);
	std::string Aeff_Sigma = "1/TMath::Sqrt( pow("+std::to_string(NoiseAmplitude[APDList[0]].Sigma)+"/fit_ampl["+APDList[0]+"], 2) + pow("+std::to_string(NoiseAmplitude[APDList[1]].Sigma)+"/fit_ampl["+APDList[1]+"], 2) )";
	h4->Draw((Aeff_Sigma+">>Aeff").c_str(), (Aeff_Sel).c_str());

	Aeff->GetXaxis()->SetRangeUser(0, 140);
	Aeff->GetXaxis()->SetTitle("Aeff/#sigma(Noise)");
	Aeff->GetYaxis()->SetTitle("events");
 
	//Drawing and saving Aeff distribution histogram
	SaveAsPngPdfRoot(Aeff, "EffectiveAmplitude/Aeff_"+APDList[0]+"_"+APDList[1]+"_"+RunStats, "");

	return Aeff;
}

std::vector<float>* XtalXtalTimeTools::AeffMeanDistribution(int NSlices, std::vector<float>* AeffMean, std::vector<float>* AeffMeanErr)
{
	if (AeffMean->size()!=0) return AeffMean;

	int i;
	TH1F* Aeff = AeffDistribution();
	int AeffNBins = Aeff->GetNbinsX();
	//Calculate mean value of Aeff Slice by Slice and storing it in a TH1
	TH1F *hAeffMean = new TH1F("hAeffMean", "", NSlices, 0, 140);

	float BinInASlice = AeffNBins/NSlices;
	//cout << endl << "AeffBins "<< AeffNBins << "  AeffMax " << Aeff->GetBinCenter(AeffNBins) << "   Bin in a Slice " << (float)AeffNBins/NSlices << endl; 
	for(i=0; i<=NSlices; i++)
	{
		Aeff->GetXaxis()->SetRange((int)i*BinInASlice, (int)(i+1)*BinInASlice);
		AeffMean->push_back(Aeff->GetMean());
		AeffMeanErr->push_back(Aeff->GetMeanError());
		hAeffMean->SetBinContent(i, Aeff->GetMean());
		hAeffMean->SetBinError(i, Aeff->GetMeanError());
		//cout << "Aeff Low " << Aeff->GetBinCenter(i*BinInASlice) << "  Aeff high " << Aeff->GetBinCenter((i+1)*BinInASlice) << " AeffMean " << Aeff->GetMean() << endl;
	}
	
	SaveAsPngPdfRoot(hAeffMean, "EffectiveAmplitude/AeffMean_"+APDList[0]+"_"+APDList[1]+"_"+RunStats, "");

	return AeffMean;
}

void XtalXtalTimeTools::TimeXTALvsXTALAeff()
{
	int NBins=35;

	//Preliminary calculus for setting selection
	for(auto& APD : APDList)
		NoiseAmplitudeDistributionFit(APD);

	cout << ">>>>> Drawing spectrum of the summed amplitude of two neighbour crystals..." << endl; 
	SummedAmplitudeFit();
	cout << ">>>>> Drawing spectrum of the effective amplitude of two neighbour crystals..." << endl;
	TH1F* Aeff = AeffDistribution();	
	AeffMeanDistribution(NBins, &SliceAmpOAeffMean[APDList[0]+"-"+APDList[1]], &SliceAmpOAeffMeanErr[APDList[0]+"-"+APDList[1]]);

	//setting selections
	std::string PosSel;
	
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "((fabs(X[0])<800 && (X[0]<-8 || X[0]>-2)) || (fabs(X[1])<800 && (X[1]<-8 || X[1]>-2))) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else if(APDList[0]=="C2" || APDList[1]=="C2") 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>3) || (fabs(Y[1])<800 && fabs(Y[1])>3))";	
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>2) || (fabs(Y[1])<800 && fabs(Y[1])>2))";	
	/*	
	
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && (fabs(Y[0])<800 || fabs(Y[1])<800)";			
	*/
	std::string SumAmpl	= "fit_ampl["+APDList[0]+"]+fit_ampl["+APDList[1]+"]";
	std::string MCPAmplSel	= "amp_max[MCP1]<1000";	

	std::string AmplSel 	= "("+SumAmpl+"-"+std::to_string(SummedAmplitude.Mean)+")<"+std::to_string(SummedAmplitude.Sigma)+" && fit_ampl["+APDList[0]+"]<fit_ampl["+APDList[1]+"] && fit_ampl["+APDList[0]+"]>0.1*"+std::to_string(SummedAmplitude.Mean)+"&& fit_ampl["+APDList[1]+"]>0.1*"+std::to_string(SummedAmplitude.Mean);
	std::string tD_XTAL_XTAL_Sel = PosSel+" && "+AmplSel;//+" && "+MCPAmplSel+" && fit_status["+APDList[0]+"]==0 && fit_status["+APDList[1]+"]==0";
	
	//define XTAL_XTAL time distribution wrt Aeff
	TH2F* tD_XTAL_XTAL_Aeff = new TH2F("tD_XTAL_XTAL_Aeff", "", NBins, 0, 140, 100, -3, 3);
	std::string Aeff_SigmaNoise = "1/TMath::Sqrt( pow("+std::to_string(NoiseAmplitude[APDList[0]].Sigma)+"/fit_ampl["+APDList[0]+"], 2) + pow("+std::to_string(NoiseAmplitude[APDList[1]].Sigma)+"/fit_ampl["+APDList[1]+"], 2) )";
	h4->Draw(("fit_time["+APDList[0]+"]-fit_time["+APDList[1]+"]:"+Aeff_SigmaNoise+">>tD_XTAL_XTAL_Aeff").c_str(), tD_XTAL_XTAL_Sel.c_str());

	tD_XTAL_XTAL_Aeff->GetXaxis()->SetTitle("Aeff/#sigma");
	tD_XTAL_XTAL_Aeff->GetYaxis()->SetTitle(("t_{"+APDList[0]+"}-t_{"+APDList[1]+"} (ns)").c_str());
	tD_XTAL_XTAL_Aeff->GetXaxis()->SetTitleSize(0.055);
	tD_XTAL_XTAL_Aeff->GetXaxis()->SetTitleOffset(0.75);
	tD_XTAL_XTAL_Aeff->GetYaxis()->SetTitleSize(0.055);
	tD_XTAL_XTAL_Aeff->GetYaxis()->SetTitleOffset(0.75);


	//plot Xtal-Xtal time distribution
	gStyle->SetOptStat(0);
	SaveAsPngPdfRoot(tD_XTAL_XTAL_Aeff, "TimeResolutionAeff/"+APDList[1]+"_"+APDList[0]+"/Time_"+APDList[0]+"-"+APDList[1]+"_"+RunStats, "COLZ");
	
	MyFitSlicesY(APDList[0], APDList[1], tD_XTAL_XTAL_Aeff, &SliceGaussFit[APDList[0]+"-"+APDList[1]]);
	FitTimeResolution(APDList[0], APDList[1]);
	TimeMeanvsAeff(APDList[0], APDList[1]);
}

TH1F* XtalXtalTimeTools::AmplitudeDistribution(std::string APD)
{
	for(auto& APD : APDList)
		NoiseAmplitudeDistributionFit(APD);

	//setting selections
	std::string PosSel;
	
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "((fabs(X[0])<800 && (X[0]<-5.5 || X[0]>0.5)) || (fabs(X[1])<800 && (X[1]<-5.5 || X[1]>0.5))) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else if(APDList[0]=="C2" || APDList[1]=="C2") 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>3) || (fabs(Y[1])<800 && fabs(Y[1])>3))";	
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>2) || (fabs(Y[1])<800 && fabs(Y[1])>2))";	
	/*
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && (fabs(Y[0])<800 || fabs(Y[1])<800)";		
	*/
	std::string MCPAmplSel	= "amp_max[MCP1]>100 && amp_max[MCP1]<1000";	
	std::string AmplSel	= "fit_ampl["+APD+"]>100";
	std::string hAmp_Sel 	= PosSel+" && "+AmplSel+" && "+MCPAmplSel;

	//Define and fill Aeff distribution histogram
	TH1F* hAmp = new TH1F("hAmp", "", 525, 0, 350);
	h4->Draw(("fit_ampl["+APD+"]/"+std::to_string(NoiseAmplitude[APD].Sigma)+">>hAmp").c_str(), (hAmp_Sel).c_str());
	hAmp->GetXaxis()->SetTitle("Amplitude/#sigma(Noise)");
	hAmp->GetYaxis()->SetTitle("events");
 
	//Drawing and saving Aeff distribution histogram
	SaveAsPngPdfRoot(hAmp, "Amp_plot/Amp_"+APD+"-MCP1_"+RunStats, "");

	return hAmp;
}

std::vector<float>* XtalXtalTimeTools::AmplitudeMeanDistribution(std::string APD, int NSlices, std::vector<float>* AmpMean, std::vector<float>* AmpMeanErr)
{
	if (AmpMean->size()!=0) return AmpMean;

	int i;
	TH1F* hAmp = AmplitudeDistribution(APD);
	int hAmpNBins = hAmp->GetNbinsX();
	//Calculate mean value of Aeff Slice by Slice and storing it in a TH1
	TH1F *hAmpMean = new TH1F("hAmpMean", "", NSlices, 0, hAmp->GetBinCenter(hAmpNBins)+0.5*hAmp->GetBinWidth(hAmpNBins));

	float BinInASlice = hAmpNBins/NSlices;
	cout << endl << "hAmpBins "<< hAmpNBins << "  hAmpMax " << hAmp->GetBinCenter(hAmpNBins) << "   Bin in a Slice " << (float)hAmpNBins/NSlices << endl; 
	for(i=0; i<=NSlices; i++)
	{
		hAmp->GetXaxis()->SetRange((int)i*BinInASlice, (int)(i+1)*BinInASlice);
		AmpMean->push_back(hAmp->GetMean());
		AmpMeanErr->push_back(hAmp->GetMeanError());
		hAmpMean->SetBinContent(i, hAmp->GetMean());
		hAmpMean->SetBinError(i, hAmp->GetMeanError());
		cout << "hAmp Low " << hAmp->GetBinCenter(i*BinInASlice) << "  hAmp high " << hAmp->GetBinCenter((i+1)*BinInASlice) << " hAmpMean " << hAmp->GetMean() << endl;
	}
	
	SaveAsPngPdfRoot(hAmpMean, "Amp_plot/AmpMean_"+APD+"-MCP1_"+RunStats, "");

	return AmpMean;
}

	
void XtalXtalTimeTools::TimeXTALvsMCP(std::string APD)
{
	int NBins=35;

	//Preliminary calculus for setting selection
	for(auto& APD : APDList)
		NoiseAmplitudeDistributionFit(APD);

	cout << ">>>>> Drawing spectrum of the amplitude of only one crystal..." << endl;
	TH1F* hAmp = AmplitudeDistribution(APD);
	AmplitudeMeanDistribution(APD, NBins, &SliceAmpOAeffMean[APD+"-MCP1"], &SliceAmpOAeffMeanErr[APD+"-MCP1"]);

	//setting selections
	std::string PosSel;
	
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "((fabs(X[0])<800 && (X[0]<-5.5 || X[0]>0.5)) || (fabs(X[1])<800 && (X[1]<-5.5 || X[1]>0.5))) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else if(APDList[0]=="C2" || APDList[1]=="C2") 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>3) || (fabs(Y[1])<800 && fabs(Y[1])>3))";	
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])>2) || (fabs(Y[1])<800 && fabs(Y[1])>2))";	

	/*
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<10 || fabs(Y[1])<10)";
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && (fabs(Y[0])<800 || fabs(Y[1])<800)";		
	*/
	std::string MCPAmplSel	= "amp_max[MCP1]>200 && amp_max[MCP1]<1000";
	std::string AmplSel 	= "fit_ampl["+APD+"]>200";
	std::string tD_XTAL_MCP_Sel = PosSel+" && "+AmplSel+" && "+MCPAmplSel+" && fit_status["+APDList[0]+"]==0 && fit_status["+APDList[1]+"]==0";
	
	//define XTAL_XTAL time distribution wrt Aeff
	TH2F* tD_XTAL_MCP = new TH2F("tD_XTAL_MCP", "", NBins, 0, hAmp->GetBinCenter(hAmp->GetNbinsX())+0.5*hAmp->GetBinWidth(1), 1000, -10, 10);
	h4->Draw(("fit_time["+APD+"]-time[MCP1]:fit_ampl["+APD+"]/"+std::to_string(NoiseAmplitude[APD].Sigma)+">>tD_XTAL_MCP").c_str(), tD_XTAL_MCP_Sel.c_str());

	tD_XTAL_MCP->GetXaxis()->SetTitle("Amplitude/#sigma");
	tD_XTAL_MCP->GetYaxis()->SetTitle(("t_{"+APD+"}-t_{MCP1} (ps)").c_str());
	tD_XTAL_MCP->GetXaxis()->SetTitleSize(0.055);
	tD_XTAL_MCP->GetXaxis()->SetTitleOffset(0.75);
	tD_XTAL_MCP->GetYaxis()->SetTitleSize(0.055);
	tD_XTAL_MCP->GetYaxis()->SetTitleOffset(0.75);

	MyFitSlicesY(APD, "MCP1", tD_XTAL_MCP, &SliceGaussFit[APD+"-MCP1"]);

	float timeMean 	= SliceGaussFit[APD+"-MCP1"][15].Mean;
	float timeWidth = SliceGaussFit[APD+"-MCP1"][15].Sigma;
	tD_XTAL_MCP->GetYaxis()->SetRangeUser(timeMean-5*timeWidth, timeMean+5*timeWidth);
	//plot Xtal-Xtal time distribution
	gStyle->SetOptStat(0);
	SaveAsPngPdfRoot(tD_XTAL_MCP, "TimeResolutionAeff/"+APDList[1]+"_"+APDList[0]+"/Time_"+APD+"-MCP1_"+RunStats, "COLZ");

	FitTimeResolution(APD, "MCP1");
	TimeMeanvsAeff(APD, "MCP1");
}

void XtalXtalTimeTools::PulseShape(std::string APD, std::string TimeRef, std::string MagOMin)
{
	int i;
	float AmpMean, AmpSigma;

	std::string PosSel;
	if(APDList[0]=="D3" || APDList[1]=="D3") 
		PosSel 		= "(fabs(X[0])<5 || fabs(X[1])<5) && (fabs(Y[0])<5 || fabs(Y[1])<5)";
	else 
		PosSel 		= "(fabs(X[0])<10 || fabs(X[1])<10) && ((fabs(Y[0])<800 && fabs(Y[0])"+MagOMin+"3) || (fabs(Y[1])<800 && fabs(Y[1])"+MagOMin+"3))";		

	std::string MCPAmplSel	= "amp_max[MCP1]>200 && amp_max[MCP1]<1000";
	std::string AmplSel 	= "fit_ampl["+APD+"]>200";
	std::string PulseShape_Sel = "WF_ch=="+APD+" && "+PosSel+" && "+AmplSel+" && "+MCPAmplSel+" && fit_status["+APDList[0]+"]==0 && fit_status["+APDList[1]+"]==0";

	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 2048, -204.8, 204.8 , 300, -0.5, 1.5, 0., 10000.);
	cout << "Drawing TH2 profile histogram of amplitude vs WF_val vs #delta t... " << endl;
        h4->Draw(("fit_ampl["+APD+"]:WF_val/amp_max["+APD+"]:WF_time-time["+TimeRef+"]>>p2D_amp_vs_time").c_str(), PulseShape_Sel.c_str());
	p2D_amp_vs_time->GetXaxis()->SetTitle("WF_time (ns)");
    	p2D_amp_vs_time->GetYaxis()->SetTitle(("WF_val/amp_max["+APD+"]").c_str());
    	p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");

	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 2048, -204.8, 204.8, 300, -0.5, 1.5);
	cout << "Drawing TH2 histogram of WF_val vs delta t... " << endl;
	h4->Draw(("WF_val/amp_max["+APD+"]:WF_time-time["+TimeRef+"]>> h2_amp_vs_time").c_str(), PulseShape_Sel.c_str());
    	h2_amp_vs_time->GetXaxis()->SetTitle("WF_time (ns)");
    	h2_amp_vs_time->GetYaxis()->SetTitle(("WF_val/amp_max["+APD+"]").c_str());
   	h2_amp_vs_time->GetZaxis()->SetTitle("amp_max");

	cout << "Drawing pulse shape... " << endl;
	TObjArray aSlices;
	h2_amp_vs_time->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	TProfile *waveForm = new TProfile(("XTAL_"+APD+"_"+RunStats+"_prof_").c_str(), "", 2048, -204.8, 204.8);
	waveForm = (TProfile*)aSlices[1];
    	waveForm->GetXaxis()->SetTitle("WF_time (ns)");
    	waveForm->GetYaxis()->SetTitle(("WF_val/amp_max["+APD+"]").c_str());

	waveForm->GetYaxis()->SetRangeUser(0, 1.05);
	waveForm->SetName(("XTAL_"+APD+"_"+RunStats+"_prof_").c_str());
    	    	
	gStyle->SetOptStat(0);

	SaveAsPngPdfRoot(h2_amp_vs_time, "PulseShapes/WFPulseShapes/WFPS_"+MagOMin+"_"+APD+"_"+RunStats+"_h2", "COLZ");
	SaveAsPngPdfRoot(p2D_amp_vs_time, "PulseShapes/WFPulseShapes/WFPS_"+MagOMin+"_"+APD+"_"+RunStats+"_Amp", "COLZ");
	SaveAsPngPdfRoot(waveForm, "PulseShapes/WFPulseShapes/WFPS_"+MagOMin+"_"+APD+"_"+RunStats+"_profile", "COLZ");
}

//Divide a TH2 in slices along Y (one slice for each bin), make a TH1 for each slcice and fit each one of them with gaus, storing each TH1 and fit	
std::vector<GaussPar>* XtalXtalTimeTools::MyFitSlicesY(std::string APD, std::string TimeRef, TH2F* Xtal_Xtal_Time, std::vector<GaussPar>* SliceGaussFit)
{      
	if (SliceGaussFit->size()!=0) return SliceGaussFit;

	GaussPar ResTemp;
	int i,j;
	int NBinsX = Xtal_Xtal_Time->GetNbinsX();
	int NBinsy = Xtal_Xtal_Time->GetNbinsY();

	gStyle->SetOptStat();
	gStyle->SetOptFit();

	gSystem->mkdir(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+APD+"-"+TimeRef+"_"+RunStats).c_str());
	gSystem->CopyFile("/afs/cern.ch/user/c/cquarant/www/Amp_plot/index.php", ("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/"+APD+"-"+TimeRef+"_"+RunStats+"/index.php").c_str());

	TH1D* hTemp;
	TH1D* hSlices[1000];

	std::string hSliceName;
	for(i=1; i<NBinsX+1; i++)
	{
		if(i<10)hSliceName	= "TimeResSlice_0"+std::to_string(i);
		else	hSliceName	= "TimeResSlice_"+std::to_string(i);
		hSlices[i]		= new TH1D(hSliceName.c_str(), "", 200, -2, 2);
		hSlices[i]		= Xtal_Xtal_Time->ProjectionY(hSliceName.c_str(), i, i, "e");
		hSlices[i]->GetXaxis()->SetTitle(("t_{"+APD+"}-t_{"+TimeRef+"} (ps)").c_str());
		hSlices[i]->GetYaxis()->SetTitle("events");	
		hSlices[i]->GetXaxis()->SetTitleSize(0.055);
		hSlices[i]->GetXaxis()->SetTitleOffset(0.75);
		hSlices[i]->GetYaxis()->SetTitleSize(0.055);
		hSlices[i]->GetYaxis()->SetTitleOffset(0.75);

		if(hSlices[i]->GetEntries()>20)
		{	
			hSlices[i]->Fit("gaus", "Q");
			ResTemp.Mean		= hSlices[i]->GetFunction("gaus")->GetParameter(1);
			ResTemp.MeanErr		= hSlices[i]->GetFunction("gaus")->GetParError(1); 
			ResTemp.Sigma		= hSlices[i]->GetFunction("gaus")->GetParameter(2);
			ResTemp.SigmaErr	= hSlices[i]->GetFunction("gaus")->GetParError(2); 
			hSlices[i]->GetXaxis()->SetRangeUser(ResTemp.Mean-6*ResTemp.Sigma, ResTemp.Mean+6*ResTemp.Sigma);
		}
		else
		{
			ResTemp.Mean		= 0;
			ResTemp.MeanErr		= 0; 
			ResTemp.Sigma		= 0;
			ResTemp.SigmaErr	= 1000;
		}
	
		SaveAsPngPdfRoot(hSlices[i], "TimeDistSlices/"+APD+"-"+TimeRef+"_"+RunStats+"/"+hSliceName+"_"+APD+"-"+TimeRef, "");
		SliceGaussFit->push_back(ResTemp);
		hSlices[i]->~TH1D();
	}

	return SliceGaussFit;	
}

void XtalXtalTimeTools::FitTimeResolution(std::string APD, std::string TimeRef)
{
	//plot Slice Resolution vs mean effective amplitude of the slices and fit it with FitFunction
	std::vector<float> SliceResolutions, SliceResolutionsErr;
	for(auto& SliceRes : SliceGaussFit[APD+"-"+TimeRef])
	{
		SliceResolutions.push_back(1000*SliceRes.Sigma);
		SliceResolutionsErr.push_back(1000*SliceRes.SigmaErr);
	}

	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + TMath::Sqrt(2)*[1]*[1] )", 15, 1000);

	fitFunc->SetParLimits(0, 10, 10000);
	fitFunc->SetParLimits(1, 0, 50);

	fitFunc->SetParameter(0, 3000);
	fitFunc->SetParameter(1, 50);

	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");

	float Xfirst = 0;
	float Xlast = 0;	
	for(unsigned int i=0; i<SliceResolutions.size(); i++)
	{	
		if(Xfirst==0 & SliceResolutions[i]>0) Xfirst = SliceAmpOAeffMean[APD+"-"+TimeRef][i]-0.5;
		if(SliceResolutions[i]>0) Xlast = SliceAmpOAeffMean[APD+"-"+TimeRef][i]+0.5;
	}
	cout << endl << Xfirst << endl;
	cout << Xlast << endl << endl;

	TGraphErrors *ResvsAeff = new TGraphErrors(SliceResolutions.size(), &SliceAmpOAeffMean[APD+"-"+TimeRef][0], &SliceResolutions[0], &SliceAmpOAeffMeanErr[APD+"-"+TimeRef][0], &SliceResolutionsErr[0]);
	ResvsAeff->GetXaxis()->SetRangeUser(Xfirst-10, Xlast+10);
	ResvsAeff->GetYaxis()->SetRangeUser(-1, *std::max_element(SliceResolutions.begin(), SliceResolutions.end())+15);


	ResvsAeff->GetXaxis()->SetTitle("Aeff/#sigma");
	ResvsAeff->GetYaxis()->SetTitle(("#sigma(t_{"+APD+"}-t_{"+TimeRef+"}) (ps)").c_str());	
	ResvsAeff->GetXaxis()->SetTitleSize(0.055);
	ResvsAeff->GetXaxis()->SetTitleOffset(0.75);
	ResvsAeff->GetYaxis()->SetTitleSize(0.055);
	ResvsAeff->GetYaxis()->SetTitleOffset(0.75);

	ResvsAeff->Fit("fitFunc", "", "", Xfirst, Xlast);

	SaveAsPngPdfRoot(ResvsAeff, "TimeResolutionAeff/"+APDList[1]+"_"+APDList[0]+"/TimeResvsAeff_"+APD+"-"+TimeRef+"_"+RunStats, "AP");
}

void XtalXtalTimeTools::TimeMeanvsAeff(std::string APD, std::string TimeRef)
{
	int NBins=35;
	std::vector<float> SliceTimeMean, SliceTimeMeanErr;
	for(auto& SliceTime : SliceGaussFit[APD+"-"+TimeRef])
	{
		SliceTimeMean.push_back(1000*SliceTime.Mean);
		SliceTimeMeanErr.push_back(1000*SliceTime.MeanErr);
	}
	
	TGraphErrors *TimevsAeff = new TGraphErrors(SliceTimeMean.size(), &SliceAmpOAeffMean[APD+"-"+TimeRef][0], &SliceTimeMean[0], &SliceAmpOAeffMeanErr[APD+"-"+TimeRef][0], &SliceTimeMeanErr[0]);
	TimevsAeff->GetXaxis()->SetRangeUser(0, 3500);
	//TimevsAeff->GetYaxis()->SetRangeUser(SliceTimeMean[7]-300, SliceTimeMean[7]+300);
	TimevsAeff->GetXaxis()->SetTitle("Aeff/#sigma");
	TimevsAeff->GetYaxis()->SetTitle(("#bar{t_{"+APD+"}-t_{"+TimeRef+"}} (ps)").c_str());	
	TimevsAeff->GetXaxis()->SetTitleSize(0.055);
	TimevsAeff->GetXaxis()->SetTitleOffset(0.75);
	TimevsAeff->GetYaxis()->SetTitleSize(0.055);
	TimevsAeff->GetYaxis()->SetTitleOffset(0.8);
	
	float Ylast=0;
	int ZeroCount=0;
	for(auto& SliceTime : SliceTimeMean)
	{
		if(SliceTime==0) ZeroCount++;
		Ylast += SliceTime;
	}
	Ylast /= SliceTimeMean.size()-ZeroCount; 

	TGaxis::SetMaxDigits(3);
	TH2F* Hset = new TH2F("Hset", "Hset", 100, 0, 350, 100, 0, 10000);
	Hset->GetYaxis()->SetRangeUser(Ylast-200, Ylast+200);
	Hset->GetXaxis()->SetTitle("Aeff/#sigma");
	Hset->GetYaxis()->SetTitle(("#bar{t_{"+APD+"}-t_{"+TimeRef+"}} (ps)").c_str());	
	Hset->GetXaxis()->SetTitleSize(0.055);
	Hset->GetXaxis()->SetTitleOffset(0.75);
	Hset->GetYaxis()->SetTitleSize(0.055);
	Hset->GetYaxis()->SetTitleOffset(0.75);

	SaveAsPngPdfRoot(Hset, TimevsAeff, "TimeMeanvsAeff/"+APDList[1]+"_"+APDList[0]+"/TimevsAeff_FarFromGap_"+APD+"_"+TimeRef+"_"+RunStats, "P");
}

void XtalXtalTimeTools::SaveAsPngPdfRoot(TObject* ToSave, std::string PathAndName, std::string DrawOpt)
{
	gStyle->SetOptFit();
	TCanvas* cOut = new TCanvas("cOut", "cOut");
	ToSave->Draw(DrawOpt.c_str());
	cOut->SaveAs(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".png").c_str());
	cOut->SaveAs(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".pdf").c_str());
	cOut->~TCanvas();

	TFile* f2 = new TFile(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".root").c_str(), "RECREATE");
	f2->cd();
	ToSave->Write();
	f2->Write();
	f2->Close();
	f2->~TFile();
}

void XtalXtalTimeTools::SaveAsPngPdfRoot(TObject* ToSave, TObject* ToSave1, std::string PathAndName, std::string DrawOpt)
{
	gStyle->SetOptFit();
	TCanvas* cOut = new TCanvas("cOut", "cOut");
	ToSave->Draw();
	ToSave1->Draw(("SAME"+DrawOpt).c_str());
	cOut->SaveAs(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".png").c_str());
	cOut->SaveAs(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".pdf").c_str());
	cOut->~TCanvas();

	TFile* f2 = new TFile(("/afs/cern.ch/user/c/cquarant/www/"+PathAndName+".root").c_str(), "RECREATE");
	f2->cd();
	ToSave->Write();
	ToSave1->Write();
	f2->Write();
	f2->Close();
	f2->~TFile();
}

