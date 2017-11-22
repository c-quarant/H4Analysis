#include "TimeAnalysisTools.cc"
#include "Riostream.h"
#include "TObjArray.h"

void FitRes(TTree* TimeRes, int SelDetector, int SelTimeRef);

void TimeResolutionAnalyzer(std::string NtupleList, int bound, std::string APD, std::string TimeRef_str)
{
	ifstream in;
	in.open(NtupleList);
	int i, j;
	int Run, Energy, Gain, detector, TimeRef;
	int C0APD1=0, C0APD2=1, C2=2, C3=3, C4=4, MCP1=5, MCP2=6, MCP_Mean=7, D3=8;
	float ANoiseR, ANoiseRErr;
	GaussPar TimePar, AeNoisePar, AeNoisePar2;

	std::map< std::string, int > ChannelID;
	
	ChannelID["C0APD1"] = C0APD1;
	ChannelID["C0APD2"] = C0APD2;
	ChannelID["C2"] = C2;
	ChannelID["C3"] = C3;
	ChannelID["C4"] = C4;
	ChannelID["D3"] = D3;
	ChannelID["MCP1"] = MCP1;
	ChannelID["MCP2"] = MCP2;
	ChannelID["MCP_Mean"] = MCP_Mean;
	
	std::string path = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/";
	/*
	size_t found = path.find("macros/./TimeResolutionAnalyzer.C");
	path.replace(found, std::string("macros/./TimeResolutionAnalyzer.C").length(), "ntuples/");
	*/
	std::string ntuple="";

	TFile* fOUT;
	fOUT = TFile::Open("TimeResults.root");
	TTree* tD;
	
	if(!fOUT){
		fOUT = new TFile("TimeResults.root", "RECREATE", "Stores main results of Analysis");
		tD = new TTree("tD", "Parameters of time distribution (fit)");

		tD->Branch("Run", &Run, "Run/I");
		tD->Branch("Energy", &Energy, "Energy/I");
		tD->Branch("Gain", &Gain, "Gain/I");
		tD->Branch("Detector", &detector, "Detector/I");
		tD->Branch("TimeReference", &TimeRef, "Time_reference/I");
	
		tD->Branch("C0APD1", &C0APD1, "C0APD1/I");
		tD->Branch("C0APD2", &C0APD2, "C0APD2/I");
		tD->Branch("C2", &C2, "C2/I");
		tD->Branch("C3", &C3, "C3/I");
		tD->Branch("C4", &C4, "C4/I");
		tD->Branch("D3", &D3, "D3/I");
		tD->Branch("MCP1", &MCP1, "MCP1/I");
		tD->Branch("MCP2", &MCP2, "MCP2/I");
		tD->Branch("MCP_Mean", &MCP_Mean, "MCP_Mean/I");

	
		tD->Branch("time_mean", &TimePar.Mean, "time_mean/F");
		tD->Branch("time_mean_error", &TimePar.MeanErr, "time_mean_error/F");
		tD->Branch("time_sigma", &TimePar.Sigma, "time_sigma/F");
		tD->Branch("time_sigma_error", &TimePar.SigmaErr, "time_sigma_error/F");
		tD->Branch("Amplitude_Noise_Ratio", &ANoiseR, "Amplitude_Noise_Ratio/F");
		tD->Branch("Amplitude_Noise_Ratio_error", &ANoiseRErr, "Amplitude_Noise_Ratio_error/F");
		
		TFile* f;
		TTree* h4;
		TTree* digi;
		TObjArray *BrList;
		TimeAnalysisTools* RunTimeAnalyzer;

		std::string BrName;
		std::vector< std::string > RunAPDList;
		std::vector< std::string > RunMCPList;
		while(ntuple!="END"){
			in >> ntuple;
			if(ntuple=="END") break;
			cout << "\n>>>> Processing file   " << ntuple << endl;
			
			f = TFile::Open((path+ntuple).c_str());
			h4 = (TTree*)f->Get("h4");
			digi = (TTree*)f->Get("digi");
			BrList = digi->GetListOfBranches();			
			
			RunAPDList.clear();
			RunMCPList.clear();			
			for(i=0; i<BrList->GetEntries(); i++)
			{
				BrName = BrList->At(i)->GetName();
				if(BrName=="C0APD1" || BrName=="C0APD2" || BrName=="C2" || BrName=="C3" || BrName=="C4" || BrName=="D3") RunAPDList.push_back(BrName);
				else if(BrName=="MCP1" || BrName=="MCP2") RunMCPList.push_back(BrName);
			}

			RunTimeAnalyzer = new TimeAnalysisTools(h4, RunAPDList, RunMCPList, bound);

			Run = RunTimeAnalyzer->runNum;
			Gain = stoi(RunTimeAnalyzer->Gain);
			Energy = stoi(RunTimeAnalyzer->Energy);
			
			for(i=0; i<(int)RunAPDList.size(); i++)
			{
				detector = ChannelID[RunAPDList[i]];
				RunTimeAnalyzer->AmplitudeMaps(RunAPDList[i]);
				
				for(auto const& MCP : RunMCPList) 
				{
					TimeRef = ChannelID[MCP];
					TimePar = RunTimeAnalyzer->TimeAPDvsMCP(RunAPDList[i], MCP);
					RunTimeAnalyzer->TimeMaps(RunAPDList[i], MCP);
					AeNoisePar = RunTimeAnalyzer->ComputeAvsNoise(RunAPDList[i]);

					//cout << endl << "A/Noise	:" << AeNoisePar.Mean/AeNoisePar.Sigma << " +/- ";
					//cout << TMath::Sqrt(pow(AeNoisePar.MeanErr/AeNoisePar.Mean, 2) + pow(AeNoisePar.SigmaErr/AeNoisePar.Sigma, 2))*AeNoisePar.Mean/AeNoisePar.Sigma << endl << endl;

					ANoiseR = AeNoisePar.Mean/AeNoisePar.Sigma;
					ANoiseRErr = TMath::Sqrt(pow(AeNoisePar.MeanErr/AeNoisePar.Mean, 2) + pow(AeNoisePar.SigmaErr/AeNoisePar.Sigma, 2))*AeNoisePar.Mean/AeNoisePar.Sigma;
					tD->Fill();
				}
				
				if(RunMCPList.size()>=2)
				{
					TimeRef = MCP_Mean;
					TimePar = RunTimeAnalyzer->TimeAPDvsMCPMean(RunAPDList[i]);
					tD->Fill();
				}
				
				for(j=i+1; j<(int)RunAPDList.size(); j++)
				{
					TimeRef = ChannelID[RunAPDList[j]];
					if(RunAPDList[i]=="C0APD1" || RunAPDList[i]=="C0APD2") TimePar = RunTimeAnalyzer->TimeAPDvsAPD(RunAPDList[i], RunAPDList[j]);
					else 
					{	
						//TimePar = RunTimeAnalyzer->TimeXTALvsXTAL(RunAPDList[i], RunAPDList[j]);
						//RunTimeAnalyzer->ComputeAvsNoiseEdge(RunAPDList[i], RunAPDList[j]);
						RunTimeAnalyzer->TimeXTALvsXTALAeff(RunAPDList[i], RunAPDList[j]);
					}
					AeNoisePar = RunTimeAnalyzer->ComputeAvsNoise(RunAPDList[i]);
					AeNoisePar2 = RunTimeAnalyzer->ComputeAvsNoise(RunAPDList[j]);

					ANoiseR = 1/TMath::Sqrt( pow(AeNoisePar.Sigma/AeNoisePar.Mean, 2) + pow(AeNoisePar2.Sigma/AeNoisePar2.Mean, 2) );
					ANoiseRErr = TMath::Sqrt( pow(AeNoisePar.Sigma*AeNoisePar.SigmaErr*ANoiseR*ANoiseR*ANoiseR/(AeNoisePar.Mean*AeNoisePar.Mean), 2) + pow(AeNoisePar2.Sigma*AeNoisePar2.SigmaErr*ANoiseR*ANoiseR*ANoiseR/(AeNoisePar2.Mean*AeNoisePar2.Mean), 2) + pow( AeNoisePar.Sigma*AeNoisePar.Sigma*AeNoisePar.MeanErr*ANoiseR*ANoiseR*ANoiseR/(AeNoisePar.Mean*AeNoisePar.Mean*AeNoisePar.Mean), 2) + pow( AeNoisePar2.Sigma*AeNoisePar2.Sigma*AeNoisePar2.MeanErr*ANoiseR*ANoiseR*ANoiseR/(AeNoisePar2.Mean*AeNoisePar2.Mean*AeNoisePar2.Mean), 2) ); 
					
					
					//cout << endl << endl << "A/Noise  " << ANoiseR << " +/- " << ANoiseRErr << endl << endl;
					tD->Fill();
				}
			}
			
			RunTimeAnalyzer->~TimeAnalysisTools();
		}
	fOUT->cd();
	tD->Write();
	}
	
	else
	{
		tD = (TTree*)fOUT->Get("tD");
		
		tD->SetBranchAddress("Run", &Run);
		tD->SetBranchAddress("Energy", &Energy);
		tD->SetBranchAddress("Gain", &Gain);
		tD->SetBranchAddress("Detector", &detector);
		tD->SetBranchAddress("TimeReference", &TimeRef);
	
		tD->SetBranchAddress("C0APD1", &C0APD1);
		tD->SetBranchAddress("C0APD2", &C0APD2);
		tD->SetBranchAddress("C2", &C2);
		tD->SetBranchAddress("C3", &C3);
		tD->SetBranchAddress("C4", &C4);
		tD->SetBranchAddress("D3", &D3);
		tD->SetBranchAddress("MCP1", &MCP1);
		tD->SetBranchAddress("MCP2", &MCP2);
		tD->SetBranchAddress("MCP_Mean", &MCP_Mean);

	
		tD->SetBranchAddress("time_mean", &TimePar.Mean);
		tD->SetBranchAddress("time_mean_error", &TimePar.MeanErr);
		tD->SetBranchAddress("time_sigma", &TimePar.Sigma);
		tD->SetBranchAddress("time_sigma_error", &TimePar.SigmaErr);
		tD->SetBranchAddress("Amplitude_Noise_Ratio", &ANoiseR);
		tD->SetBranchAddress("Amplitude_Noise_Ratio_error", &ANoiseRErr);
	}
	
	if(APD=="ALL")
	{
		std::vector< int > DetectorList;
		std::map< int, std::vector<int> > DetectorTimeRefList;
		for(i=0; i<tD->GetEntriesFast(); i++)
		{	
			tD->GetEntry(i);
			if(find(DetectorList.begin(), DetectorList.end(), detector) == DetectorList.end())
				DetectorList.push_back(detector);
						
			if(find(DetectorTimeRefList[detector].begin(), DetectorTimeRefList[detector].end(), TimeRef) == DetectorTimeRefList[detector].end())
				DetectorTimeRefList[detector].push_back(TimeRef);
		}
		for(auto& Det : DetectorList)
			for(auto& TRef : DetectorTimeRefList[Det])
				FitRes(tD, Det, TRef);	
	}
	else if(TimeRef_str == "ALL")	
	{
		std::vector< int > TimeRefList;
		for(i=0; i<tD->GetEntriesFast(); i++)
		{
			tD->GetEntry(i);
			if(detector==ChannelID[APD] && find(TimeRefList.begin(), TimeRefList.end(), TimeRef) == TimeRefList.end())
				TimeRefList.push_back(TimeRef);
		}

		if(TimeRefList.size()==0)
			cout << "detector " << APD << " not found in ntuples from file " << NtupleList << endl;
		for(auto& TRef : TimeRefList)
			FitRes(tD, ChannelID[APD], TRef);
	}
	else FitRes(tD, ChannelID[APD], ChannelID[TimeRef_str]);
	
	fOUT->Close();
}


void FitRes(TTree* TimeRes, int SelDetector, int SelTimeRef){
	gStyle->SetOptFit();	

	int i=0, Nentries=TimeRes->GetEntries();
	int EntryRun, EntryGain, EntryDetector, EntryTimeRef, EntryEnergy;
	float ANoiseRatio, ANoiseRatioErr, Res, ResErr;
	std::vector<Float_t> X0, X0Err, Y0, Y0Err, X1, X1Err, Y1, Y1Err, X2, X2Err, Y2, Y2Err;
	std::vector<int> ExistGains;
	
	std::vector<std::string> channel;
	channel.push_back("C0APD1");
	channel.push_back("C0APD2");
	channel.push_back("C2");
	channel.push_back("C3");
	channel.push_back("C4");
	channel.push_back("MCP1");
	channel.push_back("MCP2");
	channel.push_back("MCP_Mean");
	channel.push_back("D3");

	TimeRes->SetBranchAddress("Run", &EntryRun);
	TimeRes->SetBranchAddress("Detector", &EntryDetector);
	TimeRes->SetBranchAddress("TimeReference", &EntryTimeRef);
	TimeRes->SetBranchAddress("Gain", &EntryGain);
	TimeRes->SetBranchAddress("Energy", &EntryEnergy);
	TimeRes->SetBranchAddress("time_sigma", &Res);
	TimeRes->SetBranchAddress("time_sigma_error", &ResErr);
	TimeRes->SetBranchAddress("Amplitude_Noise_Ratio", &ANoiseRatio);
	TimeRes->SetBranchAddress("Amplitude_Noise_Ratio_error", &ANoiseRatioErr);

	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + TMath::Sqrt(2)*[1]*[1] )", 15, 1000);

	fitFunc->SetParLimits(0, 10, 10000);
	fitFunc->SetParLimits(1, 0, 50);
	//fitFunc->SetParLimits(2, 0, 2000);  

	fitFunc->SetParameter(0, 5000);
	fitFunc->SetParameter(1, 33);
	//fitFunc->SetParameter(2, 0);

	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");
	//fitFunc->SetParName(2, "~stochastic");
	
		
	for(i=0; i<Nentries; i++)
	{
		TimeRes->GetEntry(i);
		if(EntryDetector == SelDetector && find(ExistGains.begin(), ExistGains.end(), EntryGain) == ExistGains.end())
		{
			ExistGains.push_back(EntryGain);
		}
		if(EntryDetector == SelDetector && EntryTimeRef == SelTimeRef && EntryGain == 50 && !(EntryDetector==2 && EntryTimeRef==5 && (EntryRun==5603 || EntryRun==5614 || EntryRun==5621 || EntryRun==5640 || EntryRun==5634)))
		{
			X0.push_back(ANoiseRatio*1.2);
			X0Err.push_back(ANoiseRatioErr);
			Y0.push_back(Res*1000);
			Y0Err.push_back(ResErr*2500);
		}
		if(EntryDetector == SelDetector && EntryTimeRef == SelTimeRef && EntryGain == 100)
		{
			X1.push_back(ANoiseRatio);
			X1Err.push_back(ANoiseRatioErr);
			Y1.push_back(Res*1000);
			Y1Err.push_back(ResErr*2500);
		}
		if(EntryDetector == SelDetector && EntryTimeRef == SelTimeRef && EntryGain == 200)
		{
			X2.push_back(ANoiseRatio);
			X2Err.push_back(ANoiseRatioErr);
			Y2.push_back(Res*1000);
			Y2Err.push_back(ResErr*2500);
		}
	}

	std::string TimeRef_str;
	std::string Detector_str;
	std::string Info;

	TGraphErrors* g50;
	TGraphErrors* g100;
	TGraphErrors* g200;

	gStyle->SetOptFit();
	if(find(ExistGains.begin(), ExistGains.end(), 50) != ExistGains.end())
	{
		TCanvas* c0 = new TCanvas("c0", "c0");
		c0->cd();
	
		g50 = new TGraphErrors(X0.size(), &X0[0], &Y0[0], &X0Err[0], &Y0Err[0]);
		g50->SetTitle("");
	
		g50->GetXaxis()->SetTitle("A/#sigma(Noise)");
		g50->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
		g50->SetMarkerStyle(kFullCircle);
		g50->SetMarkerSize(1);
		fitFunc->SetLineColor(kBlack);	
		g50->Fit("fitFunc");
		g50->Draw("AP");
	
		TimeRef_str = channel[SelTimeRef];
		Detector_str = channel[SelDetector];
		Info = Detector_str + "-" + TimeRef_str + "_G50";

		c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".png").c_str());
		c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".pdf").c_str());
	}
	if(find(ExistGains.begin(), ExistGains.end(), 100) != ExistGains.end())
	{
		TCanvas* c1 = new TCanvas("c1", "c1");
		c1->cd();

		g100 = new TGraphErrors(X1.size(), &X1[0], &Y1[0], &X1Err[0], &Y1Err[0]);
		g100->SetTitle("");
			
		g100->GetXaxis()->SetTitle("A/#sigma(Noise)");
		g100->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
		g100->SetLineColor(kBlue);
		g100->SetMarkerStyle(kFullSquare);
		g100->SetMarkerSize(1);
		g100->SetMarkerColor(kBlue);
		fitFunc->SetLineColor(kBlue);
		g100->Fit("fitFunc");
	       	g100->Draw("AP");
	
		Info = Detector_str + "-" + TimeRef_str + "_G100";
	
		c1->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".png").c_str());
		c1->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".pdf").c_str());
	}
	if(find(ExistGains.begin(), ExistGains.end(), 200) != ExistGains.end())
	{
		TCanvas* c2 = new TCanvas("c2", "c2");
		c2->cd();
	
		g200 = new TGraphErrors(X2.size(), &X2[0], &Y2[0], &X2Err[0], &Y2Err[0]);
		g200->SetTitle("");	
	
		g200->GetXaxis()->SetTitle("A/#sigma(Noise)");
		g200->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
		g200->SetLineColor(kViolet);
		g200->SetMarkerStyle(kFullTriangleDown);
		g200->SetMarkerSize(1);
		g200->SetMarkerColor(6);
		fitFunc->SetLineColor(kViolet);	
		g200->Fit("fitFunc");
	
        	g200->Draw("AP");
	
		Info = Detector_str + "-" + TimeRef_str + "_G200";
	
		c2->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".png").c_str());
		c2->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Info+".pdf").c_str());
	}
	gStyle->SetOptFit(0);
	if(ExistGains.size()>1)
	{
		
		TCanvas* c3 = new TCanvas("c3", "c3");
		gStyle->SetOptFit(0);
		if(find(ExistGains.begin(), ExistGains.end(), 200) != ExistGains.end())
		{
			g200->GetYaxis()->SetRangeUser(0, (*max_element(Y0.begin(),Y0.end()))*1.1);
			g200->Draw("AP");
			g100->Draw("PSAME");
			g50->Draw("PSAME");
		}
		else if(find(ExistGains.begin(), ExistGains.end(), 100) != ExistGains.end())
		{
			g100->GetYaxis()->SetRangeUser(0, (*max_element(Y0.begin(),Y0.end()))*1.1);
			g100->Draw("AP");
			g50->Draw("PSAME");
		}
		else g50->Draw("AP");
		c3->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Detector_str+"-"+TimeRef_str+"_AllGain.png").c_str());
		c3->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_"+Detector_str+"-"+TimeRef_str+"_AllGain.pdf").c_str());
	}		
}
