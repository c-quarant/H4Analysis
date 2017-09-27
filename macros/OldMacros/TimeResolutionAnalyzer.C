#include "fitTimeVsMCPDist.C"
#include "fitTimeMCPs.C"
#include "fitTimeC0APDs.C"
#include "MyLib.h"
#include "Riostream.h"

void TimeResolutionAnalyzer(std::string NtupleList, int bound)
{
	ifstream in;
	in.open(NtupleList);
	int detector, C3=0, C0APD1=1, C0APD2=2, C0APDs=3, MCP1=0, MCP2=1, MCP_Mean=2, n_timeRef=3;
	TimeParameters TimePar;
	GaussPar TimeParMCPs, TimeParC0APDs;
	
	std::string path = (std::string)gSystem->UnixPathName(__FILE__);
	size_t found = path.find("macros/./TimeResolutionAnalyzer.C");
	path.replace(found, std::string("macros/./TimeResolutionAnalyzer.C").length(), "ntuples/");
	std::string ntuple="";
	
	TFile* fOUT;
	fOUT = TFile::Open("TimeAnalysisResults.root");
	TTree* tD;
	
	if(!fOUT){
		fOUT = new TFile("TimeAnalysisResults.root", "RECREATE", "Stores main results of Analysis");
		tD = new TTree("tD", "Parameters of time distribution (fit)");

		tD->Branch("Run", &TimePar.Run, "Run/I");
		tD->Branch("Energy", &TimePar.Energy, "Energy/I");
		tD->Branch("Gain", &TimePar.Gain, "Gain/I");
		tD->Branch("Detector", &detector, "Detector/I");
	
		tD->Branch("C3", &C3, "C3/I");
		tD->Branch("C0APD1", &C0APD1, "C0APD1/I");
		tD->Branch("C0APD2", &C0APD2, "C0APD2/I");
		tD->Branch("C0APDs", &C0APDs, "C0APDs/I");	
		
		tD->Branch("MCP1", &MCP1, "MCP1/I");
		tD->Branch("MCP2", &MCP2, "MCP2/I");
		tD->Branch("MCP_Mean", &MCP_Mean, "MCP_Mean/I");
		tD->Branch("n_timeRef", &n_timeRef, "n_timeRef/I");
	
		tD->Branch("time_mean", &TimePar.TimeMean[0], "time_mean[n_timeRef]/F");
		tD->Branch("time_mean_error", &TimePar.TimeMeanErr[0], "time_mean_error[n_timeRef]/F");
		tD->Branch("time_sigma", &TimePar.TimeSigma[0], "time_sigma[n_timeRef]/F");
		tD->Branch("time_sigma_error", &TimePar.TimeSigmaErr[0], "time_sigma_error[n_timeRef]/F");
		tD->Branch("Amplitude_Noise_Ratio", &TimePar.ANoiseR[0], "Amplitude_Noise_Ratio[n_timeRef]/F");
		tD->Branch("Amplitude_Noise_Ratio_error", &TimePar.ANoiseRErr[0], "Amplitude_Noise_Ratio_error[n_timeRef]/F");

		tD->Branch("MCPs_time_mean", &TimeParMCPs.Mean, "MCPs_time_mean/F");
		tD->Branch("MCPs_time_mean_error", &TimeParMCPs.MeanErr, "MCPs_time_mean/F");
		tD->Branch("MCPs_time_sigma", &TimeParMCPs.Sigma, "MCPs_time_sigma/F");
		tD->Branch("MCPs_time_sigma_error", &TimeParMCPs.SigmaErr, "MCPs_time_sigma/F");
		
		tD->Branch("C0APDs_time_mean", &TimeParC0APDs.Mean, "C0APDs_time_mean/F");
		tD->Branch("C0APDs_time_mean_error", &TimeParC0APDs.MeanErr, "C0APDs_time_mean/F");
		tD->Branch("C0APDs_time_sigma", &TimeParC0APDs.Sigma, "C0APDs_time_sigma/F");
		tD->Branch("C0APDs_time_sigma_error", &TimeParC0APDs.SigmaErr, "C0APDs_time_sigma/F");
	
		while(ntuple!="END"){
			in >> ntuple;
			if(ntuple=="END") break;
			cout << "\n>>>> Processing file   " << ntuple << endl;
			
			if(ntuple.find("C3")!=std::string::npos)
			{
				TimePar = fitTimeVsMCPDist((path+ntuple).c_str(), "C3", bound);
				TimeParMCPs = fitTimeMCPs((path+ntuple).c_str());
				detector = C3;

				TimeParC0APDs.Mean = -999;
				TimeParC0APDs.MeanErr = -999;
				TimeParC0APDs.Sigma = -999;
				TimeParC0APDs.SigmaErr = -999;
				
				tD->Fill();
			}			
			else{
				TimePar = fitTimeVsMCPDist((path+ntuple).c_str(), "C0APD1", bound);
				TimeParMCPs = fitTimeMCPs((path+ntuple).c_str());
				detector = C0APD1;

				TimeParC0APDs.Mean = -999;
				TimeParC0APDs.MeanErr = -999;
				TimeParC0APDs.Sigma = -999;
				TimeParC0APDs.SigmaErr = -999;

				tD->Fill();
				
				TimePar = fitTimeVsMCPDist((path+ntuple).c_str(), "C0APD2", bound);
				detector = C0APD2;

				TimeParC0APDs.Mean = -999;
				TimeParC0APDs.MeanErr = -999;
				TimeParC0APDs.Sigma = -999;
				TimeParC0APDs.SigmaErr = -999;
				
				tD->Fill();

				TimeParC0APDs = fitTimeC0APDs((path+ntuple).c_str(), &TimePar.ANoiseR[0], &TimePar.ANoiseRErr[0], bound);
				
				
				detector = C0APDs;
				
				tD->Fill();
			}
			
		}
		fOUT->cd();
		tD->Write();
		in.close();
	}
	else
	{
		cout << "Ci ENTRO" << endl;
		tD = (TTree*)fOUT->Get("tD");

		tD->SetBranchAddress("Run", &TimePar.Run);
		tD->SetBranchAddress("Energy", &TimePar.Energy);
		tD->SetBranchAddress("Gain", TimePar->Gain);
		tD->SetBranchAddress("Detector", &detector);
	
		tD->SetBranchAddress("C3", &C3);
		tD->SetBranchAddress("C0APD1", &C0APD1);
		tD->SetBranchAddress("C0APD2", &C0APD2);	
		
		tD->SetBranchAddress("MCP1", &MCP1);
		tD->SetBranchAddress("MCP2", &MCP2);
		tD->SetBranchAddress("MCP_Mean", &MCP_Mean);
		tD->SetBranchAddress("n_timeRef", &n_timeRef);
	
		tD->SetBranchAddress("time_mean", TimePar.TimeMean);
		tD->SetBranchAddress("time_mean_error", TimePar.TimeMeanErr);
		tD->SetBranchAddress("time_sigma", TimePar.TimeSigma);
		tD->SetBranchAddress("time_sigma_error", TimePar.TimeSigmaErr);
		tD->SetBranchAddress("Amplitude_Noise_Ratio", TimePar.ANoiseR);
		tD->SetBranchAddress("Amplitude_Noise_Ratio_error", TimePar.ANoiseRErr);

		tD->SetBranchAddress("MCPs_time_mean", &TimeParMCPs.Mean);
		tD->SetBranchAddress("MCPs_time_mean_error", &TimeParMCPs.MeanErr);
		tD->SetBranchAddress("MCPs_time_sigma", &TimeParMCPs.Sigma);
		tD->SetBranchAddress("MCPs_time_sigma_error", &TimeParMCPs.SigmaErr);

		tD->SetBranchAddress("C0APDs_time_mean", &TimeParC0APDs.Mean);
		tD->SetBranchAddress("C0APDs_time_mean_error", &TimeParC0APDs.MeanErr);
		tD->SetBranchAddress("C0APDs_time_sigma", &TimeParC0APDs.Sigma);
		tD->SetBranchAddress("C0APDs_time_sigma_error", &TimeParC0APDs.SigmaErr);	
	}

	FitRes(tD, C3, MCP1);
	FitRes(tD, C3, MCP2);
	FitRes(tD, C3, MCP_Mean);

	FitRes(tD, C0APD1, MCP1);
	FitRes(tD, C0APD1, MCP2);
	FitRes(tD, C0APD1, MCP_Mean);

	FitRes(tD, C0APD2, MCP1);
	FitRes(tD, C0APD2, MCP2);
	FitRes(tD, C0APD2, MCP_Mean);

	FitRes(tD, C0APDs, 0);

	
	cout << "!!!!!!!!!!! CIAO !!!!!!!!!!!!" << endl;
	std::vector<float> X0, X0Err, Y0, Y0Err, X1, X1Err, Y1, Y1Err, X2, X2Err, Y2, Y2Err;
	int SelMCP = 0, Nentries=tD->GetEntries(), i;

	for(i=0; i<Nentries; i++)
	{
	tD->GetEntry(i);
		cout << detector << "     " << TimePar.Gain << endl;
		if(detector == C0APDs && TimePar.Gain == 50)
		{
			X0.push_back(TimePar.ANoiseR[SelMCP]);
			X0Err.push_back(TimePar.ANoiseRErr[SelMCP]);
			Y0.push_back(TimeParC0APDs.Sigma*1000);
			Y0Err.push_back(TimeParC0APDs.SigmaErr*1000);
		}
		if(detector == C0APD1 && TimePar.Gain == 50)
		{
			X1.push_back(TimePar.ANoiseR[SelMCP]);
			X1Err.push_back(TimePar.ANoiseRErr[SelMCP]);
			Y1.push_back(TimePar.TimeSigma[SelMCP]*1000);
			Y1Err.push_back(TimePar.TimeSigmaErr[SelMCP]*1000);
		}
		if(detector == C0APD2 && TimePar.Gain == 50)
		{
			X2.push_back(TimePar.ANoiseR[SelMCP]);
			X2Err.push_back(TimePar.ANoiseRErr[SelMCP]);
			Y2.push_back(TimePar.TimeSigma[SelMCP]*1000);
			Y2Err.push_back(TimePar.TimeSigmaErr[SelMCP]*1000);
		}
	}	
	cout << "!!!!!!!!!!! CIAO !!!!!!!!!!!!" << endl;
	cout << endl << X1[0] << endl << endl;
	TCanvas* Can = new TCanvas("Can", "Can");
	Can->cd();

	TGraphErrors* g0 = new TGraphErrors(X0.size(), &X0[0], &Y0[0], &X0Err[0], &Y0Err[0]);
	TGraphErrors* g1 = new TGraphErrors(X1.size(), &X1[0], &Y1[0], &X1Err[0], &Y1Err[0]);
	TGraphErrors* g2 = new TGraphErrors(X2.size(), &X2[0], &Y2[0], &X2Err[0], &Y2Err[0]);
	
	g0->SetTitle("");
	g1->SetTitle("");
	g2->SetTitle("");

	g0->GetXaxis()->SetTitle("A/#sigma(Noise)");
	g0->GetYaxis()->SetTitle("#sigma(APD-MCP) (ps)");
	
	g0->SetMarkerStyle(kFullCircle);
	g0->SetMarkerSize(1);
	
	g1->SetMarkerStyle(2);
	g1->SetMarkerSize(2);
	
	g2->SetMarkerStyle(3);
	g2->SetMarkerSize(3);
	
	g0->Draw("APL");
	g1->Draw("PL");
	g2->Draw("PL");

	cout << endl << X1[0] << endl << endl;
	
	Can->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_C0APDs&normalAPD.png");
	Can->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionFit/TimeRes_vs_ANoise_C0APDs&normalAPD.pdf");

	fOUT->Close();

}
