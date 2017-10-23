#include "XtalXtalTimeTools.cc"
#include "Riostream.h"
#include "TObjArray.h"

void XtalXtalTimeResolution(std::string NtupleList, std::string OutputFile)
{
	ifstream in;
	in.open(NtupleList);
	int i, j;
	int Run, detector, TimeRes;
	float Energy, Gain;
	int C0APD1=0, C0APD2=1, C2=2, C3=3, C4=4, D3=5;
	float ANoiseR, ANoiseRErr;
	GaussPar TimePar, AeNoisePar, AeNoisePar2;

	std::map< std::string, int > ChannelID;
	
	ChannelID["C0APD1"] = C0APD1;
	ChannelID["C0APD2"] = C0APD2;
	ChannelID["C2"] = C2;
	ChannelID["C3"] = C3;
	ChannelID["C4"] = C4;
	ChannelID["D3"] = D3;
	
	std::string path = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/";
	std::string ntuple="";
	/*
	TTree* tD = new TTree("tD", "Parameters of time distribution (fit)");

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

	tD->Branch("time_mean", &TimePar.Mean, "time_mean/F");
	tD->Branch("time_mean_error", &TimePar.MeanErr, "time_mean_error/F");
	tD->Branch("time_sigma", &TimePar.Sigma, "time_sigma/F");
	tD->Branch("time_sigma_error", &TimePar.SigmaErr, "time_sigma_error/F");
	*/
	TFile* f;
	XtalXtalTimeTools* RunTimeAnalyzer;

	while(ntuple!="END")
	{
		in >> ntuple;
		if(ntuple=="END") break;
		cout << "\n>>>> Processing file   " << ntuple << endl;
			
		f = TFile::Open((path+ntuple).c_str());

		RunTimeAnalyzer = new XtalXtalTimeTools(path+ntuple);

		Run = RunTimeAnalyzer->GetRunNumber();
		Gain = RunTimeAnalyzer->GetRunGain();
		Energy = RunTimeAnalyzer->GetRunEnergy();
			
		RunTimeAnalyzer->TimeXTALvsXTALAeff();
			
		RunTimeAnalyzer->~XtalXtalTimeTools();
	}
	/*
	//Saving main result of time analysis on a root file called OutputFile
	TFile* fOUT = TFile::Open(OutputFile, "RECREATE");
	if(!fOUT) fOUT = new TFile(OutputFile, "RECREATE", "Stores main results of Time Analysis of Xtal wrt Xtal");
	fOUT->cd();
	tD->Write();	
	fOUT->Close();
	*/
}
