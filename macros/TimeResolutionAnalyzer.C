#include "TimeAnalysisTools.cc"
#include "Riostream.h"
#include "TObjArray.h"

void TimeResolutionAnalyzer(std::string NtupleList, int bound)
{
	ifstream in;
	in.open(NtupleList);
	int i, j;
	int Run, Energy, Gain, detector, TimeRef;
	int C0APD1=0, C0APD2=1, C2=2, C3=3, C4=4, MCP1=5, MCP2=6, MCP_Mean=7;
	GaussPar TimePar;

	std::map< std::string, int > ChannelID;
	
	ChannelID["C0APD1"] = C0APD1;
	ChannelID["C0APD2"] = C0APD2;
	ChannelID["C2"] = C2;
	ChannelID["C3"] = C3;
	ChannelID["C4"] = C4;
	ChannelID["MCP1"] = MCP1;
	ChannelID["MCP2"] = MCP2;
	ChannelID["MCP_Mean"] = MCP_Mean;
	
	std::string path = (std::string)gSystem->UnixPathName(__FILE__);
	size_t found = path.find("macros/./TimeResolutionAnalyzer.C");
	path.replace(found, std::string("macros/./TimeResolutionAnalyzer.C").length(), "ntuples/");
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
		tD->Branch("MCP1", &MCP1, "MCP1/I");
		tD->Branch("MCP2", &MCP2, "MCP2/I");
		tD->Branch("MCP_Mean", &MCP_Mean, "MCP_Mean/I");

	
		tD->Branch("time_mean", &TimePar.Mean, "time_mean/F");
		tD->Branch("time_mean_error", &TimePar.MeanErr, "time_mean_error/F");
		tD->Branch("time_sigma", &TimePar.Sigma, "time_sigma/F");
		tD->Branch("time_sigma_error", &TimePar.SigmaErr, "time_sigma_error/F");
		//tD->Branch("Amplitude_Noise_Ratio", &TimePar.ANoiseR[0], "Amplitude_Noise_Ratio[n_timeRef]/F");
		//tD->Branch("Amplitude_Noise_Ratio_error", &TimePar.ANoiseRErr[0], "Amplitude_Noise_Ratio_error[n_timeRef]/F");
		
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
				if(BrName=="C0APD1" || BrName=="C0APD2" || BrName=="C2" || BrName=="C3" || BrName=="C4") RunAPDList.push_back(BrName);
				else if(BrName=="MCP1" || BrName=="MCP2") RunMCPList.push_back(BrName);
			}

			RunTimeAnalyzer = new TimeAnalysisTools(h4, RunAPDList, RunMCPList, bound);

			Run = RunTimeAnalyzer->runNum;
			Gain = stoi(RunTimeAnalyzer->Gain);
			Energy = stoi(RunTimeAnalyzer->Energy);
			
			for(i=0; i<RunAPDList.size(); i++)
			{
				detector = ChannelID[RunAPDList[i]];
				
				for(auto const& MCP : RunMCPList) 
				{
					TimeRef = ChannelID[MCP];
					TimePar = RunTimeAnalyzer->TimeAPDvsMCP(RunAPDList[i], MCP);
					tD->Fill();
				}
	
				TimeRef = MCP_Mean;
				TimePar = RunTimeAnalyzer->TimeAPDvsMCPMean(RunAPDList[i]);
				tD->Fill();

				for(j=i+1; j<RunAPDList.size(); j++)
				{
					TimeRef = ChannelID[RunAPDList[j]];
					TimePar = RunTimeAnalyzer->TimeAPDvsAPD(RunAPDList[i], RunAPDList[j]);
					tD->Fill();
				}
			}
			
			RunTimeAnalyzer->~TimeAnalysisTools();
		}
	}
	fOUT->cd();
	tD->Write();	
	fOUT->Close();

}
