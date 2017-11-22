#include "TimeAnalysisTools.cc"
#include "Riostream.h"
#include "TObjArray.h"

void FourierAnalysis(std::string NtupleList, int bound)
{
	ifstream in;
	in.open(NtupleList);
	int i, j;
	int Run, Energy, Gain, detector, TimeRef;
	int C0APD1=0, C0APD2=1, C2=2, C3=3, C4=4, MCP1=5, MCP2=6, MCP_Mean=7, D3=8;
	float ANoiseR, ANoiseRErr;
	
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
	//size_t found = path.find("macros/./FourierAnalysis.C");
	//path.replace(found, std::string("macros/./FourierAnalysis.C").length(), "ntuples/");
	std::string ntuple="";
	
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

			for(auto& APD : RunAPDList)
				RunTimeAnalyzer->DrawFreqSpec(APD);

			for(auto& MCP : RunMCPList)
				RunTimeAnalyzer->DrawFreqSpecMCP(MCP);
	
			RunTimeAnalyzer->~TimeAnalysisTools();
		}
	}
