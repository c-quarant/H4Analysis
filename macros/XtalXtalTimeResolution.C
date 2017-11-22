#include "XtalXtalTimeTools.cc"
#include "Riostream.h"
#include "TObjArray.h"

void XtalXtalTimeResolution(std::string NtupleList)
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
		RunTimeAnalyzer->TimeXTALvsMCP((RunTimeAnalyzer->GetAPDList())->at(0));
		RunTimeAnalyzer->TimeXTALvsMCP((RunTimeAnalyzer->GetAPDList())->at(1));
		//RunTimeAnalyzer->PulseShape((RunTimeAnalyzer->GetAPDList())->at(0), "MCP1", ">");
		//RunTimeAnalyzer->PulseShape((RunTimeAnalyzer->GetAPDList())->at(0), "MCP1", "<");

		RunTimeAnalyzer->~XtalXtalTimeTools();
	}
}
