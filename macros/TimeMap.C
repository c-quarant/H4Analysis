#include "TimeAnalysisTools.cc"
#include "Riostream.h"
#include "TCollection.h"

void TimeMap(std::string NtupleList)
{
	PlaneCoord Center, CenterN, CenterE, CenterS;

	Center.X	=	237.5;	Center.Y	=	272.0;
	CenterN.X	=	237.5;	CenterN.Y	=	262.0;
	CenterE.X	=	222.5;	CenterE.Y	=	272.0;
	CenterS.X	=	237.5;	CenterS.Y	=	282.0; 

	Center.X_str	=	std::to_string(237.5);	Center.Y_str	=	std::to_string(272.0);
	CenterN.X_str	=	std::to_string(237.5);	CenterN.Y_str	=	std::to_string(262.0);
	CenterE.X_str	=	std::to_string(227.5);	CenterE.Y_str	=	std::to_string(272.0);
	CenterS.X_str	=	std::to_string(237.5);	CenterS.Y_str	=	std::to_string(282.0); 

	int bound = 3;
	ifstream in;
	in.open(NtupleList);

	std::string path = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/";
	//size_t found = path.find("macros/./TimeMap.C");
	//path.replace(found, std::string("macros/./TimeMap.C").length(), "ntuples/");
	std::string ntuple="";

	in >> ntuple;
	std::string Energy = ntuple;
	in >> ntuple;
	TFile* fC = TFile::Open((path+ntuple).c_str());
	in >> ntuple;
	TFile* fN = TFile::Open((path+ntuple).c_str());
	in >> ntuple;
	TFile* fE = TFile::Open((path+ntuple).c_str());
	in >> ntuple;
	TFile* fS = TFile::Open((path+ntuple).c_str());

	TTree* h4 = (TTree*)fC->Get("h4");
	TTree* h4N = (TTree*)fN->Get("h4");
	TTree* h4E = (TTree*)fE->Get("h4");
	TTree* h4S = (TTree*)fS->Get("h4");
	
	std::vector< std::string > RunAPDList, RunMCPList;
	RunMCPList.push_back("MCP1");

	RunAPDList.push_back("C3");
	TimeAnalysisTools* XTAL_C3 = new TimeAnalysisTools(h4, RunAPDList, RunMCPList, bound);  
	
	RunAPDList.push_back("C2");
	TimeAnalysisTools* XTAL_C3N = new TimeAnalysisTools(h4N, RunAPDList, RunMCPList, bound);  

	std::replace (RunAPDList.begin(), RunAPDList.end(), "C2", "D3");
	TimeAnalysisTools* XTAL_C3E = new TimeAnalysisTools(h4E, RunAPDList, RunMCPList, bound);  

	std::replace (RunAPDList.begin(), RunAPDList.end(), "D3", "C4");
	TimeAnalysisTools* XTAL_C3S = new TimeAnalysisTools(h4S, RunAPDList, RunMCPList, bound);  

	XTAL_C3->AmplitudeMaps("C3");
	XTAL_C3N->AmplitudeMapsEdge("C3", "N");
	XTAL_C3E->AmplitudeMapsEdge("C3", "E");
	XTAL_C3E->AmplitudeMapsEdge("D3", "E");
	XTAL_C3S->AmplitudeMapsEdge("C3", "S");

	GaussPar TimeC3 = XTAL_C3->TimeAPDvsMCP("C3", "MCP1");
	GaussPar TimeC3N = XTAL_C3N->TimeAPDvsMCPedge("C3", "MCP1");
	GaussPar TimeC3E = XTAL_C3E->TimeAPDvsMCPedge("C3", "MCP1");
	GaussPar TimeC3S = XTAL_C3S->TimeAPDvsMCPedge("C3", "MCP1");

	cout << endl << "Centro:" << endl << "Mean: " << TimeC3.Mean << " +/- " << TimeC3.MeanErr << endl << "Reso: " << TimeC3.Sigma << " +/- " << TimeC3.SigmaErr;
	cout << endl << "Nord:" << endl << "Mean: " << TimeC3N.Mean << " +/- " << TimeC3N.MeanErr << endl << "Reso: " << TimeC3N.Sigma << " +/- " << TimeC3N.SigmaErr;
	cout << endl << "Est:" << endl << "Mean: " << TimeC3E.Mean << " +/- " << TimeC3E.MeanErr << endl << "Reso: " << TimeC3E.Sigma << " +/- " << TimeC3E.SigmaErr;
	cout << endl << "Sud:" << endl << "Mean: " << TimeC3S.Mean << " +/- " << TimeC3S.MeanErr << endl << "Reso: " << TimeC3S.Sigma << " +/- " << TimeC3S.SigmaErr << endl;
	
	TProfile2D* TimeMap = new TProfile2D("TimeMap", "", 34, -18, 16, 34, -18, 16, "");
	
	//GaussPar temp = XTAL_C3->AmplitudeDistributionFit2Gaus	
	h4->Draw("fit_time[C3]-time[MCP1]:Y[0]:X[0]>>TimeMap", ("amp_max[MCP1]>200 && amp_max[MCP1]<2000 && fabs(time[MCP1]-("+std::to_string(XTAL_C3->MeanTimeMCP("MCP1"))+"))<7 && fabs(fit_ampl[C3]-("+std::to_string((XTAL_C3->AmplitudeDistributionFit2Gaus("C3")).Mean)+"))<"+std::to_string((XTAL_C3->AmplitudeDistributionFit2Gaus("C3")).Sigma)).c_str());
	TimeMap->GetZaxis()->SetRangeUser(TimeC3.Mean-6*TimeC3.Sigma, TimeC3.Mean+3*TimeC3.Sigma);

	//cout << endl << "fit_time[C3]:Y[0]-("+CenterN.Y_str+"):X[0]-("+CenterN.X_str+")>>+TimeMap" << endl; 
	
	h4N->Draw(("fit_time[C3]-time[MCP1]:Y[0]+("+std::to_string(CenterN.Y-Center.Y)+"):X[0]>>+TimeMap").c_str(), ("fabs(time[MCP1]-("+std::to_string(XTAL_C3N->MeanTimeMCP("MCP1"))+"))<7 && fit_ampl[C3]>400 && fabs(X[0])<5 && fabs(Y[0])<5").c_str());
	
	h4E->Draw(("fit_time[C3]-time[MCP1]:Y[0]:X[0]-("+std::to_string(CenterE.X-Center.X)+")>>+TimeMap").c_str(), ("fabs(time[MCP1]-("+std::to_string(XTAL_C3E->MeanTimeMCP("MCP1"))+"))<7 && fit_ampl[C3]>400 && X[0]>-10 && X[0]<0 && fabs(Y[0])<5").c_str());

	h4S->Draw(("fit_time[C3]-time[MCP1]:Y[0]+("+std::to_string(CenterS.Y-Center.Y)+"):X[0]>>+TimeMap").c_str(), ("fabs(time[MCP1]-("+std::to_string(XTAL_C3S->MeanTimeMCP("MCP1"))+"))<7 && fit_ampl[C3]>400 && fabs(X[0])<5 && fabs(Y[0])<5").c_str());
	
	gStyle->SetOptStat(0);
	TCanvas* cMap = new TCanvas("cMap", "cMap");
	TimeMap->Draw("COLZ");
	cMap->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeMaps/C3Edges/TimeXY_"+Energy+"Gev_G50.png").c_str());
	cMap->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeMaps/C3Edges/TimeXY_"+Energy+"Gev_G50.pdf").c_str());

	TProfile2D* AmpMap = new TProfile2D("AmpMap", "", 34, -18, 16, 34, -18, 16, "");
	
	//GaussPar temp = XTAL_C3->AmplitudeDistributionFit2Gaus	
	h4->Draw("fit_ampl[C3]:Y[0]:X[0]>>AmpMap", ("amp_max[MCP1]>200 && amp_max[MCP1]<2000 && fabs(time[MCP1]-("+std::to_string(XTAL_C3->MeanTimeMCP("MCP1"))+"))<7 && fabs(fit_ampl[C3]-("+std::to_string((XTAL_C3->AmplitudeDistributionFit2Gaus("C3")).Mean)+"))<"+std::to_string((XTAL_C3->AmplitudeDistributionFit2Gaus("C3")).Sigma)).c_str());
	TimeMap->GetZaxis()->SetRangeUser(TimeC3.Mean-6*TimeC3.Sigma, TimeC3.Mean+3*TimeC3.Sigma);

	//cout << endl << "fit_time[C3]:Y[0]-("+CenterN.Y_str+"):X[0]-("+CenterN.X_str+")>>+TimeMap" << endl; 
	
	h4N->Draw(("fit_ampl[C3]:Y[0]+("+std::to_string(CenterN.Y-Center.Y)+"):X[0]>>+AmpMap").c_str(), ("fabs(time[MCP1]-("+std::to_string(XTAL_C3N->MeanTimeMCP("MCP1"))+"))<7 && fit_ampl[C3]>400 && fabs(X[0])<5 && fabs(Y[0])<5").c_str());
	
	h4E->Draw(("fit_ampl[C3]:Y[0]:X[0]-("+std::to_string(CenterE.X-Center.X)+")>>+AmpMap").c_str(), ("fabs(time[MCP1]-("+std::to_string(XTAL_C3E->MeanTimeMCP("MCP1"))+"))<7 && fit_ampl[C3]>400 && X[0]>-10 && X[0]<0 && fabs(Y[0])<5").c_str());

	h4S->Draw(("fit_ampl[C3]:Y[0]+("+std::to_string(CenterS.Y-Center.Y)+"):X[0]>>+AmpMap").c_str(), ("fabs(time[MCP1]-("+std::to_string(XTAL_C3S->MeanTimeMCP("MCP1"))+"))<7 && fit_ampl[C3]>400 && fabs(X[0])<5 && fabs(Y[0])<5").c_str());
	
	TCanvas* cAmpMap = new TCanvas("cAmpMap", "cAmpMap");
	AmpMap->Draw("COLZ");
	cAmpMap->SaveAs(("/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/C3Edges/TimeXY_"+Energy+"Gev_G50.png").c_str());
	cAmpMap->SaveAs(("/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/C3Edges/TimeXY_"+Energy+"Gev_G50.pdf").c_str());
}
