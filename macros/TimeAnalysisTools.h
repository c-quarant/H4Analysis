#ifndef _TIME_ANALYSIS_TOOLS_H__
#define _TIME_ANALYSIS_TOOLS_H__

#include "TLeaf.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include <iostream>
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TMath.h"
#include "Math/WrappedTF1.h"
#include "Math/RootFinderAlgorithms.h"
#include "TGraphErrors.h"
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>

struct PlaneCoord{
	float X;
	float Y;

	std::string X_str;
	std::string Y_str;
};

struct GaussPar{

	float Mean;
	float MeanErr;
	float Sigma;
	float SigmaErr;

};

class TimeAnalysisTools{
public:
	int runNum;
	std::string detector;
	std::string Run;
	std::string Gain;
	std::string Energy;

	//ctor
	TimeAnalysisTools() {};
	TimeAnalysisTools(TTree* ntupleTree, std::vector<std::string> RunAPDList, std::vector< std::string > RunMCPList, float bound_);
	//dtor
	~TimeAnalysisTools(){};

	//utils
	void		GaussParInit(GaussPar* gP);
	void		SetSelections();
	float		HodoPlaneShift(std::string axis);
	PlaneCoord 	GetHodoCenter(std::string APD);
	void 		AmplitudeMaps(std::string APD);
	void		TimeMaps(std::string APD, std::string MCP);
	float 		MeanTimeMCP(std::string MCP);
	GaussPar 	AmplitudeDistribution(std::string APD);
	GaussPar	TimeAPDvsMCP(std::string APD, std::string MCP);
	GaussPar	TimeAPDvsMCPMean(std::string APD);
	GaussPar	TimeAPDvsAPD(std::string APD1, std::string APD2);
	//GaussPar	ComputeAvsNoise(){};

protected:
	TTree* h4;

	std::string RunStats;
	std::vector< std::string > APDList;
	std::vector< std::string > MCPList;	

	float bound;
	float Xshift;
	float Yshift;

	std::string bound_str;
	std::string Xshift_str;
	std::string Yshift_str;

	std::map< std::string, PlaneCoord > Center;
	float MeanTimeMCP1;
	float MeanTimeMCP2;
	GaussPar Default;
	GaussPar Amplitude;
	GaussPar TimeResults;

	//selection strings
	std::map< std::string, std::string > DeviceSelections;
};
#endif 		
