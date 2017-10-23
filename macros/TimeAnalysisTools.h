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
	std::string RunStats;

	//ctor
	TimeAnalysisTools() {};
	TimeAnalysisTools(TTree* ntupleTree, std::vector<std::string> RunAPDList, std::vector< std::string > RunMCPList, float bound_);
	TimeAnalysisTools(TTree* ntupleTree, std::vector<std::string> RunAPDList, std::vector< std::string > RunMCPList);
	//dtor
	~TimeAnalysisTools(){};

	//utils
	void		GaussParInit(GaussPar* gP);
	void		SetSelections();
	float		HodoPlaneShift(std::string axis);
	PlaneCoord 	GetHodoCenter(std::string APD);
	PlaneCoord 	GetHodoCenterEdge(std::string APD, std::string edge);
	void 		AmplitudeMaps(std::string APD);
	void 		AmplitudeMapsEdge(std::string APD, std::string edge);
	void		TimeMaps(std::string APD, std::string MCP);
	void		DrawFreqSpec(std::string APD, std::string MCP);
	void		DrawFreqSpecPedestal(std::string APD);
	float 		MeanTimeMCP(std::string MCP);
	GaussPar 	AmplitudeDistributionFit2Gaus(std::string APD);
	GaussPar 	NoiseAmplitudeDistributionFit(TTree* h4Noise, std::string APD);	
	GaussPar	TimeAPDvsMCP(std::string APD, std::string MCP);
	GaussPar	TimeAPDvsMCPedge(std::string APD, std::string MCP);
	GaussPar	TimeAPDvsMCPMean(std::string APD);
	GaussPar	TimeAPDvsAPD(std::string APD1, std::string APD2);
	GaussPar	TimeXTALvsXTAL(std::string APD1, std::string APD2);
	void		TimeXTALvsXTALAeff(std::string APD1, std::string APD2);
	GaussPar	ComputeAvsNoise(std::string APD);
	void		ComputeAvsNoiseEdge(std::string APD1, std::string APD2);

protected:
	TTree* h4;

	std::vector< std::string > APDList;
	std::vector< std::string > MCPList;	

	float bound;
	float Xshift;
	float Yshift;

	std::string bound_str;
	std::string Xshift_str;
	std::string Yshift_str;

	std::map< std::string, PlaneCoord > Center;	
	std::map< std::string, GaussPar > Amplitude, NoiseAmpl, TimeResults, NoiseAmplitude;
	float MeanTimeMCP1;
	float MeanTimeMCP2;
	GaussPar Default;
	
	//selection strings
	std::map< std::string, std::string > DeviceSelections;
};
#endif 		
