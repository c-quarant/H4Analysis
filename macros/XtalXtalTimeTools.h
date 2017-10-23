#ifndef _XTAL_XTAL_TIME_ANALYSIS_TOOLS_H__
#define _XTAL_XTAL_TIME_ANALYSIS_TOOLS_H__

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

class XtalXtalTimeTools{
public:

	//ctor
	XtalXtalTimeTools() {};
	XtalXtalTimeTools(std::string RunNtuple); //OK
	//dtor
	~XtalXtalTimeTools(){};
	
	//Getters	
	int		GetRunEntries(){return Nentries;};
	int		GetRunNumber(){return Run;};
	float		GetRunGain(){return Gain;};
	float		GetRunEnergy(){return Energy;};
	std::string	GetRunStats(){return RunStats;};

	//utils
	void		GaussParInit(GaussPar* gP);
	float		HodoPlaneShift(std::string axis); //OK
	PlaneCoord	GetXtalCenterEdge(std::string APD, std::string edge); //da rimaneggiare per adattare le selezioni
	void		AmplitudeMapsEdge(std::string APD, std::string edge); //controllare
	GaussPar	NoiseAmplitudeDistributionFit(std::string APD);
	GaussPar	SummedAmplitudeFit();	
	void		TimeXTALvsXTALAeff();
	TH1F* 		AeffDistribution();
	std::vector<float>*	AeffMeanDistribution(int NSlices);
	std::vector<GaussPar>*	MyFitSlicesY(TH2F* Xtal_Xtal_Time);

protected:
	TTree* h4;
	int Nentries;
	int Run;
	float Gain;
	float Energy;
	std::string RunStats;

	std::vector< std::string > APDList;

	PlaneCoord HodoShift;
	
	std::map< std::string, PlaneCoord > Center;	
	std::map< std::string, GaussPar > NoiseAmplitude;
	GaussPar SummedAmplitude, Default;
	std::vector< GaussPar > SliceGaussFit;
	std::vector< float >	SliceAeffMean;
	std::vector< float >	SliceAeffMeanErr;
};

#endif 		
