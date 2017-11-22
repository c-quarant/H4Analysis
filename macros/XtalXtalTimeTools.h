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
#include "TObject.h"
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
/*
float fRes(float *x, float *par)
{
    if (x[0] > 0 && x[0] < 140) {
      TF1::RejectPoint();
      return 0;
   }
   return par[0] + par[1]*x[0];
}
*/
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
	std::vector<std::string>*	GetAPDList(){return &APDList;};

	//utils
	void		GaussParInit(GaussPar* gP);
	float		HodoPlaneShift(std::string axis);
	PlaneCoord	GetXtalCenterEdge(std::string APD, std::string edge);
	void		HodoPlaneMaps();
	GaussPar	NoiseAmplitudeDistributionFit(std::string APD);
	void		AmplitudeMapsEdge(std::string APD);

	// XTAL vs XTAL time analysis
	GaussPar	SummedAmplitudeFit();	
	TH1F* 		AeffDistribution();	
	std::vector<float>*	AeffMeanDistribution(int NSlices, std::vector<float>* AeffMean, std::vector<float>* AeffMeanErr);
	void		TimeXTALvsXTALAeff();

	// XTAL vs MCP time analysis	
	TH1F* 		AmplitudeDistribution(std::string APD);
	std::vector<float>*	AmplitudeMeanDistribution(std::string APD, int NSlices, std::vector<float>* AmpMean, std::vector<float>* AmpMeanErr);
	void		TimeXTALvsMCP(std::string APD);
	void		PulseShape(std::string APD, std::string TimeRef, std::string MagOMin);

	// Plotters
	std::vector<GaussPar>*	MyFitSlicesY(std::string APD, std::string TimeRef, TH2F* Xtal_Xtal_Time, std::vector<GaussPar>* SliceGaussFit);
	void		FitTimeResolution(std::string APD, std::string TimeRef);
	void 		TimeMeanvsAeff(std::string APD, std::string TimeRef);
	void		SaveAsPngPdfRoot(TObject* ToSave, std::string PathAndName, std::string DrawOpt);
	void		SaveAsPngPdfRoot(TObject* ToSave, TObject* ToSave1, std::string PathAndName, std::string DrawOpt);

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
	std::map< std::string, std::vector<GaussPar> >	SliceGaussFit;
	std::map< std::string, std::vector<float> >	SliceAmpOAeffMean;
	std::map< std::string, std::vector<float> >	SliceAmpOAeffMeanErr;
};

#endif 		
