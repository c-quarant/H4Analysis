#include "TimeAnalysisTools.h"

TimeAnalysisTools::TimeAnalysisTools(TTree* ntupleTree, std::vector<std::string> RunAPDList, std::vector< std::string > RunMCPList, float bound_)
{
	h4	 	= ntupleTree;
	detector 	= RunAPDList.at(0);
	bound	 	= bound_;
	bound_str	= std::to_string(bound_);
	APDList	 	= RunAPDList;
	MCPList		= RunMCPList;

	h4->GetEntry(0);

	runNum = (int)h4->GetLeaf("run")->GetValue(0);
	Run = std::to_string(runNum);
	Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	RunStats = Energy+"Gev_G"+Gain;

	Xshift= -999; Yshift= -999; 
	for(auto const& APD : APDList)
	{
		Center[APD].X= -999; Center[APD].Y= -999;
	}
	MeanTimeMCP1= -999; MeanTimeMCP2= -999;

	GaussParInit(&Amplitude);
	GaussParInit(&TimeResults);
	GaussParInit(&Default);

	SetSelections();
}

void TimeAnalysisTools::GaussParInit(GaussPar* gP)
{
	gP->Mean=-999; gP->MeanErr=-999; gP->Sigma=-999; gP->SigmaErr=-999;
}

void TimeAnalysisTools::SetSelections()
{
	HodoPlaneShift("X");
	HodoPlaneShift("Y");

	for(auto& MCP : MCPList)
		DeviceSelections[MCP+"AmpSel"]	=	"amp_max["+MCP+"]>200 && amp_max["+MCP+"]<2000";

	for(auto& APD : APDList)
	{
		Center[APD]		=	GetHodoCenter(APD);
		DeviceSelections[APD+"PosSel"] 	= 	"(fabs(X[0]-("+Center[APD].X_str+"))<"+bound_str+" || fabs(X[1]-("+Center[APD].X_str+")-("+Xshift_str+"))<"+bound_str+") && (fabs(Y[0]-("+Center[APD].Y_str+"))<"+bound_str+" || fabs(Y[1]-("+Center[APD].Y_str+")-("+Yshift_str+"))<"+bound_str+") && (fabs(X[0])<5 || fabs(X[1]-"+Xshift_str+")<5) && (fabs(Y[0])<5 || fabs(Y[1]-"+Yshift_str+")<5)";
	}
	
	for(auto& MCP : MCPList)
	{
		DeviceSelections[MCP+"TimeSel"]	=	"fabs(time["+MCP+"]-("+std::to_string(MeanTimeMCP(MCP))+"))<7";
		DeviceSelections[MCP]		=	DeviceSelections[MCP+"AmpSel"] + " && " + DeviceSelections[MCP+"TimeSel"];
	}
	
	for(auto& APD : APDList)
	{
		DeviceSelections[APD+"AmpSel"] = "fabs(fit_ampl["+APD+"]-("+std::to_string(AmplitudeDistribution(APD).Mean)+"))<"+std::to_string(AmplitudeDistribution(APD).Sigma);
		DeviceSelections[APD] = DeviceSelections[APD+"AmpSel"] + " && " + DeviceSelections[APD+"PosSel"];
	}	
}
		
//returns the shift of plane 1 from plane 0 of hodoscope along selected axis
float TimeAnalysisTools::HodoPlaneShift(std::string axis)
{
	if(axis == "X" && Xshift!=-999) return Xshift;
	if(axis == "Y" && Yshift!=-999) return Yshift;
		
	auto *DXvsX = new TProfile(("D"+axis+"vs"+axis+"").c_str(), "", 128, -16, 16, -10, 10);
	
	//Filling DeltaX histogram
	h4->Draw(("("+axis+"[1]-"+axis+"[0]):"+axis+"[0]>>D"+axis+"vs"+axis+"").c_str(), (axis+"[0]>-800 && "+axis+"[1]>-800").c_str());
	
	//Fitting and Drawing DeltaX
	TCanvas *cDeltaX = new TCanvas(("cDelta"+axis+"").c_str(), ("cDelta"+axis+"").c_str());
	TH2F* Hset = new TH2F(("Hset"+axis).c_str(),"", 128, -16, 16, 100, -5, 5);     
	Hset->GetXaxis()->SetTitle((axis+"[0]").c_str());    
 	Hset->GetYaxis()->SetTitle((axis+"[1]-"+axis+"[0]").c_str());
  
	Hset->Draw();	
	DXvsX->Fit("pol1", "Q", "", -4, 4);
	DXvsX->Draw("SAME");
		
	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".png";
    	
  	cDeltaX->SaveAs(fileOutpdf.c_str(), "Q");
  	cDeltaX->SaveAs(fileOutpng.c_str(), "Q");

	if(axis == "X")
	{
		Xshift		= DXvsX->GetFunction("pol1")->GetParameter(0);
		Xshift_str	= std::to_string(Xshift);
	}
	if(axis == "Y")
	{
		Yshift		= DXvsX->GetFunction("pol1")->GetParameter(0);
		Yshift_str	= std::to_string(Yshift);
	}

	cDeltaX->~TCanvas();
	
	return DXvsX->GetFunction("pol1")->GetParameter(0);
}

//Draw amplitude profiles in X and Y calculate maximum X
PlaneCoord TimeAnalysisTools::GetHodoCenter(std::string APD)
{
	if(Center[APD].X!=-999 && Center[APD].Y!=-999) 
		return Center[APD];
 
	float fitRange=5;
	
	gStyle->SetOptStat();
	gStyle->SetOptFit();

	auto *AmpXavg = new TProfile("AmpXavg", "", 128, -16, 16, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 128, -16, 16, 0, 10000);

	//Filling Amplitude profile histograms Draw method
	std::string varexp, selection;	
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);

	std::string AmpMCPSel="1";
	for(auto const& MCP : MCPList)
		AmpMCPSel += " && " + DeviceSelections[MCP+"AmpSel"];
	
	varexp = "fit_ampl[" + APD + "]:0.5*(X[0]+X[1]-(" + Xshift_str + "))>>AmpXavg";	
	selection = AmpMCPSel+" && 0.5*(Y[0]+Y[1]-("+Yshift_str+"))>-5 && 0.5*(Y[0]+Y[1]-("+Yshift_str+"))<5";
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "fit_ampl[" + APD + "]:0.5*(Y[0]+Y[1]-(" + Yshift_str + "))>>AmpYavg";
	selection = AmpMCPSel+" && 0.5*(X[0]+X[1]-(" + Xshift_str + "))>-5 && 0.5*(X[0]+X[1]-(" + Xshift_str + "))<5";
	h4->Draw(varexp.c_str(), selection.c_str());

	
	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
 	TH2F* H1 = new TH2F("H1","", 128, -16, 16, 50, 0, (AmpXavg->GetMaximum())*1.25);
 	H1->GetXaxis()->SetTitle("Xavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpXavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	AmpXavg->Draw("SAME");

	TF1 *fitResX = AmpXavg->GetFunction("pol2");
  
	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpXAVG_"+Run+"_"+APD+"_G"+Gain+".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpXAVG_"+Run+"_"+APD+"_G"+Gain+".png";
  
	cout << "\n\nX Center Position = " << fitResX->GetMaximumX() << "\t";  
	
  	c5->SaveAs(fileOutpdf.c_str(), "Q");
  	c5->SaveAs(fileOutpng.c_str(), "Q");
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw();  
	AmpYavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	AmpYavg->Draw("SAME");
	
	TF1* fitResY = AmpYavg->GetFunction("pol2");
	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpYAVG_"+Run+"_"+APD+"_G"+Gain+".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpYAVG_"+Run+"_"+APD+"_G"+Gain+".png";

	cout << "Y Center Position = " << fitResY->GetMaximumX() << "\n" << endl;      

  	c6->SaveAs(fileOutpdf.c_str());
  	c6->SaveAs(fileOutpng.c_str());
	
	/*
	if(stof(Run)<5700){
		if(fitResX->GetMaximumX()-bound<-5) Center[APD].X = bound-5;
		else if(fitResX->GetMaximumX()+bound>5) Center[APD].X = 5-bound;
		else Center[APD].X = fitResX->GetMaximumX();
	
		if(fitResY->GetMaximumX()-bound<-5) Center[APD].Y = bound-5;
		else if(fitResY->GetMaximumX()+bound>5) Center[APD].Y = 5-bound;
		else Center[APD].Y = fitResY->GetMaximumX();
	}
	*/	
	Center[APD].X = fitResX->GetMaximumX();
	Center[APD].Y = fitResY->GetMaximumX();
	Center[APD].X_str = std::to_string(fitResX->GetMaximumX());
	Center[APD].Y_str = std::to_string(fitResY->GetMaximumX());
	

	H1->~TH2F();
	c5->~TCanvas();
	c6->~TCanvas();

	return Center[APD];
}

void TimeAnalysisTools::AmplitudeMaps(std::string APD)
{	
	gStyle->SetOptStat(0); 
    	
  	//2DHist definition 
  	auto *AmpXY0 = new TProfile2D("AmpXY0","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	auto *AmpXY1 = new TProfile2D("AmpXY1","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	auto *AmpXYM = new TProfile2D("AmpXYM","", 32, -16, 16, 32, -16, 16, 0, 10000); 
  	
  	std::string Xshift_str = std::to_string(Xshift);
  	std::string Yshift_str = std::to_string(Yshift);
	std::string Selection;
	if(stof(Run)<5700) Selection = DeviceSelections["MCP1AmpSel"];
	else Selection = DeviceSelections["MCP1AmpSel"]+" && "+DeviceSelections["MCP2AmpSel"];

   	h4->Draw(("fit_ampl["+APD+"]:Y[0]:X[0]>>AmpXY0").c_str(), (Selection + " && X[0]>-800 && Y[0]>-800").c_str());
  	h4->Draw(("fit_ampl["+APD+"]:Y[1]-("+Yshift_str+"):X[1]-("+Xshift_str+")>>AmpXY1").c_str(), (Selection + " && X[1]>-800 && Y[1]>-800").c_str());
  	h4->Draw(("fit_ampl["+APD+"]:(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>AmpXYM").c_str(), (Selection + " && X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800").c_str());

  	//Drawing p0 histogram
  	TCanvas* c1 = new TCanvas("c1","c1");
	//FPCanvasStyle(c1, "", "", 0, "", 0, 1);
  	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");
  	AmpXY0->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());
  	
  	H1->Draw();	
  	AmpXY0->Draw("COLZ SAME");

	TLine *line_left = new TLine(Center[APD].X-bound, Center[APD].Y-bound, Center[APD].X-bound, Center[APD].Y+bound);
	line_left->SetLineColor(kRed);
	line_left->Draw();

	TLine *line_right = new TLine(Center[APD].X+bound, Center[APD].Y-bound, Center[APD].X+bound, Center[APD].Y+bound);
	line_right->SetLineColor(kRed);
	line_right->Draw();

	TLine *line_up = new TLine(Center[APD].X-bound, Center[APD].Y+bound, Center[APD].X+bound, Center[APD].Y+bound);
	line_up->SetLineColor(kRed);
	line_up->Draw();

	TLine *line_down = new TLine(Center[APD].X-bound, Center[APD].Y-bound, Center[APD].X+bound, Center[APD].Y-bound);
	line_down->SetLineColor(kRed);
	line_down->Draw();
  	
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + APD + "_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + APD + "_" + RunStats + ".png";
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  

  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
  	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");

	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + APD + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + APD + "_" + RunStats + ".png";

  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	

  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");
  	
	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	

  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + APD + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + APD + "_" + RunStats + ".png";  
  	
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());
	
	H1->~TH2F();
	c1->~TCanvas();
	c2->~TCanvas();
	c3->~TCanvas();
}

//returns the mean time measured by MCP selected
float TimeAnalysisTools::MeanTimeMCP(std::string MCP)
{
	if(MCP == "MCP1" && MeanTimeMCP1!=-999) return MeanTimeMCP1;
	if(MCP == "MCP2" && MeanTimeMCP2!=-999) return MeanTimeMCP2;
	
	float MCPTimeMean;

	TH1F* MCP_time_dist = new TH1F("MCP_time_dist", "", 200, 0, 50);
	
	HodoPlaneShift("X");
	HodoPlaneShift("Y");
	GetHodoCenter(detector);

	std::string Selection = DeviceSelections[detector+"PosSel"] + " && " + DeviceSelections[MCP+"AmpSel"]; 

	h4->Draw((std::string("time["+MCP+"]>>MCP_time_dist")).c_str(), Selection.c_str());
	
	MCP_time_dist->GetXaxis()->SetTitle(("time["+MCP+"] (ns)").c_str());
	MCP_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* ca = new TCanvas("ca", "ca");
    	ca->cd();
	MCP_time_dist->Draw();

	MCPTimeMean = MCP_time_dist->GetMean();

	ca -> SaveAs(std::string("/afs/cern.ch/user/c/cquarant/www/MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".png").c_str());
    	ca -> SaveAs(std::string("/afs/cern.ch/user/c/cquarant/www/MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".pdf").c_str());

	ca->~TCanvas();
	MCP_time_dist->~TH1F();

	return MCPTimeMean;
}

//Draw Amplitude histogram and fit it with gaus
GaussPar TimeAnalysisTools::AmplitudeDistribution(std::string APD)
{
	if(Amplitude.Mean!=-999 && Amplitude.MeanErr!=-999 && Amplitude.Sigma!=-999 && Amplitude.SigmaErr!=-999)
		return Amplitude;
	
	TH1F* HAmp = new TH1F("HAmp", "", 2000, +50, 5000);

	HodoPlaneShift("X");
	HodoPlaneShift("Y");
	GetHodoCenter(detector);
	MeanTimeMCP("MCP1");
	std::string Selection = DeviceSelections[APD+"PosSel"]; 
	for(auto const& MCP : MCPList)
		Selection += " && " + DeviceSelections[MCP];	
	
	h4->Draw(("fit_ampl["+APD+"]>>HAmp").c_str(), Selection.c_str());
	
	TCanvas* c0 = new TCanvas("c0", "c0");
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Fit("gaus", "Q");	
	HAmp->Draw();
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/Amp_"+APD+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/Amp_"+APD+"_"+RunStats+".pdf").c_str());
	
	Amplitude.Mean = HAmp->GetFunction("gaus")->GetParameter(1);
	Amplitude.MeanErr = HAmp->GetFunction("gaus")->GetParError(1);
	Amplitude.Sigma = HAmp->GetFunction("gaus")->GetParameter(2);
	Amplitude.SigmaErr = HAmp->GetFunction("gaus")->GetParError(2);

	c0->~TCanvas();
	HAmp->~TH1F();

	return Amplitude;	
}


GaussPar TimeAnalysisTools::TimeAPDvsMCP(std::string APD, std::string MCP)
{
	SetSelections();

	std::string tD_APD_MCP_Sel = DeviceSelections[APD] + " && " + DeviceSelections[MCP];

	//define APD_MCP time distribution
	TH1F* tD_APD_MCP = new TH1F("tD_APD_MCP", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD_MCP->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_MCP->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+APD+"]-time["+MCP+"]>>tD_APD_MCP").c_str(), tD_APD_MCP_Sel.c_str());
	
	float Xfirst = tD_APD_MCP->GetXaxis()->GetBinCenter(tD_APD_MCP->GetMaximumBin())-1;
	float Xlast = tD_APD_MCP->GetXaxis()->GetBinCenter(tD_APD_MCP->GetMaximumBin())+1;

	//plot and fit APD_MCP1 time distribution
	gStyle->SetOptStat();
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD_MCP->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_MCP->GetXaxis()->SetTitle(("time["+APD+"]-time["+MCP+"] (ns)").c_str());
	tD_APD_MCP->GetYaxis()->SetTitle("events");
	tD_APD_MCP->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_MCP->Draw();

	TimeResults.Mean = tD_APD_MCP->GetFunction("gaus")->GetParameter(1);
	TimeResults.MeanErr = tD_APD_MCP->GetFunction("gaus")->GetParError(1);
	TimeResults.Sigma = tD_APD_MCP->GetFunction("gaus")->GetParameter(2);
	TimeResults.SigmaErr = tD_APD_MCP->GetFunction("gaus")->GetParError(2);
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-"+MCP+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-"+MCP+"_"+RunStats+".pdf").c_str());

	c0->~TCanvas();
	tD_APD_MCP->~TH1F();

	return TimeResults;
}

void TimeAnalysisTools::TimeMaps(std::string APD, std::string MCP)
{	
	TimeAPDvsMCP(APD, MCP);
	gStyle->SetOptStat(0);
  		
  	//2DHist definition 
  	auto *TimeXY0 = new TProfile2D("TimeXY0","", 32, -16, 16, 32, -16, 16, 4, 15); 
  	auto *TimeXY1 = new TProfile2D("TimeXY1","", 32, -16, 16, 32, -16, 16, 4, 15); 
  	auto *TimeXYM = new TProfile2D("TimeXYM","", 32, -16, 16, 32, -16, 16, 0, 20); 

	std::string Selection = DeviceSelections[MCP]+" && "+DeviceSelections[APD+"AmpSel"];
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);
	  	
  	h4->Draw(("fit_time["+APD+"]-time["+MCP+"]:Y[0]:X[0]>>TimeXY0").c_str(), (Selection + " && X[0]>-16 && Y[0]>-16").c_str());
  	h4->Draw(("fit_time["+APD+"]-time["+MCP+"]:Y[1]-("+Yshift_str+"):X[1]-("+Xshift_str+")>>TimeXY1").c_str(), (Selection + " && X[1]>-16 && Y[1]>-16").c_str());
  	h4->Draw(("fit_time["+APD+"]-time["+MCP+"]:(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>TimeXYM").c_str(), (Selection + " && X[0]>-16 && Y[0]>-16 && X[1]>-16 && Y[1]>-16").c_str());

  	//Drawing p0 histogram
  	TCanvas* c1 = new TCanvas("c1","c1");
	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");
  	TimeXY0->GetZaxis()->SetTitle(("fit_time["+APD+"] (ADC counts)").c_str());
  	
  	H1->Draw();
	TimeXY0->GetZaxis()->SetRangeUser(TimeResults.Mean-3*TimeResults.Sigma, TimeResults.Mean+3*TimeResults.Sigma);	
  	TimeXY0->Draw("COLZ SAME");

	TLine *line_left = new TLine(Center[APD].X-bound, Center[APD].Y-bound, Center[APD].X-bound, Center[APD].Y+bound);
	line_left->SetLineColor(kRed);
	line_left->Draw();

	TLine *line_right = new TLine(Center[APD].X+bound, Center[APD].Y-bound, Center[APD].X+bound, Center[APD].Y+bound);
	line_right->SetLineColor(kRed);
	line_right->Draw();

	TLine *line_up = new TLine(Center[APD].X-bound, Center[APD].Y+bound, Center[APD].X+bound, Center[APD].Y+bound);
	line_up->SetLineColor(kRed);
	line_up->Draw();

	TLine *line_down = new TLine(Center[APD].X-bound, Center[APD].Y-bound, Center[APD].X+bound, Center[APD].Y-bound);
	line_down->SetLineColor(kRed);
	line_down->Draw();
  	
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane0/TimeXY_" + APD + "-" + MCP + "_" + RunStats + ".pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane0/TimeXY_" + APD + "-" + MCP + "_" + RunStats + ".png";
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  
  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
  	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	
  	H1->Draw();
	TimeXY1->GetZaxis()->SetRangeUser(TimeResults.Mean-3*TimeResults.Sigma, TimeResults.Mean+3*TimeResults.Sigma);	
  	TimeXY1->Draw("COLZ SAME");

	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane1/TimeXY_" + APD + "-" + MCP + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/Plane1/TimeXY_" + APD + "-" + MCP + "_" + RunStats + ".png";
  	
	
  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	
  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");

	c3->cd();
  	
  	H1->Draw();
	TimeXYM->GetZaxis()->SetRangeUser(TimeResults.Mean-3*TimeResults.Sigma, TimeResults.Mean+3*TimeResults.Sigma);	
  	TimeXYM->Draw("COLZ SAME");
  	
	line_left->Draw();
	line_right->Draw();
	line_up->Draw();
	line_down->Draw();	
	

  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/PlaneAVG/TimeXY_" + APD + "-" + MCP + "_" + RunStats + ".pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/TimeMaps/PlaneAVG/TimeXY_" + APD + "-" + MCP + "_" + RunStats + ".png";
  	
  	
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());
	
	c1->~TCanvas();
	c2->~TCanvas();
	c3->~TCanvas();
	TimeXY0->~TProfile2D();
	TimeXY1->~TProfile2D();
	TimeXYM->~TProfile2D();
	H1->~TH2F(); 
}


GaussPar TimeAnalysisTools::TimeAPDvsMCPMean(std::string APD)
{
	if(MCPList.size()<2)
	{
		cout << "Not enough MCP for Time APD vs MCP mean!!!" << endl; 
		return Default;
	}
	SetSelections();

	std::string tD_APD_MCPMean_Sel = DeviceSelections[APD] + " && " + DeviceSelections["MCP1"] + " && " + DeviceSelections["MCP2"];

	//define APD_MCP_Mean time distribution
	TH1F* tD_APD_MCP_Mean = new TH1F("tD_APD_MCP_Mean", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD_MCP_Mean->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_MCP_Mean->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+APD+"]-0.5*(time[MCP1]+time[MCP2])>>tD_APD_MCP_Mean").c_str(), tD_APD_MCPMean_Sel.c_str());
	
	float Xfirst = tD_APD_MCP_Mean->GetXaxis()->GetBinCenter(tD_APD_MCP_Mean->GetMaximumBin())-1;
	float Xlast = tD_APD_MCP_Mean->GetXaxis()->GetBinCenter(tD_APD_MCP_Mean->GetMaximumBin())+1;

	//plot and fit APD_MCP1 time distribution
	gStyle->SetOptStat();
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD_MCP_Mean->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_MCP_Mean->GetXaxis()->SetTitle(("time["+APD+"]-MeanTimeMCP (ns)").c_str());
	tD_APD_MCP_Mean->GetYaxis()->SetTitle("events");
	tD_APD_MCP_Mean->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_MCP_Mean->Draw();

	TimeResults.Mean = tD_APD_MCP_Mean->GetFunction("gaus")->GetParameter(1);
	TimeResults.MeanErr = tD_APD_MCP_Mean->GetFunction("gaus")->GetParError(1);
	TimeResults.Sigma = tD_APD_MCP_Mean->GetFunction("gaus")->GetParameter(2);
	TimeResults.SigmaErr = tD_APD_MCP_Mean->GetFunction("gaus")->GetParError(2);
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-MCP_Mean_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-MCP_Mean_"+RunStats+".pdf").c_str());

	c0->~TCanvas();
	tD_APD_MCP_Mean->~TH1F();

	return TimeResults;
}

GaussPar TimeAnalysisTools::TimeAPDvsAPD(std::string APD1, std::string APD2)
{
	if(APDList.size()<2)
	{
		cout << "Not enough APD for Time APD vs APD mean!!!" << endl; 
		return Default;
	}

	std::string tD_APD_APD_Sel = DeviceSelections[APD1] + " && " + DeviceSelections[APD2];

	//define APD_APD time distribution
	TH1F* tD_APD_APD = new TH1F("tD_APD_APD", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD_APD->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_APD->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+APD1+"]-time["+APD2+"]>>tD_APD_APD").c_str(), tD_APD_APD_Sel.c_str());
	
	float Xfirst = tD_APD_APD->GetXaxis()->GetBinCenter(tD_APD_APD->GetMaximumBin())-1;
	float Xlast = tD_APD_APD->GetXaxis()->GetBinCenter(tD_APD_APD->GetMaximumBin())+1;

	//plot and fit APD_APD time distribution
	gStyle->SetOptStat();
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD_APD->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_APD->GetXaxis()->SetTitle(("time["+detector+"]-MeanTimeMCP (ns)").c_str());
	tD_APD_APD->GetYaxis()->SetTitle("events");
	tD_APD_APD->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_APD->Draw();

	TimeResults.Mean = tD_APD_APD->GetFunction("gaus")->GetParameter(1);
	TimeResults.MeanErr = tD_APD_APD->GetFunction("gaus")->GetParError(1);
	TimeResults.Sigma = tD_APD_APD->GetFunction("gaus")->GetParameter(2);
	TimeResults.SigmaErr = tD_APD_APD->GetFunction("gaus")->GetParError(2);
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD1+"-"+APD2+"_"+RunStats+".pdf").c_str());

	c0->~TCanvas();
	tD_APD_APD->~TH1F();

	return TimeResults;
}		
