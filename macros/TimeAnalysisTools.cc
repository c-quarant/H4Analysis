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
	RunStats = Energy+"Gev_G"+Gain+"_"+Run;

	Xshift= -999; Yshift= -999; 
	for(auto const& APD : APDList)
	{
		Center[APD].X= -999; Center[APD].Y= -999;
		GaussParInit(&Amplitude[APD]);
		GaussParInit(&TimeResults[APD]);
	}
	MeanTimeMCP1= -999; MeanTimeMCP2= -999;
	GaussParInit(&Default);

	SetSelections();
}

TimeAnalysisTools::TimeAnalysisTools(TTree* ntupleTree, std::vector<std::string> RunAPDList, std::vector< std::string > RunMCPList)
{
	h4	 	= ntupleTree;
	detector 	= RunAPDList.at(0);
	APDList	 	= RunAPDList;
	MCPList		= RunMCPList;

	h4->GetEntry(0);

	runNum = (int)h4->GetLeaf("run")->GetValue(0);
	Run = std::to_string(runNum);
	Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	RunStats = "G"+Gain+"_"+Run;

	Xshift= -999; Yshift= -999; 
	for(auto const& APD : APDList)
	{
		Center[APD].X= -999; Center[APD].Y= -999;
		GaussParInit(&Amplitude[APD]);
		GaussParInit(&TimeResults[APD]);
	}
	MeanTimeMCP1= -999; MeanTimeMCP2= -999;
	GaussParInit(&Default);
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
		AmplitudeDistributionFit2Gaus(APD);
		DeviceSelections[APD+"AmpSel"] = "fabs(fit_ampl["+APD+"]-("+std::to_string(Amplitude[APD].Mean)+"))<"+std::to_string(Amplitude[APD].Sigma);
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

//Draw amplitude profiles in X and Y calculate X Center
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

PlaneCoord TimeAnalysisTools::GetHodoCenterEdge(std::string APD, std::string edge)
{
	PlaneCoord Output;
	float fitRange=7;
	
	gStyle->SetOptStat();
	gStyle->SetOptFit();

	auto *AmpXavg = new TProfile("AmpXavg", "", 34, -18, 16, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 34, -18, 16, 0, 10000);

	//Filling Amplitude profile histograms Draw method
	std::string varexp;	
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);

	varexp = "fit_ampl[" + APD + "]:0.5*(X[0]+X[1]-(" + Xshift_str + "))>>AmpXavg";	
	std::string Selection ="0.5*(Y[0]+Y[1]-("+Yshift_str+"))>-15 && 0.5*(Y[0]+Y[1]-("+Yshift_str+"))<15";
	h4->Draw(varexp.c_str(), Selection.c_str());

	varexp = "fit_ampl[" + APD + "]:0.5*(Y[0]+Y[1]-(" + Yshift_str + "))>>AmpYavg";
	Selection = "0.5*(X[0]+X[1]-(" + Xshift_str + "))>-15 && 0.5*(X[0]+X[1]-(" + Xshift_str + "))<15";
	h4->Draw(varexp.c_str(), Selection.c_str());

	
	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
 	TH2F* H1 = new TH2F("H1","", 128, -16, 16, 50, 0, (AmpXavg->GetMaximum())*1.25);
 	H1->GetXaxis()->SetTitle("Xavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	if(edge == "E") AmpXavg->Fit("pol2", "Q", "", -16, -5);
	else AmpXavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	AmpXavg->Draw("SAME");

	TF1 *fitResX = AmpXavg->GetFunction("pol2");
  
	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpXAVG_"+Run+"_"+APD+"_G"+Gain+"_NoMCPSel.pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpXAVG_"+Run+"_"+APD+"_G"+Gain+"_NoMCPSel.png";
  
	cout << "\n\nX Center Position = " << fitResX->GetMaximumX() << "\t";  
	
  	c5->SaveAs(fileOutpdf.c_str(), "Q");
  	c5->SaveAs(fileOutpng.c_str(), "Q");
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw();  
	if(edge=="N") AmpYavg->Fit("pol2", "Q", "", 0, 15);
	else if(edge=="S") AmpYavg->Fit("pol2", "Q", "", -15, 0);
	else AmpYavg->Fit("pol2", "Q", "", -fitRange, fitRange);		
	AmpYavg->Draw("SAME");
	
	TF1* fitResY = AmpYavg->GetFunction("pol2");
	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpYAVG_"+Run+"_"+APD+"_G"+Gain+"_NoMCPSel.pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/pAVG/"+Energy+"Gev/"+"AmpYAVG_"+Run+"_"+APD+"_G"+Gain+"_NoMCPSel.png";

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
	Output.X = fitResX->GetMaximumX();
	Output.Y = fitResY->GetMaximumX();
	Output.X_str = std::to_string(fitResX->GetMaximumX());
	Output.Y_str = std::to_string(fitResY->GetMaximumX());
	

	H1->~TH2F();
	c5->~TCanvas();
	c6->~TCanvas();

	return Output;
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

void TimeAnalysisTools::AmplitudeMapsEdge(std::string APD, std::string edge)
{	
	gStyle->SetOptStat(0); 
    	
  	//2DHist definition 
  	auto *AmpXY0 = new TProfile2D("AmpXY0","", 32, -17, 15, 32, -16, 16, 0, 10000); 
  	auto *AmpXY1 = new TProfile2D("AmpXY1","", 32, -17, 15, 32, -16, 16, 0, 10000); 
  	auto *AmpXYM = new TProfile2D("AmpXYM","", 32, -17, 15, 32, -16, 16, 0, 10000); 
  	
  	std::string Xshift_str = std::to_string(Xshift);
  	std::string Yshift_str = std::to_string(Yshift);
	
   	h4->Draw(("fit_ampl["+APD+"]:Y[0]:X[0]>>AmpXY0").c_str(), "X[0]>-800 && Y[0]>-800");
  	h4->Draw(("fit_ampl["+APD+"]:Y[1]-("+Yshift_str+"):X[1]-("+Xshift_str+")>>AmpXY1").c_str(), "X[1]>-800 && Y[1]>-800");
  	h4->Draw(("fit_ampl["+APD+"]:(0.5*(Y[0]+Y[1]-("+Yshift_str+"))):(0.5*(X[0]+X[1]-("+Xshift_str+")))>>AmpXYM").c_str(), "X[0]>-800 && Y[0]>-800 && X[1]>-800 && Y[1]>-800");

	PlaneCoord TempCenter = GetHodoCenterEdge(APD, edge);  	
	//Drawing p0 histogram
  	TCanvas* c1 = new TCanvas("c1","c1");
	gStyle->SetOptStat(0);
	//FPCanvasStyle(c1, "", "", 0, "", 0, 1);
  	TH2F* H1 = new TH2F("H1","", 32, -17, 15, 32, -16, 16);     
  	H1->GetXaxis()->SetTitle("X[0]");    
  	H1->GetYaxis()->SetTitle("Y[0]");
  	AmpXY0->GetZaxis()->SetTitle(("amp_max["+APD+"] (ADC counts)").c_str());
  	
  	H1->Draw();	
  	AmpXY0->Draw("COLZ SAME");

	TLine *line_left = new TLine(TempCenter.X-bound, TempCenter.Y-bound, TempCenter.X-bound, TempCenter.Y+bound);
	line_left->SetLineColor(kRed);
	//line_left->Draw();

	TLine *line_right = new TLine(TempCenter.X+bound, TempCenter.Y-bound, TempCenter.X+bound, TempCenter.Y+bound);
	line_right->SetLineColor(kRed);
	//line_right->Draw();

	TLine *line_up = new TLine(TempCenter.X-bound, TempCenter.Y+bound, TempCenter.X+bound, TempCenter.Y+bound);
	line_up->SetLineColor(kRed);
	//line_up->Draw();

	TLine *line_down = new TLine(TempCenter.X-bound, TempCenter.Y-bound, TempCenter.X+bound, TempCenter.Y-bound);
	line_down->SetLineColor(kRed);
	//line_down->Draw();
  	
  	std::string fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + APD + "_" + RunStats + "_NoMCPSel.pdf";
  	std::string fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane0/AmpXY_" + APD + "_" + RunStats + "_NoMCPSel.png";
  	
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
  

  	//Drawing p1 histogram
  	TCanvas* c2 = new TCanvas("c2","c2");
	gStyle->SetOptStat(0);
  	H1->GetXaxis()->SetTitle("X[1]");    
  	H1->GetYaxis()->SetTitle("Y[1]");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");

	//line_left->Draw();
	//line_right->Draw();
	//line_up->Draw();
	//line_down->Draw();	
  	
  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + APD + "_" + RunStats + "_NoMCPSel.pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/Plane1/AmpXY_" + APD + "_" + RunStats + "_NoMCPSel.png";

  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());
  	

  	//Drawving pAVG histogram
  	TCanvas* c3 = new TCanvas("c3","c3");
	gStyle->SetOptStat(0);
  	H1->GetXaxis()->SetTitle("X_AVG");    
  	H1->GetYaxis()->SetTitle("Y_AVG");
  	
  	H1->Draw();	
  	AmpXY1->Draw("COLZ SAME");
  	
	//line_left->Draw();
	//line_right->Draw();
	//line_up->Draw();
	//line_down->Draw();	

  	fileOutpdf = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + APD + "_" + RunStats + "_NoMCPSel.pdf";
  	fileOutpng = "/afs/cern.ch/user/c/cquarant/www/AmplitudeMaps/PlaneAVG/AmpXY_" + APD + "_" + RunStats + "_NoMCPSel.png";  
  	
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

	if(MCP == "MCP1" && MeanTimeMCP1==-999) MeanTimeMCP1 = MCPTimeMean;
	if(MCP == "MCP2" && MeanTimeMCP2==-999) MeanTimeMCP2 = MCPTimeMean;

	return MCPTimeMean;
}

//Draw Amplitude histogram and fit it with sum of 2 gaus
GaussPar TimeAnalysisTools::AmplitudeDistributionFit2Gaus(std::string APD)
{
	HodoPlaneShift("X");
	HodoPlaneShift("Y");
	GetHodoCenter(APD);
	std::string Selection = DeviceSelections[APD+"PosSel"]; 
	for(auto const& MCP : MCPList)
	{
		MeanTimeMCP(MCP);
		Selection += " && " + DeviceSelections[MCP+"AmpSel"];
	}
	//cout << endl << endl << Selection << endl << endl;
	TH1F* HAmp = new TH1F("HAmp", "", 2000, 0, 5000);
	h4->Draw(("fit_ampl["+APD+"]>>HAmp").c_str(), Selection.c_str());
	
	float XMax = HAmp->GetBinCenter(HAmp->GetMaximumBin());
	float Xfirst = XMax*0.6;
	float Xlast = XMax*1.2;

	HAmp->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	TF1* fitFunc = new TF1("fitFunc", "[0]*TMath::Exp(-(x-[1])*(x-[1])/(2*[2]*[2])) + [3]*TMath::Exp(-(x-[4])*(x-[4])/(2*[5]*[5]))", Xfirst, Xlast);
	
	fitFunc->SetParLimits(0, HAmp->GetMaximum()*0.8, HAmp->GetMaximum()*1.39);
	fitFunc->SetParLimits(1, XMax*0.95, XMax*1.05);
	fitFunc->SetParLimits(2, 0, 150);
	fitFunc->SetParLimits(3, 0, HAmp->GetMaximum()*0.7);
	fitFunc->SetParLimits(4, XMax*0.85, XMax*0.95);
	fitFunc->SetParLimits(5, 5, 600);

	fitFunc->SetParameter(0, HAmp->GetMaximum());
	fitFunc->SetParameter(1, XMax);
	fitFunc->SetParameter(2, 40);
	fitFunc->SetParameter(3, HAmp->GetMaximum()/10);
	fitFunc->SetParameter(4, XMax*0.9);
	fitFunc->SetParameter(5, 80);	

	fitFunc->SetParNames("A0", "Mean0", "Sigma0", "A1", "Mean1", "Sigma1");

	HAmp->Fit("fitFunc", "Q");	
	

	TCanvas* c1 = new TCanvas("c1", "c1");

	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->GetXaxis()->SetTitle("Amplitude (ADC counts)");
	HAmp->GetYaxis()->SetTitle("events");
	HAmp->Draw();
	
	//draw components
	TF1* gaus0 = new TF1("gaus0", "gaus", Xfirst+0.6, Xlast-0.6);
	TF1* gaus1 = new TF1("gaus1", "gaus", Xfirst+0.6, Xlast-0.6);

	gaus0->SetParameter(0, fitFunc->GetParameter("A0"));
	gaus0->SetParameter(1, fitFunc->GetParameter("Mean0"));
	gaus0->SetParameter(2, fitFunc->GetParameter("Sigma0"));

	gaus1->SetParameter(0, fitFunc->GetParameter("A1"));
	gaus1->SetParameter(1, fitFunc->GetParameter("Mean1"));
	gaus1->SetParameter(2, fitFunc->GetParameter("Sigma1"));

	gaus0->SetLineColor(kBlue);
	gaus0->Draw("SAME");

	gaus1->SetLineColor(kGreen);
	gaus1->Draw("SAME");

	Amplitude[APD].Mean = fitFunc->GetParameter(1);
	Amplitude[APD].MeanErr = fitFunc->GetParError(1);
	Amplitude[APD].Sigma = fitFunc->GetParameter(2);
	Amplitude[APD].SigmaErr = fitFunc->GetParError(2);

	c1->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/FitAmp_"+APD+"_"+RunStats+".png").c_str());
	c1->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/FitAmp_"+APD+"_"+RunStats+".pdf").c_str());

	c1->~TCanvas();
	
	return Amplitude[APD];
}

//Draw Amplitude histogram and fit it with single gaus
GaussPar TimeAnalysisTools::NoiseAmplitudeDistributionFit(TTree* h4Noise, std::string APD) 
{
	gStyle->SetOptStat();	
	//cout << "Drawing noise amplitude distribution..." << endl;
	TH1F* HAmp = new TH1F("HAmp", "", 45, -30, 40);
	if(APD=="C2" || APD=="D3") 	h4Noise->Draw("fit_ampl[C3]>>HAmp");
	else h4Noise->Draw(("fit_ampl["+APD+"]>>HAmp").c_str());
	HAmp->GetXaxis()->SetRangeUser(HAmp->GetBinCenter(HAmp->GetMaximumBin())-30, HAmp->GetBinCenter(HAmp->GetMaximumBin())+30);
	
	TCanvas* c0 = new TCanvas("c0", "c0");
	HAmp->Fit("gaus", "Q");	
	HAmp->Draw();
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/NoiseAmp_"+APD+"_"+Gain+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/NoiseAmp_"+APD+"_"+Gain+".pdf").c_str());
	
	NoiseAmplitude[APD].Mean = HAmp->GetFunction("gaus")->GetParameter(1);
	NoiseAmplitude[APD].MeanErr = HAmp->GetFunction("gaus")->GetParError(1);
	NoiseAmplitude[APD].Sigma = HAmp->GetFunction("gaus")->GetParameter(2);
	NoiseAmplitude[APD].SigmaErr = HAmp->GetFunction("gaus")->GetParError(2);

	c0->~TCanvas();
	HAmp->~TH1F();

	return NoiseAmplitude[APD];
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
	tD_APD_MCP->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_MCP->GetXaxis()->SetTitle(("t_{"+APD+"}-t_{"+MCP+"} (ns)").c_str());
	tD_APD_MCP->GetXaxis()->SetTitleSize(0.055);
	tD_APD_MCP->GetXaxis()->SetTitleOffset(0.75);
	tD_APD_MCP->GetYaxis()->SetTitle("events");
	tD_APD_MCP->GetYaxis()->SetTitleSize(0.055);
	tD_APD_MCP->GetYaxis()->SetTitleOffset(0.75);
	tD_APD_MCP->Fit("gaus", "", "", Xfirst, Xlast);

	TimeResults[APD].Mean = tD_APD_MCP->GetFunction("gaus")->GetParameter(1);
	TimeResults[APD].MeanErr = tD_APD_MCP->GetFunction("gaus")->GetParError(1);
	TimeResults[APD].Sigma = tD_APD_MCP->GetFunction("gaus")->GetParameter(2);
	TimeResults[APD].SigmaErr = tD_APD_MCP->GetFunction("gaus")->GetParError(2);
	
	gStyle->SetOptStat();
	gStyle->SetOptFit();
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD_MCP->Draw();
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-"+MCP+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-"+MCP+"_"+RunStats+".pdf").c_str());

	TFile* f = TFile::Open(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-"+MCP+"_"+RunStats+".root").c_str(), "RECREATE");
	tD_APD_MCP->Write();
	f->Write();
	f->Close();

	c0->~TCanvas();
	tD_APD_MCP->~TH1F();
	f->~TFile();

	return TimeResults[APD];
}

GaussPar TimeAnalysisTools::TimeAPDvsMCPedge(std::string APD, std::string MCP)
{
	SetSelections();

	std::string tD_APD_MCP_Sel = "( fabs(X[0])<5 || fabs(X[1]-("+Xshift_str+"))<5 ) && ( Y[0]<0 || Y[1]-("+Yshift_str+")<0) && (Y[0]>-5 || Y[1]-("+Yshift_str+")>-5 ) && fit_ampl["+APD+"]>400 "+" && "+ DeviceSelections[MCP]; 

	//cout << endl << endl << tD_APD_MCP_Sel << endl << endl;
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

	TimeResults[APD].Mean = tD_APD_MCP->GetFunction("gaus")->GetParameter(1);
	TimeResults[APD].MeanErr = tD_APD_MCP->GetFunction("gaus")->GetParError(1);
	TimeResults[APD].Sigma = tD_APD_MCP->GetFunction("gaus")->GetParameter(2);
	TimeResults[APD].SigmaErr = tD_APD_MCP->GetFunction("gaus")->GetParError(2);
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-"+MCP+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-"+MCP+"_"+RunStats+".pdf").c_str());

	c0->~TCanvas();
	tD_APD_MCP->~TH1F();

	return TimeResults[APD];
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
	TimeXY0->GetZaxis()->SetRangeUser(TimeResults[APD].Mean-3*TimeResults[APD].Sigma, TimeResults[APD].Mean+3*TimeResults[APD].Sigma);	
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
	TimeXY1->GetZaxis()->SetRangeUser(TimeResults[APD].Mean-3*TimeResults[APD].Sigma, TimeResults[APD].Mean+3*TimeResults[APD].Sigma);	
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
	TimeXYM->GetZaxis()->SetRangeUser(TimeResults[APD].Mean-3*TimeResults[APD].Sigma, TimeResults[APD].Mean+3*TimeResults[APD].Sigma);	
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

	tD_APD_MCP_Mean->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_MCP_Mean->GetXaxis()->SetTitle(("t_{"+APD+"}-t_{MCPmean} (ns)").c_str());
	tD_APD_MCP_Mean->GetXaxis()->SetTitleSize(0.055);
	tD_APD_MCP_Mean->GetXaxis()->SetTitleOffset(0.75);
	tD_APD_MCP_Mean->GetYaxis()->SetTitle("events");
	tD_APD_MCP_Mean->GetYaxis()->SetTitleSize(0.055);
	tD_APD_MCP_Mean->GetYaxis()->SetTitleOffset(0.75);
	tD_APD_MCP_Mean->Fit("gaus", "", "", Xfirst, Xlast);

	TimeResults[APD].Mean = tD_APD_MCP_Mean->GetFunction("gaus")->GetParameter(1);
	TimeResults[APD].MeanErr = tD_APD_MCP_Mean->GetFunction("gaus")->GetParError(1);
	TimeResults[APD].Sigma = tD_APD_MCP_Mean->GetFunction("gaus")->GetParameter(2);
	TimeResults[APD].SigmaErr = tD_APD_MCP_Mean->GetFunction("gaus")->GetParError(2);

	gStyle->SetOptStat();
	gStyle->SetOptFit();
	TCanvas* c0 = new TCanvas("c0", "c0");	
	tD_APD_MCP_Mean->Draw();
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-MCP_Mean_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-MCP_Mean_"+RunStats+".pdf").c_str());

	TFile* f = TFile::Open(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD+"-MCP_MEan_"+RunStats+".root").c_str(), "RECREATE");
	tD_APD_MCP_Mean->Write();
	f->Write();
	f->Close();

	c0->~TCanvas();
	tD_APD_MCP_Mean->~TH1F();
	f->~TFile();

	return TimeResults[APD];
}

GaussPar TimeAnalysisTools::TimeAPDvsAPD(std::string APD1, std::string APD2)
{
	if(APDList.size()<2)
	{
		cout << "Not enough APD for Time APD vs APD!!!" << endl; 
		return Default;
	}

	std::string tD_APD_APD_Sel = "(("+DeviceSelections[APD1+"PosSel"]+") || ("+DeviceSelections[APD2+"PosSel"]+")) && (("+DeviceSelections[APD1+"AmpSel"]+") || ("+DeviceSelections[APD2+"AmpSel"]+"))";
	//cout << endl << tD_APD_APD_Sel << endl;

	//define APD_APD time distribution
	TH1F* tD_APD_APD = new TH1F("tD_APD_APD", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD_APD->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_APD->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+APD1+"]-fit_time["+APD2+"]>>tD_APD_APD").c_str(), tD_APD_APD_Sel.c_str());
	
	float Xfirst = tD_APD_APD->GetXaxis()->GetBinCenter(tD_APD_APD->GetMaximumBin())-1;
	float Xlast = tD_APD_APD->GetXaxis()->GetBinCenter(tD_APD_APD->GetMaximumBin())+1;

	//plot and fit APD_APD time distribution
	gStyle->SetOptStat();
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD_APD->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_APD->GetXaxis()->SetTitle(("fit_time["+APD1+"]-fit_time["+APD2+"] (ns)").c_str());
	tD_APD_APD->GetYaxis()->SetTitle("events");
	tD_APD_APD->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_APD->Draw();

	TimeResults[APD1].Mean = tD_APD_APD->GetFunction("gaus")->GetParameter(1);
	TimeResults[APD1].MeanErr = tD_APD_APD->GetFunction("gaus")->GetParError(1);
	TimeResults[APD1].Sigma = tD_APD_APD->GetFunction("gaus")->GetParameter(2);
	TimeResults[APD1].SigmaErr = tD_APD_APD->GetFunction("gaus")->GetParError(2);
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD1+"-"+APD2+"_"+RunStats+".pdf").c_str());

	c0->~TCanvas();
	tD_APD_APD->~TH1F();

	return TimeResults[APD1];
}

GaussPar TimeAnalysisTools::TimeXTALvsXTAL(std::string APD1, std::string APD2)
{
	if(APDList.size()<2)
	{
		cout << "Not enough APD for Time APD vs APD!!!" << endl; 
		return Default;
	}

	/*
	TH1F *SumAmp = new TH1F("SumAmp", "", 2000, 450, 5000);
	h4->Draw(("fit_ampl["+APD1+"]+fit_ampl["+APD2+"]>>SumAmp").c_str());

	float XMax = SumAmp->GetBinCenter(SumAmp->GetMaximumBin());
	float Xfirst = XMax*0.6;
	float Xlast = XMax*1.2;

	SumAmp->Fit("gaus", "RQ", "", Xfirst, Xlast);

	float AmplitudeSUM = SumAmp->GetFunction("gaus")->GetParameter(1);
	float AmplitudeSUMsigma = SumAmp->GetFunction("gaus")->GetParameter(2);
	std::string tD_APD_APD_Sel = "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<800 || fabs(Y[1])<800) && fabs(fit_ampl["+APD1+"]+fit_ampl["+APD2+"]-"+std::to_string(AmplitudeSUM)+")<"+std::to_string(AmplitudeSUMsigma);
	cout << endl << tD_APD_APD_Sel << endl;
	*/

	float AmplitudeAVG = 0.5*(Amplitude[APD1].Mean + Amplitude[APD2].Mean);
	float AmplitudeAVGErr = TMath::Sqrt(Amplitude[APD1].Sigma*Amplitude[APD1].Sigma + Amplitude[APD2].Sigma*Amplitude[APD2].Sigma);
	std::string tD_APD_APD_Sel = "(fabs(X[0]<"+bound_str+") || fabs(X[1]-("+Xshift_str+"))<"+bound_str+") && (fabs(Y[0]<"+bound_str+") || fabs(Y[1]-("+Xshift_str+"))<"+bound_str+") && fabs(fit_ampl["+APD1+"]+fit_ampl["+APD2+"]-"+std::to_string(AmplitudeAVG)+")<"+std::to_string(AmplitudeAVGErr) + " && fit_ampl["+APD1+"]>0.3*"+std::to_string(Amplitude[APD1].Mean)+ " && fit_ampl["+APD2+"]>0.3*"+std::to_string(Amplitude[APD2].Mean);
	//cout << endl << tD_APD_APD_Sel << endl;
	
	//define APD_APD time distribution
	TH1F* tD_APD_APD = new TH1F("tD_APD_APD", "", 2000, -20, 20);
	if(Energy == "20" && Gain!="200") tD_APD_APD->SetBins(750, -40, 40);
	if(Energy == "20" && Gain=="200") tD_APD_APD->SetBins(1000, -40, 40);
	h4->Draw(("fit_time["+APD1+"]-fit_time["+APD2+"]>>tD_APD_APD").c_str(), tD_APD_APD_Sel.c_str());
	
	float Xfirst = tD_APD_APD->GetXaxis()->GetBinCenter(tD_APD_APD->GetMaximumBin())-1;
	float Xlast = tD_APD_APD->GetXaxis()->GetBinCenter(tD_APD_APD->GetMaximumBin())+1;

	//plot and fit APD_APD time distribution
	gStyle->SetOptStat();
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_APD_APD->GetXaxis()->SetRangeUser(Xfirst, Xlast);
	tD_APD_APD->GetXaxis()->SetTitle(("fit_time["+APD1+"]-fit_time["+APD2+"] (ns)").c_str());
	tD_APD_APD->GetYaxis()->SetTitle("events");
	tD_APD_APD->Fit("gaus", "", "", Xfirst, Xlast);
	tD_APD_APD->Draw();

	TimeResults[APD1].Mean = tD_APD_APD->GetFunction("gaus")->GetParameter(1);
	TimeResults[APD1].MeanErr = tD_APD_APD->GetFunction("gaus")->GetParError(1);
	TimeResults[APD1].Sigma = tD_APD_APD->GetFunction("gaus")->GetParameter(2);
	TimeResults[APD1].SigmaErr = tD_APD_APD->GetFunction("gaus")->GetParError(2);
	
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/fitTimeDist/FinalTimeDistribution/Time_"+APD1+"-"+APD2+"_"+RunStats+".pdf").c_str());

	c0->~TCanvas();
	tD_APD_APD->~TH1F();

	return TimeResults[APD1];
}

void TimeAnalysisTools::TimeXTALvsXTALAeff(std::string APD1, std::string APD2)
{
	ComputeAvsNoise(APD1);
	ComputeAvsNoise(APD2);

	if(APDList.size()<2)
	{
		cout << "Not enough APD for Time APD vs APD!!!" << endl; 
	}

	TH1F *SumAmp = new TH1F("SumAmp", "", 2000, 450, 5000);
	if(APD1=="D3" || APD2=="D3") h4->Draw(("fit_ampl["+APD1+"]+fit_ampl["+APD2+"]>>SumAmp").c_str(), "(fabs(Y[0])<10 || fabs(Y[1])<10)");
	else h4->Draw(("fit_ampl["+APD1+"]+fit_ampl["+APD2+"]>>SumAmp").c_str(), "(fabs(X[0])<10 || fabs(X[1])<10)");

	float XMax = SumAmp->GetBinCenter(SumAmp->GetMaximumBin());
	float Xfirst = XMax*0.6;
	float Xlast = XMax*1.2;

	SumAmp->Fit("gaus", "RQ", "", Xfirst, Xlast);

	TCanvas* cAMP = new TCanvas("cAMP", "cAMP");
	cAMP->cd();
	SumAmp->GetXaxis()->SetTitle(("fit_ampl["+APD1+"]+fit_ampl["+APD2+"]").c_str());
	SumAmp->GetYaxis()->SetTitle("events");
	SumAmp->Draw();
	cAMP->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/SumAmp_"+APD1+"+"+APD2+"_"+RunStats+".png").c_str());
	cAMP->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/SumAmp_"+APD1+"+"+APD2+"_"+RunStats+".pdf").c_str());
	cAMP->~TCanvas();

	float AmplitudeSUM = SumAmp->GetFunction("gaus")->GetParameter(1);
	float AmplitudeSUMsigma = SumAmp->GetFunction("gaus")->GetParameter(2);
	std::string tD_XTAL_XTAL_Sel;
	if(APD1=="D3" || APD2=="D3") tD_XTAL_XTAL_Sel = "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<10 || fabs(Y[1])<10) && fabs(fit_ampl["+APD1+"]+fit_ampl["+APD2+"]-"+std::to_string(AmplitudeSUM)+")<"+std::to_string(AmplitudeSUMsigma);
	else tD_XTAL_XTAL_Sel = "(fabs(X[0])<10 || fabs(X[1])<10) && (fabs(Y[0])<800 || fabs(Y[1])<800) && fabs(fit_ampl["+APD1+"]+fit_ampl["+APD2+"]-"+std::to_string(AmplitudeSUM)+")<"+std::to_string(AmplitudeSUMsigma);

	//cout << endl << tD_XTAL_XTAL_Sel << endl;

	int NBins=35;
	//define XTAL_XTAL time distribution wrt Aeff
	TH2F* tD_XTAL_XTAL_Aeff = new TH2F("tD_XTAL_XTAL_Aeff", "", NBins, 0, 150, 200, -2, 2);
	h4->Draw(("(fit_time["+APD1+"]-fit_time["+APD2+"]):1/TMath::Sqrt( pow("+std::to_string(NoiseAmpl[APD1].Sigma)+"/fit_ampl["+APD1+"], 2) + pow("+std::to_string(NoiseAmpl[APD2].Sigma)+"/fit_ampl["+APD2+"], 2) )>>tD_XTAL_XTAL_Aeff").c_str(), tD_XTAL_XTAL_Sel.c_str());
	
	//plot and fit APD_APD time distribution
	gStyle->SetOptStat(0);
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_XTAL_XTAL_Aeff->GetXaxis()->SetTitleSize(0.055);
	tD_XTAL_XTAL_Aeff->GetXaxis()->SetTitleOffset(0.75);
	tD_XTAL_XTAL_Aeff->GetXaxis()->SetTitle("Aeff/#sigma");
	tD_XTAL_XTAL_Aeff->GetYaxis()->SetTitleSize(0.055);
	tD_XTAL_XTAL_Aeff->GetYaxis()->SetTitleOffset(0.75);
	tD_XTAL_XTAL_Aeff->GetYaxis()->SetTitle(("fit_time["+APD1+"]-fit_time["+APD2+"]").c_str());
	tD_XTAL_XTAL_Aeff->Draw("COLZ");

	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/Time_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/Time_"+APD1+"-"+APD2+"_"+RunStats+".pdf").c_str());

	c0->~TCanvas();

	TFile* f = new TFile(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/Time_"+APD1+"-"+APD2+"_"+RunStats+".root").c_str(), "RECREATE");
	f->cd();
	tD_XTAL_XTAL_Aeff->Write();
	f->Write();
	f->Close();

	TObjArray aSlices;
	tD_XTAL_XTAL_Aeff->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	
	TProfile *timeResAeff = new TProfile();
	timeResAeff = (TProfile*)aSlices[2];
	timeResAeff->SetTitle("");
	timeResAeff->GetXaxis()->SetTitleSize(0.055);
	timeResAeff->GetXaxis()->SetTitleOffset(0.75);
	timeResAeff->GetXaxis()->SetTitle("Aeff/#sigma");
	timeResAeff->GetYaxis()->SetTitleSize(0.055);
	timeResAeff->GetYaxis()->SetTitleOffset(0.75);
	timeResAeff->GetYaxis()->SetTitle(("#sigma(time["+APD1+"]-time["+APD2+"])").c_str());
	
	float BinValue;
	float BinError;
	for(int i=0; i<=NBins; i++)
	{
		BinValue=timeResAeff->GetBinContent(i);
		BinError=timeResAeff->GetBinError(i);
		timeResAeff->SetBinContent(i, BinValue*1000);
		timeResAeff->SetBinError(i, BinError*1000);
	}
	
	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + 2*[1]*[1] )", 15, 1000);
	
	fitFunc->SetParLimits(0, 10, 10000);
	fitFunc->SetParLimits(1, -10, 70);
	
	fitFunc->SetParameter(0, 3000);
	fitFunc->SetParameter(1, 20);
	
	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");
	
	if(Energy=="50")timeResAeff->Fit("fitFunc", "", "", 10, timeResAeff->GetBinCenter(timeResAeff->FindLastBinAbove(0))-3);
	else timeResAeff->Fit("fitFunc", "", "", 20, timeResAeff->GetBinCenter(timeResAeff->FindLastBinAbove(0))-5);
	timeResAeff->GetXaxis()->SetRangeUser(5, timeResAeff->GetBinCenter(timeResAeff->FindLastBinAbove(0))+5);

	TCanvas *cTrAeff = new TCanvas("cTrAeff", "cTrAeff");	
	cTrAeff->cd();
	timeResAeff->Draw();
	cTrAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimeRes_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());
	cTrAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimeRes_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());


	TProfile *timeMeanAeff = new TProfile();
	timeMeanAeff = (TProfile*)aSlices[1];
	timeMeanAeff->SetTitle("");
	timeMeanAeff->GetXaxis()->SetTitle("Aeff/#sigma");
	timeMeanAeff->GetYaxis()->SetTitle(("#time["+APD1+"]-time["+APD2+"]").c_str());
	timeMeanAeff->GetXaxis()->SetRangeUser(5, timeMeanAeff->GetBinCenter(timeMeanAeff->FindLastBinAbove(0))+5);

	TCanvas *cTmAeff = new TCanvas("cTmAeff", "cTmAeff");	
	cTmAeff->cd();
	timeMeanAeff->Draw();
	cTmAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimeMean_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());
	cTmAeff->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/TimeMean_"+APD1+"-"+APD2+"_"+RunStats+".png").c_str());

	f->~TFile();
	cTrAeff->~TCanvas();
	tD_XTAL_XTAL_Aeff->~TH2F();
	SumAmp->~TH1F();
	timeResAeff->~TProfile();
}

GaussPar TimeAnalysisTools::ComputeAvsNoise(std::string APD)
{
	GaussPar AmpPar, AeNoise;
	std::string NoiseNtuple; 
		
	if(APD=="C0APD1" || APD=="C0APD2")
	{
		if(Gain == "50") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBtop_Pedestal_5902.root";
		else if(Gain == "100") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBtop_Pedestal_5900.root";
		else if(Gain == "200") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBtop_Pedestal_5898.root";
	}
	else
	{
		if(Gain == "50") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5896.root";
		else if(Gain == "100") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5894.root";
		else if(Gain == "200") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5892.root";
	}
	
	TFile *fNoise = TFile::Open(NoiseNtuple.c_str());
	TTree *h4Noise = (TTree*)fNoise->Get("h4");

	h4Noise->GetEntry(0);
	std::string NoiseGain = std::to_string((int)h4Noise->GetLeaf("CHGain")->GetValue(0));

	if(NoiseGain!=Gain)
	{
		cout << endl << "!!!!!!!!!!!!!!!!!!!! Noise Gain =/= Run Gain !!!!!!!!!!!!!!!!!!!!!!11111111" << endl << endl;
		return Default;
	}

	AmpPar = AmplitudeDistributionFit2Gaus(APD);
	NoiseAmpl[APD] = NoiseAmplitudeDistributionFit(h4Noise, APD);

	AeNoise.Mean = AmpPar.Mean;
	AeNoise.MeanErr = AmpPar.MeanErr;
	AeNoise.Sigma = NoiseAmpl[APD].Sigma;
	AeNoise.SigmaErr = NoiseAmpl[APD].SigmaErr;
	/*
	*ANoiseR = AmpPar.Mean/NoisePar.Sigma;
	*ANoiseRErr = TMath::Sqrt(pow(AmpPar.MeanErr/AmpPar.Mean, 2) + pow(NoisePar.SigmaErr/NoisePar.Sigma, 2))*AmpPar.Mean/NoisePar.Sigma;
	*/
	return AeNoise;
}

void TimeAnalysisTools::ComputeAvsNoiseEdge(std::string APD1, std::string APD2)
{
	GaussPar AmpPar, NoisePar1, NoisePar2;
	float XMax, YMax;
	std::string TimeMCP, TimeShift, NoiseNtuple; 
		
	
	if(Gain == "50") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5896.root";
	else if(Gain == "100") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5894.root";
	else if(Gain == "200") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5892.root";
	
	
	TFile *fNoise = TFile::Open(NoiseNtuple.c_str());
	TTree *h4Noise = (TTree*)fNoise->Get("h4");

	h4Noise->GetEntry(0);
	std::string NoiseGain = std::to_string((int)h4Noise->GetLeaf("CHGain")->GetValue(0));

	if(NoiseGain!=Gain)
	{
		cout << endl << "!!!!!!!!!!!!!!!!!!!! Noise Gain =/= Run Gain !!!!!!!!!!!!!!!!!!!!!!11111111" << endl << endl;
	}

	NoisePar1 = NoiseAmplitudeDistributionFit(h4Noise, APD1);
	NoisePar2 = NoiseAmplitudeDistributionFit(h4Noise, APD2);
	
	TH1F* Aeff = new TH1F("Aeff", "", 100, 0, 130);
	h4->Draw(("1/TMath::Sqrt( pow("+std::to_string(NoisePar1.Sigma)+"/fit_ampl["+APD1+"], 2) + pow("+std::to_string(NoisePar2.Sigma)+"/fit_ampl["+APD2+"], 2) )>>Aeff").c_str(), "(fabs(X[0])<800 || fabs(X[1])<800) && (fabs(Y[0])<800 || fabs(Y[1])<800)"); //&& (fit_ampl["+APD1+"] + fit_ampl["+APD2+"])>"+std::to_string(0.7*(Amplitude[APD1].Mean+Amplitude[APD2].Mean))).c_str());
	Aeff->GetXaxis()->SetRangeUser(5, Aeff->GetBinCenter(Aeff->FindLastBinAbove(0))+5);
	Aeff->GetXaxis()->SetTitleSize(0.055);
	Aeff->GetXaxis()->SetTitleOffset(0.75);
	Aeff->GetXaxis()->SetTitle("Aeff/#sigma(Noise)");
	Aeff->GetYaxis()->SetTitleSize(0.055);
	Aeff->GetYaxis()->SetTitleOffset(0.75);
	Aeff->GetYaxis()->SetTitle("events");
 
	TCanvas *cH = new TCanvas("cH", "cH");
	Aeff->Draw();
	cH->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/Aeff_"+APD1+"_"+APD2+"_"+RunStats+".png").c_str());
	cH->SaveAs(("/afs/cern.ch/user/c/cquarant/www/Amp_plot/Aeff_"+APD1+"_"+APD2+"_"+RunStats+".pdf").c_str());

	cH->~TCanvas();
	Aeff->~TH1F();
}

void TimeAnalysisTools::DrawFreqSpec(std::string APD)
{
	gStyle->SetOptStat(0);
	
	TProfile* FS = new TProfile( "FS_Signal", "", 512, 0, 512);
	std::string FS_APD_MCP_Sel = DeviceSelections[APD] + " && " + DeviceSelections["MCP1"] + " && ch=="+APD;
	
	//Good event's signal frequency spectrum
	cout << "DRAWING..." << endl;
	h4->Draw("ampl:freq>>FS_Signal", FS_APD_MCP_Sel.c_str()); 
	
	FS->GetXaxis()->SetTitle("Frequency");
	FS->GetYaxis()->SetTitle("Amplitude");
	FS->SetMarkerStyle(kFullDotSmall);
	FS->SetMarkerColor(kRed);

	std::string NoiseNtuple;
	if(APD=="C0APD1" || APD=="C0APD2")
	{
		if(Gain == "50") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBtop_Pedestal_5902.root";
		else if(Gain == "100") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBtop_Pedestal_5900.root";
		else if(Gain == "200") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBtop_Pedestal_5898.root";
	}
	else
	{
		if(Gain == "50") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5896.root";
		else if(Gain == "100") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5894.root";
		else if(Gain == "200") NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5892.root";
	}	

	TFile* fNoise = TFile::Open(NoiseNtuple.c_str());
	TTree* h4Noise = (TTree*)fNoise->Get("h4");
	TProfile* FSPedestal = new TProfile("FS_Pedestal", "", 512, 0, 512);
	FSPedestal->SetMarkerStyle(kFullDotSmall);
	FSPedestal->SetMarkerColor(kBlue);
	
	cout << "DRAWING PEDESTAL..." << endl;
	h4Noise->Draw("ampl:freq>>FS_Pedestal", ("ch=="+APD).c_str());
	
	//---Ratio plot
	TProfile* Ratio = new TProfile();
	Ratio = (TProfile*)FS->Clone();	
	Ratio->SetName("Ratio");
	Ratio->Divide(FSPedestal);

	Ratio->GetXaxis()->SetRangeUser(0, 40);
	Ratio->GetXaxis()->SetTitle("Frequency");
	Ratio->GetXaxis()->SetTitleSize(0.055);
	Ratio->GetXaxis()->SetTitleOffset(0.75);

	//Ratio->GetYaxis()->SetRangeUser(0.9, 1.1);
	Ratio->GetYaxis()->SetTitle("A_{signal}/A_{pedestal}");
	Ratio->GetYaxis()->SetTitleSize(0.055);
	Ratio->GetYaxis()->SetTitleOffset(0.75);

	TLegend* legend = new TLegend(0.52, 0.7, 0.9, 0.9);
	legend->SetHeader("Signal Frequency Spectrum");
	legend->AddEntry(FS, "Physics run", "p");
	legend->AddEntry(FSPedestal, "Pedestal run", "p");

	TCanvas* c0 = new TCanvas("c0", "c0");
	TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.5,1.00,1.00);
	TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.5);
	cUp->SetBottomMargin(0.01); 	
    	cUp->Draw();
	cDown->SetTopMargin(0.01);
	cDown->Draw();

	cUp->cd();	
	cUp->SetLogy();
	FS->Draw();
	FSPedestal->Draw("SAME");
	legend->Draw("SAME");
	
	cDown->cd();
	Ratio->Draw("hist p");

	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/"+APD+"/FS_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/"+APD+"/FS_"+RunStats+".pdf").c_str());

	TFile* f = TFile::Open(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/"+APD+"/FS_"+RunStats+".root").c_str(), "RECREATE");
	f->cd();
	FS->Write();
	FSPedestal->Write();
	Ratio->Write();
	c0->Write();

	f->Write();
	f->Close();

	c0->~TCanvas();
	FS->~TProfile();
	f->~TFile();
}

void TimeAnalysisTools::DrawFreqSpecMCP(std::string MCP)
{
	gStyle->SetOptStat(0);
	
	TProfile* FS = new TProfile("FS_Signal", "", 512, 0, 512);
	std::string FS_MCP_Sel = DeviceSelections[MCP] + " && ch=="+MCP;
	
	//Good event's signal frequency spectrum
	cout << "DRAWING..." << endl;
	h4->Draw("ampl:freq>>FS_Signal", FS_MCP_Sel.c_str()); 
	
	FS->GetXaxis()->SetTitle("Frequency");
	FS->GetYaxis()->SetTitle("Amplitude");
	FS->SetMarkerStyle(kFullDotSmall);
	FS->SetMarkerColor(kRed);

	cout << "DRAWING PEDESTAL..." << endl;

	/*
	std::string NoiseNtuple = "/eos/cms/store/user/meridian/ECALTBH4/cquarant/H42016_MBbot_Pedestal_5892.root";
	TFile* fNoise = TFile::Open(NoiseNtuple.c_str());
	TTree* h4Noise = (TTree*)fNoise->Get("h4");
	*/
	
	TFile* fNoise;

	if(Gain == "50") fNoise = TFile::Open("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/C3/FS_20Gev_G50_5786.root");
	else if(Gain == "100") fNoise = TFile::Open("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/C3/FS_20Gev_G100_5785.root");
	else fNoise = TFile::Open("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/C0APD1/FS_100Gev_G200_5782.root");
	
 	

	TProfile* FSPedestal =(TProfile*)fNoise->Get("FS_Pedestal");
	FSPedestal->SetMarkerStyle(kFullDotSmall);
	FSPedestal->SetMarkerColor(kBlue);
	
	
	//h4Noise->Draw(("ampl:freq>>FS_"+MCP+"_Pedestal_"+RunStats).c_str(), "ch == C3");
	
	//---Ratio plot
	TProfile* Ratio = new TProfile();
	Ratio = (TProfile*)FS->Clone();	
	Ratio->SetName("Ratio");
	Ratio->Divide(FSPedestal);

	Ratio->GetXaxis()->SetRangeUser(0, 160);
	Ratio->GetXaxis()->SetTitle("Frequency");
	Ratio->GetXaxis()->SetTitleSize(0.055);
	Ratio->GetXaxis()->SetTitleOffset(0.75);

	//Ratio->GetYaxis()->SetRangeUser(0.9, 1.1);
	Ratio->GetYaxis()->SetTitle("A_{signal}/A_{pedestal}");
	Ratio->GetYaxis()->SetTitleSize(0.055);
	Ratio->GetYaxis()->SetTitleOffset(0.75);

	TLegend* legend = new TLegend(0.52, 0.7, 0.9, 0.9);
	legend->SetHeader("Signal Frequency Spectrum");
	legend->AddEntry(FS, "Physics run", "p");
	legend->AddEntry(FSPedestal, "Pedestal run", "p");

	TCanvas* c0 = new TCanvas("c0", "c0");
	TPad *cUp  = new TPad("pad_0","pad_0",0.00,0.5,1.00,1.00);
	TPad *cDown = new TPad("pad_1","pad_1",0.00,0.00,1.00,0.5);
	cUp->SetBottomMargin(0.01); 	
    	cUp->Draw();
	cDown->SetTopMargin(0.01);
	cDown->Draw();

	cUp->cd();	
	cUp->SetLogy();
	FS->Draw();
	FSPedestal->Draw("SAME");
	legend->Draw("SAME");
	
	cDown->cd();
	Ratio->Draw("hist p");

	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/"+MCP+"/FS_"+RunStats+".png").c_str());
	c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/"+MCP+"/FS_"+RunStats+".pdf").c_str());

	TFile* f = TFile::Open(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/"+MCP+"/FS_"+RunStats+".root").c_str(), "RECREATE");
	f->cd();
	FS->Write();
	FSPedestal->Write();
	Ratio->Write();
	c0->Write();
	f->Write();
	f->Close();

	c0->~TCanvas();
	FS->~TProfile();
	f->~TFile();
}


void TimeAnalysisTools::DrawFreqSpecPedestal(std::string APD)
{
	TFile* f = TFile::Open("FourierSpectra.root", "UPDATE");
	TCanvas* c0 = new TCanvas("c0", "c0");
	TProfile* FS = new TProfile(("FS_"+APD+"_Pedestal_"+RunStats).c_str(), "", 512, 0, 512);

	if(!f->GetListOfKeys()->Contains(("FS_"+APD+"_Pedestal_"+RunStats).c_str()))
	{
		//PedestalF frequency spectrum
		cout << "DRAWING..." << endl;
		h4->Draw(("ampl:freq>>FS_"+APD+"_Pedestal_"+RunStats).c_str());
		
		//plot and fit APD_MCP1 time distribution
		gStyle->SetOptStat();
		FS->GetXaxis()->SetTitle("Frequency");
		FS->GetYaxis()->SetTitle("Amplitude");
		FS->Draw();

		c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/FS_"+APD+"_Pedestal_"+RunStats+".png").c_str());
		c0->SaveAs(("/afs/cern.ch/user/c/cquarant/www/FourierSpectra/FS_"+APD+"_Pedestal_"+RunStats+".pdf").c_str());
		
		f->cd();
		FS->Write();
	}
	
	c0->~TCanvas();
	FS->~TProfile();
	f->Close();
}
		
