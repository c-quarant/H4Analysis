#include "XtalXtalTimeTools.cc"

void AeffMeanDistribution(TH1F* Aeff, int NSlices, std::vector<float>* SliceAeffMean, std::vector<float>* SliceAeffMeanErr);
void MyFitSlicesY(TH2F* Xtal_Xtal_Time, std::vector<GaussPar>* SliceGaussFit);
void TimeMeanvsAeff(std::vector<GaussPar>* SliceGaussFit, std::vector<float>* SliceAeffMean, std::vector<float>* SliceAeffMeanErr);

void XtalXtalTimeResolutionAllEne()
{
	int i, j;
	int Run, detector, TimeRes;
	float Energy, Gain;
	float ANoiseR, ANoiseRErr;
	GaussPar TimePar, AeNoisePar, AeNoisePar2;

	TFile* fTH2 = TFile::Open("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/Time_C3-D3_AllEne.root");
	TH2F* tD_XTAL_XTAL_Aeff = (TH2F*)fTH2->Get("tD_XTAL_XTAL_Aeff");

	//plot Xtal-Xtal time distribution
	gStyle->SetOptStat(0);
	TCanvas* c0 = new TCanvas("c0", "c0");
	tD_XTAL_XTAL_Aeff->Draw("COLZ");
	c0->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/Time_C3-D3_AllEne.png");
	c0->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/Time_C3-D3_AllEne.pdf");
	c0->~TCanvas();
	
	
	TFile* fAeff = TFile::Open("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/Aeff_C3_D3_AllEnergies.root");
	TH1F* Aeff = (TH1F*)fAeff->Get("Aeff"); 

	TCanvas* cAeff = new TCanvas("cAeff", "cAeff");
	Aeff->Draw();
	cAeff->SaveAs("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/Aeff_C3-D3_AllEne.png");
	cAeff->SaveAs("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/Aeff_C3-D3_AllEne.pdf");
	cAeff->~TCanvas();


	TFile* f200Gev = TFile::Open("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/TimeResvsAeff_C3-D3_200Gev_G50_5624.root");
	TGraphErrors* Res200Gev = (TGraphErrors*)f200Gev->Get("Graph");

	std::vector<float> SliceAeffMean, SliceAeffMeanErr;
	AeffMeanDistribution(Aeff, tD_XTAL_XTAL_Aeff->GetNbinsX(), &SliceAeffMean, &SliceAeffMeanErr);
	std::vector<GaussPar> SliceGaussFit;

	MyFitSlicesY(tD_XTAL_XTAL_Aeff, &SliceGaussFit);
	TimeMeanvsAeff(&SliceGaussFit, &SliceAeffMean, &SliceAeffMeanErr);


	std::vector<float> SliceResolutions, SliceResolutionsErr;
	for(auto& SliceRes : SliceGaussFit)
	{
		SliceResolutions.push_back(1000*SliceRes.Sigma);
		SliceResolutionsErr.push_back(1000*SliceRes.SigmaErr);
	}

	TF1 *fitFunc = new TF1("fitFunc", "TMath::Sqrt([0]*[0]/(x*x) + TMath::Sqrt(2)*[1]*[1] )", 15, 1000);
	fitFunc->SetLineColor(kBlack);	

	fitFunc->SetParLimits(0, 10, 10000);
	fitFunc->SetParLimits(1, 20, 70);

	fitFunc->SetParameter(0, 3500);
	fitFunc->SetParameter(1, 40);

	fitFunc->SetParName(0, "Noise");
	fitFunc->SetParName(1, "const");


	float Xfirst = 0;
	float Xlast = 0;	
	for(unsigned int i=0; i<SliceResolutions.size(); i++)
	{	
		if(Xfirst==0 & SliceResolutions[i]>0) Xfirst = SliceAeffMean[i]-0.5;
		if(SliceResolutions[i]>0) Xlast = SliceAeffMean[i]+0.5;
	}
	cout << endl << Xfirst << endl;
	cout << Xlast << endl << endl;

	TGraphErrors *ResvsAeff = new TGraphErrors(SliceResolutions.size(), &SliceAeffMean[0], &SliceResolutions[0], &SliceAeffMeanErr[0], &SliceResolutionsErr[0]);
	//ResvsAeff->GetXaxis()->SetRangeUser(1, 140);
	//ResvsAeff->GetYaxis()->SetRangeUser(1, 150);	

	ResvsAeff->GetXaxis()->SetTitle("Aeff/#sigma");
	ResvsAeff->GetYaxis()->SetTitle("#sigma(t_{C3}-t_{D3}) (ps)");	
	ResvsAeff->GetXaxis()->SetTitleSize(0.055);
	ResvsAeff->GetXaxis()->SetTitleOffset(0.75);
	ResvsAeff->GetYaxis()->SetTitleSize(0.055);
	ResvsAeff->GetYaxis()->SetTitleOffset(0.75);
	ResvsAeff->SetMarkerStyle(kFullSquare);
	ResvsAeff->SetMarkerSize(0.5);

	Double_t X200, Y200, XALL, YALL;
	cout << "Index\tX200Gev\tY200Gev\tXAllEne\tYAllEne" << endl;
	for(int i=0; i<ResvsAeff->GetN(); i++)
	{
		Res200Gev->GetPoint(i, X200, Y200);
		ResvsAeff->GetPoint(i, XALL, YALL);
		cout << i << "\t" << X200 << "\t" << Y200 << "\t" << XALL << "\t" << YALL << endl;
	}

	Res200Gev->SetMarkerStyle(kFullCircle);
	Res200Gev->SetMarkerSize(0.5);
	Res200Gev->GetXaxis()->SetRangeUser(1, 140);
	Res200Gev->GetYaxis()->SetRangeUser(10, 500);	
	Res200Gev->SetMarkerColor(kRed);

	//ResvsAeff->SetPointError(23, ResvsAeff->GetErrorX(23), ResvsAeff->GetErrorY(23)*100);
	ResvsAeff->Fit("fitFunc", "", "", Xfirst, Xlast);

	gStyle->SetOptFit();
	TCanvas* cResvsAeff = new TCanvas("ResvsAeff", "ResvsAeff");
	Res200Gev->Draw("AP");	
	ResvsAeff->Draw("SAMEP");
	cResvsAeff->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/TimeResvsAeff_C3-D3_AllEne.png");
	cResvsAeff->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/TimeResvsAeff_C3-D3_AllEne.pdf");
	cResvsAeff->~TCanvas();

	TFile* f1 = new TFile("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/TimeResvsAeff_C3-D3_AllEne.root", "RECREATE");
	f1->cd();
	ResvsAeff->Write();
	f1->Write();
	f1->Close();
	f1->~TFile();

	//for(i=0; i<SliceResolutions.size(); i++)
		//cout << "Aeff " << SliceAeffMean[i] << "  Res " << SliceResolutions[i] << endl;

}


void AeffMeanDistribution(TH1F* Aeff, int NSlices, std::vector<float>* SliceAeffMean, std::vector<float>* SliceAeffMeanErr)
{
	int i;
	int AeffNBins = Aeff->GetNbinsX();
	//Calculate mean value of Aeff Slice by Slice and storing it in a TH1
	TH1F *hAeffMean = new TH1F("hAeffMean", "", NSlices, 0, 140);

	float BinInASlice = AeffNBins/NSlices;
	//cout << endl << "AeffBins "<< AeffNBins << "  AeffMax " << Aeff->GetBinCenter(AeffNBins) << "   Bin in a Slice " << (float)AeffNBins/NSlices << endl; 
	for(i=0; i<=NSlices; i++)
	{
		Aeff->GetXaxis()->SetRange((int)i*BinInASlice, (int)(i+1)*BinInASlice);
		SliceAeffMean->push_back(Aeff->GetMean());
		SliceAeffMeanErr->push_back(Aeff->GetMeanError());
		hAeffMean->SetBinContent(i, Aeff->GetMean());
		hAeffMean->SetBinError(i, Aeff->GetMeanError());
		//cout << "Aeff Low " << Aeff->GetBinCenter(i*BinInASlice) << "  Aeff high " << Aeff->GetBinCenter((i+1)*BinInASlice) << " AeffMean " << Aeff->GetMean() << endl;
	}
	
	TCanvas *cAeffMean = new TCanvas("cAeffMean", "cAeffMean");
	cAeffMean->cd();
	hAeffMean->Draw();
	cAeffMean->SaveAs("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/AeffMean_C3-D3_AllEne.png");
	cAeffMean->SaveAs("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/AeffMean_C3-D3_AllEne.pdf");
	cAeffMean->~TCanvas();

	TFile* fAeffMean = new TFile("/afs/cern.ch/user/c/cquarant/www/EffectiveAmplitude/AeffMean_C3-D3_AllEne.root", "RECREATE");
	fAeffMean->cd();
	hAeffMean->Write();
	fAeffMean->Write();
	fAeffMean->Close();

	fAeffMean->~TFile();
}

//Divide a TH2 in slices along Y (one slice for each bin), make a TH1 for each slcice and fit each one of them with gaus, storing each TH1 and fit	
void MyFitSlicesY(TH2F* Xtal_Xtal_Time, std::vector<GaussPar>* SliceGaussFit)
{      
	GaussPar ResTemp;
	int i,j;
	int NBinsX = Xtal_Xtal_Time->GetNbinsX();
	int NBinsy = Xtal_Xtal_Time->GetNbinsY();

	gStyle->SetOptStat();
	gStyle->SetOptFit();

	gSystem->mkdir("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/AllEne_G50_C3-D3");
	gSystem->CopyFile("/afs/cern.ch/user/c/cquarant/www/Amp_plot/index.php", "/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/AllEne_G50_C3-D3/index.php");

	TH1D* hTemp;
	TH1D* hSlices[1000];
	TCanvas* cOUT = new TCanvas("cOUT", "cOUT");
	TObjArray* f = new TObjArray();

	std::string hSliceName;
	for(i=1; i<NBinsX+1; i++)
	{
		if(i<10)hSliceName	= "TimeResSlice_0"+std::to_string(i);
		else	hSliceName	= "TimeResSlice_"+std::to_string(i);
		hSlices[i]		= new TH1D(hSliceName.c_str(), "", 200, -2, 2);
		hSlices[i]		= Xtal_Xtal_Time->ProjectionY(hSliceName.c_str(), i, i, "e");
		hSlices[i]->GetXaxis()->SetTitle("t_{C3}-t_{D3} (ps)");
		hSlices[i]->GetYaxis()->SetTitle("events");	
		hSlices[i]->GetXaxis()->SetTitleSize(0.055);
		hSlices[i]->GetXaxis()->SetTitleOffset(0.75);
		hSlices[i]->GetYaxis()->SetTitleSize(0.055);
		hSlices[i]->GetYaxis()->SetTitleOffset(0.75);

		if(hSlices[i]->GetEntries()>20)
		{	
			hSlices[i]->Fit("gaus", "Q");
			ResTemp.Mean		= hSlices[i]->GetFunction("gaus")->GetParameter(1);
			ResTemp.MeanErr		= hSlices[i]->GetFunction("gaus")->GetParError(1); 
			ResTemp.Sigma		= hSlices[i]->GetFunction("gaus")->GetParameter(2);
			ResTemp.SigmaErr	= hSlices[i]->GetFunction("gaus")->GetParError(2); 
		}
		else
		{
			ResTemp.Mean		= 0;
			ResTemp.MeanErr		= 0; 
			ResTemp.Sigma		= 0;
			ResTemp.SigmaErr	= 0;
		}
	
		cOUT->cd();
		hSlices[i]->Draw();
		cOUT->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/AllEne_G50_C3-D3/"+hSliceName+"_C3-D3.pdf").c_str());
		cOUT->SaveAs(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/AllEne_G50_C3-D3/"+hSliceName+"_C3-D3.png").c_str());
		
		f->Add(TFile::Open(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/AllEne_G50_C3-D3/"+hSliceName+"_C3-D3.root").c_str(), "RECREATE"));	
		((TFile*)f->FindObject(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/AllEne_G50_C3-D3/"+hSliceName+"_C3-D3.root").c_str()))->cd();
		hSlices[i]->Draw();
		hSlices[i]->Write();
		((TFile*)f->FindObject(("/afs/cern.ch/user/c/cquarant/www/TimeDistSlices/AllEne_G50_C3-D3/"+hSliceName+"_C3-D3.root").c_str()))->Close();

		SliceGaussFit->push_back(ResTemp);
		hSlices[i]->~TH1D();
	}

	//hSlice->~TH1D();
	//f->~TFile();
	cOUT->~TCanvas();	
}

void TimeMeanvsAeff(std::vector<GaussPar>* SliceGaussFit, std::vector<float>* SliceAeffMean, std::vector<float>* SliceAeffMeanErr)
{
	int NBins=35;
	std::vector<float> SliceTimeMean, SliceTimeMeanErr;
	for(auto& SliceTime : *SliceGaussFit)
	{
		SliceTimeMean.push_back(1000*SliceTime.Mean);
		SliceTimeMeanErr.push_back(1000*SliceTime.MeanErr);
	}
	
	TGraphErrors *TimevsAeff = new TGraphErrors(SliceTimeMean.size(), &SliceAeffMean->at(0), &SliceTimeMean[0], &SliceAeffMeanErr->at(0), &SliceTimeMeanErr[0]);
	//TimevsAeff->GetXaxis()->SetRangeUser(SliceAeffMean[5]-1, SliceAeffMean[NBins-1]+10);
	//TimevsAeff->GetYaxis()->SetRangeUser(740, 860);
	TimevsAeff->GetXaxis()->SetTitle("Aeff/#sigma");
	TimevsAeff->GetYaxis()->SetTitle("#bar{t_{C3}-t_{D3}} (ps)");	
	TimevsAeff->GetXaxis()->SetTitleSize(0.055);
	TimevsAeff->GetXaxis()->SetTitleOffset(0.75);
	TimevsAeff->GetYaxis()->SetTitleSize(0.055);
	TimevsAeff->GetYaxis()->SetTitleOffset(0.75);
	
	TCanvas* cTimevsAeff = new TCanvas("TimevsAeff", "TimevsAeff");
	TimevsAeff->Draw("AP");
	cTimevsAeff->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/TimevsAeff_C3-D3_AllEne.png");
	cTimevsAeff->SaveAs("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/TimevsAeff_C3-D3_AllEne.pdf");
	cTimevsAeff->~TCanvas();

	TFile* f2 = new TFile("/afs/cern.ch/user/c/cquarant/www/TimeResolutionAeff/D3_C3/TimevsAeff_FarFromGap_C3-D3_AllEne.root", "RECREATE");
	f2->cd();
	TimevsAeff->Write();
	f2->Write();
	f2->Close();
	f2->~TFile();
}
