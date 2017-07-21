void AmplitudeHist(TTree *h4, std::string detector, std::string Selection, std::string pathToOut,  std::string RunStats, float* AmpMean, float* AmpSigma)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -50, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.80), (int)(HAmp->GetMaximumBin()*1.10));

	HAmp->Fit("gaus", "N");
	TCanvas* c1 = new TCanvas("c1", "c1");
	HAmp->GetXaxis()->SetRange((int)(HAmp->GetMaximumBin()*0.50), (int)(HAmp->GetMaximumBin()*1.50));
	HAmp->Draw();

	c1->SaveAs((pathToOut+"Amp_plot/Amp_"+detector+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOut+"Amp_plot/Amp_"+detector+"_"+RunStats+".pdf").c_str());
	
	*AmpMean = HAmp->GetFunction("gaus")->GetParameter(1);
	*AmpSigma = HAmp->GetFunction("gaus")->GetParameter(2);
}
