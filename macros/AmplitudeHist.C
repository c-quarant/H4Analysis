float AmplitudeHist(TTree *h4, std::string detector, std::string Selection, std::string pathToOut,  std::string RunStats)
{
	TH1F* HAmp = new TH1F("HAmp", "", 2500, -50, 5000);

	h4->Draw(("amp_max["+detector+"]>>HAmp").c_str(), Selection.c_str());
	
	HAmp->GetXaxis()->SetRange(-50, (int)(HAmp->GetMaximumBin()*1.10));
	HAmp->Fit("gaus");
	TCanvas* c1 = new TCanvas("c1", "c1");
	HAmp->Draw();

	c1->SaveAs((pathToOut+"Amp_plot/Amp_"+detector+"_"+RunStats+".png").c_str());
	c1->SaveAs((pathToOut+"Amp_plot/Amp_"+detector+"_"+RunStats+".pdf").c_str());

	return HAmp->GetFunction("gaus")->GetParameter(1) - 3*HAmp->GetFunction("gaus")->GetParameter(2);
}
