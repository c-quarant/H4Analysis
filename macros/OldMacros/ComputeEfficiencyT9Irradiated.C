#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"

#include<iostream>
#include<string>
#include<fstream>

TString AddSelection(TTree*, std::string, TString, std::string);

using namespace std;

void ComputeEfficiencyT9Irradiated() {

  gStyle->SetOptStat(0);

  // inputs
  TFile* inFile = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/TimingTB_T9_Apr2016/ntuples/v1/ntuples_v1/t92016_IRRAD_scan1_MERGED.root");
  TTree* h4 = (TTree*)inFile->Get("h4"); 
  
  TFile* inFile_noIrr = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/TimingTB_T9_Apr2016/ntuples/v1/ntuples_v1/t92016_IRRAD_scan2_MERGED.root");
  TTree* h4_noIrr = (TTree*)inFile_noIrr->Get("h4"); 

  TFile* inFile_PMT = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/TimingTB_T9_Apr2016/ntuples/v1/ntuples_v1/t92016_IRRAD_scan3_MERGED.root");
  TTree* h4_PMT = (TTree*)inFile_PMT->Get("h4"); 

  TFile* inFile_H4_scan1 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/ECALTB_H4_Fall2015/ntuples_Irradiated/v1/ntuples_v1/h42015_50GeV_IRRAD_MERGED_1.root");
  TTree* h4_H4_scan1 = (TTree*)inFile_H4_scan1->Get("h4"); 
  
  TFile* inFile_H4_scan2 = TFile::Open("root://eoscms.cern.ch//store/user/bmarzocc/ECALTB_H4_Fall2015/ntuples_Irradiated/v1/ntuples_v1/h42015_50GeV_IRRAD_MERGED_2.root");
  TTree* h4_H4_scan2 = (TTree*)inFile_H4_scan2->Get("h4"); 
 
  // histos
  TH1F* num = new TH1F("num","",10000,1800.,3100.); 
  TH1F* denum = new TH1F("denum","",10000,1800.,3100.);

  TH1F* num_noIrr = new TH1F("num_noIrr","",10000,1800.,3100.); 
  TH1F* denum_noIrr = new TH1F("denum_noIrr","",10000,1800.,3100.); 

  TH1F* num_PMT = new TH1F("num_PMT","",10000,1800.,3100.); 
  TH1F* denum_PMT = new TH1F("denum_PMT","",10000,1800.,3100.); 
  
  TH1F* num_H4_scan1 = new TH1F("num_H4_scan1","",10000,1800.,3100.); 
  TH1F* denum_H4_scan1 = new TH1F("denum_H4_scan1","",10000,1800.,3100.);  

  TH1F* num_H4_scan2 = new TH1F("num_H4_scan2","",10000,1800.,3100.); 
  TH1F* denum_H4_scan2 = new TH1F("denum_H4_scan2","",10000,1800.,3100.);  
  
  // denominators
  TString commonDen = "amp_max[Rm2]>200 && fabs(time_max[Rm2])<150";
  h4->Project("denum","HVMCP*9/11",commonDen);

  TString commonDen_noIrr = "amp_max[MCP]>30 && fabs(time_max[MCP])<150";
  h4_noIrr->Project("denum_noIrr","HVMCP*9/11",commonDen_noIrr);
  
  TString commonDen_PMT = "amp_max[Rm2]>200 && fabs(time_max[Rm2])<150";
  h4_PMT->Project("denum_PMT","HVMCP*9/11",commonDen_PMT);

  TString commonDen_H4_scan1 = "amp_max[MiB2]>200 && fabs(time_max[MiB2])<150 && X[0]>-800 && Y[0]>-800 && X[0]>-4 && X[0]<2 && Y[0]>-0 && Y[0]<5";
  h4_H4_scan1->Project("denum_H4_scan1","HVMCP*10/11",commonDen_H4_scan1);
 
  TString commonDen_H4_scan2 = "amp_max[MiB2]>200 && fabs(time_max[MiB2])<150 && X[0]>-800 && Y[0]>-800 && X[0]>-4 && X[0]<2 && Y[0]>-0 && Y[0]<5";
  h4_H4_scan2->Project("denum_H4_scan2","HVMCP*10/11",commonDen_H4_scan2);
  
  // numerators
  TString commonNum = commonDen + " && amp_max[MCP]>20.";
  commonNum = AddSelection(h4,"time[MCP]-time[Rm2]",commonNum,"1.");
  std::cout << "CommonNum = " << commonNum << std::endl;
  h4->Project("num","HVMCP*9/11",commonNum); 

  TString commonNum_noIrr = commonDen_noIrr + " && amp_max[Rm2]>18.";
  commonNum_noIrr = AddSelection(h4_noIrr,"time[MCP]-time[Rm2]",commonNum_noIrr,"1.");
  std::cout << "CommonNum_noIrr = " << commonNum_noIrr << std::endl;
  h4_noIrr->Project("num_noIrr","HVMCP*9/11",commonNum_noIrr); 

  //TString commonNum_PMT = commonDen_PMT + " && amp_max[MCP]>19.";
  TString commonNum_PMT = commonDen_PMT + " && amp_max[MCP]>21.";
  commonNum_PMT = AddSelection(h4_PMT,"time[MCP]-time[Rm2]",commonNum_PMT,"1.");
  std::cout << "CommonNum_PMT = " << commonNum_PMT << std::endl;
  h4_PMT->Project("num_PMT","HVMCP*9/11",commonNum_PMT); 

  //TString commonNum_H4_scan1 = commonDen_H4_scan1 + " && amp_max[MCP]>18.";
  TString commonNum_H4_scan1 = commonDen_H4_scan1 + " && amp_max[MCP]>14.5";
  commonNum_H4_scan1  = AddSelection(h4_H4_scan1,"time[MCP]-time[MiB2]",commonNum_H4_scan1 ,"1.");
  h4_H4_scan1->Project("num_H4_scan1","HVMCP*10/11",commonNum_H4_scan1); 

  //TString commonNum_H4_scan2 = commonDen_H4_scan2 + " && amp_max[MCP]>18.";
  TString commonNum_H4_scan2 = commonDen_H4_scan2 + " && amp_max[MCP]>14.5";
  commonNum_H4_scan2  = AddSelection(h4_H4_scan2,"time[MCP]-time[MiB2]",commonNum_H4_scan2 ,"1.");
  h4_H4_scan2->Project("num_H4_scan2","HVMCP*10/11",commonNum_H4_scan2); 
  
  std::cout << "commonNum_H4 = " << commonNum_H4_scan2  << std::endl;
 
  // efficiencies
  TGraphAsymmErrors* eff = new TGraphAsymmErrors(num,denum);
  
  eff->SetMarkerColor(kViolet+2);
  eff->SetLineColor(kViolet+2);
  eff->SetMarkerStyle(20);
  eff->SetLineStyle(8);

  TGraphAsymmErrors* eff_noIrr = new TGraphAsymmErrors(num_noIrr,denum_noIrr);
  
  eff_noIrr->SetMarkerColor(kRed+1);
  eff_noIrr->SetLineColor(kRed+1);
  eff_noIrr->SetMarkerStyle(20);
  eff_noIrr->SetLineStyle(2);

  TGraphAsymmErrors* eff_PMT = new TGraphAsymmErrors(num_PMT,denum_PMT);
  
  eff_PMT->SetMarkerColor(kMagenta+1);
  eff_PMT->SetLineColor(kMagenta+1);
  eff_PMT->SetMarkerStyle(20);
  eff_PMT->SetLineStyle(6);

  TGraphAsymmErrors* eff_H4_scan1 = new TGraphAsymmErrors(num_H4_scan1,denum_H4_scan1);
  
  eff_H4_scan1->SetMarkerColor(kBlue+1);
  eff_H4_scan1->SetLineColor(kBlue+1);
  eff_H4_scan1->SetMarkerStyle(20);
  eff_H4_scan1->SetLineStyle(3);

  TGraphAsymmErrors* eff_H4_scan2 = new TGraphAsymmErrors(num_H4_scan2,denum_H4_scan2);
  
  eff_H4_scan2->SetMarkerColor(kCyan+1);
  eff_H4_scan2->SetLineColor(kCyan+1);
  eff_H4_scan2->SetMarkerStyle(20);
  eff_H4_scan2->SetLineStyle(4);

  TFile *inputH4    = new TFile("myOutputFileH4.root"); 
  TGraphAsymmErrors *effMib3Corr = (TGraphAsymmErrors*)inputH4->Get("effMib3Corr");
  
  effMib3Corr->SetMarkerColor(kGreen+2);
  effMib3Corr->SetLineColor(kGreen+2);
  effMib3Corr->SetMarkerStyle(20);
  effMib3Corr->SetLineStyle(5);

  setStyle();

  TH2F* H2 = new TH2F("H2","",10000,1800.,3100.,100,0.,1.5);
  H2->GetXaxis()->SetTitle("HV (V)");
  H2->GetXaxis()->SetLabelSize(0.05);
  H2->GetYaxis()->SetTitle("Efficiency");
  //H2->GetYaxis()->SetLabelSize(0.04);

  TLegend* legend = new TLegend(0.19, 0.67, 0.69, 0.89);
  legend -> SetFillColor(kWhite);
  legend -> SetFillStyle(1000);
  legend -> SetLineWidth(0);
  legend -> SetLineColor(kWhite);
  legend -> SetTextFont(42);  
  legend -> SetTextSize(0.03);
  legend->AddEntry(eff_noIrr, "non-irradiated i-MCP, 5 GeV pions", "pl");
  legend->AddEntry(effMib3Corr, "non-irradiated i-MCP, 50 GeV electrons", "pl");
  legend->AddEntry(eff, "i-MCP irradiated with 10 kGy photons, 5 GeV pions", "pl");
  legend->AddEntry(eff_H4_scan1, "i-MCP irradiated with 10^{14} p/cm^{2}, 50 GeV electrons", "pl");
  legend->AddEntry(eff_PMT, "PMT-MCP irradiated with 10 kGy photons, 5 GeV pions", "pl");
  legend->AddEntry(eff_H4_scan2, "PMT-MCP irradiated with 10^{14} p/cm^{2}, 50 GeV electrons", "pl");
  
  TCanvas* c1 = new TCanvas("c1","c",1);
  FPCanvasStyle(c1);
  H2->Draw();
  eff->Draw("PC"); 
  eff_H4_scan1->Draw("PC,same"); 
  eff_H4_scan2->Draw("PC,same"); 
  eff_noIrr->Draw("PC,same"); 
  effMib3Corr->Draw("PC,same");
  eff_PMT->Draw("PC,same"); 
  legend->Draw();
  c1->SaveAs("eff_irradiated_Total.png","png");
  c1->SaveAs("eff_irradiated_Total.pdf","pdf");

  setStyle();

  TH2F* H3 = new TH2F("H3","",10000,1800.,3100.,100,0.,1.);
  H3->GetXaxis()->SetTitle("HV (V)");
  H3->GetXaxis()->SetLabelSize(0.05);
  H3->GetYaxis()->SetTitle("Efficiency");


  TLegend* legend2 = new TLegend(0.19, 0.77, 0.69, 0.89);
  legend2 -> SetFillColor(kWhite);
  legend2 -> SetFillStyle(1000);
  legend2 -> SetLineWidth(0);
  legend2 -> SetLineColor(kWhite);
  legend2 -> SetTextFont(42);  
  legend2 -> SetTextSize(0.04);
  legend2 ->AddEntry(eff_noIrr, "non-irradiated i-MCP", "pl");
  legend2 ->AddEntry(eff, "i-MCP irradiated with 10 kGy photons", "pl");
  
  TCanvas* c2 = new TCanvas("c2","c2",1);
  FPCanvasStyle(c2);
  H3->Draw();
  eff->Draw("PC"); 
  eff_noIrr->Draw("PC,same"); 
  legend2->Draw();
  TLatex latex2(0.75, 0.94,"#bf{#bf{Pions at 5 GeV}}");;
  latex2.SetTextSize(0.04);
  latex2.SetNDC(kTRUE);
  latex2.Draw();
  c2->SaveAs("eff_irradiated_T9.png","png");
  c2->SaveAs("eff_irradiated_T9.pdf","pdf");

  setStyle();

  TH2F* H4 = new TH2F("H4","",10000,1800.,3100.,100,0.,1.);
  H4->GetXaxis()->SetTitle("HV (V)");
  H4->GetXaxis()->SetLabelSize(0.05);
  H4->GetYaxis()->SetTitle("Efficiency");


  TLegend* legend3 = new TLegend(0.19, 0.77, 0.69, 0.89);
  legend3 -> SetFillColor(kWhite);
  legend3 -> SetFillStyle(1000);
  legend3 -> SetLineWidth(0);
  legend3 -> SetLineColor(kWhite);
  legend3 -> SetTextFont(42);  
  legend3 -> SetTextSize(0.04);
  legend3 -> AddEntry(effMib3Corr, "non-irradiated i-MCP", "pl");
  legend3 -> AddEntry(eff_H4_scan1, "i-MCP irradiated with 10^{14} p/cm^{2}", "pl");
  
  TCanvas* c3 = new TCanvas("c3","c3",1);
  FPCanvasStyle(c3);
  H4->Draw();
  eff_H4_scan1->Draw("PC"); 
  effMib3Corr->Draw("PC,same"); 
  legend3->Draw();
  TLatex latex3(0.65, 0.94,"#bf{#bf{Electrons at 50 GeV}}");;
  latex3.SetTextSize(0.04);
  latex3.SetNDC(kTRUE);
  latex3.Draw();
  c3->SaveAs("eff_irradiated_H4.png","png");
  c3->SaveAs("eff_irradiated_H4.pdf","pdf");
  
}

TString AddSelection(TTree* h4, std::string Var, TString Selection, std::string Cut = "1.")
{
    TH1F* h = new TH1F("h","h",8000,-40.,40.);
    TF1* g_fit = new TF1("g_fit","gaus",-40.,40.);
    h4->Draw((Var+std::string(" >> h")).c_str(),Selection);
    std::string sMean;

    char Mean [100];
    sprintf(Mean,"%f",h->GetMean());
    sMean = std::string(Mean);

    if(h->GetMean() < 0.){
       sMean.erase(sMean.begin(),sMean.begin()+1);
       Selection = Selection+" && fabs("+Var+"+"+sMean+")<"+Cut;
    }else{
       Selection = Selection+" && fabs("+Var+"-"+sMean+")<"+Cut;
    }

    delete h;
    delete g_fit;
    return Selection;
}
