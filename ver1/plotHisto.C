int plotHisto(TString filename="RADout") {
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Draw histos filled by Geant4 simulation
  //

  // Open file filled by Geant4 simulation
  TFile f( Form("%s.root",filename.Data()) );

  // Create a canvas and divide it into 2x2 pads
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  c1->Divide(2,2);

  // Draw Eabs histogram in the pad 1
  c1->cd(1);
  TH1D* hist1 = (TH1D*)f.Get("histograms/Eabs");
  hist1->Draw("HIST");

  // Draw Labs histogram in the pad 2
  c1->cd(2);
  TH1D* hist2 = (TH1D*)f.Get("histograms/Labs");
  hist2->Draw("HIST");

  // Draw Egap histogram in the pad 3
  // with logaritmic scale for y
  TH1D* hist3 = (TH1D*)f.Get("histograms/Elyso");
  c1->cd(3);
  gPad->SetLogy(1);
  hist3->Draw("HIST");

  // Draw Lgap histogram in the pad 4
  // with logaritmic scale for y
  c1->cd(4);
  gPad->SetLogy(1);
  TH1D* hist4 = (TH1D*)f.Get("histograms/Llyso");
  hist4->Draw("HIST");

  ofstream fout( Form("%s.txt",filename.Data()) );
  fout << "# lysoIdx, Edep (MeV)" << endl;
  TProfile* prof0 = (TProfile*)f.Get("histograms/LysoAbs");
  for(int i=0; i!=prof0->GetNbinsX(); ++i) {
    fout << i << ", " << prof0->GetBinContent(i+1) << endl;
  }
  fout.close();

  return 0;
}
