//-----------------------------------------------------------------------------------
// Compare the per-track cluster yield from photonic electrons between the new large
// scale simulations and the older smaller scale simulations
//
// -> CASE A: Large Scale Sims (w/ default reweighting) Sims/Cocktail090417
// -> CASE B: Small Scale Sims (w/ default reweighting) Sims/Cocktail090117

// Sept. 5 2017
//------------------------------------------------------------------------------------

//---------------------------------------------------
// Variables
//---------------------------------------------------

//Number of pT bins
const int NBINS = 3;

//Factor for rebinning
const int REBINFACTOR = 8;

int cocktailNumberA = 90417; //Large scale
int cocktailNumberB = 90117; //Small scale

//CDPHI Distributions
TH1D *h_cdphi_B0_photons_A[NBINS];
TH1D *h_cdphi_B0_pizeros_A[NBINS];
TH1D *h_cdphi_B0_etas_A[NBINS];

TH1D *h_cdphi_B1_photons_A[NBINS];
TH1D *h_cdphi_B1_pizeros_A[NBINS];
TH1D *h_cdphi_B1_etas_A[NBINS];

TH1D *h_cdphi_B0_photons_B[NBINS];
TH1D *h_cdphi_B0_pizeros_B[NBINS];
TH1D *h_cdphi_B0_etas_B[NBINS];

TH1D *h_cdphi_B1_photons_B[NBINS];
TH1D *h_cdphi_B1_pizeros_B[NBINS];
TH1D *h_cdphi_B1_etas_B[NBINS];

//Number of electron tracks
TH1D *h_ntracks_photons_A[NBINS];
TH1D *h_ntracks_photons_B[NBINS];

TH1D *h_ntracks_pizeros_A[NBINS];
TH1D *h_ntracks_pizeros_B[NBINS];

TH1D *h_ntracks_etas_A[NBINS];
TH1D *h_ntracks_etas_B[NBINS];

//---------------------------------------------------
// Functions
//---------------------------------------------------

void loadHistosA()
{
	TFile *finPhotons = new TFile(Form("Sims/Cocktail0%i/twophotons_cdphiana.root", cocktailNumberA));
	TFile *finEtas    = new TFile(Form("Sims/Cocktail0%i/twoetas_cdphiana.root", cocktailNumberA));
	TFile *finPizeros = new TFile(Form("Sims/Cocktail0%i/twopizeros_cdphiana.root", cocktailNumberA));

	for (int ibin = 0; ibin < NBINS; ibin++)
	{
		h_cdphi_B0_photons_A[ibin] = (TH1D*) finPhotons->Get(Form("h_cdphi_B0_%i", ibin));
		h_cdphi_B1_photons_A[ibin] = (TH1D*) finPhotons->Get(Form("h_cdphi_B1_%i", ibin));
		h_ntracks_photons_A[ibin]  = (TH1D*) finPhotons->Get(Form("h_ntracks_electrons_%i", ibin));

		h_cdphi_B0_pizeros_A[ibin] = (TH1D*) finPizeros->Get(Form("h_cdphi_B0_%i", ibin));
		h_cdphi_B1_pizeros_A[ibin] = (TH1D*) finPizeros->Get(Form("h_cdphi_B1_%i", ibin));
		h_ntracks_pizeros_A[ibin]  = (TH1D*) finPizeros->Get(Form("h_ntracks_electrons_%i", ibin));

		h_cdphi_B0_etas_A[ibin]    = (TH1D*) finEtas->Get(Form("h_cdphi_B0_%i", ibin));
		h_cdphi_B1_etas_A[ibin]    = (TH1D*) finEtas->Get(Form("h_cdphi_B1_%i", ibin));
		h_ntracks_etas_A[ibin]     = (TH1D*) finEtas->Get(Form("h_ntracks_electrons_%i", ibin));


		h_cdphi_B0_photons_A[ibin]->SetTitle("PHOTONS B0 A");
		h_cdphi_B1_photons_A[ibin]->SetTitle("PHOTONS B1 A");

		h_cdphi_B0_pizeros_A[ibin]->SetTitle("PIZEROS B0 A");
		h_cdphi_B1_pizeros_A[ibin]->SetTitle("PIZEROS B1 A");

		h_cdphi_B0_etas_A[ibin]->SetTitle("ETAS B0 A");
		h_cdphi_B1_etas_A[ibin]->SetTitle("ETAS B1 A");
	}
}


void loadHistosB()
{
	TFile *finPhotons = new TFile(Form("Sims/Cocktail0%i/twophotons_cdphiana.root", cocktailNumberB));
	TFile *finEtas    = new TFile(Form("Sims/Cocktail0%i/twoetas_cdphiana.root", cocktailNumberB));
	TFile *finPizeros = new TFile(Form("Sims/Cocktail0%i/twopizeros_cdphiana.root", cocktailNumberB));

	for (int ibin = 0; ibin < NBINS; ibin++)
	{
		h_cdphi_B0_photons_B[ibin] = (TH1D*) finPhotons->Get(Form("h_cdphi_B0_%i", ibin));
		h_cdphi_B1_photons_B[ibin] = (TH1D*) finPhotons->Get(Form("h_cdphi_B1_%i", ibin));
		h_ntracks_photons_B[ibin]  = (TH1D*) finPhotons->Get(Form("h_ntracks_electrons_%i", ibin));

		h_cdphi_B0_pizeros_B[ibin] = (TH1D*) finPizeros->Get(Form("h_cdphi_B0_%i", ibin));
		h_cdphi_B1_pizeros_B[ibin] = (TH1D*) finPizeros->Get(Form("h_cdphi_B1_%i", ibin));
		h_ntracks_pizeros_B[ibin]  = (TH1D*) finPizeros->Get(Form("h_ntracks_electrons_%i", ibin));

		h_cdphi_B0_etas_B[ibin]    = (TH1D*) finEtas->Get(Form("h_cdphi_B0_%i", ibin));
		h_cdphi_B1_etas_B[ibin]    = (TH1D*) finEtas->Get(Form("h_cdphi_B1_%i", ibin));
		h_ntracks_etas_B[ibin]     = (TH1D*) finEtas->Get(Form("h_ntracks_electrons_%i", ibin));

		h_cdphi_B0_photons_B[ibin]->SetTitle("PHOTONS B0 B");
		h_cdphi_B1_photons_B[ibin]->SetTitle("PHOTONS B1 B");

		h_cdphi_B0_pizeros_B[ibin]->SetTitle("PIZEROS B0 B");
		h_cdphi_B1_pizeros_B[ibin]->SetTitle("PIZEROS B1 B");

		h_cdphi_B0_etas_B[ibin]->SetTitle("ETAS B0 B");
		h_cdphi_B1_etas_B[ibin]->SetTitle("ETAS B1 B");
	}
}


void rebinHistograms()
{
	for (int ibin = 0; ibin < NBINS; ibin++)
	{
		h_cdphi_B0_photons_A[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B0_pizeros_A[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B0_etas_A[ibin]->Rebin(REBINFACTOR);

		h_cdphi_B1_photons_A[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B1_pizeros_A[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B1_etas_A[ibin]->Rebin(REBINFACTOR);

		h_cdphi_B0_photons_B[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B0_pizeros_B[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B0_etas_B[ibin]->Rebin(REBINFACTOR);

		h_cdphi_B1_photons_B[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B1_pizeros_B[ibin]->Rebin(REBINFACTOR);
		h_cdphi_B1_etas_B[ibin]->Rebin(REBINFACTOR);
	}
}


void normalizeHistograms()
{
	//Normalize all cluster distribution per track
	for (int ibin = 0; ibin < NBINS; ibin++)
	{
		float nTracksPhotonsA = h_ntracks_photons_A[ibin]->GetBinContent(1);
		float nTracksPhotonsB = h_ntracks_photons_B[ibin]->GetBinContent(1);

		float nTracksPizerosA = h_ntracks_pizeros_A[ibin]->GetBinContent(1);
		float nTracksPizerosB = h_ntracks_pizeros_B[ibin]->GetBinContent(1);

		float nTracksEtasA    = h_ntracks_etas_A[ibin]->GetBinContent(1);
		float nTracksEtasB    = h_ntracks_etas_B[ibin]->GetBinContent(1);

		h_cdphi_B0_photons_A[ibin]->Scale(1.0 / nTracksPhotonsA);
		h_cdphi_B1_photons_A[ibin]->Scale(1.0 / nTracksPhotonsA);

		h_cdphi_B0_pizeros_A[ibin]->Scale(1.0 / nTracksPizerosA);
		h_cdphi_B1_pizeros_A[ibin]->Scale(1.0 / nTracksPizerosA);

		h_cdphi_B0_etas_A[ibin]->Scale(1.0 / nTracksEtasA);
		h_cdphi_B1_etas_A[ibin]->Scale(1.0 / nTracksEtasA);

		h_cdphi_B0_photons_B[ibin]->Scale(1.0 / nTracksPhotonsB);
		h_cdphi_B1_photons_B[ibin]->Scale(1.0 / nTracksPhotonsB);

		h_cdphi_B0_pizeros_B[ibin]->Scale(1.0 / nTracksPizerosB);
		h_cdphi_B1_pizeros_B[ibin]->Scale(1.0 / nTracksPizerosB);

		h_cdphi_B0_etas_B[ibin]->Scale(1.0 / nTracksEtasB);
		h_cdphi_B1_etas_B[ibin]->Scale(1.0 / nTracksEtasB);
	}
}


void plot()
{
	gStyle->SetOptStat(0);

	TCanvas *cPizeros[NBINS];

	for (int ibin = 0; ibin < NBINS; ibin++)
	{
		cPizeros[ibin] = new TCanvas(Form("cPizeros_%i", ibin), Form("Pizeros %i", ibin), 800, 800);
		cPizeros[ibin]->Divide(2, 2);

		cPizeros[ibin]->cd(1);
		gPad->SetLogy();
		h_cdphi_B0_pizeros_A[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B0_pizeros_A[ibin]->Draw("hist");

		cPizeros[ibin]->cd(2);
		gPad->SetLogy();
		h_cdphi_B0_pizeros_B[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B0_pizeros_B[ibin]->Draw("hist");

		cPizeros[ibin]->cd(3);
		gPad->SetLogy();
		h_cdphi_B1_pizeros_A[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B1_pizeros_A[ibin]->Draw("hist");

		cPizeros[ibin]->cd(4);
		gPad->SetLogy();
		h_cdphi_B1_pizeros_B[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B1_pizeros_B[ibin]->Draw("hist");
	}

	////////////////////

	TCanvas *cEtas[NBINS];

	for (int ibin = 0; ibin < NBINS; ibin++)
	{
		cEtas[ibin] = new TCanvas(Form("cEtas_%i", ibin), Form("Etas %i", ibin), 800, 800);
		cEtas[ibin]->Divide(2, 2);

		cEtas[ibin]->cd(1);
		gPad->SetLogy();
		h_cdphi_B0_etas_A[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B0_etas_A[ibin]->Draw("hist");
		h_cdphi_B0_etas_A[ibin]->SetLineColor(kRed);

		cEtas[ibin]->cd(2);
		gPad->SetLogy();
		h_cdphi_B0_etas_B[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B0_etas_B[ibin]->Draw("hist");
		h_cdphi_B0_etas_B[ibin]->SetLineColor(kRed);

		cEtas[ibin]->cd(3);
		gPad->SetLogy();
		h_cdphi_B1_etas_A[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B1_etas_A[ibin]->Draw("hist");
		h_cdphi_B1_etas_A[ibin]->SetLineColor(kRed);

		cEtas[ibin]->cd(4);
		gPad->SetLogy();
		h_cdphi_B1_etas_B[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B1_etas_B[ibin]->Draw("hist");
		h_cdphi_B1_etas_B[ibin]->SetLineColor(kRed);
	}

	/////////////////////

	TCanvas *cPhotons[NBINS];

	for (int ibin = 0; ibin < NBINS; ibin++)
	{
		cPhotons[ibin] = new TCanvas(Form("cPhotons_%i", ibin), Form("Photons %i", ibin), 800, 800);
		cPhotons[ibin]->Divide(2, 2);

		cPhotons[ibin]->cd(1);
		gPad->SetLogy();
		h_cdphi_B0_photons_A[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B0_photons_A[ibin]->Draw("hist");
		h_cdphi_B0_photons_A[ibin]->SetLineColor(kGreen + 3);

		cPhotons[ibin]->cd(2);
		gPad->SetLogy();
		h_cdphi_B0_photons_B[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B0_photons_B[ibin]->Draw("hist");
		h_cdphi_B0_photons_B[ibin]->SetLineColor(kGreen + 3);

		cPhotons[ibin]->cd(3);
		gPad->SetLogy();
		h_cdphi_B1_photons_A[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B1_photons_A[ibin]->Draw("hist");
		h_cdphi_B1_photons_A[ibin]->SetLineColor(kGreen + 3);

		cPhotons[ibin]->cd(4);
		gPad->SetLogy();
		h_cdphi_B1_photons_B[ibin]->GetYaxis()->SetRangeUser(1E-9, 1.0);
		h_cdphi_B1_photons_B[ibin]->Draw("hist");
		h_cdphi_B1_photons_B[ibin]->SetLineColor(kGreen + 3);
	}
}


void comparePerTrackYield()
{
	loadHistosA();
	loadHistosB();
	rebinHistograms();
	normalizeHistograms();
	plot();
}