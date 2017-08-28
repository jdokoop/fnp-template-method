#include <iostream>

using namespace std;

//-----------------------
// Variables
//-----------------------

const int NBINS = 13;

//pT bins
float pT_low[NBINS]  = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
float pT_high[NBINS] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

//Save plots as PDFs
bool savePlots = false;

//Rebin factor
const int REBINF = 8;

//Reject points near peak for fit?
bool reject = !true;

//Reject sparse bins?
bool rejectSparse = !true;

//Smooth out sideband with template fit?
bool smooth = true;

//Simulated particle cocktail to use
int cocktailNumber = 82717;

//Limits for fitting data CDPHI
const float FIT_LOW = -0.15;
const float FIT_HIGH = 0.15;
const float EXCLUDE_LOW = -0.05;
const float EXCLUDE_HIGH = 0.015;

//pT bin to compute FNP
int pTBin = 2;

float pTLow = -9999;
float pTHigh = -9999;
string pTString = "";

//FNP calculated in each layer
float fnp_B0;
float fnp_B1;

float fnp_B0_err;
float fnp_B1_err;

//Functions to fit photonic electron cdphi for smoothing
TF1 *f_pos_cdphi_B0;
TF1 *f_neg_cdphi_B0;

TF1 *f_pos_cdphi_B1;
TF1 *f_neg_cdphi_B1;

//Number of inclusive electron tracks = measured - swapped electrons
long int numInclusiveElectronsB0;
long int numInclusiveElectronsB1;

//Number of photonic electrons from pizeros, etas, and direct photons
long int numPhotonicElectrons;

//Number of hadron tracks
long int numHadronsB0;
long int numHadronsB1;

//Number of clusters associated with a charged hadron track
float numClusterPerHadronTrackB0;
float numClusterPerHadronTrackB1;

//CDPHI from each source of photonic electrons
TH1D *h_cdphi_B0_pizeros;
TH1D *h_cdphi_B1_pizeros;

TH1D *h_cdphi_B0_etas;
TH1D *h_cdphi_B1_etas;

TH1D *h_cdphi_B0_photons;
TH1D *h_cdphi_B1_photons;

//CDPHI for total photonic cocktail
TH1D *h_cdphi_cocktail_B0;
TH1D *h_cdphi_cocktail_B1;
TH1D *h_cdphi_cocktail_nosmoothing_B0;
TH1D *h_cdphi_cocktail_nosmoothing_B1;

//CDPHI for photonic cocktail with multiplicity background added
TH1D *h_cdphi_cocktail_multback_B0;
TH1D *h_cdphi_cocktail_multback_B1;
TH1D *h_cdphi_cocktail_multback_beforenormalization_B0;
TH1D *h_cdphi_cocktail_multback_beforenormalization_B1;

//Photonic electron spectra
TH1D *h_elec_pT_pizeros;
TH1D *h_elec_pT_etas;
TH1D *h_elec_pT_photons;
TH1D *h_elec_pT_photonic;
TH1D *h_elec_pT_pizeros_fraction;
TH1D *h_elec_pT_etas_fraction;
TH1D *h_elec_pT_photons_fraction;

//Rebinned photonic electron spectra
TH1D *h_elec_pT_pizeros_rebinned;
TH1D *h_elec_pT_etas_rebinned;
TH1D *h_elec_pT_photons_rebinned;
TH1D *h_elec_pT_photonic_rebinned;
TH1D *h_elec_pT_pizeros_fraction_rebinned;
TH1D *h_elec_pT_etas_fraction_rebinned;
TH1D *h_elec_pT_photons_fraction_rebinned;

//CDPHI from data electrons
TH1D *h_cdphi_data_electrons_B0;
TH1D *h_cdphi_data_electrons_B1;

//CDPHI from swapped tracks
TH1D *h_cdphi_data_electrons_swapped_B0;
TH1D *h_cdphi_data_electrons_swapped_B1;

//CDPHI from charged hadrons in data
TH1D *h_cdphi_data_hadrons_B0;
TH1D *h_cdphi_data_hadrons_B1;

//CDPHI from data with swapped background subtracted
TH1D *h_cdphi_data_electrons_inclusive_B0;
TH1D *h_cdphi_data_electrons_inclusive_B1;

//Distribution of clusters per hadron track
TH1D *h_cdphi_cluspertrack_hadrons_B0;
TH1D *h_cdphi_cluspertrack_hadrons_B1;

//CDPHI distributions scaled by FNP
TH1D *h_cdphi_data_hadrons_B0_scaled;
TH1D *h_cdphi_data_hadrons_B1_scaled;

TH1D *h_cdphi_cocktail_multback_B0_scaled;
TH1D *h_cdphi_cocktail_multback_B1_scaled;

//CDPHI distributions normalized to same area in the tails
TH1D *h_cdphi_data_hadrons_B0_tailnorm;
TH1D *h_cdphi_data_hadrons_B1_tailnorm;

TH1D *h_cdphi_cocktail_multback_B0_tailnorm;
TH1D *h_cdphi_cocktail_multback_B1_tailnorm;

//Histograms to count the number of tracks in each pT bin
TH1D *h_ntracks_pizeros;
TH1D *h_ntracks_photons;
TH1D *h_ntracks_etas;

TH1D *h_ntracks_electrons;
TH1D *h_ntracks_swapped;
TH1D *h_ntracks_hadrons;

//Fit to the underlying distribution of clusters per charged hadron track
TF1 *f_multiplicitybackground_B0;
TF1 *f_multiplicitybackground_B1;

//Histograms to add errors from photonic component to data
TH1D *h_cdphi_cocktail_multback_B0_errors;
TH1D *h_cdphi_cocktail_multback_B1_errors;

//Number of data electrons (no swapped subtraction) using Tim's cuts in each pT bin
TH1D* h_elec_pT_data;
TH1D* h_elec_pT_data_rebinned;

//Fit to determine FNP
TF1 *fFitB0;
TF1 *fFitB1;

//Histogram representation of fit
TH1D *h_fit_B0;
TH1D *h_fit_B1;

//Stack for fit components + fit
THStack *hs_B0_fit;
THStack *hs_B1_fit;

//Ratio of fit to data
TH1D *h_cdphi_data_ratio_B0;
TH1D *h_cdphi_data_ratio_B1;

//FNP TGraphs
TGraphErrors *g_fnp_B0;
TGraphErrors *g_fnp_B1;
TGraphErrors *g_fnp_ratio;

//-----------------------
// Functions
//-----------------------

/*
 * Load data and simulation histograms from file
 */
void loadHistos()
{
	pTLow  = pT_low[pTBin];
	pTHigh = pT_high[pTBin];

	//Load simulation histograms
	cout << cocktailNumber << endl;
	TFile *finPhotons  = new TFile(Form("Sims/Cocktail0%i/twophotons_cdphiana.root", cocktailNumber));
	h_cdphi_B0_photons = (TH1D*) finPhotons->Get(Form("h_cdphi_B0_%i", pTBin));
	h_cdphi_B1_photons = (TH1D*) finPhotons->Get(Form("h_cdphi_B1_%i", pTBin));
	h_ntracks_photons  = (TH1D*) finPhotons->Get(Form("h_ntracks_electrons_%i", pTBin));
	h_elec_pT_photons  = (TH1D*) finPhotons->Get("h_elec_pT");

	TFile *finEtas     = new TFile(Form("Sims/Cocktail0%i/twoetas_cdphiana.root", cocktailNumber));
	h_cdphi_B0_etas    = (TH1D*) finEtas->Get(Form("h_cdphi_B0_%i", pTBin));
	h_cdphi_B1_etas    = (TH1D*) finEtas->Get(Form("h_cdphi_B1_%i", pTBin));
	h_ntracks_etas     = (TH1D*) finEtas->Get(Form("h_ntracks_electrons_%i", pTBin));
	h_elec_pT_etas     = (TH1D*) finEtas->Get("h_elec_pT");

	TFile *finPizeros   = new TFile(Form("Sims/Cocktail0%i/twopizeros_cdphiana.root", cocktailNumber));
	h_cdphi_B0_pizeros  = (TH1D*) finPizeros->Get(Form("h_cdphi_B0_%i", pTBin));
	h_cdphi_B1_pizeros  = (TH1D*) finPizeros->Get(Form("h_cdphi_B1_%i", pTBin));
	h_ntracks_pizeros   = (TH1D*) finPizeros->Get(Form("h_ntracks_electrons_%i", pTBin));
	h_elec_pT_pizeros   = (TH1D*) finPizeros->Get("h_elec_pT");

	//If pTBin == 5 or 6, add the next bin
	if (pTBin == 5 || pTBin == 6)
	{
		TH1D *h_aux_B0_photons = (TH1D*) finPhotons->Get(Form("h_cdphi_B0_%i", pTBin + 1));
		TH1D *h_aux_B0_pizeros = (TH1D*) finPizeros->Get(Form("h_cdphi_B0_%i", pTBin + 1));
		TH1D *h_aux_B0_etas    = (TH1D*) finEtas->Get(Form("h_cdphi_B0_%i", pTBin + 1));

		TH1D *h_aux_B1_photons = (TH1D*) finPhotons->Get(Form("h_cdphi_B1_%i", pTBin + 1));
		TH1D *h_aux_B1_pizeros = (TH1D*) finPizeros->Get(Form("h_cdphi_B1_%i", pTBin + 1));
		TH1D *h_aux_B1_etas    = (TH1D*) finEtas->Get(Form("h_cdphi_B1_%i", pTBin + 1));

		h_cdphi_B0_photons->Add(h_aux_B0_photons);
		h_cdphi_B0_pizeros->Add(h_aux_B0_pizeros);
		h_cdphi_B0_etas->Add(h_aux_B0_etas);

		h_cdphi_B1_photons->Add(h_aux_B1_photons);
		h_cdphi_B1_pizeros->Add(h_aux_B1_pizeros);
		h_cdphi_B1_etas->Add(h_aux_B1_etas);
	}

	//Add up combined photonic electron spectra
	h_elec_pT_photonic = (TH1D*) h_elec_pT_photons->Clone("h_elec_pT_photonic");
	h_elec_pT_photonic->Add(h_elec_pT_pizeros);
	h_elec_pT_photonic->Add(h_elec_pT_etas);

	//Load data histograms
	TFile *finData = new TFile("Data/11728.root");
	h_cdphi_data_electrons_B0                  = (TH1D*) finData->Get(Form("h_cdphi_B0_%i", pTBin));
	h_cdphi_data_electrons_B1                  = (TH1D*) finData->Get(Form("h_cdphi_B1_%i", pTBin));

	h_cdphi_data_electrons_swapped_B0          = (TH1D*) finData->Get(Form("h_cdphi_B0_sw_%i", pTBin));
	h_cdphi_data_electrons_swapped_B1          = (TH1D*) finData->Get(Form("h_cdphi_B1_sw_%i", pTBin));

	h_cdphi_data_hadrons_B0                    = (TH1D*) finData->Get(Form("h_cdphi_B0_hd_%i", pTBin));
	h_cdphi_data_hadrons_B1                    = (TH1D*) finData->Get(Form("h_cdphi_B1_hd_%i", pTBin));

	h_ntracks_electrons                        = (TH1D*) finData->Get(Form("h_ntracks_electrons_%i", pTBin));
	h_ntracks_hadrons                          = (TH1D*) finData->Get(Form("h_ntracks_hadrons_%i", pTBin));
	h_ntracks_swapped                          = (TH1D*) finData->Get(Form("h_ntracks_swapped_%i", pTBin));

	h_elec_pT_data                             = (TH1D*) finData->Get("h_elec_pT");

	if (pTBin == 5 || pTBin == 6)
	{
		TH1D *h_aux_B0                   = (TH1D*) finData->Get(Form("h_cdphi_B0_%i", pTBin + 1));
		TH1D *h_aux_swapped_B0           = (TH1D*) finData->Get(Form("h_cdphi_B0_sw_%i", pTBin + 1));
		TH1D *h_aux_hadrons_B0           = (TH1D*) finData->Get(Form("h_cdphi_B0_hd_%i", pTBin + 1));

		TH1D *h_aux_B1                   = (TH1D*) finData->Get(Form("h_cdphi_B1_%i", pTBin + 1));
		TH1D *h_aux_swapped_B1           = (TH1D*) finData->Get(Form("h_cdphi_B1_sw_%i", pTBin + 1));
		TH1D *h_aux_hadrons_B1           = (TH1D*) finData->Get(Form("h_cdphi_B1_hd_%i", pTBin + 1));

		h_cdphi_data_electrons_B0->Add(h_aux_B0);
		h_cdphi_data_electrons_swapped_B0->Add(h_aux_swapped_B0);
		h_cdphi_data_hadrons_B0->Add(h_aux_hadrons_B0);

		h_cdphi_data_electrons_B1->Add(h_aux_B1);
		h_cdphi_data_electrons_swapped_B1->Add(h_aux_swapped_B1);
		h_cdphi_data_hadrons_B1->Add(h_aux_hadrons_B1);
	}
}


void setErrors()
{
	h_cdphi_data_electrons_B0->Sumw2();
	h_cdphi_data_electrons_B1->Sumw2();

	//h_cdphi_data_electrons_swapped_B0->Sumw2();
	//h_cdphi_data_electrons_swapped_B1->Sumw2();

	//h_cdphi_data_hadrons_B0->Sumw2();
	//h_cdphi_data_hadrons_B1->Sumw2();
}

/*
 * Fit function for B0
 */
double fitFunctionB0(double *cdphi, double *par)
{
	//This is the user-defined function to fit the data with
	// ---> fit = [1-FNP] cdphi_{photonic sim} + [FNP] cdphi_{background}

	if (cdphi[0] < EXCLUDE_HIGH && cdphi[0] > EXCLUDE_LOW && reject)
	{
		TF1::RejectPoint();
		return 0.0;
	}

	double a = par[0];

	int bin = h_cdphi_cocktail_multback_B0->GetXaxis()->FindBin(cdphi[0]);

	double val = 0;

	val = (1 - a) * h_cdphi_cocktail_multback_B0->GetBinContent(bin) + a * h_cdphi_data_hadrons_B0->GetBinContent(bin);

	return val;
}


/*
 * Fit function for B1
 */
double fitFunctionB1(double *cdphi, double *par)
{
	//This is the user-defined function to fit the data with
	// ---> fit = [0] cdphi_{photonic sim} + [1] cdphi_{background}

	if (cdphi[0] < EXCLUDE_HIGH && cdphi[0] > EXCLUDE_LOW && reject)
	{
		TF1::RejectPoint();
		return 0.0;
	}

	double a = par[0];

	int bin = h_cdphi_cocktail_multback_B1->GetXaxis()->FindBin(cdphi[0]);

	double val = 0;

	val = (1 - a) * h_cdphi_cocktail_multback_B1->GetBinContent(bin) + a * h_cdphi_data_hadrons_B1->GetBinContent(bin);

	return val;
}


/*
 * Subtract swapped background from electrons in data
 */
void removeSwappedBackground()
{
	h_cdphi_data_electrons_inclusive_B0 = (TH1D*) h_cdphi_data_electrons_B0->Clone("h_cdphi_data_electrons_inclusive_B0");
	h_cdphi_data_electrons_inclusive_B1 = (TH1D*) h_cdphi_data_electrons_B1->Clone("h_cdphi_data_electrons_inclusive_B1");

	h_cdphi_data_electrons_inclusive_B0->Add(h_cdphi_data_electrons_swapped_B0, -1.0);
	h_cdphi_data_electrons_inclusive_B1->Add(h_cdphi_data_electrons_swapped_B1, -1.0);
}


/*
 * Determine the number of photonic electrons, inclusive electrons, and charged hadrons in the pT bin at hand
 */
void getNumTracks()
{
	//Number of inclusive electrons = number of measured electron tracks - number of swapped electron tracks
	numInclusiveElectronsB0 = h_ntracks_electrons->Integral() - h_ntracks_swapped->Integral();
	numInclusiveElectronsB1 = numInclusiveElectronsB0;

	//Number of photonic electrons
	long int nPhotonicEPhotons = h_ntracks_photons->Integral();
	long int nPhotonicEPizeros = h_ntracks_pizeros->Integral();
	long int nPhotonicEEtas    = h_ntracks_etas->Integral();

	numPhotonicElectrons = nPhotonicEPhotons + nPhotonicEPizeros + nPhotonicEEtas;

	//Number of hadrons
	numHadronsB0 = h_ntracks_hadrons->Integral();
	numHadronsB1 = numHadronsB0;

	//Print number of tracks
	cout << "--> Number of inclusive electrons in data B0 = " << numInclusiveElectronsB0 << endl;
	cout << "--> Number of inclusive electrons in data B1 = " << numInclusiveElectronsB1 << endl << endl;

	cout << "--> Number of photonic electrons in B0       = " << numPhotonicElectrons << endl << endl;

	cout << "--> Number of hadrons in B0                  = " << numHadronsB0 << endl;
	cout << "--> Number of hadrons in B1                  = " << numHadronsB1 << endl << endl;
}


/*
 * Rebin base histograms
 */
void rebinHistos()
{
	//Rebin CDPHI with variable binning, using wider bins in the tails, and finer bins in the peak
	/*
	const int nbins = 55;
	double binsVar[nbins + 1] = { -0.2002,
	                              -0.1902,
	                              -0.1802,
	                              -0.1702,
	                              -0.1602,
	                              -0.1502,
	                              -0.1402,
	                              -0.1302,
	                              -0.1202,
	                              -0.1102,
	                              -0.1002,
	                              -0.0902,
	                              -0.0802,
	                              -0.0702,
	                              -0.0602,
	                              -0.0502,
	                              -0.0462,
	                              -0.0422,
	                              -0.0382,
	                              -0.0342,
	                              -0.0302,
	                              -0.0262,
	                              -0.0222,
	                              -0.0182,
	                              -0.0142,
	                              -0.0102,
	                              -0.0062,
	                              -0.0022,
	                              0.0018,
	                              0.0058,
	                              0.0098,
	                              0.0138,
	                              0.0178,
	                              0.0218,
	                              0.0258,
	                              0.0298,
	                              0.0338,
	                              0.0378,
	                              0.0418,
	                              0.0458,
	                              0.0498,
	                              0.0598,
	                              0.0698,
	                              0.0798,
	                              0.0898,
	                              0.0998,
	                              0.1098,
	                              0.1198,
	                              0.1298,
	                              0.1398,
	                              0.1498,
	                              0.1598,
	                              0.1698,
	                              0.1798,
	                              0.1898,
	                              0.1998
	                            };

	h_cdphi_B0_pizeros = (TH1D*) h_cdphi_B0_pizeros->Rebin(nbins, "h_cdphi_B0_pizeros_rebinned", binsVar);
	h_cdphi_B1_pizeros = (TH1D*) h_cdphi_B1_pizeros->Rebin(nbins, "h_cdphi_B1_pizeros_rebinned", binsVar);

	h_cdphi_B0_etas = (TH1D*) h_cdphi_B0_etas->Rebin(nbins, "h_cdphi_B0_etas_rebinned", binsVar);
	h_cdphi_B1_etas = (TH1D*) h_cdphi_B1_etas->Rebin(nbins, "h_cdphi_B1_etas_rebinned", binsVar);

	h_cdphi_B0_photons = (TH1D*) h_cdphi_B0_photons->Rebin(nbins, "h_cdphi_B0_photons_rebinned", binsVar);
	h_cdphi_B1_photons = (TH1D*) h_cdphi_B1_photons->Rebin(nbins, "h_cdphi_B1_photons_rebinned", binsVar);

	h_cdphi_counter_B0_pizeros = (TH1D*) h_cdphi_counter_B0_pizeros->Rebin(nbins, "h_cdphi_counter_B0_pizeros_rebinned", binsVar);
	h_cdphi_counter_B1_pizeros = (TH1D*) h_cdphi_counter_B1_pizeros->Rebin(nbins, "h_cdphi_counter_B1_pizeros_rebinned", binsVar);

	h_cdphi_counter_B0_etas = (TH1D*) h_cdphi_counter_B0_etas->Rebin(nbins, "h_cdphi_counter_B0_etas_rebinned", binsVar);
	h_cdphi_counter_B1_etas = (TH1D*) h_cdphi_counter_B1_etas->Rebin(nbins, "h_cdphi_counter_B1_etas_rebinned", binsVar);

	h_cdphi_counter_B0_photons = (TH1D*) h_cdphi_counter_B0_photons->Rebin(nbins, "h_cdphi_counter_B0_photons_rebinned", binsVar);
	h_cdphi_counter_B1_photons = (TH1D*) h_cdphi_counter_B1_photons->Rebin(nbins, "h_cdphi_counter_B1_photons_rebinned", binsVar);

	h_cdphi_data_electrons_B0 = (TH1D*) h_cdphi_data_electrons_B0->Rebin(nbins, "h_cdphi_data_electrons_B0_rebinned", binsVar);
	h_cdphi_data_electrons_B1 = (TH1D*) h_cdphi_data_electrons_B1->Rebin(nbins, "h_cdphi_data_electrons_B1_rebinned", binsVar);

	h_cdphi_data_electrons_cr_B0 = (TH1D*) h_cdphi_data_electrons_selfcorr_B0->Rebin(nbins, "h_cdphi_data_electrons_selfcorr_B0_rebinned", binsVar);
	h_cdphi_data_electrons_selfcorr_B1 = (TH1D*) h_cdphi_data_electrons_selfcorr_B1->Rebin(nbins, "h_cdphi_data_electrons_selfcorr_B1_rebinned", binsVar);

	h_cdphi_data_electrons_swapped_B0 = (TH1D*) h_cdphi_data_electrons_swapped_B0->Rebin(nbins, "h_cdphi_data_electrons_swapped_B0_rebinned", binsVar);
	h_cdphi_data_electrons_swapped_B1 = (TH1D*) h_cdphi_data_electrons_swapped_B1->Rebin(nbins, "h_cdphi_data_electrons_swapped_B1_rebinned", binsVar);

	h_cdphi_data_selfcorr_electrons_swapped_B0 = (TH1D*) h_cdphi_data_selfcorr_electrons_swapped_B0->Rebin(nbins, "h_cdphi_data_selfcorr_electrons_swapped_B0_rebinned", binsVar);
	h_cdphi_data_selfcorr_electrons_swapped_B1 = (TH1D*) h_cdphi_data_selfcorr_electrons_swapped_B1->Rebin(nbins, "h_cdphi_data_selfcorr_electrons_swapped_B1_rebinned", binsVar);

	h_cdphi_data_hadrons_B0 = (TH1D*) h_cdphi_data_hadrons_B0->Rebin(nbins, "h_cdphi_data_hadrons_B0_rebinned", binsVar);
	h_cdphi_data_hadrons_B1 = (TH1D*) h_cdphi_data_hadrons_B1->Rebin(nbins, "h_cdphi_data_hadrons_B1_rebinned", binsVar);

	h_cdphi_data_selfcorr_hadrons_B0 = (TH1D*) h_cdphi_data_selfcorr_hadrons_B0->Rebin(nbins, "h_cdphi_data_selfcorr_hadrons_B0_rebinned", binsVar);
	h_cdphi_data_selfcorr_hadrons_B1 = (TH1D*) h_cdphi_data_selfcorr_hadrons_B1->Rebin(nbins, "h_cdphi_data_selfcorr_hadrons_B1_rebinned", binsVar);

	h_cdphi_data_electrons_inclusive_B0 = (TH1D*) h_cdphi_data_electrons_inclusive_B0->Rebin(nbins, "h_cdphi_data_electrons_inclusive_B0_rebinned", binsVar);
	h_cdphi_data_electrons_inclusive_B1 = (TH1D*) h_cdphi_data_electrons_inclusive_B1->Rebin(nbins, "h_cdphi_data_electrons_inclusive_B1_rebinned", binsVar);

	//Scale by the bin width
	for (int i = 1; i <= h_cdphi_B0_pizeros->GetNbinsX(); i++)
	{
	  if (i > 15 && i <= 40) continue;

	  h_cdphi_B0_pizeros->SetBinContent(i, h_cdphi_B0_pizeros->GetBinContent(i) / 2.5);
	  h_cdphi_B1_pizeros->SetBinContent(i, h_cdphi_B1_pizeros->GetBinContent(i) / 2.5);

	  h_cdphi_B0_etas->SetBinContent(i, h_cdphi_B0_etas->GetBinContent(i) / 2.5);
	  h_cdphi_B1_etas->SetBinContent(i, h_cdphi_B1_etas->GetBinContent(i) / 2.5);

	  h_cdphi_B0_photons->SetBinContent(i, h_cdphi_B0_photons->GetBinContent(i) / 2.5);
	  h_cdphi_B1_photons->SetBinContent(i, h_cdphi_B1_photons->GetBinContent(i) / 2.5);

	  h_cdphi_counter_B0_pizeros->SetBinContent(i, h_cdphi_counter_B0_pizeros->GetBinContent(i) / 2.5);
	  h_cdphi_counter_B1_pizeros->SetBinContent(i, h_cdphi_counter_B1_pizeros->GetBinContent(i) / 2.5);

	  h_cdphi_counter_B0_etas->SetBinContent(i, h_cdphi_counter_B0_etas->GetBinContent(i) / 2.5);
	  h_cdphi_counter_B1_etas->SetBinContent(i, h_cdphi_counter_B1_etas->GetBinContent(i) / 2.5);

	  h_cdphi_counter_B0_photons->SetBinContent(i, h_cdphi_counter_B0_photons->GetBinContent(i) / 2.5);
	  h_cdphi_counter_B1_photons->SetBinContent(i, h_cdphi_counter_B1_photons->GetBinContent(i) / 2.5);

	  h_cdphi_data_electrons_B0->SetBinContent(i, h_cdphi_data_electrons_B0->GetBinContent(i) / 2.5);
	  h_cdphi_data_electrons_B1->SetBinContent(i, h_cdphi_data_electrons_B1->GetBinContent(i) / 2.5);

	  h_cdphi_data_electrons_selfcorr_B0->SetBinContent(i, h_cdphi_data_electrons_selfcorr_B0->GetBinContent(i) / 2.5);
	  h_cdphi_data_electrons_selfcorr_B1->SetBinContent(i, h_cdphi_data_electrons_selfcorr_B1->GetBinContent(i) / 2.5);

	  h_cdphi_data_electrons_swapped_B0->SetBinContent(i, h_cdphi_data_electrons_swapped_B0->GetBinContent(i) / 2.5);
	  h_cdphi_data_electrons_swapped_B1->SetBinContent(i, h_cdphi_data_electrons_swapped_B1->GetBinContent(i) / 2.5);

	  h_cdphi_data_selfcorr_electrons_swapped_B0->SetBinContent(i, h_cdphi_data_selfcorr_electrons_swapped_B0->GetBinContent(i) / 2.5);
	  h_cdphi_data_selfcorr_electrons_swapped_B1->SetBinContent(i, h_cdphi_data_selfcorr_electrons_swapped_B1->GetBinContent(i) / 2.5);

	  h_cdphi_data_hadrons_B0->SetBinContent(i, h_cdphi_data_hadrons_B0->GetBinContent(i) / 2.5);
	  h_cdphi_data_hadrons_B1->SetBinContent(i, h_cdphi_data_hadrons_B1->GetBinContent(i) / 2.5);

	  h_cdphi_data_selfcorr_hadrons_B0->SetBinContent(i, h_cdphi_data_selfcorr_hadrons_B0->GetBinContent(i) / 2.5);
	  h_cdphi_data_selfcorr_hadrons_B1->SetBinContent(i, h_cdphi_data_selfcorr_hadrons_B1->GetBinContent(i) / 2.5);

	  h_cdphi_data_electrons_inclusive_B0->SetBinContent(i, h_cdphi_data_electrons_inclusive_B0->GetBinContent(i) / 2.5);
	  h_cdphi_data_electrons_inclusive_B1->SetBinContent(i, h_cdphi_data_electrons_inclusive_B1->GetBinContent(i) / 2.5);
	}
	*/

	//Rebin CDPHI distributions by constant factor
	h_cdphi_B0_pizeros->Rebin(REBINF);

	h_cdphi_B1_pizeros->Rebin(REBINF);

	h_cdphi_B0_etas->Rebin(REBINF);
	h_cdphi_B1_etas->Rebin(REBINF);

	h_cdphi_B0_photons->Rebin(REBINF);
	h_cdphi_B1_photons->Rebin(REBINF);

	h_cdphi_data_electrons_B0->Rebin(REBINF);
	h_cdphi_data_electrons_B1->Rebin(REBINF);

	h_cdphi_data_electrons_swapped_B0->Rebin(REBINF);
	h_cdphi_data_electrons_swapped_B1->Rebin(REBINF);

	h_cdphi_data_hadrons_B0->Rebin(REBINF);
	h_cdphi_data_hadrons_B1->Rebin(REBINF);

	h_cdphi_data_electrons_inclusive_B0->Rebin(REBINF);
	h_cdphi_data_electrons_inclusive_B1->Rebin(REBINF);

	//Rebin photonic electron spectra with variable bin width
	double bins[13] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0};
	h_elec_pT_pizeros_rebinned = (TH1D*) h_elec_pT_pizeros->Rebin(12, "h_elec_pT_pizeros_rebinned", bins);
	h_elec_pT_etas_rebinned = (TH1D*) h_elec_pT_etas->Rebin(12, "h_elec_pT_etas_rebinned", bins);
	h_elec_pT_photons_rebinned = (TH1D*) h_elec_pT_photons->Rebin(12, "h_elec_pT_photons_rebinned", bins);

	h_elec_pT_photonic_rebinned = (TH1D*) h_elec_pT_photons_rebinned->Clone("h_elec_pT_photonic_rebinned");
	h_elec_pT_photonic_rebinned->Add(h_elec_pT_pizeros_rebinned);
	h_elec_pT_photonic_rebinned->Add(h_elec_pT_etas_rebinned);

	//Rebin data electron spectrum with variable bin width
	double binsData[12] = {1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
	h_elec_pT_data_rebinned = (TH1D*) h_elec_pT_data->Rebin(11, "h_elec_pT_data_rebinned", binsData);

}


/*
 * Set a histogram's title and format for axis titles and labels
 */
void formatHistogram(TH1D *& h, string xaxis, string yaxis, string title)
{
	h->SetTitle(title.c_str());
	h->GetXaxis()->SetTitle(xaxis.c_str());
	h->GetYaxis()->SetTitle(yaxis.c_str());
	h->GetXaxis()->SetTitleFont(62);
	h->GetXaxis()->SetLabelFont(62);
	h->GetYaxis()->SetTitleFont(62);
	h->GetYaxis()->SetLabelFont(62);
}


void calculatePhotonicRatios()
{
	h_elec_pT_photons_fraction_rebinned = (TH1D*) h_elec_pT_photons_rebinned->Clone("h_elec_pT_photons_fraction_rebinned");
	h_elec_pT_pizeros_fraction_rebinned = (TH1D*) h_elec_pT_pizeros_rebinned->Clone("h_elec_pT_pizeros_fraction_rebinned");
	h_elec_pT_etas_fraction_rebinned = (TH1D*) h_elec_pT_etas_rebinned->Clone("h_elec_pT_etas_fraction_rebinned");

	h_elec_pT_photons_fraction_rebinned->Divide(h_elec_pT_photonic_rebinned);
	h_elec_pT_pizeros_fraction_rebinned->Divide(h_elec_pT_photonic_rebinned);
	h_elec_pT_etas_fraction_rebinned->Divide(h_elec_pT_photonic_rebinned);

	formatHistogram(h_elec_pT_pizeros_rebinned, "p_{T} [GeV/c]", "Counts", "");
	formatHistogram(h_elec_pT_photons_rebinned, "p_{T} [GeV/c]", "Counts", "");
	formatHistogram(h_elec_pT_etas_rebinned, "p_{T} [GeV/c]", "Counts", "");

	h_elec_pT_pizeros_fraction_rebinned->SetLineWidth(2);
	h_elec_pT_photons_fraction_rebinned->SetLineWidth(2);
	h_elec_pT_etas_fraction_rebinned->SetLineWidth(2);

	h_elec_pT_pizeros_rebinned->SetLineWidth(2);
	h_elec_pT_photons_rebinned->SetLineWidth(2);
	h_elec_pT_etas_rebinned->SetLineWidth(2);

	h_elec_pT_pizeros_fraction_rebinned->SetLineColor(kAzure - 3);
	h_elec_pT_photons_fraction_rebinned->SetLineColor(kOrange - 2);
	h_elec_pT_etas_fraction_rebinned->SetLineColor(kPink);

	h_elec_pT_pizeros_rebinned->SetLineColor(kAzure - 3);
	h_elec_pT_photons_rebinned->SetLineColor(kOrange - 2);
	h_elec_pT_etas_rebinned->SetLineColor(kPink);
}


/*
 * Define histogram titles, labels, fonts, and colors
 */
void setFormat()
{
	//Set range for x axis
	h_cdphi_B0_pizeros->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);
	h_cdphi_B1_pizeros->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);

	h_cdphi_B0_etas->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);
	h_cdphi_B1_etas->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);

	h_cdphi_B0_photons->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);
	h_cdphi_B1_photons->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);

	h_cdphi_data_electrons_B0->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);
	h_cdphi_data_electrons_B1->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);

	h_cdphi_data_electrons_swapped_B0->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);
	h_cdphi_data_electrons_swapped_B1->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);

	//h_cdphi_data_hadrons_B0->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);
	//h_cdphi_data_hadrons_B1->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);

	h_cdphi_data_electrons_inclusive_B0->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);
	h_cdphi_data_electrons_inclusive_B1->GetXaxis()->SetRangeUser(FIT_LOW, FIT_HIGH);

	//Set titles and fonts
	formatHistogram(h_cdphi_B0_pizeros, "charge #times #Delta#phi", "Counts", "B0");
	formatHistogram(h_cdphi_B1_pizeros, "charge #times #Delta#phi", "Counts", "B1");

	formatHistogram(h_cdphi_B0_etas, "charge #times #Delta#phi", "Counts", "B0");
	formatHistogram(h_cdphi_B1_etas, "charge #times #Delta#phi", "Counts", "B1");

	formatHistogram(h_cdphi_B0_photons, "charge #times #Delta#phi", "Counts", "B0");
	formatHistogram(h_cdphi_B1_photons, "charge #times #Delta#phi", "Counts", "B1");

	formatHistogram(h_cdphi_data_electrons_B0, "charge #times #Delta#phi", "Counts", "B0");
	formatHistogram(h_cdphi_data_electrons_B1, "charge #times #Delta#phi", "Counts", "B1");

	formatHistogram(h_cdphi_data_electrons_swapped_B0, "charge #times #Delta#phi", "Counts", "B0");
	formatHistogram(h_cdphi_data_electrons_swapped_B1, "charge #times #Delta#phi", "Counts", "B1");

	formatHistogram(h_cdphi_data_hadrons_B0, "charge #times #Delta#phi", "Counts", "B0");
	formatHistogram(h_cdphi_data_hadrons_B1, "charge #times #Delta#phi", "Counts", "B1");

	formatHistogram(h_cdphi_data_electrons_inclusive_B0, "charge #times #Delta#phi", "Counts", "B0");
	formatHistogram(h_cdphi_data_electrons_inclusive_B1, "charge #times #Delta#phi", "Counts", "B1");
}


/*
 * Add up the contributions from each cocktail particle
 */
void addCocktail()
{
	h_cdphi_cocktail_B0 = (TH1D*) h_cdphi_B0_photons->Clone("h_cdphi_cocktail_B0");
	h_cdphi_cocktail_B0->Add(h_cdphi_B0_pizeros);
	h_cdphi_cocktail_B0->Add(h_cdphi_B0_etas);

	h_cdphi_cocktail_B1 = (TH1D*) h_cdphi_B1_photons->Clone("h_cdphi_cocktail_B1");
	h_cdphi_cocktail_B1->Add(h_cdphi_B1_pizeros);
	h_cdphi_cocktail_B1->Add(h_cdphi_B1_etas);
}


/*
 * Fit the negative tail and positive sideband and tail of the photonic electron cdphi distribution
 */
void fitPhotonicSideband()
{
	//Make a copy of the cdphi distribution from the simulate cocktail, but set
	//the errors to a constant value in order to fit the tails of the distribution
	TH1D *h_cdphi_cocktail_B1_copy = (TH1D*) h_cdphi_cocktail_B1->Clone("h_cdphi_cocktail_B1_copy");
	for (int i = 1; i <= h_cdphi_cocktail_B1_copy->GetNbinsX(); i++)
	{
		h_cdphi_cocktail_B1_copy->SetBinError(i, 1.0);
	}

	TH1D *h_cdphi_cocktail_B0_copy = (TH1D*) h_cdphi_cocktail_B0->Clone("h_cdphi_cocktail_B0_copy");
	for (int i = 1; i <= h_cdphi_cocktail_B0_copy->GetNbinsX(); i++)
	{
		h_cdphi_cocktail_B0_copy->SetBinError(i, 1.0);
	}

	float fitHighPos = 0.15;
	float fitLowPos = 0.02;

	float fitHighNeg = -0.01;
	float fitLowNeg = -0.15;

	//----------B0-------------

	if (pTBin == 2)
	{
		f_pos_cdphi_B0 = new TF1("f_pos_cdphi_B0", "[0]*TMath::Exp([1]*x) + [2]*TMath::Exp([3]*x)", fitLowPos, fitHighPos);

		f_pos_cdphi_B0->SetParameter(0, 8.62982e+09);
		f_pos_cdphi_B0->SetParameter(1, -4.65486e+01);
		f_pos_cdphi_B0->SetParameter(2, 2.33832e+09);
		f_pos_cdphi_B0->SetParameter(3, -2.04280e+01);
	}
	else if (pTBin == 3)
	{
		f_pos_cdphi_B0 = new TF1("f_pos_cdphi_B0", "[0]*TMath::Exp([1]*x) + [2]*TMath::Exp([3]*x)", fitLowPos, fitHighPos);

		f_pos_cdphi_B0->SetParameter(0, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_pos_cdphi_B0->SetParameter(1, -4.33578e+01);
		f_pos_cdphi_B0->SetParameter(2, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_pos_cdphi_B0->SetParameter(3, -3.81166e+00);
	}
	else if (pTBin == 4)
	{
		f_pos_cdphi_B0 = new TF1("f_pos_cdphi_B0", "[0]*TMath::Exp([1]*x) + [2]*TMath::Exp([3]*x)", fitLowPos, fitHighPos);

		f_pos_cdphi_B0->SetParameter(0, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_pos_cdphi_B0->SetParameter(1, -4.33578e+01);
		f_pos_cdphi_B0->SetParameter(2, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_pos_cdphi_B0->SetParameter(3, -3.81166e+00);
	}
	else if (pTBin == 5)
	{
		f_pos_cdphi_B0 = new TF1("f_pos_cdphi_B0", "1.2*[0]*TMath::Exp([1]*x) + 1.2*[2]*TMath::Exp([3]*x)", fitLowPos, fitHighPos);
		f_pos_cdphi_B0->FixParameter(0, 2.24691e+07);
		f_pos_cdphi_B0->FixParameter(1, -7.68560e+01);
		f_pos_cdphi_B0->FixParameter(2, 5.12410e+05);
		f_pos_cdphi_B0->FixParameter(3, -1.64956e+01);
	}
	else if (pTBin == 6)
	{
		f_pos_cdphi_B0 = new TF1("f_pos_cdphi_B0", "1.2*[0]*TMath::Exp([1]*x) + 1.2*[2]*TMath::Exp([3]*x)", fitLowPos, fitHighPos);

		f_pos_cdphi_B0->FixParameter(0, 2.54E06);
		f_pos_cdphi_B0->FixParameter(1, -44.3022);
		f_pos_cdphi_B0->FixParameter(2, 1.09225E06);
		f_pos_cdphi_B0->FixParameter(3, -55.046);
	}

	h_cdphi_cocktail_B0_copy->Fit(f_pos_cdphi_B0, "Q0R");

	//Fit negative tail
	f_neg_cdphi_B0 = new TF1("f_neg_cdphi_B0", "[0]*TMath::Exp([1]*x) + [2]*TMath::Exp([3]*x)", fitLowNeg, fitHighNeg);

	if (pTBin == 2)
	{
		f_neg_cdphi_B0->FixParameter(0, 0.0);
		f_neg_cdphi_B0->FixParameter(1, 75.53);
		f_neg_cdphi_B0->FixParameter(2, 4.22E8);
		f_neg_cdphi_B0->FixParameter(3, 48.45);
		/*
		f_neg_cdphi_B0->SetParameter(0, 2.94E09);
		f_neg_cdphi_B0->SetParameter(1, 345.8);
		f_neg_cdphi_B0->SetParameter(2, 1E09);
		f_neg_cdphi_B0->SetParameter(3, 33.4);
		*/
	}
	else if (pTBin == 3)
	{
		f_neg_cdphi_B0->SetParameter(0, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(1, 1.54722e+02);
		f_neg_cdphi_B0->SetParameter(2, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(3, 2.80063e+01);
	}
	else if (pTBin == 4)
	{
		f_neg_cdphi_B0->SetParameter(0, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(1, 1.54722e+02);
		f_neg_cdphi_B0->SetParameter(2, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(3, 2.80063e+01);
	}
	else if (pTBin == 5)
	{
		f_neg_cdphi_B0->SetParameter(0, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(1, 1.54722e+02);
		f_neg_cdphi_B0->SetParameter(2, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(3, 2.80063e+01);
	}
	else if (pTBin == 6)
	{
		f_neg_cdphi_B0->SetParameter(0, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(1, 1.54722e+02);
		f_neg_cdphi_B0->SetParameter(2, h_cdphi_cocktail_B0_copy->GetMaximum());
		f_neg_cdphi_B0->SetParameter(3, 2.80063e+01);
	}
	h_cdphi_cocktail_B0_copy->Fit(f_neg_cdphi_B0, "R");


	//----------B1-------------
	f_pos_cdphi_B1 = new TF1("f_pos_cdphi_B1", "[0]*TMath::Exp([1]*x) + [2]*TMath::Exp([3]*x)", fitLowPos, fitHighPos);

	if (pTBin == 2)
	{
		f_pos_cdphi_B1->SetParameter(0, 1.05E08);
		f_pos_cdphi_B1->SetParameter(1, 3.01968);
		f_pos_cdphi_B1->SetParameter(2, 2.0308E10);
		f_pos_cdphi_B1->SetParameter(3, -40.3869);
	}
	else if (pTBin == 3)
	{
		f_pos_cdphi_B1->FixParameter(0, 6.89994e+08);
		f_pos_cdphi_B1->FixParameter(1, -23.0516);
		f_pos_cdphi_B1->FixParameter(2, 3.23243e+09);
		f_pos_cdphi_B1->FixParameter(3, -56.2);
	}
	else if (pTBin == 4)
	{
		f_pos_cdphi_B1->SetParameter(0, 6.89994e+08);
		f_pos_cdphi_B1->SetParameter(1, -23.0516);
		f_pos_cdphi_B1->SetParameter(2, 3.23243e+09);
		f_pos_cdphi_B1->SetParameter(3, -56.2);
	}
	else if (pTBin == 5)
	{
		f_pos_cdphi_B1->SetParameter(0, 6.89994e+08);
		f_pos_cdphi_B1->SetParameter(1, -23.0516);
		f_pos_cdphi_B1->SetParameter(2, 3.23243e+09);
		f_pos_cdphi_B1->SetParameter(3, -56.2);
	}
	else if (pTBin == 6)
	{
		f_pos_cdphi_B1 = new TF1("f_pos_cdphi_B1", "1.35*[0]*TMath::Exp([1]*x) + 1.35*[2]*TMath::Exp([3]*x)", fitLowPos, fitHighPos);

		f_pos_cdphi_B1->FixParameter(0, 42641.3);
		f_pos_cdphi_B1->FixParameter(1, -2.23027);
		f_pos_cdphi_B1->FixParameter(2, 5.91443E06);
		f_pos_cdphi_B1->FixParameter(3, -48.1219);
	}

	h_cdphi_cocktail_B1_copy->Fit(f_pos_cdphi_B1, "Q0R");

	//Fit negative tail
	f_neg_cdphi_B1 = new TF1("f_neg_cdphi_B1", "[0]*TMath::Exp([1]*x) + [2]*TMath::Exp([3]*x)", fitLowNeg, fitHighNeg);

	if (pTBin == 2)
	{
		f_neg_cdphi_B1->SetParameter(0, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(1, 7.85768e+01);
		f_neg_cdphi_B1->SetParameter(2, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(3, 3.49565e+00);
	}
	else if (pTBin == 3)
	{
		f_neg_cdphi_B1->SetParameter(0, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(1, 7.85768e+01);
		f_neg_cdphi_B1->SetParameter(2, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(3, 3.49565e+00);
	}
	else if (pTBin == 4)
	{
		f_neg_cdphi_B1->SetParameter(0, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(1, 7.85768e+01);
		f_neg_cdphi_B1->SetParameter(2, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(3, 3.49565e+00);
	}
	else if (pTBin == 5)
	{
		f_neg_cdphi_B1->SetParameter(0, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(1, 7.85768e+01);
		f_neg_cdphi_B1->SetParameter(2, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(3, 3.49565e+00);
	}
	else if (pTBin == 6)
	{
		f_neg_cdphi_B1->SetParameter(0, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(1, 7.85768e+01);
		f_neg_cdphi_B1->SetParameter(2, h_cdphi_cocktail_B1_copy->GetMaximum());
		f_neg_cdphi_B1->SetParameter(3, 3.49565e+00);
	}

	h_cdphi_cocktail_B1_copy->Fit(f_neg_cdphi_B1, "Q0R");


	//Make a copy of photonic cdphi before smoothing out the tails
	h_cdphi_cocktail_nosmoothing_B0 = (TH1D*) h_cdphi_cocktail_B0->Clone("h_cdphi_cocktail_nosmoothing_B0");
	h_cdphi_cocktail_nosmoothing_B1 = (TH1D*) h_cdphi_cocktail_B1->Clone("h_cdphi_cocktail_nosmoothing_B1");

	//Now that we've smoothed out the tails by fitting them with a function, let's replace
	//that region in the simulated cocktai histogram with the template shape
	if (!smooth) return;

	for (int i = 1; i <= h_cdphi_cocktail_B1->GetNbinsX(); i++)
	{
		float cdphi = h_cdphi_cocktail_B1->GetBinCenter(i);

		if (cdphi < fitHighNeg)
		{
			h_cdphi_cocktail_B1->SetBinContent(i, f_neg_cdphi_B1->Eval(cdphi));
		}
		if (cdphi > fitLowPos)
		{
			h_cdphi_cocktail_B1->SetBinContent(i, f_pos_cdphi_B1->Eval(cdphi));
		}
	}

	for (int i = 1; i <= h_cdphi_cocktail_B0->GetNbinsX(); i++)
	{
		float cdphi = h_cdphi_cocktail_B0->GetBinCenter(i);

		if (cdphi < fitHighNeg)
		{
			h_cdphi_cocktail_B0->SetBinContent(i, f_neg_cdphi_B0->Eval(cdphi));
		}
		if (cdphi > fitLowPos)
		{
			h_cdphi_cocktail_B0->SetBinContent(i, f_pos_cdphi_B0->Eval(cdphi));
		}
	}
}


/*
 * Get the distribution of clusters per charged hadron track
 * and use that to model the effect of the underlying event
 * multiplicity on the photonic electron cdphi distribution
 */
/*
void constructMultiplicityBackground()
{
	//Get the distributions of clusters per charged hadron track
	h_cdphi_cluspertrack_hadrons_B0 = (TH1D*) h_cdphi_data_hadrons_B0->Clone("h_cdphi_cluspertrack_hadrons_B0");
	h_cdphi_cluspertrack_hadrons_B1 = (TH1D*) h_cdphi_data_hadrons_B1->Clone("h_cdphi_cluspertrack_hadrons_B1");

	h_cdphi_cluspertrack_hadrons_B0->Scale(1.0 / numHadronsB0);
	h_cdphi_cluspertrack_hadrons_B1->Scale(1.0 / numHadronsB1);

	//Fit with a triangular function excluding the peak (peak defined as |cdphi| < 0.04)
	TF1 *f_multiplicitybackground_aux_B0 = new TF1("f_multiplicitybackground_aux_B0", "(x<-0.04)*([0]*x) + (x>0.04)*([1]*x) + (x<-0.04 || x>0.04)*[2]", -0.15, 0.15);
	h_cdphi_cluspertrack_hadrons_B0->Fit(f_multiplicitybackground_aux_B0, "Q0R");

	TF1 *f_multiplicitybackground_aux_B1 = new TF1("f_multiplicitybackground_aux_B1", "(x<-0.04)*([0]*x) + (x>0.04)*([1]*x) + (x<-0.04 || x>0.04)*[2]", -0.15, 0.15);
	h_cdphi_cluspertrack_hadrons_B1->Fit(f_multiplicitybackground_aux_B1, "Q0R");

	//Extrapolate the fit functions to the peak region
	f_multiplicitybackground_B0 = new TF1("f_multiplicitybackground_B0", "(x<0)*([0]*x + [2]) + (x>=0)*([1]*x + [2])", -0.15, 0.15);
	f_multiplicitybackground_B0->SetParameter(0, f_multiplicitybackground_aux_B0->GetParameter(0));
	f_multiplicitybackground_B0->SetParameter(1, f_multiplicitybackground_aux_B0->GetParameter(1));
	f_multiplicitybackground_B0->SetParameter(2, f_multiplicitybackground_aux_B0->GetParameter(2));

	f_multiplicitybackground_B1 = new TF1("f_multiplicitybackground_B1", "(x<0)*([0]*x + [2]) + (x>=0)*([1]*x + [2])", -0.15, 0.15);
	f_multiplicitybackground_B1->SetParameter(0, f_multiplicitybackground_aux_B1->GetParameter(0));
	f_multiplicitybackground_B1->SetParameter(1, f_multiplicitybackground_aux_B1->GetParameter(1));
	f_multiplicitybackground_B1->SetParameter(2, f_multiplicitybackground_aux_B1->GetParameter(2));

	//Determine the number of clusters per hadron track within the cdphi region of interest
	numClusterPerHadronTrackB0 = 0;
	int binlow  = h_cdphi_cluspertrack_hadrons_B0->GetXaxis()->FindBin(FIT_LOW);
	int binhigh = h_cdphi_cluspertrack_hadrons_B0->GetXaxis()->FindBin(FIT_HIGH);

	for (int i = binlow; i < binhigh; i++)
	{
		numClusterPerHadronTrackB0 += f_multiplicitybackground_B0->Eval(h_cdphi_cluspertrack_hadrons_B0->GetBinCenter(i));
	}

	binlow  = h_cdphi_cluspertrack_hadrons_B1->GetXaxis()->FindBin(FIT_LOW);
	binhigh = h_cdphi_cluspertrack_hadrons_B1->GetXaxis()->FindBin(FIT_HIGH);
	numClusterPerHadronTrackB1 = 0;

	for (int i = binlow; i < binhigh; i++)
	{
		numClusterPerHadronTrackB1 += f_multiplicitybackground_B1->Eval(h_cdphi_cluspertrack_hadrons_B1->GetBinCenter(i));
	}

	cout << " Clusters per Hadron Track B0 = " << numClusterPerHadronTrackB0 << endl;
	cout << " Clusters per Hadron Track B1 = " << numClusterPerHadronTrackB1 << endl << endl;

	//Add clusters to photonic electron cocktail to account for underlying event multiplicity background
	h_cdphi_cocktail_multback_B0 = (TH1D*) h_cdphi_cocktail_B0->Clone("h_cdphi_cocktail_multback_B0");
	h_cdphi_cocktail_multback_B1 = (TH1D*) h_cdphi_cocktail_B1->Clone("h_cdphi_cocktail_multback_B1");

	long int numExtraClusB0 = 0;
	long int numExtraClusB1 = 0;

	if (pTBin < 3)
	{
		numExtraClusB0 = TMath::Floor(numClusterPerHadronTrackB0 * numPhotonicElectrons / 100);
		numExtraClusB1 = TMath::Floor(numClusterPerHadronTrackB1 * numPhotonicElectrons / 100);
	}
	else
	{
		numExtraClusB0 = numClusterPerHadronTrackB0 * numPhotonicElectrons;
		numExtraClusB1 = numClusterPerHadronTrackB1 * numPhotonicElectrons;
	}

	cout << "--> Adding " << numExtraClusB0 << " clusters to photonic electron CDPHI in B0" << endl;
	cout << "--> Adding " << numExtraClusB1 << " clusters to photonic electron CDPHI in B1" << endl << endl;

	for (int i = 0; i < numExtraClusB0; i++)
	{
		if (i % 5000000 == 0)
		{
			cout << "... Processing " << (float) i * 100 / numExtraClusB0 << "%" << endl;
		}

		float cdphi = f_multiplicitybackground_B0->GetRandom(FIT_LOW, FIT_HIGH);
		int bin = h_cdphi_cocktail_multback_B0->GetXaxis()->FindBin(cdphi);

		if (pTBin < 3)
		{
			h_cdphi_cocktail_multback_B0->SetBinContent(bin, h_cdphi_cocktail_multback_B0->GetBinContent(bin) + 100.0);
		}
		else
		{
			h_cdphi_cocktail_multback_B0->SetBinContent(bin, h_cdphi_cocktail_multback_B0->GetBinContent(bin) + 1.0);
		}
	}

	for (int i = 0; i < numExtraClusB1; i++)
	{
		if (i % 5000000 == 0)
		{
			cout << "... Processing " << (float) i * 100 / numExtraClusB1 << "%" << endl;
		}

		float cdphi = f_multiplicitybackground_B1->GetRandom(FIT_LOW, FIT_HIGH);
		int bin = h_cdphi_cocktail_multback_B1->GetXaxis()->FindBin(cdphi);

		if (pTBin < 3)
		{
			h_cdphi_cocktail_multback_B1->SetBinContent(bin, h_cdphi_cocktail_multback_B1->GetBinContent(bin) + 100.0);
		}
		else
		{
			h_cdphi_cocktail_multback_B1->SetBinContent(bin, h_cdphi_cocktail_multback_B1->GetBinContent(bin) + 1.0);
		}
	}
}
*/

void constructMultiplicityBackground()
{
	//Get the distributions of clusters per charged hadron track
	h_cdphi_cluspertrack_hadrons_B0 = (TH1D*) h_cdphi_data_hadrons_B0->Clone("h_cdphi_cluspertrack_hadrons_B0");
	h_cdphi_cluspertrack_hadrons_B1 = (TH1D*) h_cdphi_data_hadrons_B1->Clone("h_cdphi_cluspertrack_hadrons_B1");

	h_cdphi_cluspertrack_hadrons_B0->Scale(1.0 / numHadronsB0);
	h_cdphi_cluspertrack_hadrons_B1->Scale(1.0 / numHadronsB1);

	//Fit with a triangular function excluding the peak (peak defined as |cdphi| < 0.04)
	TF1 *f_multiplicitybackground_aux_B0 = new TF1("f_multiplicitybackground_aux_B0", "(x<-0.06)*(pol5) + (x>0.06)*(pol5) + (x<-0.06 || x>0.06)*[6]", -0.15, 0.15);
	h_cdphi_cluspertrack_hadrons_B0->Fit(f_multiplicitybackground_aux_B0, "Q0R");

	TF1 *f_multiplicitybackground_aux_B1 = new TF1("f_multiplicitybackground_aux_B1", "(x<-0.06)*(pol5) + (x>0.06)*(pol5) + (x<-0.06 || x>0.06)*[6]", -0.15, 0.15);
	h_cdphi_cluspertrack_hadrons_B1->Fit(f_multiplicitybackground_aux_B1, "Q0R");

	//Extrapolate the fit functions to the peak region
	f_multiplicitybackground_B0 = new TF1("f_multiplicitybackground_B0",  "(x<0)*(pol5 + [6]) + (x>=0)*(pol5 + [6])", -0.15, 0.15);
	for (int ipar = 0; ipar < 7; ipar++)
	{
		f_multiplicitybackground_B0->SetParameter(ipar, f_multiplicitybackground_aux_B0->GetParameter(ipar));

	}

	f_multiplicitybackground_B1 = new TF1("f_multiplicitybackground_B1",  "(x<0)*(pol5 + [6]) + (x>=0)*(pol5 + [6])", -0.15, 0.15);
	for (int ipar = 0; ipar < 7; ipar++)
	{
		f_multiplicitybackground_B1->SetParameter(ipar, f_multiplicitybackground_aux_B1->GetParameter(ipar));

	}

	//Determine the number of clusters per hadron track within the cdphi region of interest
	numClusterPerHadronTrackB0 = 0;
	int binlow  = h_cdphi_cluspertrack_hadrons_B0->GetXaxis()->FindBin(FIT_LOW);
	int binhigh = h_cdphi_cluspertrack_hadrons_B0->GetXaxis()->FindBin(FIT_HIGH);

	for (int i = binlow; i < binhigh; i++)
	{
		numClusterPerHadronTrackB0 += f_multiplicitybackground_B0->Eval(h_cdphi_cluspertrack_hadrons_B0->GetBinCenter(i));
	}

	binlow  = h_cdphi_cluspertrack_hadrons_B1->GetXaxis()->FindBin(FIT_LOW);
	binhigh = h_cdphi_cluspertrack_hadrons_B1->GetXaxis()->FindBin(FIT_HIGH);
	numClusterPerHadronTrackB1 = 0;

	for (int i = binlow; i < binhigh; i++)
	{
		numClusterPerHadronTrackB1 += f_multiplicitybackground_B1->Eval(h_cdphi_cluspertrack_hadrons_B1->GetBinCenter(i));
	}

	cout << " Clusters per Hadron Track B0 = " << numClusterPerHadronTrackB0 << endl;
	cout << " Clusters per Hadron Track B1 = " << numClusterPerHadronTrackB1 << endl << endl;

	//Add clusters to photonic electron cocktail to account for underlying event multiplicity background
	h_cdphi_cocktail_multback_B0 = (TH1D*) h_cdphi_cocktail_B0->Clone("h_cdphi_cocktail_multback_B0");
	h_cdphi_cocktail_multback_B1 = (TH1D*) h_cdphi_cocktail_B1->Clone("h_cdphi_cocktail_multback_B1");

	long int numExtraClusB0 = 0;
	long int numExtraClusB1 = 0;

	if (pTBin < 3)
	{
		numExtraClusB0 = TMath::Floor(numClusterPerHadronTrackB0 * numPhotonicElectrons / 100);
		numExtraClusB1 = TMath::Floor(numClusterPerHadronTrackB1 * numPhotonicElectrons / 100);
	}
	else
	{
		numExtraClusB0 = numClusterPerHadronTrackB0 * numPhotonicElectrons;
		numExtraClusB1 = numClusterPerHadronTrackB1 * numPhotonicElectrons;
	}

	cout << "--> Adding " << numExtraClusB0 << " clusters to photonic electron CDPHI in B0" << endl;
	cout << "--> Adding " << numExtraClusB1 << " clusters to photonic electron CDPHI in B1" << endl << endl;

	for (int i = 0; i < numExtraClusB0; i++)
	{
		if (i % 5000000 == 0)
		{
			cout << "... Processing " << (float) i * 100 / numExtraClusB0 << "%" << endl;
		}

		float cdphi = f_multiplicitybackground_B0->GetRandom(FIT_LOW, FIT_HIGH);//h_cdphi_cluspertrack_hadrons_B0->GetRandom();
		int bin = h_cdphi_cocktail_multback_B0->GetXaxis()->FindBin(cdphi);

		if (pTBin < 3)
		{
			h_cdphi_cocktail_multback_B0->SetBinContent(bin, h_cdphi_cocktail_multback_B0->GetBinContent(bin) + 100.0);
		}
		else
		{
			h_cdphi_cocktail_multback_B0->SetBinContent(bin, h_cdphi_cocktail_multback_B0->GetBinContent(bin) + 1.0);
		}
	}

	for (int i = 0; i < numExtraClusB1; i++)
	{
		if (i % 5000000 == 0)
		{
			cout << "... Processing " << (float) i * 100 / numExtraClusB1 << "%" << endl;
		}

		float cdphi = f_multiplicitybackground_B1->GetRandom(FIT_LOW, FIT_HIGH);//h_cdphi_cluspertrack_hadrons_B1->GetRandom();
		int bin = h_cdphi_cocktail_multback_B1->GetXaxis()->FindBin(cdphi);

		if (pTBin < 3)
		{
			h_cdphi_cocktail_multback_B1->SetBinContent(bin, h_cdphi_cocktail_multback_B1->GetBinContent(bin) + 100.0);
		}
		else
		{
			h_cdphi_cocktail_multback_B1->SetBinContent(bin, h_cdphi_cocktail_multback_B1->GetBinContent(bin) + 1.0);
		}
	}
}


void plotClustersPerHadronTrack()
{
	TCanvas *cClusPerHadronTrackB0 = new TCanvas("cClusPerHadronTrackB0", "Clusters per Hadron Track in B0", 500, 500);
	cClusPerHadronTrackB0->SetLogy();
	h_cdphi_cluspertrack_hadrons_B0->SetTitle("");
	h_cdphi_cluspertrack_hadrons_B0->GetYaxis()->SetTitle("Clusters per Track");
	h_cdphi_cluspertrack_hadrons_B0->GetYaxis()->SetTitleOffset(1.3);
	h_cdphi_cluspertrack_hadrons_B0->Draw();
	f_multiplicitybackground_B0->Draw("same");

	TLatex *tlB0 = new TLatex(0.15, 0.76, "B0");
	tlB0->SetNDC();
	tlB0->SetTextSize(0.12);
	tlB0->DrawClone("same");

	TLatex *tlNclusHistoB0 = new TLatex(0.2, 0.22, Form("Clusters per Hadron Track [histo] = %.3g", h_cdphi_cluspertrack_hadrons_B0->Integral(h_cdphi_cluspertrack_hadrons_B0->GetXaxis()->FindBin(FIT_LOW), h_cdphi_cluspertrack_hadrons_B0->GetXaxis()->FindBin(FIT_HIGH))));
	tlNclusHistoB0->SetNDC();
	tlNclusHistoB0->SetTextSize(0.03);
	tlNclusHistoB0->DrawClone("same");

	TLatex *tlNclusFitB0 = new TLatex(0.2, 0.18, Form("Clusters per Hadron Track [fit]     = %.3g", numClusterPerHadronTrackB0));
	tlNclusFitB0->SetNDC();
	tlNclusFitB0->SetTextSize(0.03);
	tlNclusFitB0->DrawClone("same");

	TCanvas *cClusPerHadronTrackB1 = new TCanvas("cClusPerHadronTrackB1", "Clusters per Hadron Track in B1", 500, 500);
	cClusPerHadronTrackB1->SetLogy();
	h_cdphi_cluspertrack_hadrons_B1->SetTitle("");
	h_cdphi_cluspertrack_hadrons_B1->GetYaxis()->SetTitle("Clusters per Track");
	h_cdphi_cluspertrack_hadrons_B1->GetYaxis()->SetTitleOffset(1.3);
	h_cdphi_cluspertrack_hadrons_B1->Draw();
	f_multiplicitybackground_B1->Draw("same");

	TLatex *tlB1 = new TLatex(0.15, 0.76, "B1");
	tlB1->SetNDC();
	tlB1->SetTextSize(0.12);
	tlB1->DrawClone("same");

	TLatex *tlNclusHistoB1 = new TLatex(0.2, 0.22, Form("Clusters per Hadron Track [histo] = %.3g", h_cdphi_cluspertrack_hadrons_B1->Integral(h_cdphi_cluspertrack_hadrons_B1->GetXaxis()->FindBin(FIT_LOW), h_cdphi_cluspertrack_hadrons_B1->GetXaxis()->FindBin(FIT_HIGH))));
	tlNclusHistoB1->SetNDC();
	tlNclusHistoB1->SetTextSize(0.03);
	tlNclusHistoB1->DrawClone("same");

	TLatex *tlNclusFitB1 = new TLatex(0.2, 0.18, Form("Clusters per Hadron Track [fit]     = %.3g", numClusterPerHadronTrackB1));
	tlNclusFitB1->SetNDC();
	tlNclusFitB1->SetTextSize(0.03);
	tlNclusFitB1->DrawClone("same");
}


void normalizeHistograms()
{
	//Calculate errors on photonic background before normalizing, then add to data
	h_cdphi_cocktail_multback_B0_errors = (TH1D*) h_cdphi_cocktail_multback_B0->Clone("h_cdphi_cocktail_multback_B0_errors");
	h_cdphi_cocktail_multback_B1_errors = (TH1D*) h_cdphi_cocktail_multback_B1->Clone("h_cdphi_cocktail_multback_B1_errors");

	//Normalize everything per track
	h_cdphi_cocktail_multback_B0->Scale(1.0 / numPhotonicElectrons);
	h_cdphi_cocktail_multback_B1->Scale(1.0 / numPhotonicElectrons);
	h_cdphi_data_electrons_inclusive_B0->Scale(1.0 / numInclusiveElectronsB0);
	h_cdphi_data_electrons_inclusive_B1->Scale(1.0 / numInclusiveElectronsB1);
	h_cdphi_data_hadrons_B0->Scale(1.0 / numHadronsB0);
	h_cdphi_data_hadrons_B1->Scale(1.0 / numHadronsB1);

	for (int i = 1; i <= h_cdphi_cocktail_multback_B0_errors->GetNbinsX(); i++)
	{
		double error = h_cdphi_cocktail_multback_B0->GetBinError(i);
		h_cdphi_cocktail_multback_B0_errors->SetBinContent(i, error);
	}

	for (int i = 1; i <= h_cdphi_cocktail_multback_B1_errors->GetNbinsX(); i++)
	{
		double error = h_cdphi_cocktail_multback_B1->GetBinError(i);
		h_cdphi_cocktail_multback_B1_errors->SetBinContent(i, error);
	}

	//for (int i = 1; i <= h_cdphi_cocktail_multback_B0->GetNbinsX(); i++)
	//{
	//  cout << h_cdphi_cocktail_multback_B0->GetBinError(i) << endl;
	//}

	//Normalize to the same area in the tails
	float cocktail_multback_B0_tailarea = h_cdphi_cocktail_multback_B0->Integral(h_cdphi_cocktail_multback_B0->GetXaxis()->FindBin(FIT_LOW), h_cdphi_cocktail_multback_B0->GetXaxis()->FindBin(-0.1));
	//cocktail_multback_B0_tailarea += h_cdphi_cocktail_multback_B0->Integral(h_cdphi_cocktail_multback_B0->GetXaxis()->FindBin(0.1), h_cdphi_cocktail_multback_B0->GetXaxis()->FindBin(FIT_HIGH));

	float cocktail_multback_B1_tailarea = h_cdphi_cocktail_multback_B1->Integral(h_cdphi_cocktail_multback_B1->GetXaxis()->FindBin(FIT_LOW), h_cdphi_cocktail_multback_B1->GetXaxis()->FindBin(-0.1));
	//cocktail_multback_B1_tailarea += h_cdphi_cocktail_multback_B1->Integral(h_cdphi_cocktail_multback_B1->GetXaxis()->FindBin(0.1), h_cdphi_cocktail_multback_B1->GetXaxis()->FindBin(FIT_HIGH));

	float data_hadrons_B0_tailarea = h_cdphi_data_hadrons_B0->Integral(h_cdphi_data_hadrons_B0->GetXaxis()->FindBin(FIT_LOW), h_cdphi_data_hadrons_B0->GetXaxis()->FindBin(-0.1));
	//data_hadrons_B0_tailarea += h_cdphi_data_hadrons_B0->Integral(h_cdphi_data_hadrons_B0->GetXaxis()->FindBin(0.1), h_cdphi_data_hadrons_B0->GetXaxis()->FindBin(FIT_HIGH));

	float data_hadrons_B1_tailarea = h_cdphi_data_hadrons_B1->Integral(h_cdphi_data_hadrons_B1->GetXaxis()->FindBin(FIT_LOW), h_cdphi_data_hadrons_B1->GetXaxis()->FindBin(-0.1));
	//data_hadrons_B1_tailarea += h_cdphi_data_hadrons_B1->Integral(h_cdphi_data_hadrons_B1->GetXaxis()->FindBin(0.1), h_cdphi_data_hadrons_B1->GetXaxis()->FindBin(FIT_HIGH));

	h_cdphi_cocktail_multback_B0_tailnorm = (TH1D*) h_cdphi_cocktail_multback_B0->Clone("h_cdphi_cocktail_multback_B0_tailnorm");
	h_cdphi_cocktail_multback_B1_tailnorm = (TH1D*) h_cdphi_cocktail_multback_B1->Clone("h_cdphi_cocktail_multback_B1_tailnorm");

	h_cdphi_data_hadrons_B0_tailnorm = (TH1D*) h_cdphi_data_hadrons_B0->Clone("h_cdphi_data_hadrons_B0_tailnorm");
	h_cdphi_data_hadrons_B1_tailnorm = (TH1D*) h_cdphi_data_hadrons_B1->Clone("h_cdphi_data_hadrons_B1_tailnorm");

	h_cdphi_cocktail_multback_B0_tailnorm->Scale(1.0 / cocktail_multback_B0_tailarea);
	h_cdphi_cocktail_multback_B1_tailnorm->Scale(1.0 / cocktail_multback_B1_tailarea);

	h_cdphi_data_hadrons_B0_tailnorm->Scale(1.0 / data_hadrons_B0_tailarea);
	h_cdphi_data_hadrons_B1_tailnorm->Scale(1.0 / data_hadrons_B1_tailarea);
}


void addPhotonicErrorToData()
{
	for (int i = 1; i <= h_cdphi_cocktail_multback_B0_errors->GetNbinsX(); i++)
	{
		double photonicError = h_cdphi_cocktail_multback_B0_errors->GetBinContent(i);
		double dataError = h_cdphi_data_electrons_inclusive_B0->GetBinError(i);
		double totalError = TMath::Sqrt(photonicError * photonicError + dataError * dataError);

		h_cdphi_data_electrons_inclusive_B0->SetBinError(i, totalError);
	}

	for (int i = 1; i <= h_cdphi_cocktail_multback_B1_errors->GetNbinsX(); i++)
	{
		double photonicError = h_cdphi_cocktail_multback_B1_errors->GetBinContent(i);
		double dataError = h_cdphi_data_electrons_inclusive_B1->GetBinError(i);
		double totalError = TMath::Sqrt(photonicError * photonicError + dataError * dataError);

		h_cdphi_data_electrons_inclusive_B1->SetBinError(i, totalError);
	}
}


void plotFitComponents()
{
	TCanvas *cFitComponentsB0 = new TCanvas("cFitComponentsB0", "cFitComponentsB0", 600, 600);
	cFitComponentsB0->SetLogy();

	h_cdphi_data_electrons_inclusive_B0->SetMarkerStyle(20);
	h_cdphi_data_electrons_inclusive_B0->SetMarkerSize(0.5);
	h_cdphi_data_electrons_inclusive_B0->SetLineColor(kBlack);
	h_cdphi_cocktail_multback_B0->SetLineColor(kBlue);
	h_cdphi_data_hadrons_B0->SetLineColor(kOrange + 7);

	h_cdphi_data_electrons_inclusive_B0->GetYaxis()->SetTitle("Clusters per Track");
	h_cdphi_data_electrons_inclusive_B0->Draw();
	h_cdphi_cocktail_multback_B0->Draw("hist,same");
	h_cdphi_data_hadrons_B0->Draw("hist,same");

	TLegend *legB0 = new TLegend(FIT_HIGH, 0.6, 0.65, 0.8);
	legB0->SetLineColor(kWhite);
	legB0->AddEntry(h_cdphi_data_electrons_inclusive_B0, "Inclusive Electrons - Swapped Electrons", "P");
	legB0->AddEntry(h_cdphi_cocktail_multback_B0, "Photonic Component with HM Background", "L");
	legB0->AddEntry(h_cdphi_data_hadrons_B0, "Non-Photonic Background", "L");
	legB0->Draw("same");

	TLatex *tlB0 = new TLatex(0.68, 0.73, "B0");
	tlB0->SetNDC();
	tlB0->SetTextSize(0.15);
	tlB0->Draw("same");

	TCanvas *cFitComponentsB1 = new TCanvas("cFitComponentsB1", "cFitComponentsB1", 600, 600);
	cFitComponentsB1->SetLogy();

	h_cdphi_data_electrons_inclusive_B1->SetMarkerStyle(20);
	h_cdphi_data_electrons_inclusive_B1->SetMarkerSize(0.5);
	h_cdphi_data_electrons_inclusive_B1->SetLineColor(kBlack);
	h_cdphi_cocktail_multback_B1->SetLineColor(kBlue);
	h_cdphi_data_hadrons_B1->SetLineColor(kOrange + 7);

	h_cdphi_data_electrons_inclusive_B1->GetYaxis()->SetTitle("Clusters per Track");
	h_cdphi_data_electrons_inclusive_B1->Draw();
	h_cdphi_cocktail_multback_B1->Draw("hist,same");
	h_cdphi_data_hadrons_B1->Draw("hist,same");

	TLegend *legB1 = new TLegend(FIT_HIGH, 0.6, 0.65, 0.8);
	legB1->SetLineColor(kWhite);
	legB1->AddEntry(h_cdphi_data_electrons_inclusive_B1, "Inclusive Electrons - Swapped Electrons", "P");
	legB1->AddEntry(h_cdphi_cocktail_multback_B1, "Photonic Component with HM Background", "L");
	legB1->AddEntry(h_cdphi_data_hadrons_B1, "Non-Photonic Background", "L");
	legB1->Draw("same");

	TLatex *tlB1 = new TLatex(0.68, 0.73, "B1");
	tlB1->SetNDC();
	tlB1->SetTextSize(0.15);
	tlB1->Draw("same");

	//Verify integrals
	int binLow = h_cdphi_data_electrons_inclusive_B0->GetXaxis()->FindBin(-0.14);
	int binHigh = h_cdphi_data_electrons_inclusive_B0->GetXaxis()->FindBin(0.14);
	cout << "DATA INTEGRAL B0 = " << h_cdphi_data_electrons_inclusive_B0->Integral(binLow, binHigh) << endl;
	cout << "PHOT INTEGRAL B0 = " << h_cdphi_cocktail_multback_B0->Integral(binLow, binHigh) << endl;
	cout << "NPHO INTEGRAL B0 = " << h_cdphi_data_hadrons_B0->Integral(binLow, binHigh) << endl;
}


void plotPhotonicComponentWithoutBackground()
{
	TCanvas *cPhotonicNoBackgB0 = new TCanvas("cPhotonicNoBackgB0", "cPhotonicNoBackgB0", 600, 600);
	cPhotonicNoBackgB0->SetLogy();
	h_cdphi_cocktail_nosmoothing_B0->GetYaxis()->SetRangeUser(1E4, 1E10);
	h_cdphi_cocktail_nosmoothing_B0->Draw();
	h_cdphi_cocktail_nosmoothing_B0->SetMarkerStyle(20);
	h_cdphi_cocktail_nosmoothing_B0->SetMarkerSize(0.5);
	h_cdphi_cocktail_nosmoothing_B0->SetMarkerColor(kBlue);
	h_cdphi_cocktail_nosmoothing_B0->SetLineColor(kBlue);
	f_pos_cdphi_B0->Draw("same");
	f_neg_cdphi_B0->Draw("same");

	TCanvas *cPhotonicNoBackgB1 = new TCanvas("cPhotonicNoBackgB1", "cPhotonicNoBackgB1", 600, 600);
	cPhotonicNoBackgB1->SetLogy();
	h_cdphi_cocktail_nosmoothing_B1->GetYaxis()->SetRangeUser(1E4, 1E10);
	h_cdphi_cocktail_nosmoothing_B1->Draw();
	h_cdphi_cocktail_nosmoothing_B1->SetMarkerStyle(20);
	h_cdphi_cocktail_nosmoothing_B1->SetMarkerSize(0.5);
	h_cdphi_cocktail_nosmoothing_B1->SetMarkerColor(kBlue);
	h_cdphi_cocktail_nosmoothing_B1->SetLineColor(kBlue);
	f_pos_cdphi_B1->Draw("same");
	f_neg_cdphi_B1->Draw("same");
}


/*
 * Plot the cdphi distribution with reweighting for every cocktail component
 */
void plotPhotonicReweighted()
{
	TCanvas *cPhotonicReweightedB1 = new TCanvas("cPhotonicReweightedB1", "Weighted CDPHI B1", 600, 800);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.66, 1, 1);
	pad1->SetLogy();
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	h_cdphi_B1_pizeros->SetTitle("");
	h_cdphi_B1_pizeros->SetLineWidth(2);
	h_cdphi_B1_pizeros->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_B1_pizeros->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_B1_pizeros->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_B1_pizeros->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_B1_pizeros->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_B1_pizeros->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_B1_pizeros->GetXaxis()->SetRangeUser(-0.15, 0.15);
	h_cdphi_B1_pizeros->Draw("hist");

	TLatex *tPizeroB1 = new TLatex(FIT_HIGH, 0.4, "Pizero B1");
	tPizeroB1->SetNDC();
	tPizeroB1->SetTextSize(0.07);
	tPizeroB1->Draw("same");

	cPhotonicReweightedB1->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.33, 1, 0.66);
	pad2->SetLogy();
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();
	h_cdphi_B1_etas->SetTitle("");
	h_cdphi_B1_etas->SetLineWidth(2);
	h_cdphi_B1_etas->GetYaxis()->CenterTitle();
	h_cdphi_B1_etas->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_B1_etas->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_B1_etas->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_B1_etas->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_B1_etas->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_B1_etas->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_B1_etas->GetXaxis()->SetRangeUser(-0.15, 0.15);
	h_cdphi_B1_etas->Draw("hist");

	TLatex *tEtaB1 = new TLatex(FIT_HIGH, 0.4, "Eta B1");
	tEtaB1->SetNDC();
	tEtaB1->SetTextSize(0.07);
	tEtaB1->Draw("same");

	cPhotonicReweightedB1->cd();
	TPad *pad3 = new TPad("pad3", "pad3", 0, 0, 1, 0.33);
	pad3->SetLogy();
	pad3->SetTopMargin(0);
	pad3->SetBottomMargin(0.15);
	pad3->Draw();
	pad3->cd();
	pad3->SetTickx();
	pad3->SetTicky();
	h_cdphi_B1_photons->SetTitle("");
	h_cdphi_B1_photons->SetLineWidth(2);
	h_cdphi_B1_photons->GetYaxis()->CenterTitle();
	h_cdphi_B1_photons->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_B1_photons->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_B1_photons->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_B1_photons->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_B1_photons->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_B1_photons->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_B1_photons->GetXaxis()->SetRangeUser(-0.15, 0.15);
	h_cdphi_B1_photons->Draw("hist");

	TLatex *tPhotonB1 = new TLatex(FIT_HIGH, 0.4, "Photon B1");
	tPhotonB1->SetNDC();
	tPhotonB1->SetTextSize(0.07);
	tPhotonB1->Draw("same");

	cPhotonicReweightedB1->cd();

	//////////////////////////////

	TCanvas *cPhotonicReweightedB0 = new TCanvas("cPhotonicReweightedB0", "Weighted CDPHI B0", 600, 800);
	TPad *pad1B0 = new TPad("pad1B0", "pad1B0", 0, 0.66, 1, 1);
	pad1B0->SetLogy();
	pad1B0->SetTickx();
	pad1B0->SetTicky();
	pad1B0->SetBottomMargin(0);
	pad1B0->Draw();
	pad1B0->cd();
	h_cdphi_B0_pizeros->SetTitle("");
	h_cdphi_B0_pizeros->SetLineWidth(2);
	h_cdphi_B0_pizeros->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_B0_pizeros->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_B0_pizeros->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_B0_pizeros->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_B0_pizeros->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_B0_pizeros->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_B0_pizeros->GetXaxis()->SetRangeUser(-0.15, 0.15);
	h_cdphi_B0_pizeros->Draw("hist");

	TLatex *tPizeroB0 = new TLatex(FIT_HIGH, 0.4, "Pizero B0");
	tPizeroB0->SetNDC();
	tPizeroB0->SetTextSize(0.07);
	tPizeroB0->Draw("same");

	cPhotonicReweightedB0->cd();
	TPad *pad2B0 = new TPad("pad2B0", "pad2B0", 0, 0.33, 1, 0.66);
	pad2B0->SetLogy();
	pad2B0->SetTopMargin(0);
	pad2B0->SetBottomMargin(0);
	pad2B0->Draw();
	pad2B0->cd();
	pad2B0->SetTickx();
	pad2B0->SetTicky();
	h_cdphi_B0_etas->SetTitle("");
	h_cdphi_B0_etas->SetLineWidth(2);
	h_cdphi_B0_etas->GetYaxis()->CenterTitle();
	h_cdphi_B0_etas->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_B0_etas->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_B0_etas->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_B0_etas->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_B0_etas->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_B0_etas->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_B0_etas->GetXaxis()->SetRangeUser(-0.15, 0.15);
	h_cdphi_B0_etas->Draw("hist");

	TLatex *tEtaB0 = new TLatex(FIT_HIGH, 0.4, "Eta B0");
	tEtaB0->SetNDC();
	tEtaB0->SetTextSize(0.07);
	tEtaB0->Draw("same");

	cPhotonicReweightedB0->cd();
	TPad *pad3B0 = new TPad("pad3B0", "pad3B0", 0, 0, 1, 0.33);
	pad3B0->SetLogy();
	pad3B0->SetTopMargin(0);
	pad3B0->SetBottomMargin(0.15);
	pad3B0->Draw();
	pad3B0->cd();
	pad3B0->SetTickx();
	pad3B0->SetTicky();
	h_cdphi_B0_photons->SetTitle("");
	h_cdphi_B0_photons->SetLineWidth(2);
	h_cdphi_B0_photons->GetYaxis()->CenterTitle();
	h_cdphi_B0_photons->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_B0_photons->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_B0_photons->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_B0_photons->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_B0_photons->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_B0_photons->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_B0_photons->GetXaxis()->SetRangeUser(-0.15, 0.15);
	h_cdphi_B0_photons->Draw("hist");

	TLatex *tPhotonB0 = new TLatex(FIT_HIGH, 0.4, "Photon B0");
	tPhotonB0->SetNDC();
	tPhotonB0->SetTextSize(0.07);
	tPhotonB0->Draw("same");

	cPhotonicReweightedB0->cd();
}


/*
 * Determine FMP by integrating the CDPHI distributions as follows:
 * Clusters per track data = (1-FNP) Clusters per track photonic + FNP Clusters per track non-photonic
 */
void integrateFNP()
{
	double nd_b0 = 0.0;
	float nd_b1 = 0.0;
	float np_b0 = 0.0;
	float np_b1 = 0.0;
	float nnp_b0 = 0.0;
	float nnp_b1 = 0.0;

	for (int i = 1; i <= h_cdphi_data_electrons_inclusive_B0->GetNbinsX(); i++)
	{
		float binCenter = h_cdphi_data_electrons_inclusive_B0->GetBinCenter(i);
		if ((binCenter < FIT_LOW || binCenter > FIT_HIGH || (binCenter > EXCLUDE_LOW && binCenter < EXCLUDE_HIGH)) && reject) continue;
		if ((binCenter < FIT_LOW || binCenter > FIT_HIGH) && !reject) continue;

		nd_b0 += h_cdphi_data_electrons_inclusive_B0->GetBinContent(i);
		nd_b1 += h_cdphi_data_electrons_inclusive_B1->GetBinContent(i);

		np_b0 += h_cdphi_cocktail_multback_B0->GetBinContent(i);
		np_b1 += h_cdphi_cocktail_multback_B1->GetBinContent(i);

		nnp_b0 += h_cdphi_data_hadrons_B0->GetBinContent(i);
		nnp_b1 += h_cdphi_data_hadrons_B1->GetBinContent(i);
	}

	/*
		//Get the bins for integration
		int bin1 = h_cdphi_data_electrons_inclusive_B0->GetXaxis()->FindBin(-0.15);
		int bin2 = h_cdphi_data_electrons_inclusive_B1->GetXaxis()->FindBin(0.15);

		//First, integrate the measured CDPHI distribution to get the total number of clusters per track in data
		//Make sure that the histogram is properly normalized by the total number of tracks
		double nd_b0 = h_cdphi_data_electrons_inclusive_B0->Integral(bin1, bin2);
		double nd_b1 = h_cdphi_data_electrons_inclusive_B1->Integral(bin1, bin2);

		//Now, integrate the CDPHI distribution for the photonic component with the underlying event background added
		double np_b0 = h_cdphi_cocktail_multback_B0->Integral(bin1, bin2);
		double np_b1 = h_cdphi_cocktail_multback_B1->Integral(bin1, bin2);

		//Finally, integrate CDPHI distribution for charged hadrons in data (non-photonic contribution)
		double nnp_b0 = h_cdphi_data_hadrons_B0->Integral(bin1, bin2);
		double nnp_b1 = h_cdphi_data_hadrons_B1->Integral(bin1, bin2);
		*/

	double fnp_b0 = (nd_b0 - np_b0) / (nnp_b0 - np_b0);
	double fnp_b1 = (nd_b1 - np_b1) / (nnp_b1 - np_b1);

	cout << "---------------------------------------" << endl;
	cout << "FNP Through Integration" << endl << endl;
	cout << "FNP B0 = " << fnp_b0 << endl;
	cout << "FNP B1 = " << fnp_b1 << endl;
	cout << "---------------------------------------" << endl << endl;
}



/*
 * Fit the cdphi distribution of inclusive electrons with a combination of photonic and non-photonic components
 * using FNP as the only free parameter
 */
void fitFNP()
{
	fFitB0 = new TF1("fFitB0", fitFunctionB0, FIT_LOW, FIT_HIGH, 1);
	fFitB1 = new TF1("fFitB1", fitFunctionB1, FIT_LOW, FIT_HIGH, 1);

	fFitB0->SetParLimits(0, 0.0, 1.0);
	fFitB1->SetParLimits(0, 0.0, 1.0);

	h_cdphi_data_electrons_inclusive_B0->Fit(fFitB0, "Q0R");
	h_cdphi_data_electrons_inclusive_B1->Fit(fFitB1, "Q0R");

	fnp_B0 = fFitB0->GetParameter(0);
	fnp_B1 = fFitB1->GetParameter(0);

	fnp_B0_err = fFitB0->GetParError(0);
	fnp_B1_err = fFitB1->GetParError(0);

	cout << "--> FNP_B0 = " << fnp_B0 << " ± " << fnp_B0_err << endl;
	cout << "--> FNP_B1 = " << fnp_B1 << " ± " << fnp_B1_err << endl;

	//Take ratio of data to fit
	h_cdphi_data_ratio_B0 = (TH1D*) h_cdphi_data_electrons_inclusive_B0->Clone("h_cdphi_data_ratio_B0");
	h_cdphi_data_ratio_B1 = (TH1D*) h_cdphi_data_electrons_inclusive_B1->Clone("h_cdphi_data_ratio_B1");
	h_cdphi_data_ratio_B0->Reset();
	h_cdphi_data_ratio_B1->Reset();

	for (int i = 1; i < h_cdphi_data_ratio_B0->GetNbinsX(); i++)
	{
		float cdphi   = h_cdphi_data_ratio_B0->GetBinCenter(i);

		float content;
		if (cdphi >= EXCLUDE_LOW && cdphi <= EXCLUDE_HIGH && reject)
		{
			content = -9999;
			h_cdphi_data_ratio_B0->SetBinContent(i, content);
			h_cdphi_data_ratio_B0->SetBinError(i, 0.0);
		}
		else
		{
			content = h_cdphi_data_electrons_inclusive_B0->GetBinContent(i) / fFitB0->Eval(cdphi);
			h_cdphi_data_ratio_B0->SetBinContent(i, content);
			h_cdphi_data_ratio_B0->SetBinError(i, h_cdphi_data_electrons_inclusive_B0->GetBinError(i) / fFitB0->Eval(cdphi));
		}
	}

	for (int i = 1; i < h_cdphi_data_ratio_B1->GetNbinsX(); i++)
	{
		float cdphi   = h_cdphi_data_ratio_B1->GetBinCenter(i);

		float content;
		if (cdphi >= EXCLUDE_LOW && cdphi <= EXCLUDE_HIGH && reject)
		{
			content = -9999;
			h_cdphi_data_ratio_B1->SetBinContent(i, content);
			h_cdphi_data_ratio_B1->SetBinError(i, 0.0);
		}
		else
		{
			content = h_cdphi_data_electrons_inclusive_B1->GetBinContent(i) / fFitB1->Eval(cdphi);
			h_cdphi_data_ratio_B1->SetBinContent(i, content);
			h_cdphi_data_ratio_B1->SetBinError(i, h_cdphi_data_electrons_inclusive_B1->GetBinError(i) / fFitB1->Eval(cdphi));
		}
	}
}


/*
 *
 */
void constructFitVisualization()
{
	//Scale the fit components by FNP
	h_cdphi_data_hadrons_B0_scaled = (TH1D*) h_cdphi_data_hadrons_B0->Clone("h_cdphi_data_hadrons_B0_scaled");
	h_cdphi_data_hadrons_B1_scaled = (TH1D*) h_cdphi_data_hadrons_B1->Clone("h_cdphi_data_hadrons_B1_scaled");

	h_cdphi_cocktail_multback_B0_scaled = (TH1D*) h_cdphi_cocktail_multback_B0->Clone("h_cdphi_cocktail_multback_B0_scaled");
	h_cdphi_cocktail_multback_B1_scaled = (TH1D*) h_cdphi_cocktail_multback_B1->Clone("h_cdphi_cocktail_multback_B1_scaled");

	h_cdphi_data_hadrons_B0_scaled->Scale(fnp_B0);
	h_cdphi_data_hadrons_B1_scaled->Scale(fnp_B1);

	h_cdphi_cocktail_multback_B0_scaled->Scale(1.0 - fnp_B0);
	h_cdphi_cocktail_multback_B1_scaled->Scale(1.0 - fnp_B1);

	//Define aesthetics
	h_cdphi_cocktail_multback_B0_scaled->SetFillColorAlpha(kAzure, 0.2);
	//h_cdphi_cocktail_multback_B0_scaled->SetFillStyle(3002);
	h_cdphi_cocktail_multback_B0_scaled->SetLineColor(kAzure);

	h_cdphi_data_hadrons_B0_scaled->SetFillColorAlpha(kOrange + 1, 0.2);
	//h_cdphi_data_hadrons_B0_scaled->SetFillStyle(3002);
	h_cdphi_data_hadrons_B0_scaled->SetLineColor(kOrange + 1);

	h_cdphi_cocktail_multback_B1_scaled->SetFillColorAlpha(kAzure, 0.2);
	//h_cdphi_cocktail_multback_B1_scaled->SetFillStyle(3002);
	h_cdphi_cocktail_multback_B1_scaled->SetLineColor(kAzure);

	h_cdphi_data_hadrons_B1_scaled->SetFillColorAlpha(kOrange + 1, 0.2);
	//h_cdphi_data_hadrons_B1_scaled->SetFillStyle(3002);
	h_cdphi_data_hadrons_B1_scaled->SetLineColor(kOrange + 1);

	h_cdphi_cocktail_multback_B0_scaled->SetTitle("");
	h_cdphi_cocktail_multback_B1_scaled->SetTitle("");

	h_cdphi_data_hadrons_B0_scaled->SetTitle("");
	h_cdphi_data_hadrons_B1_scaled->SetTitle("");

	//Create histogram version of fits
	h_fit_B0 = (TH1D*) h_cdphi_data_electrons_inclusive_B0->Clone("h_fit_B0");
	h_fit_B1 = (TH1D*) h_cdphi_data_electrons_inclusive_B1->Clone("h_fit_B1");
	h_fit_B0->Reset();
	h_fit_B1->Reset();

	for (int i = 1; i <= h_fit_B0->GetNbinsX(); i++)
	{
		float cdphi = h_fit_B0->GetBinCenter(i);
		float val = fFitB0->Eval(cdphi);
		h_fit_B0->SetBinContent(i, val);
	}

	for (int i = 1; i <= h_fit_B1->GetNbinsX(); i++)
	{
		float cdphi = h_fit_B1->GetBinCenter(i);
		float val = fFitB1->Eval(cdphi);
		h_fit_B1->SetBinContent(i, val);
	}

	//Aesthetics for histogram version of fit
	h_fit_B0->SetLineColor(kRed);
	h_fit_B0->SetMarkerColor(kRed);
	h_fit_B0->SetLineWidth(2);

	h_fit_B1->SetLineColor(kRed);
	h_fit_B1->SetMarkerColor(kRed);
	h_fit_B1->SetLineWidth(2);

	//Define stacks to show fit components and total fit
	hs_B0_fit = new THStack("hs_B0_fit", "hs_B0_fit");
	hs_B0_fit->Add(h_cdphi_data_hadrons_B0_scaled);
	hs_B0_fit->Add(h_cdphi_cocktail_multback_B0_scaled);

	hs_B1_fit = new THStack("hs_B1_fit", "hs_B1_fit");
	hs_B1_fit->Add(h_cdphi_data_hadrons_B1_scaled);
	hs_B1_fit->Add(h_cdphi_cocktail_multback_B1_scaled);
}


void plotFitStacked()
{
	//------------- B0

	TCanvas *cFitStackedB0 = new TCanvas("cFitStackedB0", "Stacked Representation of Fit", 600, 600);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
	pad1->SetLogy();
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	h_cdphi_data_electrons_inclusive_B0->GetYaxis()->SetRangeUser(8E-5, 12);
	h_cdphi_data_electrons_inclusive_B0->SetTitle("");
	h_cdphi_data_electrons_inclusive_B0->Draw("p");
	h_cdphi_data_electrons_inclusive_B0->SetMarkerStyle(20);
	h_cdphi_data_electrons_inclusive_B0->SetMarkerSize(0.5);
	h_cdphi_data_electrons_inclusive_B0->SetLineColor(kBlack);
	hs_B0_fit->Draw("hist,same");
	h_fit_B0->Draw("same");

	TLatex *tlB0 = new TLatex(0.68, 0.7, "B0");
	tlB0->SetNDC();
	tlB0->SetTextSize(0.2);
	tlB0->Draw("same");

	TLatex *tlFNPB0 = new TLatex(FIT_HIGH, 0.4, Form("F_{NP} = %.3g #pm %.3g", fnp_B0, fnp_B0_err));
	tlFNPB0->SetNDC();
	tlFNPB0->SetTextSize(0.04);
	tlFNPB0->Draw("same");

	TLatex *tlFNPB0pT = new TLatex(FIT_HIGH, 0.45, Form("%g < p_{T} [GeV/c] < %g", pTLow, pTHigh));
	tlFNPB0pT->SetNDC();
	tlFNPB0pT->SetTextSize(0.04);
	tlFNPB0pT->Draw("same");

	TLegend *legB0Fit = new TLegend(FIT_HIGH, 0.6, 0.65, 0.8);
	legB0Fit->SetLineColor(kWhite);
	legB0Fit->AddEntry(h_cdphi_data_electrons_inclusive_B0, "Inclusive Electrons - Swapped Electrons", "P");
	legB0Fit->AddEntry(h_fit_B0, "Fit", "L");
	legB0Fit->AddEntry(h_cdphi_cocktail_multback_B0_scaled, "Photonic Component with HM Background #times (1-F_{NP})", "F");
	legB0Fit->AddEntry(h_cdphi_data_hadrons_B0_scaled, "Non-Photonic Background #times F_{NP}", "F");
	legB0Fit->Draw("same");

	cFitStackedB0->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();
	h_cdphi_data_ratio_B0->SetMarkerStyle(7);
	h_cdphi_data_ratio_B0->SetTitle("");
	h_cdphi_data_ratio_B0->GetYaxis()->CenterTitle();
	h_cdphi_data_ratio_B0->GetYaxis()->SetRangeUser(0.2, 1.67);
	h_cdphi_data_ratio_B0->GetYaxis()->SetTitle("Data / Fit");
	h_cdphi_data_ratio_B0->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B0->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_data_ratio_B0->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B0->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B0->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_data_ratio_B0->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B0->Draw("p");

	TBox *b20 = new TBox(FIT_LOW, 0.8, FIT_HIGH, 1.2);
	b20->SetFillColorAlpha(kGray, 0.3);
	b20->Draw("same");

	cFitStackedB0->cd();

	//------------- B1

	TCanvas *cFitStackedB1 = new TCanvas("cFitStackedB1", "Stacked Representation of Fit", 600, 600);
	TPad *pad1B1 = new TPad("pad1B1", "pad1B1", 0, 0.3, 1, 1);
	pad1B1->SetLogy();
	pad1B1->SetTickx();
	pad1B1->SetTicky();
	pad1B1->SetBottomMargin(0);
	pad1B1->Draw();
	pad1B1->cd();
	h_cdphi_data_electrons_inclusive_B1->SetTitle("");
	h_cdphi_data_electrons_inclusive_B1->GetYaxis()->SetRangeUser(8E-5, 12);
	h_cdphi_data_electrons_inclusive_B1->Draw("p");
	h_cdphi_data_electrons_inclusive_B1->SetMarkerStyle(20);
	h_cdphi_data_electrons_inclusive_B1->SetMarkerSize(0.5);
	h_cdphi_data_electrons_inclusive_B1->SetLineColor(kBlack);
	hs_B1_fit->Draw("hist,same");
	h_fit_B1->Draw("same");

	TLatex *tlB1 = new TLatex(0.68, 0.7, "B1");
	tlB1->SetNDC();
	tlB1->SetTextSize(0.2);
	tlB1->Draw("same");

	TLatex *tlFNPB1 = new TLatex(FIT_HIGH, 0.4, Form("F_{NP} = %.3g #pm %.3g", fnp_B1, fnp_B1_err));
	tlFNPB1->SetNDC();
	tlFNPB1->SetTextSize(0.04);
	tlFNPB1->Draw("same");

	TLatex *tlFNPB1pT = new TLatex(FIT_HIGH, 0.45, Form("%g < p_{T} [GeV/c] < %g", pTLow, pTHigh));
	tlFNPB1pT->SetNDC();
	tlFNPB1pT->SetTextSize(0.04);
	tlFNPB1pT->Draw("same");

	TLegend *legB1Fit = new TLegend(FIT_HIGH, 0.6, 0.65, 0.8);
	legB1Fit->SetLineColor(kWhite);
	legB1Fit->AddEntry(h_cdphi_data_electrons_inclusive_B1, "Inclusive Electrons - Swapped Electrons", "P");
	legB1Fit->AddEntry(h_fit_B1, "Fit", "L");
	legB1Fit->AddEntry(h_cdphi_cocktail_multback_B1_scaled, "Photonic Component with HM Background #times (1-F_{NP})", "F");
	legB1Fit->AddEntry(h_cdphi_data_hadrons_B1_scaled, "Non-Photonic Background #times F_{NP}", "F");
	legB1Fit->Draw("same");

	cFitStackedB1->cd();
	TPad *pad2B1 = new TPad("pad2B1", "pad2B1", 0, 0, 1, 0.3);
	pad2B1->SetTopMargin(0);
	pad2B1->SetBottomMargin(0.3);
	pad2B1->Draw();
	pad2B1->cd();
	pad2B1->SetTickx();
	pad2B1->SetTicky();
	h_cdphi_data_ratio_B1->SetTitle("");
	h_cdphi_data_ratio_B1->SetMarkerStyle(7);
	h_cdphi_data_ratio_B1->GetYaxis()->CenterTitle();
	h_cdphi_data_ratio_B1->GetYaxis()->SetRangeUser(0.2, 1.67);
	h_cdphi_data_ratio_B1->GetYaxis()->SetTitle("Data / Fit");
	h_cdphi_data_ratio_B1->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B1->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_data_ratio_B1->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B1->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B1->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_data_ratio_B1->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B1->Draw("p");

	b20->Draw("same");

	cFitStackedB1->cd();

	if (savePlots)
	{
		if (pTBin == 5)
		{
			cFitStackedB0->SaveAs(Form("FNPPlots/fit_stacked_B0_%s.pdf", "4_6"));
			cFitStackedB1->SaveAs(Form("FNPPlots/fit_stacked_B1_%s.pdf", "4_6"));
		}
		else if (pTBin == 6)
		{
			cFitStackedB0->SaveAs(Form("FNPPlots/fit_stacked_B0_%s.pdf", "6_8"));
			cFitStackedB1->SaveAs(Form("FNPPlots/fit_stacked_B1_%s.pdf", "6_8"));
		}
		else
		{
			cFitStackedB0->SaveAs(Form("FNPPlots/fit_stacked_B0_%s.pdf", pTString.c_str()));
			cFitStackedB1->SaveAs(Form("FNPPlots/fit_stacked_B1_%s.pdf", pTString.c_str()));
		}
	}
}



void plotFitDirect()
{
	//------------- B0

	TCanvas *cFitDirectB0 = new TCanvas("cFitDirectB0", "Direct Representation of Fit", 600, 600);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
	pad1->SetLogy();
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	h_cdphi_data_electrons_inclusive_B0->GetYaxis()->SetRangeUser(8E-5, 12);
	h_cdphi_data_electrons_inclusive_B0->Draw("p");
	h_cdphi_data_electrons_inclusive_B0->SetMarkerStyle(20);
	h_cdphi_data_electrons_inclusive_B0->SetMarkerSize(0.5);
	h_cdphi_cocktail_multback_B0_scaled->Draw("same");
	h_cdphi_data_hadrons_B0_scaled->Draw("same");
	h_fit_B0->Draw("same");

	TLatex *tlB0 = new TLatex(0.68, 0.7, "B0");
	tlB0->SetNDC();
	tlB0->SetTextSize(0.2);
	tlB0->Draw("same");

	TLatex *tlFNPB0 = new TLatex(FIT_HIGH, 0.4, Form("F_{NP} = %.3g #pm %.3g", fnp_B0, fnp_B0_err));
	tlFNPB0->SetNDC();
	tlFNPB0->SetTextSize(0.04);
	tlFNPB0->Draw("same");

	TLatex *tlFNPB0pT = new TLatex(FIT_HIGH, 0.45, Form("%g < p_{T} [GeV/c] < %g", pTLow, pTHigh));
	tlFNPB0pT->SetNDC();
	tlFNPB0pT->SetTextSize(0.04);
	tlFNPB0pT->Draw("same");

	TLegend *legB0Fit = new TLegend(FIT_HIGH, 0.6, 0.65, 0.8);
	legB0Fit->SetLineColor(kWhite);
	legB0Fit->AddEntry(h_cdphi_data_electrons_inclusive_B0, "Inclusive Electrons - Swapped Electrons", "P");
	legB0Fit->AddEntry(h_fit_B0, "Fit", "L");
	legB0Fit->AddEntry(h_cdphi_cocktail_multback_B0_scaled, "Photonic Component with HM Background", "F");
	legB0Fit->AddEntry(h_cdphi_data_hadrons_B0_scaled, "Non-Photonic Background", "F");
	legB0Fit->Draw("same");

	cFitDirectB0->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();
	h_cdphi_data_ratio_B0->SetMarkerStyle(7);
	h_cdphi_data_ratio_B0->GetYaxis()->CenterTitle();
	h_cdphi_data_ratio_B0->GetYaxis()->SetRangeUser(0.2, 1.67);
	h_cdphi_data_ratio_B0->GetYaxis()->SetTitle("Data / Fit");
	h_cdphi_data_ratio_B0->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B0->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_data_ratio_B0->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B0->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B0->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_data_ratio_B0->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B0->Draw("p");

	TBox *b20 = new TBox(FIT_LOW, 0.8, FIT_HIGH, 1.2);
	b20->SetFillColorAlpha(kGray, 0.3);
	b20->Draw("same");

	cFitDirectB0->cd();

	//------------- B1

	TCanvas *cFitDirectB1 = new TCanvas("cFitDirectB1", "Direct Representation of Fit", 600, 600);
	TPad *pad1B1 = new TPad("pad1B1", "pad1B1", 0, 0.3, 1, 1);
	pad1B1->SetLogy();
	pad1B1->SetTickx();
	pad1B1->SetTicky();
	pad1B1->SetBottomMargin(0);
	pad1B1->Draw();
	pad1B1->cd();
	h_cdphi_data_electrons_inclusive_B1->GetYaxis()->SetRangeUser(8E-5, 12);
	h_cdphi_data_electrons_inclusive_B1->Draw("p");
	h_cdphi_data_electrons_inclusive_B1->SetMarkerStyle(20);
	h_cdphi_data_electrons_inclusive_B1->SetMarkerSize(0.5);
	h_cdphi_cocktail_multback_B1_scaled->Draw("same");
	h_cdphi_data_hadrons_B1_scaled->Draw("same");
	h_fit_B1->Draw("same");

	TLatex *tlB1 = new TLatex(0.68, 0.7, "B1");
	tlB1->SetNDC();
	tlB1->SetTextSize(0.2);
	tlB1->Draw("same");

	TLatex *tlFNPB1 = new TLatex(FIT_HIGH, 0.4, Form("F_{NP} = %.3g #pm %.3g", fnp_B1, fnp_B1_err));
	tlFNPB1->SetNDC();
	tlFNPB1->SetTextSize(0.04);
	tlFNPB1->Draw("same");

	TLatex *tlFNPB1pT = new TLatex(FIT_HIGH, 0.45, Form("%g < p_{T} [GeV/c] < %g", pTLow, pTHigh));
	tlFNPB1pT->SetNDC();
	tlFNPB1pT->SetTextSize(0.04);
	tlFNPB1pT->Draw("same");

	TLegend *legB1Fit = new TLegend(FIT_HIGH, 0.6, 0.65, 0.8);
	legB1Fit->SetLineColor(kWhite);
	legB1Fit->AddEntry(h_cdphi_data_electrons_inclusive_B1, "Inclusive Electrons - Swapped Electrons", "P");
	legB1Fit->AddEntry(h_fit_B1, "Fit", "L");
	legB1Fit->AddEntry(h_cdphi_cocktail_multback_B1_scaled, "Photonic Component with HM Background", "F");
	legB1Fit->AddEntry(h_cdphi_data_hadrons_B1_scaled, "Non-Photonic Background", "F");
	legB1Fit->Draw("same");

	cFitDirectB1->cd();
	TPad *pad2B1 = new TPad("pad2B1", "pad2B1", 0, 0, 1, 0.3);
	pad2B1->SetTopMargin(0);
	pad2B1->SetBottomMargin(0.3);
	pad2B1->Draw();
	pad2B1->cd();
	pad2B1->SetTickx();
	pad2B1->SetTicky();
	h_cdphi_data_ratio_B1->SetMarkerStyle(7);
	h_cdphi_data_ratio_B1->GetYaxis()->CenterTitle();
	h_cdphi_data_ratio_B1->GetYaxis()->SetRangeUser(0.2, 1.67);
	h_cdphi_data_ratio_B1->GetYaxis()->SetTitle("Data / Fit");
	h_cdphi_data_ratio_B1->GetYaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B1->GetYaxis()->SetTitleOffset(0.4);
	h_cdphi_data_ratio_B1->GetYaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B1->GetXaxis()->SetTitleSize(0.085);
	h_cdphi_data_ratio_B1->GetXaxis()->SetTitleOffset(1.2);
	h_cdphi_data_ratio_B1->GetXaxis()->SetLabelSize(0.085);
	h_cdphi_data_ratio_B1->Draw("p");

	b20->Draw("same");

	cFitDirectB1->cd();
}

void plotTailNorm()
{
	h_cdphi_data_hadrons_B0_tailnorm->SetLineColor(kOrange + 7);
	h_cdphi_data_hadrons_B0_tailnorm->SetLineWidth(2);
	h_cdphi_cocktail_multback_B0_tailnorm->SetLineColor(kBlue);
	h_cdphi_cocktail_multback_B0_tailnorm->SetLineWidth(2);

	h_cdphi_data_hadrons_B1_tailnorm->SetLineColor(kOrange + 7);
	h_cdphi_data_hadrons_B1_tailnorm->SetLineWidth(2);
	h_cdphi_cocktail_multback_B1_tailnorm->SetLineColor(kBlue);
	h_cdphi_cocktail_multback_B1_tailnorm->SetLineWidth(2);

	TCanvas *cTailNormalizationB0 = new TCanvas("cTailNormalizationB0", "cTailNormalizationB0", 600, 600);
	cTailNormalizationB0->SetLogy();
	h_cdphi_data_hadrons_B0_tailnorm->SetTitle("");
	h_cdphi_data_hadrons_B0_tailnorm->GetYaxis()->SetTitle("A.U.");
	h_cdphi_data_hadrons_B0_tailnorm->GetYaxis()->SetRangeUser(1E-3, 1);
	h_cdphi_data_hadrons_B0_tailnorm->Draw();
	h_cdphi_cocktail_multback_B0_tailnorm->Draw("same");

	TLatex *tlB0 = new TLatex(0.15, 0.76, "B0");
	tlB0->SetNDC();
	tlB0->SetTextSize(0.12);
	tlB0->DrawClone("same");

	TCanvas *cTailNormalizationB1 = new TCanvas("cTailNormalizationB1", "cTailNormalizationB1", 600, 600);
	cTailNormalizationB1->SetLogy();
	h_cdphi_data_hadrons_B1_tailnorm->SetTitle("");
	h_cdphi_data_hadrons_B1_tailnorm->GetYaxis()->SetTitle("A.U.");
	h_cdphi_data_hadrons_B1_tailnorm->GetYaxis()->SetRangeUser(1E-3, 1);
	h_cdphi_data_hadrons_B1_tailnorm->Draw();
	h_cdphi_cocktail_multback_B1_tailnorm->Draw("same");

	TLatex *tlB1 = new TLatex(0.15, 0.76, "B1");
	tlB1->SetNDC();
	tlB1->SetTextSize(0.12);
	tlB1->DrawClone("same");
}


void constructFNP()
{
	float pT[7] = {1.25, 1.75, 2.25, 2.75, 3.5, 5.0, 7.0};

	//FNP with deadmap and |dz| < 0.14 cut
	float fnp_B0[7] = {8.90288e-13, 0.00262489, 0.10425, 0.185711, 0.367047, 0.597861, 1.40443e-14};
	float fnp_err_B0[7] = {0.120909, 0.22654, 0.0574432, 0.0622212, 0.0539654, 0.0595762, 2.17788e-05};

	float fnp_B1[7] = {0.0784, 0.280456, 0.279993, 0.410501, 0.492351, 1, 1};
	float fnp_err_B1[7] = {0.0315309, 0.0351332, 0.0388977, 0.0390457, 0.0339418, 4.56275e-07, 9.81743e-05};

	float err_x[7] = {0.0};

	//Construct ratio
	float fnp_ratio_B0overB1[7] = {0.0};
	float fnp_ratio_B0overB1_error[7] = {0.0};

	for (int i = 0; i < 7; i++)
	{
		fnp_ratio_B0overB1[i] = fnp_B0[i] / fnp_B1[i];
		fnp_ratio_B0overB1_error[i] = fnp_ratio_B0overB1[i] * TMath::Sqrt(pow(fnp_err_B0[i] / fnp_B0[i], 2.0) + pow(fnp_err_B1[i] / fnp_B1[i], 2.0));
	}

	g_fnp_B0 = new TGraphErrors(7, pT, fnp_B0, err_x, fnp_err_B0);
	g_fnp_B1 = new TGraphErrors(7, pT, fnp_B1, err_x, fnp_err_B1);

	g_fnp_B0->SetTitle("");
	g_fnp_B1->SetTitle("");

	g_fnp_B0->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_fnp_B0->GetXaxis()->SetLabelFont(62);
	g_fnp_B0->GetXaxis()->SetTitleFont(62);

	g_fnp_B0->GetYaxis()->SetTitle("F_{NP}");
	g_fnp_B0->GetYaxis()->SetLabelFont(62);
	g_fnp_B0->GetYaxis()->SetTitleFont(62);

	g_fnp_B1->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_fnp_B1->GetXaxis()->SetLabelFont(62);
	g_fnp_B1->GetXaxis()->SetTitleFont(62);

	g_fnp_B1->GetYaxis()->SetTitle("F_{NP}");
	g_fnp_B1->GetYaxis()->SetLabelFont(62);
	g_fnp_B1->GetYaxis()->SetTitleFont(62);

	g_fnp_ratio = new TGraphErrors(7, pT, fnp_ratio_B0overB1, err_x, fnp_ratio_B0overB1_error);
}

void plotFNP()
{
	//FNP in B0
	TCanvas *cFNPB0 = new TCanvas("cFNPB0", "cFNPB0", 600, 600);
	g_fnp_B0->SetMarkerStyle(20);
	g_fnp_B0->GetYaxis()->SetRangeUser(0, 1.0);
	g_fnp_B0->GetXaxis()->SetTitleOffset(1.2);
	g_fnp_B0->DrawClone("AP");

	TLatex *tlB0 = new TLatex(0.68, 0.76, "B0");
	tlB0->SetNDC();
	tlB0->SetTextSize(0.16);
	tlB0->DrawClone("same");

	//FNP in B1
	TCanvas *cFNPB1 = new TCanvas("cFNPB1", "cFNPB1", 600, 600);
	g_fnp_B1->SetMarkerStyle(20);
	g_fnp_B1->GetYaxis()->SetRangeUser(0, 0.49);
	g_fnp_B1->GetXaxis()->SetTitleOffset(1.2);
	g_fnp_B1->Draw("AP");

	TLatex *tlB1 = new TLatex(0.68, 0.76, "B1");
	tlB1->SetNDC();
	tlB1->SetTextSize(0.16);
	tlB1->Draw("same");

	TCanvas *cFNP = new TCanvas("cFNP", "cFNP", 600, 600);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	g_fnp_B0->GetYaxis()->SetRangeUser(0.0, 1.0);
	g_fnp_B0->GetXaxis()->SetTitleOffset(1.2);
	g_fnp_B0->SetMarkerColor(kRed);
	g_fnp_B0->SetLineColor(kRed);
	g_fnp_B0->Draw("AP");
	g_fnp_B1->Draw("P, same");

	TLegend *legFNP = new TLegend(0.18, 0.55, 0.4, 0.8);
	legFNP->SetLineColor(kWhite);
	legFNP->AddEntry(g_fnp_B0, "B0", "P");
	legFNP->AddEntry(g_fnp_B1, "B1", "P");
	legFNP->Draw("same");

	cFNP->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();
	g_fnp_ratio->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	g_fnp_ratio->GetYaxis()->SetTitle("Ratio B0/B1");
	g_fnp_ratio->SetTitle("");
	g_fnp_ratio->SetMarkerStyle(20);
	g_fnp_ratio->GetYaxis()->SetRangeUser(0.4, 1.5);
	g_fnp_ratio->GetYaxis()->SetTitleSize(0.085);
	g_fnp_ratio->GetYaxis()->SetTitleOffset(0.4);
	g_fnp_ratio->GetYaxis()->SetTitleFont(62);
	g_fnp_ratio->GetYaxis()->SetLabelFont(62);
	g_fnp_ratio->GetYaxis()->SetLabelSize(0.085);
	g_fnp_ratio->GetXaxis()->SetTitleSize(0.085);
	g_fnp_ratio->GetXaxis()->SetTitleOffset(1.2);
	g_fnp_ratio->GetXaxis()->SetLabelSize(0.085);
	g_fnp_ratio->GetXaxis()->SetTitleFont(62);
	g_fnp_ratio->GetXaxis()->SetLabelFont(62);
	g_fnp_ratio->Draw("AP");

	TLine *l = new TLine(1.0, 1.0, 7.5, 1.0);
	l->SetLineStyle(7);
	l->Draw("same");

	cFNP->cd();
}

void plotElectronsData()
{
	TCanvas *cElectronsData = new TCanvas("cElectronsData", "cElectronsData", 600, 600);
	cElectronsData->SetLogy();
	h_elec_pT_data_rebinned->SetLineWidth(2);
	h_elec_pT_data_rebinned->GetXaxis()->SetTitleOffset(1.2);
	h_elec_pT_data_rebinned->GetYaxis()->SetTitleOffset(1.4);
	formatHistogram(h_elec_pT_data_rebinned, "p_{T} [GeV/c]", "Number of Electrons", "");
	h_elec_pT_data_rebinned->Draw();
}

void plotPhotonicRatios()
{
	TCanvas *cPhotonicRatios = new TCanvas("cPhotonicRatios", "cPhotonicRatios", 600, 600);
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1);
	pad1->SetTickx();
	pad1->SetTicky();
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->cd();
	pad1->SetLogy();
	h_elec_pT_pizeros_rebinned->GetXaxis()->SetTitleOffset(1.2);
	h_elec_pT_pizeros_rebinned->Draw();
	h_elec_pT_etas_rebinned->Draw("same");
	h_elec_pT_photons_rebinned->Draw("same");

	TLegend *legPhotonicElectrons = new TLegend(0.4, 0.55, 0.8, 0.8);
	legPhotonicElectrons->SetLineColor(kWhite);
	legPhotonicElectrons->AddEntry(h_elec_pT_pizeros_rebinned, "Electrons from #pi^{0} (Dalitz + Conversion)", "L");
	legPhotonicElectrons->AddEntry(h_elec_pT_etas_rebinned, "Electrons from #eta (Dalitz + Conversion)", "L");
	legPhotonicElectrons->AddEntry(h_elec_pT_photons_rebinned, "Electrons from Direct #gamma (Conversion)", "L");
	legPhotonicElectrons->Draw("same");

	cPhotonicRatios->cd();
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.3);
	pad2->Draw();
	pad2->cd();
	pad2->SetTickx();
	pad2->SetTicky();
	h_elec_pT_pizeros_fraction_rebinned->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->SetTitle("Ratio to Photonic");
	h_elec_pT_pizeros_fraction_rebinned->SetTitle("");
	h_elec_pT_pizeros_fraction_rebinned->SetMarkerStyle(20);
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->CenterTitle();
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->SetRangeUser(-0.19, 1.11);
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->SetTitleSize(0.085);
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->SetTitleOffset(0.4);
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->SetTitleFont(62);
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->SetLabelFont(62);
	h_elec_pT_pizeros_fraction_rebinned->GetYaxis()->SetLabelSize(0.085);
	h_elec_pT_pizeros_fraction_rebinned->GetXaxis()->SetTitleSize(0.085);
	h_elec_pT_pizeros_fraction_rebinned->GetXaxis()->SetTitleOffset(1.2);
	h_elec_pT_pizeros_fraction_rebinned->GetXaxis()->SetLabelSize(0.085);
	h_elec_pT_pizeros_fraction_rebinned->GetXaxis()->SetTitleFont(62);
	h_elec_pT_pizeros_fraction_rebinned->GetXaxis()->SetLabelFont(62);
	h_elec_pT_pizeros_fraction_rebinned->Draw();
	h_elec_pT_etas_fraction_rebinned->Draw("same");
	h_elec_pT_photons_fraction_rebinned->Draw("same");

	TLine *l = new TLine(0.0, 0.5, 10.0, 0.5);
	l->SetLineStyle(7);
	l->Draw("same");

	cPhotonicRatios->cd();
}

void calculateFNP()
{
	//Do analysis
	loadHistos();
	setErrors();
	removeSwappedBackground();
	getNumTracks();
	rebinHistos();

	calculatePhotonicRatios();
	setFormat();
	addCocktail();
	fitPhotonicSideband();

	constructMultiplicityBackground();
	normalizeHistograms();
	//addPhotonicErrorToData();

	integrateFNP();
	fitFNP();
	constructFitVisualization();
	constructFNP();

	//Do all of the plotting
	gStyle->SetOptStat(0);
	//plotElectronsData();
	plotClustersPerHadronTrack();
	plotFitComponents();
	plotFitStacked();
	plotPhotonicComponentWithoutBackground();
	//plotPhotonicReweighted();
	//plotFitDirect();
	//plotTailNorm();
	//plotFNP();
	//plotPhotonicRatios();
}
