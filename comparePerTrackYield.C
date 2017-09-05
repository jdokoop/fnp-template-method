//-----------------------------------------------------------------------------------
// Compare the per-track cluster yield from photonic electrons between the new large
// scale simulations and the older smaller scale simulations
//
// -> Large Scale Sims (w/ default reweighting) Sims/Cocktail090417
// -> Small Scale Sims (w/ default reweighting) Sims/Cocktail090117

// Sept. 5 2017
//------------------------------------------------------------------------------------

//---------------------------------------------------
// Variables
//---------------------------------------------------

//Number of pT bins
const int NBINS = 13;

//pT bins
float pT_low[NBINS]  = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
float pT_high[NBINS] = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};

//---------------------------------------------------
// Functions
//---------------------------------------------------

void comparePerTrackYield()
{
	
}