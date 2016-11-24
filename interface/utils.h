//****************************************************************************************
//*** collection of variuos utilities used by different plugins
//***
//****************************************************************************************

#include <string>

#include "TH1.h"
#include "TH2.h"

using namespace std;

//----------Smart Profile-----------------------------------------------------------------
//---generate profile from TH2 istograms using an adaptive window instead of the full
//   Y range as in TProfile.
//   For the moment the parameters are not configurable...
void GetIterativeProfile(TH2* h2d, TH1* prof);

//---------*------------------------------------------------------------------------------
