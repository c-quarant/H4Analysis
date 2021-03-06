#ifndef __WFCLASS_H__
#define __WFCLASS_H__

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "Math/Interpolator.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TH1F.h"
#include "TF1.h"

using namespace std;

#define MAX_INTERPOLATOR_POINTS 10000

struct WFBaseline
{
    float baseline;
    float rms;
    float slope;
    float k;
    float chi2;
};

struct WFFitResults
{
    double ampl;
    double time;
    double chi2;
    int    fitStatus_;
};      

class WFClass : public TObject
{
public:
    //---ctors---
    WFClass() {};
    WFClass(int polarity, float tUnit);
    //---dtor---
    ~WFClass() {};

    //---getters---
    inline const vector<double>* GetSamples() {return &samples_;};
    inline const vector<double>* GetNoiseFiltSamples() {return &noiseFiltSamples_;};
    inline float           GetBaseline() {return baseline_;}
    inline int             GetNSample() {return samples_.size();};
    inline float           GetTUnit() {return tUnit_;};
    float                  GetAmpMax(int min=-1, int max=-1);
    WFFitResults           GetInterpolatedAmpMax(int min=-1, int max=-1, int nFitSamples=7);
    pair< pair<float, float>, float>     GetTime(string method, vector<float>& params); 
    pair< pair<float, float>, float>     GetTimeCF(float frac, int nFitSamples=5, int min=-1, int max=-1);
    pair< pair<float, float>, float>     GetTimeLE(float thr, int nmFitSamples=1, int npFitSamples=3, int min=-1, int max=-1);
    float                  GetIntegral(int min=-1, int max=-1);
    float                  GetModIntegral(int min=-1, int max=-1);
    virtual float          GetSignalIntegral(int riseWin, int fallWin);
    //---setters---
    inline void            SetTrigRef(float trigRef){trigRef_ = trigRef;};
    void                   SetSignalWindow(int min, int max);
    void                   SetBaselineWindow(int min, int max);
    void                   SetTemplate(TH1* templateWF=NULL);
    //---utils---
    void                   Reset();
    void                   AddSample(float sample) {samples_.push_back(polarity_*sample);};
    void                   AddFiltSample(float sample) {noiseFiltSamples_.push_back(polarity_*sample);};
    WFBaseline             SubtractBaseline(int min=-1, int max=-1);
    WFFitResults           TemplateFit(bool NoiseCut, float offset=0., int lW=0, int hW=0);//---if NoiseCut==1 use noise filtered samples
    void                   EmulatedWF(WFClass& wf, float rms, float amplitude, float time);
    void                   FFT(WFClass& wf, float tau, int cut);
    void                   BWFilter(WFClass& wf, float tau, int cut);
    void                   Print();
    float                  GetTriggerRef();
    //---operators---
    WFClass&               operator=(const WFClass& origin);
    WFClass                operator-(const WFClass& sub);
    WFClass                operator+(const WFClass& add);
    WFClass&               operator-=(const WFClass& sub);
    WFClass&               operator+=(const WFClass& add);
protected:
    //---utils---
    float                 BaselineRMS();
    float                 LinearInterpolation(float& A, float& B, const int& min, const int& max);
    double                TemplateChi2(const double* par=NULL);//---if NoiseCut==1 use noise filtered samples
    
protected:
    vector<double> samples_;
    vector<double> noiseFiltSamples_;
    vector<double> chi2samples_;
    float          tUnit_;
    int            polarity_;
    float          trigRef_;
    int            sWinMin_;
    int            sWinMax_;
    int            bWinMin_;
    int            bWinMax_;
    int            maxSample_;
    float          fitAmpMax_;
    float          fitTimeMax_;
    float          fitChi2Max_;
    float          baseline_;
    float          bRMS_;
    int            cfSample_;
    int            cfSampleMirror_;
    float          cfFrac_;
    float          cfTime_;
    float          cfTimeMirror_;
    int            leSample_;
    int            leSampleMirror_;
    float          leThr_;
    float          leTime_;
    float          leTimeMirror_;
    float          chi2cf_;
    float          chi2le_;
    int            fWinMin_;
    int            fWinMax_;
    float          tempFitTime_;
    float          tempFitAmp_;
    int            fStatus_;
    ROOT::Math::Interpolator* interpolator_;
};

#endif
