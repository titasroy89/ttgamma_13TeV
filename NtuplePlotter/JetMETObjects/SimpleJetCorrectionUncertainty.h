#ifndef SimpleJetCorrectionUncertainty_h
#define SimpleJetCorrectionUncertainty_h

//#include "/uscms_data/d3/troy2012/cmssw/CondFormats/Serialization/interface/Serializable.h"

#include <string>
#include <vector>
class JetCorrectorParameters;

class SimpleJetCorrectionUncertainty 
{
 public:
  SimpleJetCorrectionUncertainty();
  SimpleJetCorrectionUncertainty(const std::string& fDataFile);
  SimpleJetCorrectionUncertainty(const JetCorrectorParameters& fParameters);
  ~SimpleJetCorrectionUncertainty();
  const JetCorrectorParameters& parameters() const {return *mParameters;}
  float uncertainty(const std::vector<float>& fX, float fY, bool fDirection) const;

 private:
  SimpleJetCorrectionUncertainty(const SimpleJetCorrectionUncertainty&);
  SimpleJetCorrectionUncertainty& operator= (const SimpleJetCorrectionUncertainty&);
  int findBin(const std::vector<float>& v, float x) const;
  float uncertaintyBin(unsigned fBin, float fY, bool fDirection) const;
  float linearInterpolation (float fZ, const float fX[2], const float fY[2]) const;
  JetCorrectorParameters* mParameters;
};

#endif


