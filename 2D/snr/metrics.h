#ifndef __METRICS_H__
#define __METRICS_H__

#include <iostream>
#include <cmath>


class Metrics 
{

public:

  Metrics(bool text=false, double snrThresh=200, double relErrThresh=1e-3, double absErrThresh=1e-6):
    isText(text), SNR_THRESH(snrThresh), REL_ERR_THRESH(relErrThresh), ABS_ERR_THRESH(absErrThresh),
    sigPow(0), errPow(0), nRelErr(0), nAbsErr(0), nErr(0), nKept(0), nDiscard(0), nInf(0), nNaN(0)
  { }

  void update(std::istream & in1, std::istream & in2);
  void summary(void);
  bool text_next(std::istream & in, double & t);
  bool binary_next(std::istream & in, double & t);
  void text_update(std::istream & in1, std::istream & in2);
  void binary_update(std::istream & in1, std::istream & in2);
  void update_snr(double t1, double t2);

  double getSNR() {
    return 20.0 * log10(sigPow / errPow);
  }

  size_t getSampleCount() const {
    return nKept + nDiscard;
  }

  size_t getRelErrCount() const {
    return nRelErr;
  }

  size_t getAbsErrCount() const {
    return nAbsErr;
  }

  size_t getErrCount() const {
    return nErr;
  }

  size_t getKeptCount() const {
    return nKept;
  }

  size_t getDiscardedCount() const {
    return nDiscard;
  }

  size_t getInfSampleCount() const {
    return nInf;
  }
  
  size_t getNaNSampleCount() const {
    return nNaN;
  }

  bool good() {
    return (getSNR() >= SNR_THRESH);
  }

private:

  const double SNR_THRESH;
  const double REL_ERR_THRESH;
  const double ABS_ERR_THRESH;

  bool isText;
  long double sigPow;
  long double errPow;
  size_t nRelErr;
  size_t nAbsErr;
  size_t nErr;
  size_t nKept;
  size_t nDiscard;
  size_t nInf;
  size_t nNaN;

};

#endif
