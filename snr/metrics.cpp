#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "metrics.h"

using namespace std;

void Metrics::update(istream & in1, istream & in2) 
{
  if (in1.good() && in2.good()) {
    if (isText) {
      text_update(in1, in2);
    } else {
      binary_update(in1, in2);
    }
  }
}

void Metrics::summary(void)
{
  cout << "SNR THRESHOLD: " << SNR_THRESH << " dB" << endl;
  cout << "ABS THRESHOLD: " << ABS_ERR_THRESH << endl;
  cout << "REL THRESHOLD: " << REL_ERR_THRESH << endl;
  cout << "SNR: " << getSNR() << " dB" << endl;

  double percent = (double)nErr / (double)nKept * 100;
  cout << nErr << " / " << nKept << " samples above error thresholds (" << percent << "%)" << endl;

  if (nDiscard) {
    cout << "NOTE: " << nDiscard << " samples were discarded (e.g. not numeric)." << endl;
  }
  if (nInf) {
    cout << "WARNING: " << nInf << " samples were inf." << endl;
  }
  if (nNaN) {
    cout << "WARNING: " << nNaN << " samples were NaN." << endl;
  }
}

bool Metrics::text_next(istream & in, double & t) 
{
  in >> t;
  if (in.fail() && !in.eof()) {
    in.clear();
    
    // Is it Inf or NaN?
    string s;
    in >> s;
    transform(s.begin(), s.end() ,s.begin(), ::toupper);
    if (s == "INF") {
      t = INFINITY;
    } else if (s == "NAN") {
      t = NAN;
    } else {
      return false;
    }
  }
  return true;
}

void Metrics::text_update(istream & in1, istream & in2) 
{
  double t1, t2;
  bool r1 = text_next(in1, t1);
  bool r2 = text_next(in2, t2);

  if (!(r1 && r2)) {
    ++nDiscard;
  } else if (isnan(t1) || isnan(t2)) {
    ++nNaN;
    ++nDiscard;
  } else if (isinf(t1) || isinf(t2)) {
    ++nInf;
    ++nDiscard;
  } else {
    update_snr(t1, t2);
    ++nKept;
  }
}

bool Metrics::binary_next(std::istream & in, int & t) 
{
  t = in.get();
  return in.good();
}

void Metrics::binary_update(std::istream & in1, std::istream & in2) 
{
  int t1, t2;
  bool r1 = binary_next(in1, t1);
  bool r2 = binary_next(in2, t2);
  if (!(r1 && r2)) {
    ++nDiscard;
  } else {
    update_snr(t1, t2);
    ++nKept;
  }
}

void Metrics::update_snr(double t1, double t2) 
{
  double mean = (t1 + t2) / 2.0;
  double sigSq = t1 * t1;
  double errSq = (t1 - t2) * (t1 - t2);
  double absErr = sqrt(errSq);
  double relErr = absErr / mean;

  sigPow += sigSq;
  errPow += errSq;
  nRelErr += (relErr > REL_ERR_THRESH);
  nAbsErr += (absErr > ABS_ERR_THRESH);

  if (relErr > REL_ERR_THRESH && absErr > ABS_ERR_THRESH) {
    ++nErr;
    printf("%e != %e    (abs=%e, rel=%e)\n", t1, t2, absErr, relErr);
  }
}

