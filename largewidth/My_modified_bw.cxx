#include "Riostream.h"

#include "My_modified_bw.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include "Math/PdfFuncMathCore.h"
#include <math.h>
#include "TMath.h"
using namespace ROOT::Math;

ClassImp(My_modified_bw)

My_modified_bw::My_modified_bw(const char *name, const char *title,
			       RooAbsReal& _x,
			       RooAbsReal& _width,
			       RooAbsReal& _m0,
             RooAbsReal& _k,
             RooAbsReal& _mass) :
RooAbsPdf(name,title),
  x("x","x",this,_x),
  width("width","width",this,_width),
  m0("m0","m0",this,_m0),
  k("k","k",this,_k),
  mass("mass","mass",this,_mass)
{
}

My_modified_bw::My_modified_bw(const My_modified_bw& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  width("width",this,other.width),
  m0("m0",this,other.m0),
  k("k",this,other.k),
  mass("mass",this,other.mass)
{
}


Double_t My_modified_bw::evaluate() const
{
  Double_t ln_k = TMath::Abs(TMath::Log(k));
  Double_t ln_m0 = TMath::Log(m0);
  // Double_t ln_k = TMath::Abs(k);
  // Double_t ln_m0 = TMath::Abs(m0);
  Double_t xvalue = TMath::Power(x, 1.);
  Double_t widthvalue = TMath::Power(width, 1.);
  Double_t massvalue = TMath::Power(mass, 1.);
  // if (x + massvalue < 0) {
  //   // return -(x + massvalue) * 100;
  //   return 0;
  // }
  Double_t bw = 1/(xvalue * xvalue + 0.25 * widthvalue * widthvalue);
  Double_t ret = ROOT::Math::lognormal_pdf(xvalue, ln_m0, ln_k, massvalue) * bw;
  return ret;

}
