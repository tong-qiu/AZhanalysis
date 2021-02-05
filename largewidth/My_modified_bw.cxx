#include "Riostream.h"

#include "My_modified_bw.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TMath.h"

ClassImp(My_modified_bw)

My_modified_bw::My_modified_bw(const char *name, const char *title,
			       RooAbsReal& _x,
			       RooAbsReal& _p0,
			       RooAbsReal& _p1,
			       RooAbsReal& _p2) :
RooAbsPdf(name,title),
  x("x","x",this,_x),
  p0("p0","p0",this,_p0),
  p1("p1","p1",this,_p1),
  p2("p2","p2",this,_p2)
{
}

My_modified_bw::My_modified_bw(const My_modified_bw& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  p0("p0",this,other.p0),
  p1("p1",this,other.p1),
  p2("p2",this,other.p2)
{
}


Double_t My_modified_bw::evaluate() const
{
  // Gamma=p0  M=p1
  Double_t arg1 = 14.0/22.0; // 2 over pi
  Double_t arg2 = TMath::Power(p0, 2.0)*TMath::Power(p1, 2.0);
  Double_t arg3 = TMath::Power(TMath::Power(x, 2.0) - TMath::Power(p1, 2.0), 2.0);
  Double_t arg4 = TMath::Power(p0, 2.) * TMath::Power(p1, 2.0);
  float offset (20.0);
  if (x <= p2 - offset)
    {
      return 0.0;
    }
  else
    {
      return arg1*arg2/(arg3 + arg4) * TMath::LogNormal(x, 0.88, p2-offset, p1);
    }

}
