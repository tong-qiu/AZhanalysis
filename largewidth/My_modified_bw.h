#ifndef MY_MODIFIED_BW
#define MY_MODIFIED_BW

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class My_modified_bw : public RooAbsPdf {
 public:
  My_modified_bw() {} ;
  My_modified_bw(const char *name, const char *title,
		 RooAbsReal& _x,
		 RooAbsReal& _width,
		 RooAbsReal& _m0,
     RooAbsReal& _k,
     RooAbsReal& _mass);
  My_modified_bw(const My_modified_bw& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new My_modified_bw(*this,newname); }
  inline virtual ~My_modified_bw() { }

 protected:

  RooRealProxy x ;
  RooRealProxy width ;
  RooRealProxy m0 ;
  RooRealProxy k ;
  RooRealProxy mass ;

  Double_t evaluate() const ;

 private:

  ClassDef(My_modified_bw,1) // Your description goes here...
    };

#endif
