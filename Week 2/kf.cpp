#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  PARAMETER(logSdRw);
  PARAMETER(logSdObs);
  PARAMETER(lam0);
  int N=y.size();
  vector<Type> lam(N), lamVar(N), lamPred(N), lamPredVar(N), lamSmooth(N), lamSmoothVar(N), yPredVar(N), w(N);
  Type varLam, varY, varFrac;
  Type ans=Type(0.0);
  varLam=exp(2.0*logSdRw);
  varY=exp(2.0*logSdObs);
  lamPred(0)=lam0;
  lamPredVar(0)=varLam+Type(1000.0); 
  yPredVar(0)=lamPredVar(0)+varY;
  w(0)=y(0)-lamPred(0); 
  lam(0)=lamPred(0)+lamPredVar(0)/yPredVar(0)*w(0); 
  lamVar(0)=lamPredVar(0)-lamPredVar(0)/yPredVar(0)*lamPredVar(0); 
  ans+=Type(0.5)*log(yPredVar(0))+Type(0.5)*w(0)/yPredVar(0)*w(0);
  for(int i=1; i<N; ++i){
    lamPred(i)=lam(i-1);
    lamPredVar(i)=lamVar(i-1)+varLam; 
    yPredVar(i)=lamPredVar(i)+varY;
    w(i)=y(i)-lamPred(i); 
    lam(i)=lamPred(i)+lamPredVar(i)/yPredVar(i)*w(i); 
    lamVar(i)=lamPredVar(i)-lamPredVar(i)/yPredVar(i)*lamPredVar(i); 
    ans+=0.5*log(yPredVar(i))+0.5*w(i)/yPredVar(i)*w(i);
  }
  
  // ----------------- smoothing part 
  lamSmooth(N-1)=lam(N-1); 
  lamSmoothVar(N-1)=lamVar(N-1); 
  for(int i=N-2; i>=0; --i){
    varFrac=lamVar(i)/lamPredVar(i+1); 
    lamSmooth(i)=lam(i)+varFrac*(lamSmooth(i+1)-lamPred(i+1)); 
    lamSmoothVar(i)=lamVar(i)+varFrac*(lamSmoothVar(i+1)-lamPredVar(i+1))*varFrac; 
  }

  REPORT(lamSmooth);
  REPORT(lamSmoothVar);
  return ans;
}



