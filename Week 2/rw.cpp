#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y);
  PARAMETER(logSdRw);
  PARAMETER(logSdObs);
  PARAMETER(lam0);
  PARAMETER_VECTOR(lam);

  int timeSteps=y.size();
  Type sdRw=exp(logSdRw);
  Type sdObs=exp(logSdObs);

  Type ans=-dnorm(lam(0),lam0,sdRw,true);

  for(int i=1;i<timeSteps;i++){    
    ans+=-dnorm(lam(i),lam(i-1),sdRw,true);
  }

  for(int i=0;i<timeSteps;i++){
    ans+=-dnorm(y(i),lam(i),sdObs,true); 
  } 

  return ans;
}
