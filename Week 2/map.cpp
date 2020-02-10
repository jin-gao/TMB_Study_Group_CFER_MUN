#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator()(){
DATA_VECTOR(count); 
DATA_FACTOR(spray); 
PARAMETER_VECTOR(logAlpha);
Type nll = 0;
for(int i=0; i<count.size(); ++i){
Type lambda = exp(logAlpha(spray(i))); 
nll += -dpois(count(i),lambda,true);
}
return nll;
}