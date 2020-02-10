# include <TMB.hpp> 

template<class Type>
vector<Type> trans(vector<Type> alpha){
	int dim=alpha.size(); 
	vector<Type> p(dim+1); 
	vector<Type> expa=exp(alpha); 
	Type s=sum(expa);
    Type lastp=1;
    for(int i=0; i<dim; ++i){
	    p(i)=expa(i)/(Type(1)+s); 
	    lastp-=p(i);
	} 
	p(dim)=lastp; 
	return p;
	}

template<class Type>
Type objective_function<Type>::operator()(){
	DATA_VECTOR(X); 
	PARAMETER_VECTOR(alpha); 
	vector<Type> p=trans(alpha); 
	Type nll = -dmultinom(X,p,true); 
	ADREPORT(p);
    return nll;
}