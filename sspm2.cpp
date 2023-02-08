#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x)); 
}

template<class Type>
Type posfun(Type x, Type eps, Type &pen){
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(y1);
  DATA_VECTOR(y2);
  DATA_VECTOR(Catch); 
  DATA_SCALAR(eps);
  
  PARAMETER(logit_D0);
  PARAMETER(log_r);
  PARAMETER(log_K);
  PARAMETER(log_z);
  PARAMETER(log_q1);
  PARAMETER(log_q2);
  PARAMETER(log_sig1);
  PARAMETER(log_sig2);
  PARAMETER(log_tau);
  PARAMETER_VECTOR(log_Dep);
  
  //int n=y.size();
  int n = Catch.size();
  Type r = exp(log_r);
  Type K = exp(log_K);
  Type z = exp(log_z);
  Type q1 = exp(log_q1);
  Type q2 = exp(log_q2);
  Type sig1 = exp(log_sig1);
  Type sig2 = exp(log_sig2);
  Type tau = exp(log_tau);
  Type MSY = r*K/(pow((z+1),((z+1)/z)));
  Type Bmsy = K*pow((z+1), -1/z);
  
  vector<Type> Dep(n); 
  vector<Type> B(n); 
  vector<Type> mu(n); 
  //Type MSY;
  //Fill in the formula below
  //MSY = 0;
  
  Type pen;
  pen=0;
  Type nll;
  nll=0;
  
  for(int i=0;i<n;i++){ 
    Dep(i) = exp(log_Dep(i));
    B(i) = K*Dep(i);
  }
  
  mu(0) = 1/(1+exp(-logit_D0));
  nll += (-1.0)*dnorm(log(Dep(0)),log(mu(0)),tau,true);
  //Dep(0) = mu(0);
  
  for(int i=1;i<n;i++){   
    mu(i) = posfun(Dep(i-1) + (r/z)*Dep(i-1)*(1-pow(Dep(i-1),z)) - Catch(i-1)/K, eps, pen);
    nll += (-1.0)*dnorm(log_Dep(i),log(mu(i)),tau,true);
    nll += pen;
  }
  
  for(int i=0;i<n;i++){
    if(!isNA(y1(i))){
      nll += (-1.0) * dnorm(log(y1(i)), log(q1*B(i)), sig1, true); 
    }
  } 
  for(int i=0;i<n;i++){
    if(!isNA(y2(i))){
      nll += (-1.0) * dnorm(log(y2(i)), log(q2*B(i)), sig2, true); 
    }
  } 
  
  ADREPORT(r);
  ADREPORT(K);
  ADREPORT(MSY);
  ADREPORT(Bmsy);
  ADREPORT(Dep);
  ADREPORT(B);
  ADREPORT(q1);
  ADREPORT(q2);
  ADREPORT(sig1);
  ADREPORT(sig2);
  ADREPORT(tau);
  
  return nll;
}
