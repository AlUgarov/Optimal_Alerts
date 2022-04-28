 clear all
 cd C:\Tornado_warnings\Experiment\Alerts_Experiment\Code
 version 14
 mata:
   
   function bayes_adj(p,p11, p10, alpha, beta)
   {
     pw=((p11^alpha)*(p^beta))/((p11^alpha)*(p^beta)+(p10^alpha)*((1-p)^beta))
     return(pw)
   }
   mata mosave bayes_adj(), replace
   
end
