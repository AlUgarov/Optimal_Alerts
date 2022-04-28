 clear all
 cd C:\Tornado_warnings\Experiment\Alerts_Experiment\Code
 version 14
 mata:
   //the function calculates the difference in expected utilities for buying a signal and not buying a signal
   //>0 - it is beneficial to buy a signal
   //Y0- endowment, prot_cost- protection cost, L - potential loss, 
   //V - price of the signal
   //p - prior prob of the bad state, phintWW - probability to receive white hint conditional on white/good state.
   //phintBB -prob to receive a black hint conditional on black/bad state, theta - relative risk aversion coefficient
   function infoval_diff_cara(V, Y0,prot_cost, L, p, phintWW, phintBB, theta)
   {
      phintBW=1-phintWW
	  phintWB=1-phintBB
	  if ((p<0)||(p>1)){
	    diff=-10000
	  }
	  else {
	    pblackhint=p*phintBB+(1-p)*phintBW
	    ut_signal=pblackhint*cara(Y0-prot_cost-V,theta)+(1-p)*phintWW*cara(Y0-V,theta)+p*phintWB*cara(Y0-V-L,theta)
		ut_nosignal=max((cara(Y0-prot_cost,theta),p*cara(Y0-L,theta)+(1-p)*cara(Y0,theta)))
	    diff=ut_signal-ut_nosignal
	  }
	  return(diff)
   }
   
   //function myfunc(V) return(infoval_diff(V, 30,5,20,0.5,1,1,0.01))
   //mm_root(V=., &myfunc(), 0.001, 5, 0.0001, 1000)
   //V
   //infoval_diff(30,5,20,0,0.5,1,1,0.01)
   //infoval_diff(30,5,20,1,0.5,1,1,0.01)
   //infoval_diff(30,5,20,2,0.5,1,1,0.01)
   //infoval_diff(30,5,20,5,0.5,1,1,0.01)
     
   mata mosave infoval_diff_cara(), replace
end
