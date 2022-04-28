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
   function infoval_diff_probtransform(V, Y0,prot_cost, L, p, phintWW, phintBB, theta, gam)
   {
      phintBW=1-phintWW
	  phintWB=1-phintBB
	  if ((p<0)||(p>1)){
	    diff=-10000
	  }
	  else {
	    pblackhint=kahn_dist(p*phintBB+(1-p)*phintBW, gam)
		ptrue_neg=kahn_dist((1-p)*phintWW, gam)
		pfalse_neg=kahn_dist(p*phintWB, gam)
		pw=kahn_dist(p,gam)
	    ut_signal=pblackhint*crra(Y0-prot_cost-V,theta)+ptrue_neg*crra(Y0-V,theta)+pfalse_neg*crra(Y0-V-L,theta)
		ut_nosignal=max((crra(Y0-prot_cost,theta),pw*crra(Y0-L,theta)+(1-pw)*crra(Y0,theta)))
	    diff=ut_signal-ut_nosignal
	  }
	  return(diff)
   }
   
    
   mata mosave infoval_diff_probtransform(), replace
end
