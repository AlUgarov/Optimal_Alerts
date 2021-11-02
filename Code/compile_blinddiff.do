 clear all
 cd C:\Tornado_warnings\Experiment\Alerts_Experiment\Code
 version 14
 mata:
   crra(2,0.5)
   
   function blind_diff(theta, Y0,prot_cost, L, p)
   {
      if ((p<0)||(p>1)){
	    diff=-10000
	  }
	  else {
	    diff=crra(Y0-prot_cost,theta)-p*crra(Y0-L,theta)-(1-p)*crra(Y0,theta)
	  }
	  return(diff)
   
   }
   
   blind_diff(0.05, 30,5,20,0.5)
   mata mosave blind_diff(), replace
end
