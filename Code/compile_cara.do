 clear all
 cd C:\Tornado_warnings\Experiment\Alerts_Experiment\Code
 version 14
 mata:
   
   function cara(c,theta)
   {
     if (theta==0){
        u=c
     }
     else {
        u = -(1/theta)*exp(-theta*c)
     }
     return(u)
   }
   mata mosave cara(), replace
   
end
