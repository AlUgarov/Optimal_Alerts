 clear all
 cd C:\Tornado_warnings\Experiment\Alerts_Experiment\Code
 version 14
 mata:
   
   function hara(c,theta)
   {
     if (theta==1){
        u=c
     }
     else {
        u = -(1/theta)*exp(-theta*c)
     }
     return(u)
   }
   mata mosave hara(), replace
   
end
