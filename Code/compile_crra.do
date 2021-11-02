 clear all
 cd C:\Tornado_warnings\Experiment\Alerts_Experiment\Code
 version 14
 mata:
   
   function crra(c,theta)
   {
     if (theta==1){
       u=log(c)
     }
     else {
        u = (c^(1-theta)-1)/(1-theta)
     }
     return(u)
   }
   mata mosave crra(), replace
   
end
