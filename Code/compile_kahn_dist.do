 clear all
 cd C:\Tornado_warnings\Experiment\Alerts_Experiment\Code
 version 14
 mata:
   
   function kahn_dist(p,gamma)
   {
     pw=p^gamma/(p^gamma+(1-p)^gamma)^(1/gamma)
     return(pw)
   }
   mata mosave kahn_dist(), replace
   
end
