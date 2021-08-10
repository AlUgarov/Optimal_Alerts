# Calculate the risk aversion coefficient in the CRRA fn based on blind protection choices 
program define riskaversion_crra, eclass
   version 14
 
   syntax varlist(max=1)
   
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
   
   end
   ereturn u
end
