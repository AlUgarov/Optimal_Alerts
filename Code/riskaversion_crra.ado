# Calculate the risk aversion coefficient in the CRRA fn based on blind protection choices 
program define _griskaversion_crra, eclass
   version 14
 
   syntax varlist(max=1)
   
   syntax newvarname=/exp
   
   function myfunc(V) return(infoval_diff(30,5,20,V,0.5,1,1,0.01))
   mm_root(V=., &myfunc(), 0.001, 5, 0.0001, 1000)
   tempname z g
   gen `z' = (`exp' - `mean')/`stdev'
   gen `g' = gammap(0.604, exp(-1.896 * (`z' + 0.805)))
   gen `typlist' `varlist' = crra 
end 
   
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
   
   
   err=1.0
   min_error=10^-6
   theta_max
   theta_min
   while (err<min_error)&(iter<1000){
   err=crra()
   }
   //here the mata ends
   end
   ereturn u
end
