function [PSimbas,upper_CI, lower_CI] = jointband_simbas_Michelle_sdcorrect(MCMC_P,alpha)
  % function to calculate simbas
  
   % Inputs
   %    'MCMC_P' - a B by T matrix containing MCMC samples. 
   %        B - number of MCMC iterations (MCMCspecs.B)
   %        T - number of function samples (number of columns in Y)
      
   % Outputs
  
   % MCMC_P=X_BetaDataSpace1;
   [B,T] = size(MCMC_P); 
   sd_P=std(MCMC_P,1,1);
  
   index_sd=sd_P~=0;
   index_sdn=sd_P==0;
   MCMC_P(:,index_sdn)=[];
   mean_P = mean(MCMC_P);
   sd_P(:,index_sdn)=[];
   
   z_P = NaN(1,B);
   for j=1:B
        z_P(j) = max(abs((MCMC_P(j,:)-mean_P)./sd_P));
   end
   a_T=abs(mean_P./sd_P);
   T2=size(MCMC_P,2);
   PSimbas0=NaN(1,T2);
   for t=1:T2
    PSimbas0(t)=sum(a_T(t)<=z_P)./B;
   end
  
PSimbas=ones(1,T);
PSimbas(index_sd)=PSimbas0;


cb2 = quantile(z_P,1-alpha);
upper_CI_0 = mean_P + cb2*sd_P;
lower_CI_0= mean_P - cb2*sd_P;  
  
upper_CI = zeros(1,T);
lower_CI = zeros(1,T);
upper_CI(index_sd)=upper_CI_0;
lower_CI(index_sd)=lower_CI_0;

