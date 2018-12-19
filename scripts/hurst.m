addpath(genpath('./'));

%% David Study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load david.mat;
c = @cmu.colors;

%% Human A %%%%%%%%%%% 
[M_A,N_A] = size(stoolA);
data_A = stoolA ./ repmat(sum(stoolA),M_A,1);
days_A = stoolA_days;

%Find most prevalent species;
prevs_A = sum(data_A>0,2);
most_prev_A = N_A/2;
most_prev_inds_A = find(prevs_A >= most_prev_A);

%Find most abundant species;
abunds_A = mean(data_A,2);

%% Do analysis on abundant and prevalent OTUs
most_abd_A = -3;
abd_inds_A = find(log10(abunds_A) > most_abd_A);
abd_inds_A = intersect(abd_inds_A,most_prev_inds_A);

%% Update data
data_A = data_A(abd_inds_A,:);
[M_A,N_A] = size(data_A);

%Find most abundant species;
abunds_A = mean(data_A,2);

%Update tax info
tax_A = tax(abd_inds_A);
ptax_A = ptax(abd_inds_A);
ctax_A = ctax(abd_inds_A);
otax_A = otax(abd_inds_A);
ftax_A = ftax(abd_inds_A);
gtax_A = gtax(abd_inds_A);


%% Hurst exponent (averaged over all OTUs)  
max_lag_A = 100;
lags_A = 1:max_lag_A; 
msd_A = [];

for i = 1:length(lags_A)
   lag = lags_A(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_A)
       for k = (j+1):length(days_A)
           d_i = days_A(k)-days_A(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end           
   tX0_i = data_A(:,inds0);
   tX1_i = data_A(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_A(i) = var(tdx_i);
end
[Gamma_A,s_A] = polyfit(log10(lags_A),log10(msd_A),1); 
ste_Gamma_A = sqrt(diag(inv(s_A.R)*inv(s_A.R')).*s_A.normr.^2./s_A.df);


%% Hurst exponent of individual OTUs
 msd_k_A = [];
for i = 1:length(lags_A)
   lag = lags_A(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_A)
       for k = (j+1):length(days_A)
           d_i = days_A(k)-days_A(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end 
   
   tX0_i = data_A(:,inds0);
   tX1_i = data_A(:,inds1);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdxs_i(find(abs(tdxs_i) == Inf)) = NaN; %Ignore counts of zero
   msd_k_A(:,i) = nanvar(tdxs_i')';
end

gammas_A = [];
for i = 1:size(msd_k_A,1)
   h = polyfit(log10(lags_A),log10(msd_k_A(i,:)),1);
   gammas_A(i) = h(1);
end



%% Human B %%%%%%%%%%% 
[M_B,N_B] = size(stoolB);
data_B = stoolB ./ repmat(sum(stoolB),M_B,1);
days_B = stoolB_days;

%Find most prevalent species;
prevs_B = sum(data_B>0,2);
most_prev_B = N_B/2;
most_prev_inds_B = find(prevs_B >= most_prev_B);

%Find most abundant species;
abunds_B = mean(data_B,2);


%% Do analysis on abundant and prevalent OTUs
most_abd_B = -3;
abd_inds_B = find(log10(abunds_B) > most_abd_B);
abd_inds_B = intersect(abd_inds_B,most_prev_inds_B);

%% Update data
data_B = data_B(abd_inds_B,:);
[M_B,N_B] = size(data_B);

%Find most abundant species;
abunds_B = mean(data_B,2);

%Update tax info
tax_B = tax(abd_inds_B);
ptax_B = ptax(abd_inds_B);
ctax_B = ctax(abd_inds_B);
otax_B = otax(abd_inds_B);
ftax_B = ftax(abd_inds_B);
gtax_B = gtax(abd_inds_B);



%% Hurst exponent (averaged over all OTUs)  
max_lag_B = 100;
lags_B = 1:max_lag_B; 
msd_B = [];

for i = 1:length(lags_B)
   lag = lags_B(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_B)
       for k = (j+1):length(days_B)
           d_i = days_B(k)-days_B(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end           
   tX0_i = data_B(:,inds0);
   tX1_i = data_B(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_B(i) = var(tdx_i);
end
[Gamma_B,s_B] = polyfit(log10(lags_B),log10(msd_B),1); 
ste_Gamma_B = sqrt(diag(inv(s_B.R)*inv(s_B.R')).*s_B.normr.^2./s_B.df);


%% Hurst exponent of individual OTUs
 msd_k_B = [];
for i = 1:length(lags_B)
   lag = lags_B(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_B)
       for k = (j+1):length(days_B)
           d_i = days_B(k)-days_B(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end 
   
   tX0_i = data_B(:,inds0);
   tX1_i = data_B(:,inds1);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdxs_i(find(abs(tdxs_i) == Inf)) = NaN; %Ignore counts of zero
   msd_k_B(:,i) = nanvar(tdxs_i')';
end

gammas_B = [];
for i = 1:size(msd_k_B,1)
   h = polyfit(log10(lags_B),log10(msd_k_B(i,:)),1);
   gammas_B(i) = h(1);
end

























%% Caporaso study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load caporaso.mat;



%% Human M3 %%%%%%%%%%% 
[M_m3,N_m3] = size(m3);
data_m3 = m3 ./ repmat(sum(m3),M_m3,1);
days_m3 = m3_days;

%Find most prevalent species;
prevs_m3 = sum(data_m3>0,2);
most_prev_m3 = N_m3/2;
most_prev_inds_m3 = find(prevs_m3 >= most_prev_m3);

%Find most abundant species;
abunds_m3 = mean(data_m3,2);


%% Do analysis on only abundant and prevalent OTUs
most_abd_m3 = -3;
abd_inds_m3 = find(log10(abunds_m3) > most_abd_m3);
abd_inds_m3 = intersect(abd_inds_m3,most_prev_inds_m3);

%% Update data
data_m3 = data_m3(abd_inds_m3,:);
[M_m3,N_m3] = size(data_m3);

%Find most abundant species;
abunds_m3 = mean(data_m3,2);

%Update tax info
tax_m3 = tax(abd_inds_m3);
ptax_m3 = ptax(abd_inds_m3);
ctax_m3 = ctax(abd_inds_m3);
otax_m3 = otax(abd_inds_m3);
ftax_m3 = ftax(abd_inds_m3);
gtax_m3 = gtax(abd_inds_m3);




%% Hurst exponent (averaged over all OTUs)  
max_lag_m3 = 100;
lags_m3 = 1:max_lag_m3; 
msd_m3 = [];

for i = 1:length(lags_m3)
   lag = lags_m3(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_m3)
       for k = (j+1):length(days_m3)
           d_i = days_m3(k)-days_m3(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end           
   tX0_i = data_m3(:,inds0);
   tX1_i = data_m3(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_m3(i) = var(tdx_i);
end
[Gamma_m3,s_m3] = polyfit(log10(lags_m3),log10(msd_m3),1); 
ste_Gamma_m3 = sqrt(diag(inv(s_m3.R)*inv(s_m3.R')).*s_m3.normr.^2./s_m3.df);


%% Hurst exponent of individual OTUs
 msd_k_m3 = [];
for i = 1:length(lags_m3)
   lag = lags_m3(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_m3)
       for k = (j+1):length(days_m3)
           d_i = days_m3(k)-days_m3(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end 
   
   tX0_i = data_m3(:,inds0);
   tX1_i = data_m3(:,inds1);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdxs_i(find(abs(tdxs_i) == Inf)) = NaN; %Ignore counts of zero
   msd_k_m3(:,i) = nanvar(tdxs_i')';
end

gammas_m3 = [];
for i = 1:size(msd_k_m3,1)
   h = polyfit(log10(lags_m3),log10(msd_k_m3(i,:)),1);
   gammas_m3(i) = h(1);
end









%% Human F4 %%%%%%%%%%% 
[M_f4,N_f4] = size(f4);
data_f4 = f4 ./ repmat(sum(f4),M_f4,1);
days_f4 = f4_days;

%Find most prevalent species;
prevs_f4 = sum(data_f4>0,2);
most_prev_f4 = N_f4/2;
most_prev_inds_f4 = find(prevs_f4 >= most_prev_f4);

%Find most abundant species;
abunds_f4 = mean(data_f4,2);


%% Do analysis on only abundant and prevalent OTUS (cutoff made prior to 16S normalization)
most_abd_f4 = -3;
abd_inds_f4 = find(log10(abunds_f4) > most_abd_f4);
abd_inds_f4 = intersect(abd_inds_f4,most_prev_inds_f4);

%% Update data
data_f4 = data_f4(abd_inds_f4,:);
[M_f4,N_f4] = size(data_f4);

%Find most abundant species;
abunds_f4 = mean(data_f4,2);

%Update tax info
tax_f4 = tax(abd_inds_f4);
ptax_f4 = ptax(abd_inds_f4);
ctax_f4 = ctax(abd_inds_f4);
otax_f4 = otax(abd_inds_f4);
ftax_f4 = ftax(abd_inds_f4);
gtax_f4 = gtax(abd_inds_f4);



%% Hurst exponent (averaged over all OTUs)  
max_lag_f4 = 100;
lags_f4 = 1:max_lag_f4; 
msd_f4 = [];

for i = 1:length(lags_f4)
   lag = lags_f4(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_f4)
       for k = (j+1):length(days_f4)
           d_i = days_f4(k)-days_f4(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end           
   tX0_i = data_f4(:,inds0);
   tX1_i = data_f4(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_f4(i) = var(tdx_i);
end
[Gamma_f4,s_f4] = polyfit(log10(lags_f4),log10(msd_f4),1); 
ste_Gamma_f4 = sqrt(diag(inv(s_f4.R)*inv(s_f4.R')).*s_f4.normr.^2./s_f4.df);


%% Hurst exponent of individual OTUs
 msd_k_f4 = [];
for i = 1:length(lags_f4)
   lag = lags_f4(i);
   
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_f4)
       for k = (j+1):length(days_f4)
           d_i = days_f4(k)-days_f4(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end 
   
   tX0_i = data_f4(:,inds0);
   tX1_i = data_f4(:,inds1);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdxs_i(find(abs(tdxs_i) == Inf)) = NaN; %Ignore counts of zero
   msd_k_f4(:,i) = nanvar(tdxs_i')';
end

gammas_f4 = [];
for i = 1:size(msd_k_f4,1)
   h = polyfit(log10(lags_f4),log10(msd_k_f4(i,:)),1);
   gammas_f4(i) = h(1);
end




















%% Carmody study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load carmody.mat;


%% Set up mice data %%%%%%%%%%%%%%%%%%%%% 
[M,N] = size(lf1);
data_hf1 = hf1 ./ repmat(sum(hf1),M,1);
data_hf2 = hf2 ./ repmat(sum(hf2),M,1);
data_hf3 = hf3 ./ repmat(sum(hf3),M,1);
data_lf1 = lf1 ./ repmat(sum(lf1),M,1);
data_lf2 = lf2 ./ repmat(sum(lf2),M,1);
data_lf3 = lf3 ./ repmat(sum(lf3),M,1);
days_hf1 = hf1_days;
days_hf2 = hf2_days;
days_hf3 = hf3_days;
days_lf1 = lf1_days;
days_lf2 = lf2_days;
days_lf3 = lf3_days;

%Find most prevalent species;
[M_hf1,N_hf1] = size(data_hf1);
prevs_hf1 = sum(data_hf1>0,2);
most_prev_hf1 = N_hf1/2;
most_prev_inds_hf1 = find(prevs_hf1 >= most_prev_hf1);

[M_hf2,N_hf2] = size(data_hf2);
prevs_hf2 = sum(data_hf2>0,2);
most_prev_hf2 = N_hf2/2;
most_prev_inds_hf2 = find(prevs_hf2 >= most_prev_hf2);

[M_hf3,N_hf3] = size(data_hf3);
prevs_hf3 = sum(data_hf3>0,2);
most_prev_hf3 = N_hf3/2;
most_prev_inds_hf3 = find(prevs_hf3 >= most_prev_hf3);

[M_lf1,N_lf1] = size(data_lf1);
prevs_lf1 = sum(data_lf1>0,2);
most_prev_lf1 = N_lf1/2;
most_prev_inds_lf1 = find(prevs_lf1 >= most_prev_lf1);

[M_lf2,N_lf2] = size(data_lf2);
prevs_lf2 = sum(data_lf2>0,2);
most_prev_lf2 = N_lf2/2;
most_prev_inds_lf2 = find(prevs_lf2 >= most_prev_lf2);

[M_lf3,N_lf3] = size(data_lf3);
prevs_lf3 = sum(data_lf3>0,2);
most_prev_lf3 = N_lf3/2;
most_prev_inds_lf3 = find(prevs_lf3 >= most_prev_lf3);


%Find most abundant species;
abunds_hf1 = mean(data_hf1,2);
abunds_hf2 = mean(data_hf2,2);
abunds_hf3 = mean(data_hf3,2);
abunds_lf1 = mean(data_lf1,2);
abunds_lf2 = mean(data_lf2,2);
abunds_lf3 = mean(data_lf3,2);


%% Do analysis on abundant and prevalent OTUs
most_abd_hf1 = -3;
abd_inds_hf1 = find(log10(abunds_hf1) > most_abd_hf1);
abd_inds_hf1 = intersect(abd_inds_hf1,most_prev_inds_hf1);
most_abd_hf2 = -3;
abd_inds_hf2 = find(log10(abunds_hf2) > most_abd_hf2);
abd_inds_hf2 = intersect(abd_inds_hf2,most_prev_inds_hf2);
most_abd_hf3 = -3;
abd_inds_hf3 = find(log10(abunds_hf3) > most_abd_hf3);
abd_inds_hf3 = intersect(abd_inds_hf3,most_prev_inds_hf3);
most_abd_lf1 = -3;
abd_inds_lf1 = find(log10(abunds_lf1) > most_abd_lf1);
abd_inds_lf1 = intersect(abd_inds_lf1,most_prev_inds_lf1);
most_abd_lf2 = -3;
abd_inds_lf2 = find(log10(abunds_lf2) > most_abd_lf2);
abd_inds_lf2 = intersect(abd_inds_lf2,most_prev_inds_lf2);
most_abd_lf3 = -3;
abd_inds_lf3 = find(log10(abunds_lf3) > most_abd_lf3);
abd_inds_lf3 = intersect(abd_inds_lf3,most_prev_inds_lf3);



%% Update data
data_hf1 = data_hf1(abd_inds_hf1,:);
data_hf2 = data_hf2(abd_inds_hf2,:);
data_hf3 = data_hf3(abd_inds_hf3,:);
data_lf1 = data_lf1(abd_inds_lf1,:);
data_lf2 = data_lf2(abd_inds_lf2,:);
data_lf3 = data_lf3(abd_inds_lf3,:);
[M,N] = size(data_hf1);

%Find Means;
abunds_hf1 = mean(data_hf1,2);
abunds_hf2 = mean(data_hf2,2);
abunds_hf3 = mean(data_hf3,2);
abunds_lf1 = mean(data_lf1,2);
abunds_lf2 = mean(data_lf2,2);
abunds_lf3 = mean(data_lf3,2);

%Find variances
vars_hf1 = var(data_hf1')';
vars_hf2 = var(data_hf2')';
vars_hf3 = var(data_hf3')';
vars_lf1 = var(data_lf1')';
vars_lf2 = var(data_lf2')';
vars_lf3 = var(data_lf3')';


%Update tax info
tax_hf1 = tax(abd_inds_hf1);
tax_hf2 = tax(abd_inds_hf2);
tax_hf3 = tax(abd_inds_hf3);
tax_lf1 = tax(abd_inds_lf1);
tax_lf2 = tax(abd_inds_lf2);
tax_lf3 = tax(abd_inds_lf3);

ptax_hf1 = ptax(abd_inds_hf1);
ptax_hf2 = ptax(abd_inds_hf2);
ptax_hf3 = ptax(abd_inds_hf3);
ptax_lf1 = ptax(abd_inds_lf1);
ptax_lf2 = ptax(abd_inds_lf2);
ptax_lf3 = ptax(abd_inds_lf3);

ctax_hf1 = ctax(abd_inds_hf1);
ctax_hf2 = ctax(abd_inds_hf2);
ctax_hf3 = ctax(abd_inds_hf3);
ctax_lf1 = ctax(abd_inds_lf1);
ctax_lf2 = ctax(abd_inds_lf2);
ctax_lf3 = ctax(abd_inds_lf3);

otax_hf1 = otax(abd_inds_hf1);
otax_hf2 = otax(abd_inds_hf2);
otax_hf3 = otax(abd_inds_hf3);
otax_lf1 = otax(abd_inds_lf1);
otax_lf2 = otax(abd_inds_lf2);
otax_lf3 = otax(abd_inds_lf3);

ftax_hf1 = ftax(abd_inds_hf1);
ftax_hf2 = ftax(abd_inds_hf2);
ftax_hf3 = ftax(abd_inds_hf3);
ftax_lf1 = ftax(abd_inds_lf1);
ftax_lf2 = ftax(abd_inds_lf2);
ftax_lf3 = ftax(abd_inds_lf3);

gtax_hf1 = gtax(abd_inds_hf1);
gtax_hf2 = gtax(abd_inds_hf2);
gtax_hf3 = gtax(abd_inds_hf3);
gtax_lf1 = gtax(abd_inds_lf1);
gtax_lf2 = gtax(abd_inds_lf2);
gtax_lf3 = gtax(abd_inds_lf3);

otus_hf1 = otu_ids(abd_inds_hf1);
otus_hf2 = otu_ids(abd_inds_hf2);
otus_hf3 = otu_ids(abd_inds_hf3);
otus_lf1 = otu_ids(abd_inds_lf1);
otus_lf2 = otu_ids(abd_inds_lf2);
otus_lf3 = otu_ids(abd_inds_lf3);


%% Hurst 


%% HFHS 1
max_lag_hf1 = 15;
lags_hf1 = 1:max_lag_hf1; 


% Hurst (Averaged over all OTUs)  
msd_hf1 = [];

for i = 1:length(lags_hf1)
   lag = lags_hf1(i);
   
   %Loop over all pairs of days
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_hf1)
       for k = (j+1):length(days_hf1)
           d_i = days_hf1(k)-days_hf1(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end  
   
   tX0_i = data_hf1(:,inds0);
   tX1_i = data_hf1(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_hf1(i) = var(tdx_i);
end
[Gamma_hf1,s_hf1] = polyfit(log10(lags_hf1),log10(msd_hf1),1); 
ste_Gamma_hf1 = sqrt(diag(inv(s_hf1.R)*inv(s_hf1.R')).*s_hf1.normr.^2./s_hf1.df);



%% HFHS 2
max_lag_hf2 = 15;
lags_hf2 = 1:max_lag_hf2; 

% Hurst (Averaged over all OTUs)  
msd_hf2 = [];

for i = 1:length(lags_hf2)
   lag = lags_hf2(i);
   
   %Loop over all pairs of days
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_hf2)
       for k = (j+1):length(days_hf2)
           d_i = days_hf2(k)-days_hf2(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end  
   
   tX0_i = data_hf2(:,inds0);
   tX1_i = data_hf2(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_hf2(i) = var(tdx_i);
end
[Gamma_hf2,s_hf2] = polyfit(log10(lags_hf2),log10(msd_hf2),1); 
ste_Gamma_hf2 = sqrt(diag(inv(s_hf2.R)*inv(s_hf2.R')).*s_hf2.normr.^2./s_hf2.df);


%% HFHS 3
max_lag_hf3 = 15;
lags_hf3 = 1:max_lag_hf3; 


% Hurst (Averaged over all OTUs)  
msd_hf3 = [];

for i = 1:length(lags_hf3)
   lag = lags_hf3(i);
      
   %Loop over all pairs of days
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_hf3)
       for k = (j+1):length(days_hf3)
           d_i = days_hf3(k)-days_hf3(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end  
   
   tX0_i = data_hf3(:,inds0);
   tX1_i = data_hf3(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_hf3(i) = var(tdx_i);
end
[Gamma_hf3,s_hf3] = polyfit(log10(lags_hf3),log10(msd_hf3),1); 
ste_Gamma_hf3 = sqrt(diag(inv(s_hf3.R)*inv(s_hf3.R')).*s_hf3.normr.^2./s_hf3.df);



%% LFPP 1
max_lag_lf1 = 15;
lags_lf1 = 1:max_lag_lf1; 


% Hurst (Averaged over all OTUs)  
msd_lf1 = [];

for i = 1:length(lags_lf1)
   lag = lags_lf1(i);
   
   %Loop over all pairs of days
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_lf1)
       for k = (j+1):length(days_lf1)
           d_i = days_lf1(k)-days_lf1(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end  
   
   tX0_i = data_lf1(:,inds0);
   tX1_i = data_lf1(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_lf1(i) = var(tdx_i);
end
[Gamma_lf1,s_lf1] = polyfit(log10(lags_lf1),log10(msd_lf1),1); 
ste_Gamma_lf1 = sqrt(diag(inv(s_lf1.R)*inv(s_lf1.R')).*s_lf1.normr.^2./s_lf1.df);



%% LFPP 2
max_lag_lf2 = 15;
lags_lf2 = 1:max_lag_lf2; 


% Hurst (Averaged over all OTUs)  
msd_lf2 = [];

for i = 1:length(lags_lf2)
   lag = lags_lf2(i);
   
   %Loop over all pairs of days
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_lf2)
       for k = (j+1):length(days_lf2)
           d_i = days_lf2(k)-days_lf2(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end  
   
   tX0_i = data_lf2(:,inds0);
   tX1_i = data_lf2(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_lf2(i) = var(tdx_i);
end
[Gamma_lf2,s_lf2] = polyfit(log10(lags_lf2),log10(msd_lf2),1); 
ste_Gamma_lf2 = sqrt(diag(inv(s_lf2.R)*inv(s_lf2.R')).*s_lf2.normr.^2./s_lf2.df);



%% LFPP 3
max_lag_lf3 = 15;
lags_lf3 = 1:max_lag_lf3; 


% Hurst (Averaged over all OTUs)  
msd_lf3 = [];

for i = 1:length(lags_lf3)
   lag = lags_lf3(i);
      
   %Loop over all pairs of days
   inds0 = [];
   inds1 = [];
   for j = 1:length(days_lf3)
       for k = (j+1):length(days_lf3)
           d_i = days_lf3(k)-days_lf3(j);
           if d_i == lag
               inds0 = [inds0 j];
               inds1 = [inds1 k];
           end
       end       
   end  
   
   tX0_i = data_lf3(:,inds0);
   tX1_i = data_lf3(:,inds1);
   tX0_inds_i = find(tX0_i);
   tX1_inds_i = find(tX1_i);
   tdXs_i = tX1_i ./ tX0_i;
   tdxs_i = log10(tdXs_i); 
   tdx_inds_i = intersect(tX0_inds_i,tX1_inds_i);
   tdx_i = tdxs_i(tdx_inds_i);
   msd_lf3(i) = var(tdx_i);
end
[Gamma_lf3,s_lf3] = polyfit(log10(lags_lf3),log10(msd_lf3),1); 
ste_Gamma_lf3 = sqrt(diag(inv(s_lf3.R)*inv(s_lf3.R')).*s_lf3.normr.^2./s_lf3.df);




%% Combine diets
lags_hf = lags_hf1;
lags_lf = lags_lf1;

msd_hf = mean([msd_hf1; msd_hf2; msd_hf3],1); %Averaged over three mice
[Gamma_hf,s_hf] = polyfit(log10(lags_hf),log10(msd_hf),1); 

msd_lf = mean([msd_lf1; msd_lf2; msd_lf3],1); %Averaged over three mie
[Gamma_lf,s_lf] = polyfit(log10(lags_lf),log10(msd_lf),1); 

Gammas_lf = [Gamma_lf1(1) Gamma_lf2(1) Gamma_lf3(1)];
Gammas_hf = [Gamma_hf1(1) Gamma_hf2(1) Gamma_hf3(1)];





























%% Plotting
f1 = figure;


%Human A
subplot(1,3,1);

x_A = -.2:.1:2.4;
p1a = plot(x_A,polyval(Gamma_A,x_A),'--');
p1a.Color = c('matlab maroon');
p1a.LineWidth = 1.5;
hold on
p1b = plot(log10(lags_A),log10(msd_A),'.');
p1b.Color = c('black');
p1b.MarkerSize = 22;
hold off
set(gca,'xlim',[-.2,2.4])
set(gca,'XTick',[0 0.5 1 1.5 2]);
Xlabel = {};
Xlabel{1} = '10^{0}';
Xlabel{2} = '';
Xlabel{3} = '10^{1}';
Xlabel{4} = '';
Xlabel{5} = '10^{2}';
set(gca,'XTickLabel',Xlabel);
set(gca,'ylim',[-1,.1])
set(gca,'YTick',[-1:.25:0]);
set(gca,'YTickLabel',{'10^{-1}','','10^{-0.5}','','10^{0}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Time, \Delta t (days)');
ylabel('Mean-squared displacement');
title('Human A','FontSize',16,'FontWeight','Normal');
box off


%Human M3
subplot(1,3,2);

x_m3 = -.2:.1:2.4;
p2a = plot(x_m3,polyval(Gamma_m3,x_m3),'--');
p2a.Color = c('matlab blue');
p2a.LineWidth = 1.5;
hold on
p2b = plot(log10(lags_m3),log10(msd_m3),'.');
p2b.Color = c('black');
p2b.MarkerSize = 22;
hold off
set(gca,'xlim',[-.2,2.4])
set(gca,'XTick',[0 0.5 1 1.5 2]);
Xlabel = {};
Xlabel{1} = '10^{0}';
Xlabel{2} = '';
Xlabel{3} = '10^{1}';
Xlabel{4} = '';
Xlabel{5} = '10^{2}';
set(gca,'XTickLabel',Xlabel);
set(gca,'ylim',[-1,.1])
set(gca,'YTick',[-1:.25:0]);
set(gca,'YTickLabel',{'10^{-1}','','10^{-0.5}','','10^{0}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Time, \Delta t (days)');
title('Human M3','FontSize',16,'FontWeight','Normal');
box off

%Mice LF
subplot(1,3,3);

x_lf = -.2:.1:2.2;
p3a = plot(x_lf,polyval(Gamma_lf,x_lf),'--');
p3a.Color = c('matlab green');
p3a.LineWidth = 1.5;
hold on
p3b = plot(log10(lags_lf),log10(msd_lf),'.');
p3b.Color = c('black');
p3b.MarkerSize = 22;
hold off
set(gca,'xlim',[-.2,1.5])
set(gca,'XTick',[0 0.5 1 1.5]);
Xlabel = {};
Xlabel{1} = '10^{0}';
Xlabel{2} = '';
Xlabel{3} = '10^{1}';
Xlabel{4} = '';
set(gca,'XTickLabel',Xlabel);
set(gca,'ylim',[-1,.1])
set(gca,'YTick',[-1:.25:0]);
set(gca,'YTickLabel',{'10^{-1}','','10^{-0.5}','','10^{0}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Time, \Delta t (days)');
title('Mice (LFPP Diet)','FontSize',16,'FontWeight','Normal');
box off




