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


%% Daily growth rates

%Get samples that are 1 day apart
inds0_A = [];
inds1_A = [];
for j = 1:N_A
    for k = (j+1):N_A
        d_i = days_A(k)-days_A(j);
        if d_i == 1
            inds0_A = [inds0_A j];
            inds1_A = [inds1_A k];
        end
    end
end
   
%Growth rates in log10 base
tX0_A = data_A(:,inds0_A);
tX1_A = data_A(:,inds1_A);
tX0_inds_A = find(tX0_A);
tX1_inds_A = find(tX1_A);
tdXs_A = tX1_A ./ tX0_A;
tdxs_A = log10(tdXs_A);
tdx_inds_A = intersect(tX0_inds_A,tX1_inds_A);
tdx_A = tdxs_A(tdx_inds_A);
txo_A = log10(tX0_A(tdx_inds_A));
tx1_A = log10(tX1_A(tdx_inds_A));

%MLE of laplace exponent (in the natural log space)
params_A = fit_laplace(log(10.^(tdx_A)));
b_lp_A = params_A.b;



%% Daily growth rate variability versus mean daily abundance
pc1_A = (txo_A + tx1_A)/2;
pc2_A = (tx1_A - txo_A);
bins_xm_A = [-3:.4:-.6 0];
stds_A = [];
for i = 1:length(bins_xm_A)-1
    inds_i = find(pc1_A > bins_xm_A(i) & pc1_A < bins_xm_A(i+1));
    pc2_i = pc2_A(inds_i); %Growth rates
    stds_A(i) = std(pc2_i);
end
[r_A,s_A] = polyfit(bins_xm_A(1:end-1),stds_A,1);
ste_r_A = sqrt(diag(inv(s_A.R)*inv(s_A.R')).*s_A.normr.^2./s_A.df);



%% Variability of Laplace parameter within human A
W_A = floor(N_A / 6); %window size
bs_lp_A = [];
rs_A = [];

for i = 1:floor(N_A/W_A)
    
   %Calculate growth rates within each window
   
   %Take 1/6 of samples
   inds = (i-1)*W_A + (1:W_A); 
   data_W = data_A(:,inds);
   days_i = days_A(inds);
   
   inds0_i = [];
   inds1_i = [];
   for j = 1:W_A
       for k = (j+1):W_A
           d_i = days_i(k)-days_i(j);
           if d_i == 1
               inds0_i = [inds0_i j];
               inds1_i = [inds1_i k];
           end
       end
   end
   
   tX0_W = data_W(:,inds0_i);
   tX1_W = data_W(:,inds1_i);
   tX0_inds_W = find(tX0_W);
   tX1_inds_W = find(tX1_W);
   tdXs_W = tX1_W ./ tX0_W;
   tdxs_W = log10(tdXs_W);
   tdx_inds_W = intersect(tX0_inds_W,tX1_inds_W);
   tdx_W = tdxs_W(tdx_inds_W);
   txo_W = log10(tX0_W(tdx_inds_W));
   tx1_W = log10(tX1_W(tdx_inds_W));
   
   %Parameter fitting
   params_W = fit_laplace(log(10.^(tdx_W))); %fit in natural log space
   b_lp_W = params_W.b;
   bs_lp_A(i) = b_lp_W;
   
   %Growth rate variability versus mean daily abundance
   pc1_W = (txo_W + tx1_W)/2;
   pc2_W = (tx1_W - txo_W);
   bins_xm_W = [-3:.4:-.6 0];
   stds_W = [];
   for j = 1:length(bins_xm_W)-1
        inds_j = find(pc1_W > bins_xm_W(j) & pc1_W < bins_xm_W(j+1));
        pc2_j = pc2_W(inds_j); 
        stds_W(j) = std(pc2_j);
   end
   r_inds = find(~isnan(stds_W));
   lb_bins_xm_W = bins_xm_W(1:end-1);
   [r_W,s_w] = polyfit(lb_bins_xm_W(r_inds),stds_W(r_inds),1);
   rs_A(i) = r_W(1);
   
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


%% Daily growth rates

%Get samples that are 1 day apart
inds0_B = [];
inds1_B = [];
for j = 1:N_B
    for k = (j+1):N_B
        d_i = days_B(k)-days_B(j);
        if d_i == 1
            inds0_B = [inds0_B j];
            inds1_B = [inds1_B k];
        end
    end
end
   
%Growth rates in log10 base
tX0_B = data_B(:,inds0_B);
tX1_B = data_B(:,inds1_B);
tX0_inds_B = find(tX0_B);
tX1_inds_B = find(tX1_B);
tdXs_B = tX1_B ./ tX0_B;
tdxs_B = log10(tdXs_B);
tdx_inds_B = intersect(tX0_inds_B,tX1_inds_B);
tdx_B = tdxs_B(tdx_inds_B);
txo_B = log10(tX0_B(tdx_inds_B));
tx1_B = log10(tX1_B(tdx_inds_B));

%MLE of laplace exponent (in the natural log space)
params_B = fit_laplace(log(10.^(tdx_B)));
b_lp_B = params_B.b;



%% Daily growth rate variability versus mean daily abundance
pc1_B = (txo_B + tx1_B)/2;
pc2_B = (tx1_B - txo_B);
bins_xm_B = [-3:.4:-.6 0];
stds_B = [];
for i = 1:length(bins_xm_B)-1
    inds_i = find(pc1_B > bins_xm_B(i) & pc1_B < bins_xm_B(i+1));
    pc2_i = pc2_B(inds_i); %Growth rates
    stds_B(i) = std(pc2_i);
end
[r_B,s_B] = polyfit(bins_xm_B(1:end-1),stds_B,1);
ste_r_B = sqrt(diag(inv(s_B.R)*inv(s_B.R')).*s_B.normr.^2./s_B.df);



%% Variability of Laplace parameter within human B
W_B = floor(N_B / 6); %window size
bs_lp_B = [];
rs_B = [];

for i = 1:floor(N_B/W_B)
    
   %Calculate growth rates within each window
   
   %Take 1/6 of samples
   inds = (i-1)*W_B + (1:W_B); 
   data_W = data_B(:,inds);
   days_i = days_B(inds);
   
   inds0_i = [];
   inds1_i = [];
   for j = 1:W_B
       for k = (j+1):W_B
           d_i = days_i(k)-days_i(j);
           if d_i == 1
               inds0_i = [inds0_i j];
               inds1_i = [inds1_i k];
           end
       end
   end
   
   tX0_W = data_W(:,inds0_i);
   tX1_W = data_W(:,inds1_i);
   tX0_inds_W = find(tX0_W);
   tX1_inds_W = find(tX1_W);
   tdXs_W = tX1_W ./ tX0_W;
   tdxs_W = log10(tdXs_W);
   tdx_inds_W = intersect(tX0_inds_W,tX1_inds_W);
   tdx_W = tdxs_W(tdx_inds_W);
   txo_W = log10(tX0_W(tdx_inds_W));
   tx1_W = log10(tX1_W(tdx_inds_W));
   
   %Parameter fitting
   params_W = fit_laplace(log(10.^(tdx_W))); %fit based on natural log
   b_lp_W = params_W.b;
   bs_lp_B(i) = b_lp_W;
   
   %Growth rate variability versus abundance
   pc1_W = (txo_W + tx1_W)/2;
   pc2_W = (tx1_W - txo_W);
   bins_xm_W = [-3:.4:-.6 0];
   stds_W = [];
   for j = 1:length(bins_xm_W)-1
        inds_j = find(pc1_W > bins_xm_W(j) & pc1_W < bins_xm_W(j+1));
        pc2_j = pc2_W(inds_j); %Growth rates
        stds_W(j) = std(pc2_j);
   end
   r_inds = find(~isnan(stds_W));
   lb_bins_xm_W = bins_xm_W(1:end-1);
   [r_W,s_w] = polyfit(lb_bins_xm_W(r_inds),stds_W(r_inds),1);
   rs_B(i) = r_W(1);
   
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


%% Daily growth rates

%Get samples that are 1 day apart
inds0_m3 = [];
inds1_m3 = [];
for j = 1:N_m3
    for k = (j+1):N_m3
        d_i = days_m3(k)-days_m3(j);
        if d_i == 1
            inds0_m3 = [inds0_m3 j];
            inds1_m3 = [inds1_m3 k];
        end
    end
end
   
%Growth rates in log10 base
tX0_m3 = data_m3(:,inds0_m3);
tX1_m3 = data_m3(:,inds1_m3);
tX0_inds_m3 = find(tX0_m3);
tX1_inds_m3 = find(tX1_m3);
tdXs_m3 = tX1_m3 ./ tX0_m3;
tdxs_m3 = log10(tdXs_m3);
tdx_inds_m3 = intersect(tX0_inds_m3,tX1_inds_m3);
tdx_m3 = tdxs_m3(tdx_inds_m3);
txo_m3 = log10(tX0_m3(tdx_inds_m3));
tx1_m3 = log10(tX1_m3(tdx_inds_m3));

%MLE of laplace exponent (in the natural log space)
params_m3 = fit_laplace(log(10.^(tdx_m3)));
b_lp_m3 = params_m3.b;



%% Daily growth rate variability versus mean daily abundance
pc1_m3 = (txo_m3 + tx1_m3)/2;
pc2_m3 = (tx1_m3 - txo_m3);
bins_xm_m3 = [-3:.4:-.6 0];
stds_m3 = [];
for i = 1:length(bins_xm_m3)-1
    inds_i = find(pc1_m3 > bins_xm_m3(i) & pc1_m3 < bins_xm_m3(i+1));
    pc2_i = pc2_m3(inds_i); %Growth rates
    stds_m3(i) = std(pc2_i);
end
[r_m3,s_m3] = polyfit(bins_xm_m3(1:end-1),stds_m3,1);
ste_r_m3 = sqrt(diag(inv(s_m3.R)*inv(s_m3.R')).*s_m3.normr.^2./s_m3.df);



%% Variability of Laplace parameter within human M3
W_m3 = floor(N_m3 / 6); %window size
bs_lp_m3 = [];
rs_m3 = [];

for i = 1:floor(N_m3/W_m3)
    
   %Calculate growth rates within each window
   
   %Take 1/6 of samples
   inds = (i-1)*W_m3 + (1:W_m3); 
   data_W = data_m3(:,inds);
   days_i = days_m3(inds);
   
   inds0_i = [];
   inds1_i = [];
   for j = 1:W_m3
       for k = (j+1):W_m3
           d_i = days_i(k)-days_i(j);
           if d_i == 1
               inds0_i = [inds0_i j];
               inds1_i = [inds1_i k];
           end
       end
   end
   
   tX0_W = data_W(:,inds0_i);
   tX1_W = data_W(:,inds1_i);
   tX0_inds_W = find(tX0_W);
   tX1_inds_W = find(tX1_W);
   tdXs_W = tX1_W ./ tX0_W;
   tdxs_W = log10(tdXs_W);
   tdx_inds_W = intersect(tX0_inds_W,tX1_inds_W);
   tdx_W = tdxs_W(tdx_inds_W);
   txo_W = log10(tX0_W(tdx_inds_W));
   tx1_W = log10(tX1_W(tdx_inds_W));
   
   %Parameter fitting
   params_W = fit_laplace(log(10.^(tdx_W))); %fit based on natural log
   b_lp_W = params_W.b;
   bs_lp_m3(i) = b_lp_W;
   
   %Growth rate variability versus abundance
   pc1_W = (txo_W + tx1_W)/2;
   pc2_W = (tx1_W - txo_W);
   bins_xm_W = [-3:.4:-.6 0];
   stds_W = [];
   for j = 1:length(bins_xm_W)-1
        inds_j = find(pc1_W > bins_xm_W(j) & pc1_W < bins_xm_W(j+1));
        pc2_j = pc2_W(inds_j); %Growth rates
        stds_W(j) = std(pc2_j);
   end
   r_inds = find(~isnan(stds_W));
   lb_bins_xm_W = bins_xm_W(1:end-1);
   [r_W,s_w] = polyfit(lb_bins_xm_W(r_inds),stds_W(r_inds),1);
   rs_m3(i) = r_W(1);
   
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


%% Daily growth rates

%Get samples that are 1 day apart
inds0_f4 = [];
inds1_f4 = [];
for j = 1:N_f4
    for k = (j+1):N_f4
        d_i = days_f4(k)-days_f4(j);
        if d_i == 1
            inds0_f4 = [inds0_f4 j];
            inds1_f4 = [inds1_f4 k];
        end
    end
end
   
%Growth rates in log10 base
tX0_f4 = data_f4(:,inds0_f4);
tX1_f4 = data_f4(:,inds1_f4);
tX0_inds_f4 = find(tX0_f4);
tX1_inds_f4 = find(tX1_f4);
tdXs_f4 = tX1_f4 ./ tX0_f4;
tdxs_f4 = log10(tdXs_f4);
tdx_inds_f4 = intersect(tX0_inds_f4,tX1_inds_f4);
tdx_f4 = tdxs_f4(tdx_inds_f4);
txo_f4 = log10(tX0_f4(tdx_inds_f4));
tx1_f4 = log10(tX1_f4(tdx_inds_f4));

%MLE of laplace exponent (in the natural log space)
params_f4 = fit_laplace(log(10.^(tdx_f4)));
b_lp_f4 = params_f4.b;


%% Daily growth rate variability versus mean daily abundance
pc1_f4 = (txo_f4 + tx1_f4)/2;
pc2_f4 = (tx1_f4 - txo_f4);
bins_xm_f4 = [-3:.4:-.6 0];
stds_f4 = [];
for i = 1:length(bins_xm_f4)-1
    inds_i = find(pc1_f4 > bins_xm_f4(i) & pc1_f4 < bins_xm_f4(i+1));
    pc2_i = pc2_f4(inds_i); %Growth rates
    stds_f4(i) = std(pc2_i);
end
[r_f4,s_f4] = polyfit(bins_xm_f4(1:end-1),stds_f4,1);
ste_r_f4 = sqrt(diag(inv(s_f4.R)*inv(s_f4.R')).*s_f4.normr.^2./s_f4.df);



%% Variability of Laplace parameter within human F4
W_f4 = floor(N_f4 / 6); %window size
bs_lp_f4 = [];
rs_f4 = [];

for i = 1:floor(N_f4/W_f4)
    
   %Calculate growth rates within each window
   
   %Take 1/6 of samples
   inds = (i-1)*W_f4 + (1:W_f4); 
   data_W = data_f4(:,inds);
   days_i = days_f4(inds);
   
   inds0_i = [];
   inds1_i = [];
   for j = 1:W_f4
       for k = (j+1):W_f4
           d_i = days_i(k)-days_i(j);
           if d_i == 1
               inds0_i = [inds0_i j];
               inds1_i = [inds1_i k];
           end
       end
   end
   
   tX0_W = data_W(:,inds0_i);
   tX1_W = data_W(:,inds1_i);
   tX0_inds_W = find(tX0_W);
   tX1_inds_W = find(tX1_W);
   tdXs_W = tX1_W ./ tX0_W;
   tdxs_W = log10(tdXs_W);
   tdx_inds_W = intersect(tX0_inds_W,tX1_inds_W);
   tdx_W = tdxs_W(tdx_inds_W);
   txo_W = log10(tX0_W(tdx_inds_W));
   tx1_W = log10(tX1_W(tdx_inds_W));
   
   %Parameter fitting
   params_W = fit_laplace(log(10.^(tdx_W))); %fit based on natural log
   b_lp_W = params_W.b;
   bs_lp_f4(i) = b_lp_W;
   
   %Growth rate variability versus abundance
   pc1_W = (txo_W + tx1_W)/2;
   pc2_W = (tx1_W - txo_W);
   bins_xm_W = [-3:.4:-.6 0];
   stds_W = [];
   for j = 1:length(bins_xm_W)-1
        inds_j = find(pc1_W > bins_xm_W(j) & pc1_W < bins_xm_W(j+1));
        pc2_j = pc2_W(inds_j); %Growth rates
        stds_W(j) = std(pc2_j);
   end
   r_inds = find(~isnan(stds_W));
   lb_bins_xm_W = bins_xm_W(1:end-1);
   [r_W,s_w] = polyfit(lb_bins_xm_W(r_inds),stds_W(r_inds),1);
   rs_f4(i) = r_W(1);
   
end























%% Carmody study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


%% Do analysis on only abundant and prevalent OTUs
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


%% Daiy growth rates

%% HFHS 1
inds0_hf1 = [];
inds1_hf1 = [];
for j = 1:N_hf1
    for k = (j+1):N_hf1
        d_i = days_hf1(k)-days_hf1(j);
        if d_i == 1
            inds0_hf1 = [inds0_hf1 j];
            inds1_hf1 = [inds1_hf1 k];
        end
    end
end
tX0_hf1 = data_hf1(:,inds0_hf1);
tX1_hf1 = data_hf1(:,inds1_hf1);
tX0_inds_hf1 = find(tX0_hf1);
tX1_inds_hf1 = find(tX1_hf1);
tdXs_hf1 = tX1_hf1 ./ tX0_hf1;
tdxs_hf1 = log10(tdXs_hf1);
tdx_inds_hf1 = intersect(tX0_inds_hf1,tX1_inds_hf1);
tdx_hf1 = tdxs_hf1(tdx_inds_hf1);
txo_hf1 = log10(tX0_hf1(tdx_inds_hf1));
tx1_hf1 = log10(tX1_hf1(tdx_inds_hf1));

%MLE of laplace exponent
params_hf1 = fit_laplace(log(10.^(tdx_hf1)));
b_lp_hf1 = params_hf1.b;

% Daily growth rate variability versus mean daily abundance
pc1_hf1 = (txo_hf1 + tx1_hf1)/2;
pc2_hf1 = (tx1_hf1 - txo_hf1);
bins_xm_hf1 = [-3:.4:-.6 0];
stds_hf1 = [];
for i = 1:length(bins_xm_hf1)-1
    inds_i = find(pc1_hf1 > bins_xm_hf1(i) & pc1_hf1 < bins_xm_hf1(i+1));
    pc2_i = pc2_hf1(inds_i); %Growth rates
    stds_hf1(i) = std(pc2_i);
end
[r_hf1,s_hf1] = polyfit(bins_xm_hf1(1:end-1),stds_hf1,1);
ste_r_hf1 = sqrt(diag(inv(s_hf1.R)*inv(s_hf1.R')).*s_hf1.normr.^2./s_hf1.df);



%% HFHS 2
inds0_hf2 = [];
inds1_hf2 = [];
for j = 1:N_hf2
    for k = (j+1):N_hf2
        d_i = days_hf2(k)-days_hf2(j);
        if d_i == 1
            inds0_hf2 = [inds0_hf2 j];
            inds1_hf2 = [inds1_hf2 k];
        end
    end
end
tX0_hf2 = data_hf2(:,inds0_hf2);
tX1_hf2 = data_hf2(:,inds1_hf2);
tX0_inds_hf2 = find(tX0_hf2);
tX1_inds_hf2 = find(tX1_hf2);
tdXs_hf2 = tX1_hf2 ./ tX0_hf2;
tdxs_hf2 = log10(tdXs_hf2);
tdx_inds_hf2 = intersect(tX0_inds_hf2,tX1_inds_hf2);
tdx_hf2 = tdxs_hf2(tdx_inds_hf2);
txo_hf2 = log10(tX0_hf2(tdx_inds_hf2));
tx1_hf2 = log10(tX1_hf2(tdx_inds_hf2));

%MLE of laplace exponent
params_hf2 = fit_laplace(log(10.^(tdx_hf2)));
b_lp_hf2 = params_hf2.b;

% Daily growth rate variability versus mean daily abundance
pc1_hf2 = (txo_hf2 + tx1_hf2)/2;
pc2_hf2 = (tx1_hf2 - txo_hf2);
bins_xm_hf2 = [-3:.4:-.6 0];
stds_hf2 = [];
for i = 1:length(bins_xm_hf2)-1
    inds_i = find(pc1_hf2 > bins_xm_hf2(i) & pc1_hf2 < bins_xm_hf2(i+1));
    pc2_i = pc2_hf2(inds_i); %Growth rates
    stds_hf2(i) = std(pc2_i);
end
[r_hf2,s_hf2] = polyfit(bins_xm_hf2(1:end-1),stds_hf2,1);
ste_r_hf2 = sqrt(diag(inv(s_hf2.R)*inv(s_hf2.R')).*s_hf2.normr.^2./s_hf2.df);


%% HFHS 3
inds0_hf3 = [];
inds1_hf3 = [];
for j = 1:N_hf3
    for k = (j+1):N_hf3
        d_i = days_hf3(k)-days_hf3(j);
        if d_i == 1
            inds0_hf3 = [inds0_hf3 j];
            inds1_hf3 = [inds1_hf3 k];
        end
    end
end
tX0_hf3 = data_hf3(:,inds0_hf3);
tX1_hf3 = data_hf3(:,inds1_hf3);
tX0_inds_hf3 = find(tX0_hf3);
tX1_inds_hf3 = find(tX1_hf3);
tdXs_hf3 = tX1_hf3 ./ tX0_hf3;
tdxs_hf3 = log10(tdXs_hf3);
tdx_inds_hf3 = intersect(tX0_inds_hf3,tX1_inds_hf3);
tdx_hf3 = tdxs_hf3(tdx_inds_hf3);
txo_hf3 = log10(tX0_hf3(tdx_inds_hf3));
tx1_hf3 = log10(tX1_hf3(tdx_inds_hf3));

%MLE of laplace exponent
params_hf3 = fit_laplace(log(10.^(tdx_hf3)));
b_lp_hf3 = params_hf3.b;

% Daily growth rate variability versus mean daily abundance
pc1_hf3 = (txo_hf3 + tx1_hf3)/2;
pc2_hf3 = (tx1_hf3 - txo_hf3);
bins_xm_hf3 = [-3:.4:-.6 0];
stds_hf3 = [];
for i = 1:length(bins_xm_hf3)-1
    inds_i = find(pc1_hf3 > bins_xm_hf3(i) & pc1_hf3 < bins_xm_hf3(i+1));
    pc2_i = pc2_hf3(inds_i); %Growth rates
    stds_hf3(i) = std(pc2_i);
end
[r_hf3,s_hf3] = polyfit(bins_xm_hf3(1:end-1),stds_hf3,1);
ste_r_hf3 = sqrt(diag(inv(s_hf3.R)*inv(s_hf3.R')).*s_hf3.normr.^2./s_hf3.df);



%% LFPP 1
inds0_lf1 = [];
inds1_lf1 = [];
for j = 1:N_lf1
    for k = (j+1):N_lf1
        d_i = days_lf1(k)-days_lf1(j);
        if d_i == 1
            inds0_lf1 = [inds0_lf1 j];
            inds1_lf1 = [inds1_lf1 k];
        end
    end
end
tX0_lf1 = data_lf1(:,inds0_lf1);
tX1_lf1 = data_lf1(:,inds1_lf1);
tX0_inds_lf1 = find(tX0_lf1);
tX1_inds_lf1 = find(tX1_lf1);
tdXs_lf1 = tX1_lf1 ./ tX0_lf1;
tdxs_lf1 = log10(tdXs_lf1);
tdx_inds_lf1 = intersect(tX0_inds_lf1,tX1_inds_lf1);
tdx_lf1 = tdxs_lf1(tdx_inds_lf1);
txo_lf1 = log10(tX0_lf1(tdx_inds_lf1));
tx1_lf1 = log10(tX1_lf1(tdx_inds_lf1));

%MLE of laplace exponent
params_lf1 = fit_laplace(log(10.^(tdx_lf1)));
b_lp_lf1 = params_lf1.b;

% Daily growth rate variability versus mean daily abundance
pc1_lf1 = (txo_lf1 + tx1_lf1)/2;
pc2_lf1 = (tx1_lf1 - txo_lf1);
bins_xm_lf1 = [-3:.4:-.6 0];
stds_lf1 = [];
for i = 1:length(bins_xm_lf1)-1
    inds_i = find(pc1_lf1 > bins_xm_lf1(i) & pc1_lf1 < bins_xm_lf1(i+1));
    pc2_i = pc2_lf1(inds_i); %Growth rates
    stds_lf1(i) = std(pc2_i);
end
[r_lf1,s_lf1] = polyfit(bins_xm_lf1(1:end-1),stds_lf1,1);
ste_r_lf1 = sqrt(diag(inv(s_lf1.R)*inv(s_lf1.R')).*s_lf1.normr.^2./s_lf1.df);



%% LFPP 2
inds0_lf2 = [];
inds1_lf2 = [];
for j = 1:N_lf2
    for k = (j+1):N_lf2
        d_i = days_lf2(k)-days_lf2(j);
        if d_i == 1
            inds0_lf2 = [inds0_lf2 j];
            inds1_lf2 = [inds1_lf2 k];
        end
    end
end
tX0_lf2 = data_lf2(:,inds0_lf2);
tX1_lf2 = data_lf2(:,inds1_lf2);
tX0_inds_lf2 = find(tX0_lf2);
tX1_inds_lf2 = find(tX1_lf2);
tdXs_lf2 = tX1_lf2 ./ tX0_lf2;
tdxs_lf2 = log10(tdXs_lf2);
tdx_inds_lf2 = intersect(tX0_inds_lf2,tX1_inds_lf2);
tdx_lf2 = tdxs_lf2(tdx_inds_lf2);
txo_lf2 = log10(tX0_lf2(tdx_inds_lf2));
tx1_lf2 = log10(tX1_lf2(tdx_inds_lf2));

%MLE of laplace exponent
params_lf2 = fit_laplace(log(10.^(tdx_lf2)));
b_lp_lf2 = params_lf2.b;

% Daily growth rate variability versus mean daily abundance
pc1_lf2 = (txo_lf2 + tx1_lf2)/2;
pc2_lf2 = (tx1_lf2 - txo_lf2);
bins_xm_lf2 = [-3:.4:-.6 0];
stds_lf2 = [];
for i = 1:length(bins_xm_lf2)-1
    inds_i = find(pc1_lf2 > bins_xm_lf2(i) & pc1_lf2 < bins_xm_lf2(i+1));
    pc2_i = pc2_lf2(inds_i); %Growth rates
    stds_lf2(i) = std(pc2_i);
end
[r_lf2,s_lf2] = polyfit(bins_xm_lf2(1:end-1),stds_lf2,1);
ste_r_lf2 = sqrt(diag(inv(s_lf2.R)*inv(s_lf2.R')).*s_lf2.normr.^2./s_lf2.df);


%% LFPP 3
inds0_lf3 = [];
inds1_lf3 = [];
for j = 1:N_lf3
    for k = (j+1):N_lf3
        d_i = days_lf3(k)-days_lf3(j);
        if d_i == 1
            inds0_lf3 = [inds0_lf3 j];
            inds1_lf3 = [inds1_lf3 k];
        end
    end
end
tX0_lf3 = data_lf3(:,inds0_lf3);
tX1_lf3 = data_lf3(:,inds1_lf3);
tX0_inds_lf3 = find(tX0_lf3);
tX1_inds_lf3 = find(tX1_lf3);
tdXs_lf3 = tX1_lf3 ./ tX0_lf3;
tdxs_lf3 = log10(tdXs_lf3);
tdx_inds_lf3 = intersect(tX0_inds_lf3,tX1_inds_lf3);
tdx_lf3 = tdxs_lf3(tdx_inds_lf3);
txo_lf3 = log10(tX0_lf3(tdx_inds_lf3));
tx1_lf3 = log10(tX1_lf3(tdx_inds_lf3));


%MLE of laplace exponent
params_lf3 = fit_laplace(log(10.^(tdx_lf3)));
b_lp_lf3 = params_lf3.b;

% Daily growth rate variability versus mean daily abundance
pc1_lf3 = (txo_lf3 + tx1_lf3)/2;
pc2_lf3 = (tx1_lf3 - txo_lf3);
bins_xm_lf3 = [-3:.4:-.6 0];
stds_lf3 = [];
for i = 1:length(bins_xm_lf3)-1
    inds_i = find(pc1_lf3 > bins_xm_lf3(i) & pc1_lf3 < bins_xm_lf3(i+1));
    pc2_i = pc2_lf3(inds_i); %Growth rates
    stds_lf3(i) = std(pc2_i);
end
[r_lf3,s_lf3] = polyfit(bins_xm_lf3(1:end-1),stds_lf3,1);
ste_r_lf3 = sqrt(diag(inv(s_lf3.R)*inv(s_lf3.R')).*s_lf3.normr.^2./s_lf3.df);




%% Combine HFHS and LFPP data across mice on each diet

tdx_hf = [tdx_hf1; tdx_hf2; tdx_hf3];
txo_hf = [txo_hf1; txo_hf2; txo_hf3];
tx1_hf = [tx1_hf1; tx1_hf2; tx1_hf3];

tdx_lf = [tdx_lf1; tdx_lf2; tdx_lf3];
txo_lf = [txo_lf1; txo_lf2; txo_lf3];
tx1_lf = [tx1_lf1; tx1_lf2; tx1_lf3];

%Laplace exponent of combined HFHS growth rates
params_hf = fit_laplace(log(10.^(tdx_hf)));
b_lp_hf = params_hf.b;
bs_lp_hf = [b_lp_hf1 b_lp_hf2 b_lp_hf3];

%Laplace exponent of combined LFPP growth rates
params_lf = fit_laplace(log(10.^(tdx_lf)));
b_lp_lf = params_lf.b;
bs_lp_lf = [b_lp_lf1 b_lp_lf2 b_lp_lf3];



%% Daily growth rate variability versus mean daily abundance

%% HFHS
pc1_hf = (txo_hf + tx1_hf)/2;
pc2_hf = (tx1_hf - txo_hf);
bins_xm_hf = [-3:.4:-.6 0];
stds_hf = [];
for i = 1:length(bins_xm_hf)-1
    inds_i = find(pc1_hf > bins_xm_hf(i) & pc1_hf < bins_xm_hf(i+1));
    pc2_i = pc2_hf(inds_i); %Growth rates
    stds_hf(i) = std(pc2_i);
end
[r_hf,s_hf] = polyfit(bins_xm_hf(1:end-1),stds_hf,1);
ste_r_hf = sqrt(diag(inv(s_hf.R)*inv(s_hf.R')).*s_hf.normr.^2./s_hf.df);
rs_hf = [r_hf1(1) r_hf2(1) r_hf3(1)];

%% LFPP
pc1_lf = (txo_lf + tx1_lf)/2;
pc2_lf = (tx1_lf - txo_lf);
bins_xm_lf = [-3:.4:-.6 0];
stds_lf = [];
for i = 1:length(bins_xm_lf)-1
    inds_i = find(pc1_lf > bins_xm_lf(i) & pc1_lf < bins_xm_lf(i+1));
    pc2_i = pc2_lf(inds_i); %Growth rates
    stds_lf(i) = std(pc2_i);
end
[r_lf,s_lf] = polyfit(bins_xm_lf(1:end-1),stds_lf,1);
ste_r_lf = sqrt(diag(inv(s_lf.R)*inv(s_lf.R')).*s_lf.normr.^2./s_lf.df);
rs_lf = [r_lf1(1) r_lf2(1) r_lf3(1)];








%% Plotting 
f1 = figure;



%% Alm study %%%%%%%%%%


%Plot pdf of temporal growth rates
subplot(2,3,1);

bins_A = linspace(-4,4,40);
[tcount_A,tval_A] = hist(tdx_A,bins_A);
params_log10_A = fit_laplace(tdx_A);
u_A = params_log10_A.u;
b_A = params_log10_A.b;
lp_A = @(x) (1./(2*b_A)) * exp(-abs(x-u_A) ./ b_A);
x_A = -5:.1:5;
dx_A = bins_A(2)-bins_A(1);

p1a = plot(x_A,lp_A(x_A)*dx_A,'-');
p1a.Color = c('matlab maroon');
p1a.LineWidth = 1.5;
hold on
p1b = plot(tval_A,tcount_A/sum(tcount_A),'.');
p1b.Color = c('black');
p1b.MarkerSize = 22;
hold off
set(gca,'xlim',[-5,5]);
set(gca,'XTick',[-4,-2,0,2,4])
set(gca,'XTickLabel',{'-4','-2','0','2','4'})
set(gca,'ylim',[10^(-4),1]);
set(gca,'yscale','log')
set(gca,'YTick',[1e-4,1e-3,1e-2,1e-1,1])
set(gca,'YTickLabel',{'10^{-4}','','10^{-2}','','10^{0}'})
set(gca,'YMinorTick','off');
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Daily growth rate, \mu');
ylabel('Probability');
title('Human A','FontWeight','Normal','FontSize',17);
box off




% Growth rate variability versus abundance
subplot(2,3,4);

x_A = -3.4:.1:0;
p4a = plot(x_A,polyval(r_A,x_A),'--');
p4a.LineWidth = 1.5;
p4a.Color = c('matlab maroon');
hold on
p4b = plot(bins_xm_A(1:end-1),stds_A,'.');
p4b.Color = c('black');
p4b.MarkerSize = 22;
hold off
set(gca,'xlim',[-3.4,0]);
set(gca,'ylim',[0,.65])
set(gca,'XTick',[-3,-2,-1,0]);
set(gca,'XTickLabel',{'-3','-2','-1','0'});
set(gca,'YTick',[0,.2,.4,.6]);
set(gca,'YTickLabel',{'0','0.2','0.4','0.6'});
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Mean daily abundance, x_{m}');
ylabel('S.d. of daily growth rates \sigma_{\mu}')
box off



%% Caparaso study %%%%%%%%%%


%Plot pdf of temporal growth rates
subplot(2,3,2);


bins_m3 = linspace(-4,4,40);
[tcount_m3,tval_m3] = hist(tdx_m3,bins_m3);
params_log10_m3 = fit_laplace(tdx_m3);
u_m3 = params_log10_m3.u;
b_m3 = params_log10_m3.b;
lp_m3 = @(x) (1./(2*b_m3)) * exp(-abs(x-u_m3) ./ b_m3);
x_m3 = -5:.1:5;
dx_m3 = bins_m3(2)-bins_m3(1);

p2a = plot(x_m3,lp_m3(x_m3)*dx_m3,'-');
p2a.Color = c('matlab blue');
p2a.LineWidth = 1.5;
hold on
p2b = plot(tval_m3,tcount_m3/sum(tcount_m3),'.');
p2b.Color = c('black');
p2b.MarkerSize = 22;
hold off
set(gca,'xlim',[-5,5]);
set(gca,'XTick',[-4,-2,0,2,4])
set(gca,'XTickLabel',{'-4','-2','0','2','4'})
set(gca,'ylim',[10^(-4),1]);
set(gca,'yscale','log')
set(gca,'YTick',[1e-4,1e-3,1e-2,1e-1,1])
set(gca,'YTickLabel',{'10^{-4}','','10^{-2}','','10^{0}'})
set(gca,'YMinorTick','off');
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Daily growth rate, \mu');
title('Human M3','FontWeight','Normal','FontSize',17);
box off




% Growth rate deviations as function of geometric mean
subplot(2,3,5);

x_m3 = -3.4:.1:0;
p5a = plot(x_m3,polyval(r_m3,x_m3),'--');
p5a.LineWidth = 1.5;
p5a.Color = c('matlab blue');
hold on
p5b = plot(bins_xm_m3(1:end-1),stds_m3,'.');
p5b.Color = c('black');
p5b.MarkerSize = 22;
hold off
set(gca,'xlim',[-3.4,0]);
set(gca,'ylim',[0,.65])
set(gca,'XTick',[-3,-2,-1,0]);
set(gca,'XTickLabel',{'-3','-2','-1','0'});
set(gca,'YTick',[0,.2,.4,.6]);
set(gca,'YTickLabel',{'0','0.2','0.4','0.6'});
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Mean daily abundance, x_{m}');
box off


%% Carmody study %%%%%%%%%%

%Plot pdf of temporal growth rates
subplot(2,3,3);


bins_lf = linspace(-4,4,40);
[tcount_lf,tval_lf] = hist(tdx_lf,bins_lf);
params_log10_lf = fit_laplace(tdx_lf);
u_lf = params_log10_lf.u;
b_lf = params_log10_lf.b;
lp_lf = @(x) (1./(2*b_lf)) * exp(-abs(x-u_lf) ./ b_lf);
x_lf = -5:.1:5;
dx_lf = bins_lf(2)-bins_lf(1);

p3a = plot(x_lf,lp_lf(x_lf)*dx_lf,'-');
p3a.Color = c('matlab green');
p3a.LineWidth = 1.5;
hold on
p3b = plot(tval_lf,tcount_lf/sum(tcount_lf),'.');
p3b.Color = c('black');
p3b.MarkerSize = 22;
hold off
set(gca,'xlim',[-3.5,3.5]);
set(gca,'XTick',[-3,-1.5,0,1.5,3])
set(gca,'XTickLabel',{'-3','-1.5','0','1.5','3'})
set(gca,'ylim',[10^(-3),10^(0)]);
set(gca,'yscale','log')
set(gca,'YTick',[1e-3,1e-2,1e-1,1])
set(gca,'YTickLabel',{'10^{-3}','10^{-2}','10^{-1}','10^{0}'})
set(gca,'YMinorTick','off');
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
title('Mice (LFPP Diet)','FontWeight','Normal','FontSize',17);
xlabel('Daily growth rate, \mu');
box off



% Growth rate deviations as function of geometric mean
subplot(2,3,6);

x_lf = -3.4:.1:0;
p6a = plot(x_lf,polyval(r_lf,x_lf),'--');
p6a.LineWidth = 1.5;
p6a.Color = c('matlab green');
hold on
p6b = plot(bins_xm_lf(1:end-1),stds_lf,'.');
p6b.Color = c('black');
p6b.MarkerSize = 22;
hold off
set(gca,'xlim',[-3.4,0]);
set(gca,'ylim',[0,.65])
set(gca,'XTick',[-3,-2,-1,0]);
set(gca,'XTickLabel',{'-3','-2','-1','0'});
set(gca,'YTick',[0,.2,.4,.6]);
set(gca,'YTickLabel',{'0','0.2','0.4','0.6'});
set(gca,'FontSize',14);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Mean daily abundance, x_{m}');
box off








