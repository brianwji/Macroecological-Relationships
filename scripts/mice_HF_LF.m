addpath(genpath('./'));

%% Carmody study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load carmody.mat;
c = @cmu.colors;

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
means_hf1 = mean(data_hf1,2);
means_hf2 = mean(data_hf2,2);
means_hf3 = mean(data_hf3,2);
means_lf1 = mean(data_lf1,2);
means_lf2 = mean(data_lf2,2);
means_lf3 = mean(data_lf3,2);


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



%% Daily growth rates

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




%% Combine HFHS and LFPP diets

tdx_hf = [tdx_hf1; tdx_hf2; tdx_hf3];
txo_hf = [txo_hf1; txo_hf2; txo_hf3];
tx1_hf = [tx1_hf1; tx1_hf2; tx1_hf3];

tdx_lf = [tdx_lf1; tdx_lf2; tdx_lf3];
txo_lf = [txo_lf1; txo_lf2; txo_lf3];
tx1_lf = [tx1_lf1; tx1_lf2; tx1_lf3];

%HFHS MLE of laplace exponent
params_hf = fit_laplace(log(10.^(tdx_hf)));
b_lp_hf = params_hf.b;
bs_lp_hf = [b_lp_hf1 b_lp_hf2 b_lp_hf3];

%LFPP MLE of laplace exponent
params_lf = fit_laplace(log(10.^(tdx_lf)));
b_lp_lf = params_lf.b;
bs_lp_lf = [b_lp_lf1 b_lp_lf2 b_lp_lf3];



%% Growth rate variability for HFHS and LFPP diets

%HFHS
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
rs_hf = [r_hf1 r_hf2 r_hf3];

%LFPP
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
rs_lf = [r_lf1 r_lf2 r_lf3];





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



%% LF1
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
ste_Gamma_hf = sqrt(diag(inv(s_hf.R)*inv(s_hf.R')).*s_hf.normr.^2./s_hf.df);

msd_lf = mean([msd_lf1; msd_lf2; msd_lf3],1); %Averaged over three mie
[Gamma_lf,s_lf] = polyfit(log10(lags_lf),log10(msd_lf),1); 
ste_Gamma_lf = sqrt(diag(inv(s_lf.R)*inv(s_lf.R')).*s_lf.normr.^2./s_lf.df);


Gammas_lf = [Gamma_lf1(1) Gamma_lf2(1) Gamma_lf3(1)];
Gammas_hf = [Gamma_hf1(1) Gamma_hf2(1) Gamma_hf3(1)];












%% Taylors Law

%HFHS 1
[coefs_hf1,s_hf1] = polyfit(log10(means_hf1),log10(vars_hf1),1);
ste_coefs_hf1 = sqrt(diag(inv(s_hf1.R)*inv(s_hf1.R')).*s_hf1.normr.^2./s_hf1.df);
beta_hf1 = coefs_hf1(1);
ste_beta_hf1 = ste_coefs_hf1(1);

%HFHS 2
[coefs_hf2,s_hf2] = polyfit(log10(means_hf2),log10(vars_hf2),1);
ste_coefs_hf2 = sqrt(diag(inv(s_hf2.R)*inv(s_hf2.R')).*s_hf2.normr.^2./s_hf2.df);
beta_hf2 = coefs_hf2(1);
ste_beta_hf2 = ste_coefs_hf2(1);

%HFHS 3
[coefs_hf3,s_hf3] = polyfit(log10(means_hf3),log10(vars_hf3),1);
ste_coefs_hf3 = sqrt(diag(inv(s_hf3.R)*inv(s_hf3.R')).*s_hf3.normr.^2./s_hf3.df);
beta_hf3 = coefs_hf3(1);
ste_beta_hf3 = ste_coefs_hf3(1);

%LFPP 1
[coefs_lf1,s_lf1] = polyfit(log10(means_lf1),log10(vars_lf1),1);
ste_coefs_lf1 = sqrt(diag(inv(s_lf1.R)*inv(s_lf1.R')).*s_lf1.normr.^2./s_lf1.df);
beta_lf1 = coefs_lf1(1);
ste_beta_lf1 = ste_coefs_lf1(1);

%LFFP 2
[coefs_lf2,s_lf2] = polyfit(log10(means_lf2),log10(vars_lf2),1);
ste_coefs_lf2 = sqrt(diag(inv(s_lf2.R)*inv(s_lf2.R')).*s_lf2.normr.^2./s_lf2.df);
beta_lf2 = coefs_lf2(1);
ste_beta_lf2 = ste_coefs_lf2(1);

%LFPP 3
[coefs_lf3,s_lf3] = polyfit(log10(means_lf3),log10(vars_lf3),1);
ste_coefs_lf3 = sqrt(diag(inv(s_lf3.R)*inv(s_lf3.R')).*s_lf3.normr.^2./s_lf3.df);
beta_lf3 = coefs_lf3(1);
ste_beta_lf3 = ste_coefs_lf3(1);




%% Combine diets

%HFHS
means_hf = [means_hf1;means_hf2;means_hf3];
vars_hf = [vars_hf1;vars_hf2;vars_hf3];
[coefs_hf,s_hf] = polyfit(log10(means_hf),log10(vars_hf),1);
ste_coefs_hf = sqrt(diag(inv(s_hf.R)*inv(s_hf.R')).*s_hf.normr.^2./s_hf.df);
beta_hf = coefs_hf(1);
ste_beta_hf = ste_coefs_hf(1);

tax_hf = [tax_hf1;tax_hf2;tax_hf3];
ptax_hf = [ptax_hf1;ptax_hf2;ptax_hf3];
ctax_hf = [ctax_hf1;ctax_hf2;ctax_hf3];
otax_hf = [otax_hf1;otax_hf2;otax_hf3];
ftax_hf = [ftax_hf1;ftax_hf2;ftax_hf3];
gtax_hf = [gtax_hf1;gtax_hf2;gtax_hf3];
otus_hf = [otus_hf1;otus_hf2;otus_hf3];

%LFPP
means_lf = [means_lf1;means_lf2;means_lf3];
vars_lf = [vars_lf1;vars_lf2;vars_lf3];
[coefs_lf,s_lf] = polyfit(log10(means_lf),log10(vars_lf),1);
ste_coefs_lf = sqrt(diag(inv(s_lf.R)*inv(s_lf.R')).*s_lf.normr.^2./s_lf.df);
beta_lf = coefs_lf(1);
ste_beta_lf = ste_coefs_lf(1);

tax_lf = [tax_lf1;tax_lf2;tax_lf3];
ptax_lf = [ptax_lf1;ptax_lf2;ptax_lf3];
ctax_lf = [ctax_lf1;ctax_lf2;ctax_lf3];
otax_lf = [otax_lf1;otax_lf2;otax_lf3];
ftax_lf = [ftax_lf1;ftax_lf2;ftax_lf3];
gtax_lf = [gtax_lf1;gtax_lf2;gtax_lf3];
otus_lf = [otus_lf1;otus_lf2;otus_lf3];

betas_hf = [beta_hf1 beta_hf2 beta_hf3];
betas_lf = [beta_lf1 beta_lf2 beta_lf3];


%% Outliers 

%% HFHS
%Likelihood appraoch
all_coefs_hf = [];
ps_hf = [];

for i = 1:length(means_hf)
    x = log10(means_hf(i));
    y = log10(vars_hf(i));
    means_hf_i = means_hf;
    means_hf_i(i) = [];
    vars_hf_i = vars_hf;
    vars_hf_i(i) = [];
    xs = log10(means_hf_i);
    ys = log10(vars_hf_i);
    coefs_i = polyfit(xs,ys,1);
    rs = ys - polyval(coefs_i,xs);
    sigma = std(rs); %Estimate residual error
    test_stat = y - polyval(coefs_i,x);
    ps_hf(i) = normcdf(test_stat,0,sigma);
    all_coefs_hf(i,:) = coefs_i; 
end

out_inds_hf = find(ps_hf < .025 | ps_hf > .975); 
nonout_inds_hf = setdiff(1:length(ps_hf),out_inds_hf);

%by phyla
bacter_inds_hf = find(strcmp('Bacteroidetes',ptax_hf));
nonbacter_inds_hf = setdiff(1:length(ptax_hf),bacter_inds_hf);
firmi_inds_hf = find(strcmp('Firmicutes',ptax_hf));
tener_inds_hf = find(strcmp('Tenericutes',ptax_hf));
other_inds_hf = tener_inds_hf;

%by phyla and outliers
bacter_out_inds_hf = intersect(bacter_inds_hf,out_inds_hf);
firmi_out_inds_hf = intersect(firmi_inds_hf,out_inds_hf);
tener_out_inds_hf = intersect(tener_inds_hf,out_inds_hf);

%No Bacteroidetes
[coefs_nonbacter_hf,s_nonbacter_hf] = polyfit(log10(means_hf(nonbacter_inds_hf)),log10(vars_hf(nonbacter_inds_hf)),1);
ste_coefs_nonbacter_hf = sqrt(diag(inv(s_nonbacter_hf.R)*inv(s_nonbacter_hf.R')).*s_nonbacter_hf.normr.^2./s_nonbacter_hf.df);
beta_nonbacter_hf = coefs_nonbacter_hf(1);
ste_beta_nonbacter_hf = ste_coefs_nonbacter_hf(1);



%% LFPP

all_coefs_lf = [];
ps_lf = [];

for i = 1:length(means_lf)
    x = log10(means_lf(i));
    y = log10(vars_lf(i));
    means_lf_i = means_lf;
    means_lf_i(i) = [];
    vars_lf_i = vars_lf;
    vars_lf_i(i) = [];
    xs = log10(means_lf_i);
    ys = log10(vars_lf_i);
    coefs_i = polyfit(xs,ys,1);
    rs = ys - polyval(coefs_i,xs);
    sigma = std(rs); %Estimate variance of residuals
    test_stat = y - polyval(coefs_i,x);
    ps_lf(i) = normcdf(test_stat,0,sigma);
    all_coefs_lf(i,:) = coefs_i; 
end

out_inds_lf = find(ps_lf < .025 | ps_lf > .975); 
nonout_inds_lf = setdiff(1:length(ps_lf),out_inds_lf);

%by phyla
actino_inds_lf = find(strcmp('Actinobacteria',ptax_lf));
aq_inds_lf = find(strcmp('Aquificae',ptax_lf));
bacter_inds_lf = find(strcmp('Bacteroidetes',ptax_lf));
nonbacter_inds_lf = setdiff(1:length(ptax_lf),bacter_inds_lf);
firmi_inds_lf = find(strcmp('Firmicutes',ptax_lf));
proteo_inds_lf = find(strcmp('Proteobacteria',ptax_lf));
tener_inds_lf = find(strcmp('Tenericutes',ptax_lf));
other_inds_lf = [aq_inds_lf;tener_inds_lf];

%by phyla and outliers
actino_out_inds_lf = intersect(actino_inds_lf,out_inds_lf);
aq_out_inds_lf = intersect(aq_inds_lf,out_inds_lf);
bacter_out_inds_lf = intersect(bacter_inds_lf,out_inds_lf);
firmi_out_inds_lf = intersect(firmi_inds_lf,out_inds_lf);
proteo_out_inds_lf = intersect(proteo_inds_lf,out_inds_lf);
tener_out_inds_lf = intersect(tener_inds_lf,out_inds_lf);

%No Bacteroidetes
[coefs_nonbacter_lf,s_nonbacter_lf] = polyfit(log10(means_lf(nonbacter_inds_lf)),log10(vars_lf(nonbacter_inds_lf)),1);
ste_coefs_nonbacter_lf = sqrt(diag(inv(s_nonbacter_lf.R)*inv(s_nonbacter_lf.R')).*s_nonbacter_lf.normr.^2./s_nonbacter_lf.df);
beta_nonbacter_lf = coefs_nonbacter_lf(1);
ste_beta_nonbacter_lf = ste_coefs_nonbacter_lf(1);


%% Plotting
f1 = figure;



%% Daily growth rate variability
subplot(2,2,1);

%High fat
p1a = plot(bins_xm_hf(1:end-1),stds_hf,'.');
p1a.Color = c('matlab purple');
p1a.MarkerSize = 23;
p1a.LineWidth = 1.5;
hold on
x_hf = -3.4:.1:0;
p1b = plot(x_hf,polyval(r_hf,x_hf),'--');
p1b.Color = c('matlab purple');
p1b.LineWidth = 1;

%Low fat
p1c = plot(bins_xm_lf(1:end-1),stds_lf,'.');
p1c.Color = c('matlab green');
p1c.MarkerSize = 23;
p1c.LineWidth = 1.5;
hold on
x_lf = -3.4:.1:0;
p1d = plot(x_lf,polyval(r_lf,x_lf),'--');
p1d.Color = c('matlab green');
p1d.LineWidth = 1;
hold off
xlabel('Mean daily abundance, x_{m}')
ylabel('S.d. of daily growth rates, \sigma_{\mu}')
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'xlim',[-3.4,0]);
set(gca,'ylim',[0,.6]);
set(gca,'YTick',[0,.2,.4,.6]);
box off
ll = legend([p1a,p1c],'HFHS','LFPP');
ll.FontSize = 12;
ll.Location = 'southwest';
legend boxoff;



%% Hurst averaged over mice in each group
subplot(2,2,2);
p2a = plot(log10(lags_hf),log10(msd_hf),'k.');
p2a.MarkerSize = 22;
p2a.Color = c('matlab purple');
hold on
x_hf = -.2:.1:2.2;
p2b = plot(x_hf,polyval(Gamma_hf,x_hf),'--k');
p2b.LineWidth = 1;
p2b.Color = c('matlab purple');
hold on
p2c = plot(log10(lags_lf),log10(msd_lf),'k.');
p2c.MarkerSize = 22;
p2c.Color = c('matlab green');
hold on
x_lf = -.2:.1:2.2;
p2d = plot(x_lf,polyval(Gamma_lf,x_lf),'--k');
p2d.LineWidth = 1;
p2d.Color = c('matlab green');
hold off
xlabel('Time, \Delta t (days)');
ylabel('Mean-squared displacement');
set(gca,'xlim',[-.2,1.5])
set(gca,'ylim',[-.9,-.2])
set(gca,'XTick',[0 .5 1 1.5]);
Xlabel = {};
Xlabel{1} = '10^{0}';
Xlabel{2} = '';
Xlabel{3} = '10^{1}';
Xlabel{4} = '';
set(gca,'XTickLabel',Xlabel);
set(gca,'ylim',[-1,0])
set(gca,'YTick',[-1:.25:0]);
set(gca,'YTickLabel',{'10^{-1}','','10^{-0.5}','','10^{0}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1); set(gca,'FontName','Arial');
ll = legend([p2a,p2c],{'HFHS','LFPP'});
ll.FontSize = 12;
ll.Location = 'northwest';
legend boxoff
box off



%% Low Fat
cc = colormap(pink(100));
ccb = colorGradient(c('matlab yellow'),c('matlab purple'),4);
ccb(5,:) = c('matlab blue');


subplot(2,2,3);

%Plot regression line
p3a = plot(-3.6:.1:-.2,polyval(coefs_lf,-3.6:.1:-.2),'--');
p3a.LineWidth = 1.2;
p3a.Color = c('matlab green');
hold on

%Plot all OTUs
p3b = plot(log10(means_lf(firmi_inds_lf)),log10(vars_lf(firmi_inds_lf)),'o');
p3b.MarkerSize = 4;
p3b.Color = cc(60,:);
p3b.MarkerFaceColor = cc(60,:);
hold on
p3c = plot(log10(means_lf(bacter_inds_lf)),log10(vars_lf(bacter_inds_lf)),'o');
p3c.MarkerSize = 4;
p3c.Color = ccb(5,:);
p3c.MarkerFaceColor = ccb(5,:);
hold on
p3d = plot(log10(means_lf(actino_inds_lf)),log10(vars_lf(actino_inds_lf)),'o');
p3d.MarkerSize = 4;
p3d.Color = ccb(1,:);
p3d.MarkerFaceColor = ccb(1,:);
hold on
p3e = plot(log10(means_lf(aq_inds_lf)),log10(vars_lf(aq_inds_lf)),'o');
p3e.MarkerSize = 4;
p3e.Color = ccb(4,:);
p3e.MarkerFaceColor = ccb(4,:);
hold on
p3f = plot(log10(means_lf(proteo_inds_lf)),log10(vars_lf(proteo_inds_lf)),'o');
p3f.MarkerSize = 4;
p3f.Color = c('matlab maroon');
p3f.MarkerFaceColor = c('matlab maroon');
hold on
p3g = plot(log10(means_lf(tener_inds_lf)),log10(vars_lf(tener_inds_lf)),'o');
p3g.MarkerSize = 4;
p3g.Color = ccb(4,:);
p3g.MarkerFaceColor = ccb(4,:);


hold off
set(gca,'xlim',[-3.7,.1]);
set(gca,'XTick',[-3,-2,-1,0]);
set(gca,'XTickLabel',{'10^{-3}','10^{-2}','10^{-1}','10^{0}'});
set(gca,'ylim',[-8,-1]);
set(gca,'YTick',[-7,-5,-3,-1]);
set(gca,'YTickLabel',{'10^{-7}','10^{-5}','10^{-3}','10^{-1}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1); 
set(gca,'FontName','Arial');
xlabel('Average abundance, <X>');
ylabel('Temporal variance, \sigma^2_{X}')  
title('Mice (LFPP Diet)','FontSize',16,'FontWeight','Normal');
box off



%% HFHS
subplot(2,2,4);

%Plot regression line
p4a = plot(-3.6:.1:-.2,polyval(coefs_hf,-3.6:.1:-.2),'--');
p4a.LineWidth = 1.2;
p4a.Color = c('matlab purple');
hold on

%All OTUs
p4b = plot(log10(means_hf(firmi_inds_hf)),log10(vars_hf(firmi_inds_hf)),'o');
p4b.MarkerSize = 4;
p4b.Color = cc(60,:);
p4b.MarkerFaceColor = cc(60,:);
hold on
p4c = plot(log10(means_hf(bacter_inds_hf)),log10(vars_hf(bacter_inds_hf)),'o');
p4c.MarkerSize = 4;
p4c.Color = ccb(5,:);
p4c.MarkerFaceColor = ccb(5,:);
hold on
p4d = plot(log10(means_hf(tener_inds_hf)),log10(vars_hf(tener_inds_hf)),'o');
p4d.MarkerSize = 4;
p4d.Color = ccb(4,:);
p4d.MarkerFaceColor = ccb(4,:);
hold off

set(gca,'xlim',[-3.7,.1]);
set(gca,'XTick',[-3,-2,-1,0]);
set(gca,'XTickLabel',{'10^{-3}','10^{-2}','10^{-1}','10^{0}'});
set(gca,'ylim',[-8,-1]);
set(gca,'YTick',[-7,-5,-3,-1]);
set(gca,'YTickLabel',{'10^{-7}','10^{-5}','10^{-3}','10^{-1}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Average abundance, <X>');
ylabel('Temporal variance, \sigma^2_{X}')  
title('Mice (HFHS Diet)','FontSize',16,'FontWeight','Normal')
box off




%% Legend
f2 = figure;

p1 = plot(1,1,'o');
p1.MarkerSize = 8;
p1.Color = ccb(1,:);
p1.MarkerFaceColor = ccb(1,:);
hold on
p2 = plot(1,1,'o');
p2.MarkerSize = 8;
p2.Color = ccb(2,:);
p2.MarkerFaceColor = ccb(2,:);
hold on
p3 = plot(1,1,'o');
p3.MarkerSize = 8;
p3.Color = c('matlab maroon');
p3.MarkerFaceColor = c('matlab maroon');
hold on
p4 = plot(1,1,'o');
p4.MarkerSize = 8;
p4.Color = ccb(4,:);
p4.MarkerFaceColor = ccb(4,:);
hold on
p5 = plot(1,1,'o');
p5.MarkerSize = 8;
p5.Color = ccb(5,:);
p5.MarkerFaceColor = ccb(5,:);
hold on
p6 = plot(1,1,'o');
p6.MarkerSize = 8;
p6.Color = cc(60,:);
p6.MarkerFaceColor = cc(60,:);
hold off
ll = legend([p6,p1,p3,p5,p4],{'Firmicutes','Actinobacteria','Proteobacteria','Bacteroidetes','Other'});
ll.FontSize = 12;
legend boxoff;


    
