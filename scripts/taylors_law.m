%% Alm Study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load alm_uclust_25K.mat;
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

proteo_inds_A = find(strcmp('Proteobacteria',ptax_A));

%% Taylor's Law
means_A = mean(data_A,2);
vars_A = var(data_A')';
[coefs_A,s_A] = polyfit(log10(means_A),log10(vars_A),1);
ste_coefs_A = sqrt(diag(inv(s_A.R)*inv(s_A.R')).*s_A.normr.^2./s_A.df);
beta_A = coefs_A(1);
ste_beta_A = ste_coefs_A(1);

%Variability of Taylor's law exponent within human A
W_A = floor(N_A / 6);
betas_A = [];
for i = 1:floor(N_A/W_A)
   inds = (i-1)*W_A + (1:W_A); 
   data_W = data_A(:,inds);
   means_W = mean(data_W,2);
   vars_W = var(data_W')';
   keep_inds = find(log10(means_W) > -100);
   [beta_W,s_W] = polyfit(log10(means_W(keep_inds)),log10(vars_W(keep_inds)),1);
   betas_A(i) = beta_W(1);
end

%Find spiking strains on a single day
fc_cut_A = 25;
SW_A = 1;
fcs_A = [];
for i = 1:N_A-SW_A+1;
   inds = i + (1:SW_A) - 1; 
   noninds = setdiff(1:N_A,inds);
   data_W = data_A(:,inds);
   data_NW = data_A(:,noninds);
   means_W = mean(data_W,2);
   means_NW = mean(data_NW,2);
   fc_A = means_W./means_NW;
   fcs_A(:,i) = fc_A;
end
max_fcs_A = max(fcs_A')';
spike_inds_A = find(max_fcs_A > fc_cut_A); % > 25x fold-change in single day  
nonspike_inds_A = setdiff(1:length(abunds_A),spike_inds_A);
proteo_spike_inds_A = intersect(proteo_inds_A,spike_inds_A);

%Find spiking OTUs during travel period
travel_inds_A = find(stoolA_days >= 71 & stoolA_days <= 122);
nontravel_inds_A = setdiff(1:length(stoolA_days),travel_inds_A);
travel_fcs_A = fcs_A(:,travel_inds_A);
max_travel_fcs_A = max(travel_fcs_A')';
travel_spike_inds_A = find(max_travel_fcs_A > fc_cut_A); % > 25x fold-change in single day  
nontravel_spike_inds_A = setdiff(1:length(abunds_A),travel_spike_inds_A);


%Outliers
out_inds_A = [];
nonout_inds_A = [];

for i = 1:length(abunds_A)
    x = log10(abunds_A(i));
    y = log10(vars_A(i));
    yhat = polyval(coefs_A,x);
    if abs(y-yhat) > log10(6) %6 fold change
        out_inds_A = [out_inds_A i];
    else
        nonout_inds_A = [nonout_inds_A i];
    end  
end


%Deviations of OTUs from Taylor's law
deviations_A = [];
for i = 1:length(abunds_A)
    x = log10(abunds_A(i));
    y = log10(vars_A(i));
    yhat = polyval(coefs_A,x);
    deviations_A(i) = 10^(y-yhat);       
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

proteo_inds_B = find(strcmp('Proteobacteria',ptax_B));

%% Taylor's Law
means_B = mean(data_B,2);
vars_B = var(data_B')';
[coefs_B,s_B] = polyfit(log10(means_B),log10(vars_B),1);
ste_coefs_B = sqrt(diag(inv(s_B.R)*inv(s_B.R')).*s_B.normr.^2./s_B.df);
beta_B = coefs_B(1);
ste_beta_B = ste_coefs_B(1);

%Variability of Taylor's law exponent within human B
W_B = floor(N_B / 6);
betas_B = [];
for i = 1:floor(N_B/W_B)
   inds = (i-1)*W_B + (1:W_B); 
   data_W = data_B(:,inds);
   means_W = mean(data_W,2);
   vars_W = var(data_W')';
   keep_inds = find(log10(means_W) > -100);
   [beta_W,s_W] = polyfit(log10(means_W(keep_inds)),log10(vars_W(keep_inds)),1);
   betas_B(i) = beta_W(1);
end

%Find spiking strains on a single day
fc_cut_B = 25;
SW_B = 1;
fcs_B = [];
for i = 1:N_B-SW_B+1;
   inds = i + (1:SW_B) - 1; 
   noninds = setdiff(1:N_B,inds);
   data_W = data_B(:,inds);
   data_NW = data_B(:,noninds);
   means_W = mean(data_W,2);
   means_NW = mean(data_NW,2);
   fc_B = means_W./means_NW;
   fcs_B(:,i) = fc_B;
end
max_fcs_B = max(fcs_B')';
spike_inds_B = find(max_fcs_B > fc_cut_B); % > 25x fold-change in single day  
nonspike_inds_B = setdiff(1:length(abunds_B),spike_inds_B);
proteo_spike_inds_B = intersect(proteo_inds_B,spike_inds_B);

%Find spiking OTUs during illness period
ill_inds_B = find(stoolB_days >= 151 & stoolB_days <= 159);
nonill_inds_B = setdiff(1:length(stoolB_days),ill_inds_B);
ill_fcs_B = fcs_B(:,ill_inds_B);
max_ill_fcs_B = max(ill_fcs_B')';
ill_spike_inds_B = find(max_ill_fcs_B > fc_cut_B); % > 25x fold-change in single day  
nonill_spike_inds_B = setdiff(1:length(abunds_B),ill_spike_inds_B);

%Outliers
out_inds_B = [];
nonout_inds_B = [];

for i = 1:length(abunds_B)
    x = log10(abunds_B(i));
    y = log10(vars_B(i));
    yhat = polyval(coefs_B,x);
    if abs(y-yhat) > log10(6) %6 fold change
        out_inds_B = [out_inds_B i];
    else
        nonout_inds_B = [nonout_inds_B i];
    end  
end

%Deviations of OTUs from Taylor's law
deviations_B = [];
for i = 1:length(abunds_B)
    x = log10(abunds_B(i));
    y = log10(vars_B(i));
    yhat = polyval(coefs_B,x);
    deviations_B(i) = 10^(y-yhat);       
end







%% Caporaso study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load caporaso_uclust_17K.mat;



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


%% Taylor's Law
means_m3 = mean(data_m3,2);
vars_m3 = var(data_m3')';
[coefs_m3,s_m3] = polyfit(log10(means_m3),log10(vars_m3),1);
ste_coefs_m3 = sqrt(diag(inv(s_m3.R)*inv(s_m3.R')).*s_m3.normr.^2./s_m3.df);
beta_m3 = coefs_m3(1);
ste_beta_m3 = ste_coefs_m3(1);

%Variability of Taylor's law exponent within human M3
W_m3 = floor(N_m3 / 6);
betas_m3 = [];
for i = 1:floor(N_m3/W_m3)
   inds = (i-1)*W_m3 + (1:W_m3); 
   data_W = data_m3(:,inds);
   means_W = mean(data_W,2);
   vars_W = var(data_W')';
   keep_inds = find(log10(means_W) > -100);
   [beta_W,s_W] = polyfit(log10(means_W(keep_inds)),log10(vars_W(keep_inds)),1);
   betas_m3(i) = beta_W(1);
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




%% Taylor's Law
means_f4 = mean(data_f4,2);
vars_f4 = var(data_f4')';
[coefs_f4,s_f4] = polyfit(log10(means_f4),log10(vars_f4),1);
ste_coefs_f4 = sqrt(diag(inv(s_f4.R)*inv(s_f4.R')).*s_f4.normr.^2./s_f4.df);
beta_f4 = coefs_f4(1);
ste_beta_f4 = ste_coefs_f4(1);

%Variability of Taylor's law exponent within human F4
W_f4 = floor(N_f4 / 6);
betas_f4 = [];
for i = 1:floor(N_f4/W_f4)
   inds = (i-1)*W_f4 + (1:W_f4); 
   data_W = data_f4(:,inds);
   means_W = mean(data_W,2);
   vars_W = var(data_W')';
   keep_inds = find(log10(means_W) > -100);
   [beta_W,s_W] = polyfit(log10(means_W(keep_inds)),log10(vars_W(keep_inds)),1);
   betas_f4(i) = beta_W(1);
end














%% Carmody study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load mice_uclust_25K.mat;


%% Set up mice data %%%%%%%%%%%%%%%%%%%%% 
[M,N] = size(lf1);
data_hf1 = hf1 ./ repmat(sum(hf1),M,1);
data_hf2 = hf2 ./ repmat(sum(hf2),M,1);
data_hf3 = hf3 ./ repmat(sum(hf3),M,1);
data_lf1 = lf1 ./ repmat(sum(lf1),M,1);
data_lf2 = lf2 ./ repmat(sum(lf2),M,1);
data_lf3 = lf3 ./ repmat(sum(lf3),M,1);

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


%% Do analysis on only abundant and prevalent OTUS
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






%% Estimate coefficients

%HF1
[coefs_hf1,s_hf1] = polyfit(log10(means_hf1),log10(vars_hf1),1);
ste_coefs_hf1 = sqrt(diag(inv(s_hf1.R)*inv(s_hf1.R')).*s_hf1.normr.^2./s_hf1.df);
beta_hf1 = coefs_hf1(1);
ste_beta_hf1 = ste_coefs_hf1(1);

%HF2
[coefs_hf2,s_hf2] = polyfit(log10(means_hf2),log10(vars_hf2),1);
ste_coefs_hf2 = sqrt(diag(inv(s_hf2.R)*inv(s_hf2.R')).*s_hf2.normr.^2./s_hf2.df);
beta_hf2 = coefs_hf2(1);
ste_beta_hf2 = ste_coefs_hf2(1);

%HF3
[coefs_hf3,s_hf3] = polyfit(log10(means_hf3),log10(vars_hf3),1);
ste_coefs_hf3 = sqrt(diag(inv(s_hf3.R)*inv(s_hf3.R')).*s_hf3.normr.^2./s_hf3.df);
beta_hf3 = coefs_hf3(1);
ste_beta_hf3 = ste_coefs_hf3(1);

%LF1
[coefs_lf1,s_lf1] = polyfit(log10(means_lf1),log10(vars_lf1),1);
ste_coefs_lf1 = sqrt(diag(inv(s_lf1.R)*inv(s_lf1.R')).*s_lf1.normr.^2./s_lf1.df);
beta_lf1 = coefs_lf1(1);
ste_beta_lf1 = ste_coefs_lf1(1);

%LF2
[coefs_lf2,s_lf2] = polyfit(log10(means_lf2),log10(vars_lf2),1);
ste_coefs_lf2 = sqrt(diag(inv(s_lf2.R)*inv(s_lf2.R')).*s_lf2.normr.^2./s_lf2.df);
beta_lf2 = coefs_lf2(1);
ste_beta_lf2 = ste_coefs_lf2(1);

%LF3
[coefs_lf3,s_lf3] = polyfit(log10(means_lf3),log10(vars_lf3),1);
ste_coefs_lf3 = sqrt(diag(inv(s_lf3.R)*inv(s_lf3.R')).*s_lf3.normr.^2./s_lf3.df);
beta_lf3 = coefs_lf3(1);
ste_beta_lf3 = ste_coefs_lf3(1);



%% Combine high fat and low fat diets

%High Fat
means_hf = [means_hf1;means_hf2;means_hf3];
vars_hf = [vars_hf1;vars_hf2;vars_hf3];
[coefs_hf,s_hf] = polyfit(log10(means_hf),log10(vars_hf),1);
ste_coefs_hf = sqrt(diag(inv(s_hf.R)*inv(s_hf.R')).*s_hf.normr.^2./s_hf.df);
beta_hf = coefs_hf(1);
ste_beta_hf = ste_coefs_hf(1);

%Low Fat
means_lf = [means_lf1;means_lf2;means_lf3];
vars_lf = [vars_lf1;vars_lf2;vars_lf3];
[coefs_lf,s_lf] = polyfit(log10(means_lf),log10(vars_lf),1);
ste_coefs_lf = sqrt(diag(inv(s_lf.R)*inv(s_lf.R')).*s_lf.normr.^2./s_lf.df);
beta_lf = coefs_lf(1);
ste_beta_lf = ste_coefs_lf(1);


betas_hf = [beta_hf1 beta_hf2 beta_hf3];
betas_lf = [beta_lf1 beta_lf2 beta_lf3];














%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;



%Human A
subplot(1,3,1);
p1a = plot(-3.5:.1:0,polyval(coefs_A,-3.5:.1:0),'--');
p1a.LineWidth = 1.5;
p1a.Color = c('matlab maroon');
hold on
p1b = plot(log10(means_A(nonspike_inds_A)),log10(vars_A(nonspike_inds_A)),'.');
p1b.Color = c('black');
p1b.MarkerSize = 22;
hold on
p1c = plot(log10(means_A(spike_inds_A)),log10(vars_A(spike_inds_A)),'.');
p1c.Color = c('salmon');
p1c.MarkerSize = 22;
hold on
p1d = plot(log10(means_A(travel_spike_inds_A)),log10(vars_A(travel_spike_inds_A)),'.');
p1d.Color = c('deep sky blue');
p1d.MarkerSize = 22;
hold off
set(gca,'xlim',[-3.5,-.5]);
set(gca,'XTick',[-3,-2,-1]);
set(gca,'XTickLabel',{'10^{-3}','10^{-2}','10^{-1}'});
set(gca,'ylim',[-8,-1]);
set(gca,'YTick',[-7,-5,-3,-1]);
set(gca,'YTickLabel',{'10^{-7}','10^{-5}','10^{-3}','10^{-1}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Average abundance, <X>');
ylabel('Temporal variance, \sigma^2_{X}');
title('Human A','FontSize',16,'FontWeight','Normal');
box off
ll1 = legend([p1d,p1c],{'Travel-related OTUs','Other spiking OTUs'});
ll1.FontSize = 11;
ll1.Location = 'SouthEast';
legend boxoff;

%Human B
subplot(1,3,2);
p2a = plot(-3.5:.1:0,polyval(coefs_B,-3.5:.1:0),'--');
p2a.LineWidth = 1.5;
p2a.Color = c('matlab blue');
hold on
p2b = plot(log10(means_B(nonspike_inds_B)),log10(vars_B(nonspike_inds_B)),'.');
p2b.Color = c('black');
p2b.MarkerSize = 22;
hold on
p2c = plot(log10(means_B(spike_inds_B)),log10(vars_B(spike_inds_B)),'.');
p2c.Color = c('salmon');
p2c.MarkerSize = 22;
hold on
p2d = plot(log10(means_B(ill_spike_inds_B)),log10(vars_B(ill_spike_inds_B)),'.');
p2d.Color = c('deep sky blue');
p2d.MarkerSize = 22;
hold off
set(gca,'xlim',[-3.5,-.5]);
set(gca,'XTick',[-3,-2,-1]);
set(gca,'XTickLabel',{'10^{-3}','10^{-2}','10^{-1}'});
set(gca,'ylim',[-8,-1]);
set(gca,'YTick',[-7,-5,-3,-1]);
set(gca,'YTickLabel',{'10^{-7}','10^{-5}','10^{-3}','10^{-1}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Average abundance, <X>');
%ylabel('log_{10} Temporal variance (\sigma^2_{X})');
title('Human B','FontSize',16,'FontWeight','Normal');
box off
ll2 = legend([p2d,p2c],{'Infection-related OTUs','Other spiking OTUs'});
ll2.FontSize = 11;
ll2.Location = 'SouthEast';
legend boxoff;

%Mice LFPP
subplot(1,3,3);
p3a = plot(-3.5:.1:0,polyval(coefs_lf,-3.5:.1:0),'--');
p3a.LineWidth = 1.5;
p3a.Color = c('matlab green');
hold on
p3b = plot(log10(means_lf),log10(vars_lf),'.');
p3b.Color = c('black');
p3b.MarkerSize = 22;
hold off
set(gca,'xlim',[-3.5,0]);
set(gca,'XTick',[-3,-2,-1,0]);
set(gca,'XTickLabel',{'10^{-3}','10^{-2}','10^{-1}','10^{0}'});
set(gca,'ylim',[-8,-1]);
set(gca,'YTick',[-7,-5,-3,-1]);
set(gca,'YTickLabel',{'10^{-7}','10^{-5}','10^{-3}','10^{-1}'});
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Average abundance, <X>');
title('Mice (LFPP Diet)','FontSize',16,'FontWeight','Normal');
box off






%% Numbers

%A
beta_A;
betas_A;

%B
beta_B;
betas_B;

%M3
betas_m3;
betas_m3;

%F4
coefs_f4;
betas_f4;

%HF
beta_hf;
beta_hf1;
beta_hf2;
beta_hf3;
mean(betas_hf);
std(betas_hf);

%LF
beta_lf;
beta_lf1;
beta_lf2;
beta_lf3;
mean(betas_lf);
std(betas_lf);











