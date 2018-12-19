addpath(genpath('./');

%% David study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load david.mat;
c = @cmu.colors;


%% Human A %%%%%%%%%%% 
[M_A,N_A] = size(stoolA);
data_A = stoolA ./ repmat(sum(stoolA),M_A,1);

%Find most prevalent species;
prevs_A = sum(data_A>0,2);
most_prev_A = N_A/2;
most_prev_inds_A = find(prevs_A >= most_prev_A);

%Find most abundant species;
abunds_A = mean(data_A,2);


%% Waiting time distributions
ldata_A = logical(data_A);
pcell_A = {};
acell_A = {};
ptimes_A = [];
atimes_A = [];


%Absence times (return times)
for i = 1:M_A
    sig = ldata_A(i,:);
    first_inds = find(sig); %present days
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1; %fill in all ones before first appearance 
        last_ind = first_inds(end);
        sig(last_ind:end) = 1; %fill in all ones after last appearance
    else
        sig(1:end) = 1; %Ignore entire row of absence
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0); %First day of absence
    end_ind = find(dsig > 0) - 1; %Last day of absence
    atime = end_ind - start_ind + 1;
    acell_A{i} = atime;
    if atime ~= N_A
    atimes_A = [atimes_A atime];
    end
end


%Presence times (residence times) (residence times)
ldata_not_A = ~ldata_A;
for i = 1:M_A
    sig = ldata_not_A(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_A{i} = ptime;
    ptimes_A = [ptimes_A ptime];
end




%% Mean residence/return time for each OTU

mean_ptimes_A = [];
for i = 1:length(pcell_A)
    if length(pcell_A{i}) > 0
        mean_ptime = mean(pcell_A{i});
        mean_ptimes_A(i) = mean_ptime;
    else
        mean_ptimes_A(i) = NaN;
    end
end


mean_atimes_A = [];
for i = 1:length(acell_A) 
    if length(acell_A{i}) > 0
        mean_atime = mean(acell_A{i});
        mean_atimes_A(i) = mean_atime;
    else
        mean_atimes_A(i) = NaN;
    end
end





%% Human B %%%%%%%%%%% 
[M_B,N_B] = size(stoolB);
data_B = stoolB ./ repmat(sum(stoolB),M_B,1);

%Find most prevalent species;
prevs_B = sum(data_B>0,2);
most_prev_B = N_B/2;
most_prev_inds_B = find(prevs_B >= most_prev_B);

%Find most abundant species;
abunds_B = mean(data_B,2);




%% Waiting time distributions
ldata_B = logical(data_B);
pcell_B = {};
acell_B = {};
ptimes_B = [];
atimes_B = [];

%Absence times (return times)
for i = 1:M_B
    sig = ldata_B(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1;        
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_B{i} = atime;
    if atime ~= N_B
    atimes_B = [atimes_B atime];
    end
end

%Presence times (residence times) (residence times)
ldata_not_B = ~ldata_B;
for i = 1:M_B
    sig = ldata_not_B(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_B{i} = ptime;
    ptimes_B = [ptimes_B ptime];
end

%% Mean residence/return time for each OTU

mean_ptimes_B = [];
for i = 1:length(pcell_B)
    if length(pcell_B{i}) > 0
        mean_ptime = mean(pcell_B{i});
        mean_ptimes_B(i) = mean_ptime;
    else
        mean_ptimes_B(i) = NaN;
    end
end

mean_atimes_B = [];
for i = 1:length(acell_B) 
    if length(acell_B{i}) > 0
        mean_atime = mean(acell_B{i});
        mean_atimes_B(i) = mean_atime;
    else
        mean_atimes_B(i) = NaN;
    end
end



















%% Caparaso study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load caporaso.mat;


%% Human M3 %%%%%%%%%%% 
[M_m3,N_m3] = size(m3);
data_m3 = m3 ./ repmat(sum(m3),M_m3,1);

%Find most prevalent species;
prevs_m3 = sum(data_m3>0,2);
most_prev_m3 = N_m3/2;
most_prev_inds_m3 = find(prevs_m3 >= most_prev_m3);

%Find most abundant species;
abunds_m3 = mean(data_m3,2);



%% Waiting time distributions
ldata_m3 = logical(data_m3);
pcell_m3 = {};
acell_m3 = {};
ptimes_m3 = [];
atimes_m3 = [];

%Absence times (return times)
for i = 1:M_m3
    sig = ldata_m3(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_m3{i} = atime;
    if atime ~= N_m3
    atimes_m3 = [atimes_m3 atime];
    end
end

%Presence times (residence times) (residence times)
ldata_not_m3 = ~ldata_m3;
for i = 1:M_m3
    sig = ldata_not_m3(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_m3{i} = ptime;
    ptimes_m3 = [ptimes_m3 ptime];
end

%% Mean residence/return time for each OTU

mean_ptimes_m3 = [];
for i = 1:length(pcell_m3)
    if length(pcell_m3{i}) > 0
        mean_ptime = mean(pcell_m3{i});
        mean_ptimes_m3(i) = mean_ptime;
    else
        mean_ptimes_m3(i) = NaN;
    end
end

mean_atimes_m3 = [];
for i = 1:length(acell_m3) 
    if length(acell_m3{i}) > 0
        mean_atime = mean(acell_m3{i});
        mean_atimes_m3(i) = mean_atime;
    else
        mean_atimes_m3(i) = NaN;
    end
end




%% Human F4 %%%%%%%%%%% 
[M_f4,N_f4] = size(f4);
data_f4 = f4 ./ repmat(sum(f4),M_f4,1);

%Find most prevalent species;
prevs_f4 = sum(data_f4>0,2);
most_prev_f4 = N_f4/2;
most_prev_inds_f4 = find(prevs_f4 >= most_prev_f4);

%Find most abundant species;
abunds_f4 = mean(data_f4,2);


%% Waiting time distributions
ldata_f4 = logical(data_f4);
pcell_f4 = {};
acell_f4 = {};
ptimes_f4 = [];
atimes_f4 = [];


%Absence times (return times)
for i = 1:M_f4
    sig = ldata_f4(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_f4{i} = atime;
    if atime ~= N_f4
    atimes_f4 = [atimes_f4 atime];
    end
end

%Presence times (residence times) (residence times)
ldata_not_f4 = ~ldata_f4;
for i = 1:M_f4
    sig = ldata_not_f4(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_f4{i} = ptime;
    ptimes_f4 = [ptimes_f4 ptime];
end

%% Mean residence/return time for each OTU

mean_ptimes_f4 = [];
for i = 1:length(pcell_f4)
    if ~isempty(pcell_f4{i})
        mean_ptime = mean(pcell_f4{i});
        mean_ptimes_f4(i) = mean_ptime;
    else
        mean_ptimes_f4(i) = NaN;
    end
end

mean_atimes_f4 = [];
for i = 1:length(acell_f4) 
    if length(acell_f4{i}) > 0
        mean_atime = mean(acell_f4{i});
        mean_atimes_f4(i) = mean_atime;
    else
        mean_atimes_f4(i) = NaN;
    end
end









%% Carmody study %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load carmody.mat;


%% Set up mice data %%%%%%%%%%%%%%%%%%%%% 
[M,N] = size(lf1);

data_hf1 = hf1 ./ repmat(sum(hf1),M,1);
data_hf2 = hf2 ./ repmat(sum(hf2),M,1);
data_hf3 = hf3 ./ repmat(sum(hf3),M,1);
data_lf1 = lf1 ./ repmat(sum(lf1),M,1);
data_lf2 = lf2 ./ repmat(sum(lf2),M,1);
data_lf3 = lf3 ./ repmat(sum(lf3),M,1);


%Find most abundant species;
abunds_hf1 = mean(hf1,2);
abunds_hf2 = mean(hf2,2);
abunds_hf3 = mean(hf3,2);
abunds_lf1 = mean(lf1,2);
abunds_lf2 = mean(lf2,2);
abunds_lf3 = mean(lf3,2);



%% Waiting time distributions

%% HFHS 1
ldata_hf1 = logical(hf1);
[m,n] = size(hf1);
pcell_hf1 = {};
acell_hf1 = {};
ptimes_hf1 = [];
atimes_hf1 = [];

%Absence times (return times)
for i = 1:m
    sig = ldata_hf1(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_hf1{i} = atime;
    if atime ~= n
    atimes_hf1 = [atimes_hf1 atime];
    end
end

%Presence times (residence times)
ldata_not_hf1 = ~ldata_hf1;
for i = 1:m
    sig = ldata_not_hf1(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_hf1{i} = ptime;
    ptimes_hf1 = [ptimes_hf1 ptime];
end

% Mean residence/return time for each OTU
mean_ptimes_hf1 = [];
for i = 1:length(pcell_hf1)
   mean_ptime_hf1 = mean(pcell_hf1{i});
   mean_ptimes_hf1(i) = mean_ptime_hf1; 
end

mean_atimes_hf1 = [];
for i = 1:length(acell_hf1)   
   mean_atime_hf1 = mean(acell_hf1{i});
   mean_atimes_hf1(i) = mean_atime_hf1;
end

%% HFHS 2

ldata_hf2 = logical(hf2);
[m,n] = size(hf2);
pcell_hf2 = {};
acell_hf2 = {};
ptimes_hf2 = [];
atimes_hf2 = [];

%Absence times (return times)
for i = 1:m
    sig = ldata_hf2(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_hf2{i} = atime;
    if atime ~= n
    atimes_hf2 = [atimes_hf2 atime];
    end
end

%Presence times (residence times)
ldata_not_hf2 = ~ldata_hf2;
for i = 1:m
    sig = ldata_not_hf2(i,:);
    first_inds = find(sig);
    if ~isempty(first_inds)
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_hf2{i} = ptime;
    ptimes_hf2 = [ptimes_hf2 ptime];
end

% Mean residence/return time for each OTU
mean_ptimes_hf2 = [];
for i = 1:length(pcell_hf2)
   mean_ptime_hf2 = mean(pcell_hf2{i});
   mean_ptimes_hf2(i) = mean_ptime_hf2; 
end

mean_atimes_hf2 = [];
for i = 1:length(acell_hf2)   
   mean_atime_hf2 = mean(acell_hf2{i});
   mean_atimes_hf2(i) = mean_atime_hf2;
end

%%  HFHS 3
ldata_hf3 = logical(hf3);
[m,n] = size(hf3);
pcell_hf3 = {};
acell_hf3 = {};
ptimes_hf3 = [];
atimes_hf3 = [];

%Absence times (return times)
for i = 1:m
    sig = ldata_hf3(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_hf3{i} = atime;
    if atime ~= n
    atimes_hf3 = [atimes_hf3 atime];
    end
end

%Presence times (residence times)
ldata_not_hf3 = ~ldata_hf3;
for i = 1:m
    sig = ldata_not_hf3(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_hf3{i} = ptime;
    ptimes_hf3 = [ptimes_hf3 ptime];
end

% Mean residence/return time for each OTU
mean_ptimes_hf3 = [];
for i = 1:length(pcell_hf3)
   mean_ptime_hf3 = mean(pcell_hf3{i});
   mean_ptimes_hf3(i) = mean_ptime_hf3; 
end

mean_atimes_hf3 = [];
for i = 1:length(acell_hf3)   
   mean_atime_hf3 = mean(acell_hf3{i});
   mean_atimes_hf3(i) = mean_atime_hf3;
end


%%  LFPP 1
ldata_lf1 = logical(lf1);
[m,n] = size(lf1);
pcell_lf1 = {};
acell_lf1 = {};
ptimes_lf1 = [];
atimes_lf1 = [];

%Absence times (return times)
for i = 1:m
    sig = ldata_lf1(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_lf1{i} = atime;
    if atime ~= n
    atimes_lf1 = [atimes_lf1 atime];
    end
end

%Presence times (residence times)
ldata_not_lf1 = ~ldata_lf1;
for i = 1:m
    sig = ldata_not_lf1(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_lf1{i} = ptime;
    ptimes_lf1 = [ptimes_lf1 ptime];
end

% Mean residence/return time for each OTU
mean_ptimes_lf1 = [];
for i = 1:length(pcell_lf1)
   mean_ptime_lf1 = mean(pcell_lf1{i});
   mean_ptimes_lf1(i) = mean_ptime_lf1; 
end

mean_atimes_lf1 = [];
for i = 1:length(acell_lf1)   
   mean_atime_lf1 = mean(acell_lf1{i});
   mean_atimes_lf1(i) = mean_atime_lf1;
end


%%  LFPP 2
ldata_lf2 = logical(lf2);
[m,n] = size(lf2);
pcell_lf2 = {};
acell_lf2 = {};
ptimes_lf2 = [];
atimes_lf2 = [];

%Absence times (return times)
for i = 1:m
    sig = ldata_lf2(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_lf2{i} = atime;
    if atime ~= n
    atimes_lf2 = [atimes_lf2 atime];
    end
end

%Presence times (residence times)
ldata_not_lf2 = ~ldata_lf2;
for i = 1:m
    sig = ldata_not_lf2(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_lf2{i} = ptime;
    ptimes_lf2 = [ptimes_lf2 ptime];
end

% Mean residence/return time for each OTU
mean_ptimes_lf2 = [];
for i = 1:length(pcell_lf2)
   mean_ptime_lf2 = mean(pcell_lf2{i});
   mean_ptimes_lf2(i) = mean_ptime_lf2; 
end

mean_atimes_lf2 = [];
for i = 1:length(acell_lf2)   
   mean_atime_lf2 = mean(acell_lf2{i});
   mean_atimes_lf2(i) = mean_atime_lf2;
end



%%  LFPP 3
ldata_lf3 = logical(lf3);
[m,n] = size(lf3);
pcell_lf3 = {};
acell_lf3 = {};
ptimes_lf3 = [];
atimes_lf3 = [];

%Absence times (return times)
for i = 1:m
    sig = ldata_lf3(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    atime = end_ind - start_ind + 1;
    acell_lf3{i} = atime;
    if atime ~= n
    atimes_lf3 = [atimes_lf3 atime];
    end
end

%Presence times (residence times)
ldata_not_lf3 = ~ldata_lf3;
for i = 1:m
    sig = ldata_not_lf3(i,:);
    first_inds = find(sig);
    if length(first_inds) > 0
        first_ind = first_inds(1);
        sig(1:first_ind) = 1;
        last_ind = first_inds(end);
        sig(last_ind:end) = 1;
    else
        sig(1:end) = 1; 
    end
    dsig = diff([1 sig 1]);
    start_ind = find(dsig < 0);
    end_ind = find(dsig > 0) - 1;
    ptime = end_ind - start_ind + 1;
    pcell_lf3{i} = ptime;
    ptimes_lf3 = [ptimes_lf3 ptime];
end

% Mean residence/return time for each OTU
mean_ptimes_lf3 = [];
for i = 1:length(pcell_lf3)
   mean_ptime_lf3 = mean(pcell_lf3{i});
   mean_ptimes_lf3(i) = mean_ptime_lf3; 
end

mean_atimes_lf3 = [];
for i = 1:length(acell_lf3)   
   mean_atime_lf3 = mean(acell_lf3{i});
   mean_atimes_lf3(i) = mean_atime_lf3;
end




%% Combine diets

ptimes_hf = [ptimes_hf1 ptimes_hf2 ptimes_hf3];
atimes_hf = [atimes_hf1 atimes_hf2 atimes_hf3];

ptimes_lf = [ptimes_lf1 ptimes_lf2 ptimes_lf3];
atimes_lf = [atimes_lf1 atimes_lf2 atimes_lf3];


%% Parameter fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Human A

%params_a_A = fit_plcut_discrete(atimes_A,2);
%alpha_a_A = params_a_A.alpha_est;
alpha_a_A = 1.1419;
%lambda_a_A = params_a_A.lambda_est;
lambda_a_A = .0202;
%params_p_A = fit_plcut_discrete(ptimes_A,2);
%alpha_p_A = params_p_A.alpha_est;
alpha_p_A = 2.328;
%lambda_p_A = params_p_A.lambda_est;
lambda_p_A = .0148;


%Human B

%params_a_B = fit_plcut_discrete(atimes_B,2);
%alpha_a_B = params_a_B.alpha_est;
alpha_a_B = 1.1544;
%lambda_a_B = params_a_B.lambda_est;
lambda_a_B = .0225;
%params_p_B = fit_plcut_discrete(ptimes_B,2);
%alpha_p_B = params_p_B.alpha_est;
alpha_p_B = 2.2059;
%lambda_p_B = params_p_B.lambda_est;
lambda_p_B = .0057;


%Human M3

%params_a_m3 = fit_plcut_discrete(atimes_m3,2);
%alpha_a_m3 = params_a_m3.alpha_est;
alpha_a_m3 = 1.1910;
%lambda_a_m3 = params_a_m3.lambda_est;
lambda_a_m3 = .0247;
%params_p_m3 = fit_plcut_discrete(ptimes_m3,2);
%alpha_p_m3 = params_p_m3.alpha_est;
alpha_p_m3 = 2.2248;
%lambda_p_m3 = params_p_m3.lambda_est;
lambda_p_m3 = .0057;


%Human F4

%params_a_f4 = fit_plcut_discrete(atimes_f4,2);
%alpha_a_f4 = params_a_f4.alpha_est;
alpha_a_f4 = 1.0913;
%lambda_a_f4 = params_a_f4.lambda_est;
lambda_a_f4 = .0471;
%params_p_f4 = fit_plcut_discrete(ptimes_f4,2);
%alpha_p_f4 = params_p_f4.alpha_est;
alpha_p_f4 = 2.1472;
%lambda_p_f4 = params_p_f4.lambda_est;
lambda_p_f4 = .0327;


%HFHS Mice

%params_a_hf = fit_plcut_discrete(atimes_hf,2);
%alpha_a_hf = params_a_hf.alpha_est;
alpha_a_hf = .6799;
%lambda_a_hf = params_a_hf.lambda_est;
lambda_a_hf = .1817;
%params_p_hf = fit_plcut_discrete(ptimes_hf,2);
%alpha_p_hf = params_p_hf.alpha_est;
alpha_p_hf = 2.1944;
%lambda_p_hf = params_p_hf.lambda_est;
lambda_p_hf = .1234;


%LFPP Mice

%params_a_lf = fit_plcut_discrete(atimes_lf,2);
%alpha_a_lf = params_a_lf.alpha_est;
alpha_a_lf = .7243;
%lambda_a_lf = params_a_lf.lambda_est;
lambda_a_lf = .1596;
%params_p_lf = fit_plcut_discrete(ptimes_lf,2);
%alpha_p_lf = params_p_lf.alpha_est;
alpha_p_lf = 2.2033;
%lambda_p_lf = params_p_lf.lambda_est;
lambda_p_lf = .1139;


%HFHS 1

%params_a_hf1 = fit_plcut_discrete(atimes_hf1,2);
%alpha_a_hf1 = params_a_hf1.alpha_est;
alpha_a_hf1 = .7142;
%lambda_a_hf1 = params_a_hf1.lambda_est;
lambda_a_hf1 = .1767;
%params_p_hf1 = fit_plcut_discrete(ptimes_hf1,2);
%alpha_p_hf1 = params_p_hf1.alpha_est;
alpha_p_hf1 = 2.2023;
%lambda_p_hf1 = params_p_hf1.lambda_est;
lambda_p_hf1 = .1136;


%HFHS 2

%params_a_hf2 = fit_plcut_discrete(atimes_hf2,2);
%alpha_a_hf2 = params_a_hf2.alpha_est;
alpha_a_hf2 = .7057;
%lambda_a_hf2 = params_a_hf2.lambda_est;
lambda_a_hf2 = .1670;
%params_p_hf2 = fit_plcut_discrete(ptimes_hf2,2);
%alpha_p_hf2 = params_p_hf2.alpha_est;
alpha_p_hf2 = 2.2077;
%lambda_p_hf2 = params_p_hf2.lambda_est;
lambda_p_hf2 = .1514;


%HFHS 3

%params_a_hf3 = fit_plcut_discrete(atimes_hf3,2);
%alpha_a_hf3 = params_a_hf3.alpha_est;
alpha_a_hf3 = .6006;
%lambda_a_hf3 = params_a_hf3.lambda_est;
lambda_a_hf3 = .2076;
%params_p_hf3 = fit_plcut_discrete(ptimes_hf3,2);
%alpha_p_hf3 = params_p_hf3.alpha_est;
alpha_p_hf3 = 2.1518;
%lambda_p_hf3 = params_p_hf3.lambda_est;
lambda_p_hf3 = .1165;


%LFPP 1

%params_a_lf1 = fit_plcut_discrete(atimes_lf1,2);
%alpha_a_lf1 = params_a_lf1.alpha_est;
alpha_a_lf1 = .7569;
%lambda_a_lf1 = params_a_lf1.lambda_est;
lambda_a_lf1 = .1630;
%params_p_lf1 = fit_plcut_discrete(ptimes_lf1,2);
%alpha_p_lf1 = params_p_lf1.alpha_est;
alpha_p_lf1 = 2.2487;
%lambda_p_lf1 = params_p_lf1.lambda_est;
lambda_p_lf1 = .1078;


%LFPP 2

%params_a_lf2 = fit_plcut_discrete(atimes_lf2,2);
%alpha_a_lf2 = params_a_lf2.alpha_est;
alpha_a_lf2 = .6897;
%lambda_a_lf2 = params_a_lf2.lambda_est;
lambda_a_lf2 = .1615;
%params_p_lf2 = fit_plcut_discrete(ptimes_lf2,2);
%alpha_p_lf2 = params_p_lf2.alpha_est;
alpha_p_lf2 = 2.1698;
%lambda_p_lf2 = params_p_lf2.lambda_est;
lambda_p_lf2 = .1329;


%LFPP 3

%params_a_lf3 = fit_plcut_discrete(atimes_lf3,2);
%alpha_a_lf3 = params_a_lf3.alpha_est;
alpha_a_lf3 = .7237;
%lambda_a_lf3 = params_a_lf3.lambda_est;
lambda_a_lf3 = .1550;
%params_p_lf3 = fit_plcut_discrete(ptimes_lf3,2);
%alpha_p_lf3 = params_p_lf3.alpha_est;
alpha_p_lf3 = 2.1827;
%lambda_p_lf3 = params_p_lf3.lambda_est;
lambda_p_lf3 = .1041;




%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = figure;

%Human A
subplot(1,3,1);

plcut_atimes_A = @(x) 1/polyg(alpha_a_A,exp(-lambda_a_A)) * x.^(-alpha_a_A) .* exp(-lambda_a_A .* x);
plcut_ptimes_A = @(x) 1/polyg(alpha_p_A,exp(-lambda_p_A)) * x.^(-alpha_p_A) .* exp(-lambda_p_A .* x);
[counts_a_A,bins_a_A] = hist(atimes_A,min(atimes_A):max(atimes_A));
[counts_p_A,bins_p_A] = hist(ptimes_A,min(ptimes_A):max(ptimes_A));

p1a = plot(bins_a_A,plcut_atimes_A(bins_a_A),'-');
p1a.Color = c('matlab maroon');
p1a.LineWidth = 1.5;
hold on
p1b = plot(bins_a_A,counts_a_A/sum(counts_a_A),'.');
p1b.MarkerSize = 22;
p1b.Color = c('gray');
p1b.MarkerFaceColor = c('gray');
hold on
p1c = plot(min(ptimes_A):max(ptimes_A),plcut_ptimes_A(min(ptimes_A):max(ptimes_A)),'-');
p1c.Color = c('matlab maroon');
p1c.LineWidth = 1.5;
hold on
p1d = plot(bins_p_A,counts_p_A/sum(counts_p_A),'.');
p1d.MarkerSize = 22;
p1d.Color = c('black');
p1d.MarkerFaceColor = c('black');
hold off
set(gca,'xlim',[.8,500])
set(gca,'XTick',[10^(0),10^(.5),10^(1),10^(1.5),10^(2)])
set(gca,'XTickLabel',{'10^{0}','','10^{1}','','10^{2}'})
set(gca,'xscale','log');
set(gca,'XMinorTick','off');
set(gca,'ylim',[10^(-5.25),5])
set(gca,'YTick',[1e-5,1e-4,1e-3,1e-2,1e-1,1]);
set(gca,'YTickLabel',{'','10^{-4}','','10^{-2}','','10^{0}'});
set(gca,'YMinorTick','off');
set(gca,'yscale','log');
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Residence or return time');
ylabel('Probability');
title('Human A','FontSize',16,'FontWeight','Normal');
box off




% Human M3
subplot(1,3,2);

plcut_atimes_m3 = @(x) 1/polyg(alpha_a_m3,exp(-lambda_a_m3)) * x.^(-alpha_a_m3) .* exp(-lambda_a_m3 .* x);
plcut_ptimes_m3 = @(x) 1/polyg(alpha_p_m3,exp(-lambda_p_m3)) * x.^(-alpha_p_m3) .* exp(-lambda_p_m3 .* x);
[counts_a_m3,bins_a_m3] = hist(atimes_m3,min(atimes_m3):max(atimes_m3));
[counts_p_m3,bins_p_m3] = hist(ptimes_m3,min(ptimes_m3):max(ptimes_m3));

p2a = plot(bins_a_m3,plcut_atimes_m3(bins_a_m3),'-');
p2a.Color = c('matlab blue');
p2a.LineWidth = 1.5;
hold on
p2b = plot(bins_a_m3,counts_a_m3/sum(counts_a_m3),'.');
p2b.MarkerSize = 22;
p2b.Color = c('gray');
p2b.MarkerFaceColor = c('gray');
hold on
p2c = plot(min(ptimes_m3):max(ptimes_m3),plcut_ptimes_m3(min(ptimes_m3):max(ptimes_m3)),'-');
p2c.Color = c('matlab blue');
p2c.LineWidth = 1.5;
hold on
p2d = plot(bins_p_m3,counts_p_m3/sum(counts_p_m3),'.');
p2d.MarkerSize = 22;
p2d.Color = c('black');
p2d.MarkerFaceColor = c('black');
hold off
set(gca,'xlim',[.8,500])
set(gca,'XTick',[10^(0),10^(.5),10^(1),10^(1.5),10^(2)])
set(gca,'XTickLabel',{'10^{0}','','10^{1}','','10^{2}'})
set(gca,'xscale','log');
set(gca,'XMinorTick','off');
set(gca,'ylim',[10^(-4.95),5])
set(gca,'YTick',[1e-5,1e-4,1e-3,1e-2,1e-1,1]);
set(gca,'YTickLabel',{'','10^{-4}','','10^{-2}','','10^{0}'});
set(gca,'YMinorTick','off');
set(gca,'yscale','log');
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Residence or return time');
title('Human M3','FontSize',16,'FontWeight','Normal');
box off



% Residence/return times for LFPP mice
subplot(1,3,3);


plcut_atimes_lf = @(x) 1/polyg(alpha_a_lf,exp(-lambda_a_lf)) * x.^(-alpha_a_lf) .* exp(-lambda_a_lf .* x);
plcut_ptimes_lf = @(x) 1/polyg(alpha_p_lf,exp(-lambda_p_lf)) * x.^(-alpha_p_lf) .* exp(-lambda_p_lf .* x);
[counts_a_lf,bins_a_lf] = hist(atimes_lf,min(atimes_lf):max(atimes_lf));
[counts_p_lf,bins_p_lf] = hist(ptimes_lf,min(ptimes_lf):max(ptimes_lf));

p3a = plot(bins_a_lf,plcut_atimes_lf(bins_a_lf),'-');
p3a.Color = c('matlab green');
p3a.LineWidth = 1.5;
hold on
p3b = plot(bins_a_lf,counts_a_lf/sum(counts_a_lf),'.');
p3b.MarkerSize = 22;
p3b.Color = c('gray');
p3b.MarkerFaceColor = c('gray');
hold on
p3c = plot(min(ptimes_lf):max(ptimes_lf),plcut_ptimes_lf(min(ptimes_lf):max(ptimes_lf)),'-');
p3c.Color = c('matlab green');
p3c.LineWidth = 1.5;
hold on
p3d = plot(bins_p_lf,counts_p_lf/sum(counts_p_lf),'.');
p3d.MarkerSize = 22;
p3d.Color = c('black');
p3d.MarkerFaceColor = c('black');
hold off
set(gca,'xlim',[.8,50])
set(gca,'XTick',[10^(0),10^(.5),10^(1),10^(1.5),10^(2)])
set(gca,'XTickLabel',{'10^{0}','','10^{1}','','10^{2}'})
set(gca,'xscale','log');
set(gca,'XMinorTick','off');
set(gca,'ylim',[10^(-4.2),5])
set(gca,'YTick',[1e-4,1e-3,1e-2,1e-1,1]);
set(gca,'YTickLabel',{'10^{-4}','','10^{-2}','','10^{0}'});
set(gca,'YMinorTick','off');
set(gca,'yscale','log');
set(gca,'FontSize',13);
set(gca,'LineWidth',1);
set(gca,'FontName','Arial');
xlabel('Residence or return time');
title('Mice (LFPP Diet)','FontSize',16,'FontWeight','Normal');
box off
ll3 = legend([p3d,p3b],{'Residence times (t_{res })','Return times (t_{ret })'});
ll3.FontSize = 12;
ll3.LineWidth = 1;



