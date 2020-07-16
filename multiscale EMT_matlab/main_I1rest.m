%%
clear all
%4 steady states for TGF0=0.2,k0O =0.35
% k0O =0.35;
% i1=[0.1061650568352754, 1.035058126788898, 0.0029161627276987684, 0.0026592443269308663, 0.009362293186243556, 0.03582472960609682, 0.2596290968433698, 0.0007542737619013734, 0.0001530683818322264, 3.106273674315562e-5, 6.303626688642198e-6, 1.2791933126240401e-6, 0.06082067886128277, 1.2422605158994386, 0.04879770015972695, 201.36829994112483, 135.25103547317727, 1.548101672206943];
% epi = [0.023794182494844874, 0.04626328687111362, 0.060599959272124235, 0.019167853807733514, 0.008382951493345918, 0.028700893177463348, 0.28676680164986645, 0.0006848546180541147, 0.00015751185985241303, 3.6226646269333935e-5, 8.331879905126656e-6, 1.916275170418355e-6, 0.059574747835966384, 1.239953236733271, 0.04893613912933707, 446.2040639679741, 100.81667306882437, 1.5487808823635656];
% i2 =[0.13304383688086044, 1.305491548438854, 0.002685773433924973, 0.0024946820369750457, 0.029309071317134866, 0.2700702703164508, 0.04254131960560801, 0.00025323718470898396, 2.288706320469354e-6, 2.068486358895298e-8, 1.8694560235486306e-10, 1.6895764426729076e-12, 0.21028155749434824, 1.51903992128583, 0.03219093805618353, 156.98572754269918, 154.71450280831382, 1.450870110973855];
% mes = [0.2041827148056699, 2.0185991737933886, 0.002437867196322382, 0.002322797426331047, 0.16579211291162238, 1.5961125944685055, 0.002684357302410029, 3.906780888546334e-5, 9.216934033112352e-9, 2.1744724210100256e-12, 5.130046816809946e-16, 1.210288072102617e-19, 0.47972378834333135, 2.0180070154506136, 0.0022529124062965174, 104.14666510719147, 209.78313874068004, 0.6394946144420348];

TGF0=0.18;k0O =0.35;
i1 = [0.104505, 1.01831, 0.00293619, 0.00267363, 0.00932251, 0.0354959, 0.260837, 0.000751749, 0.000153455, 3.13247e-5, 6.39432e-6, 1.30527e-6, 0.0607589, 1.24215, 0.0488046, 201.485, 135.094, 1.54814];
mes = [0.203673, 2.0135, 0.00243876, 0.00232337, 0.165027, 1.58873, 0.00269959, 3.91302e-5, 9.28936e-9, 2.20526e-12, 5.23519e-16, 1.24281e-19, 0.479603, 2.01778, 0.00226636, 104.181, 209.681, 0.641536];
i2 = [0.130952, 1.28448, 0.00269865, 0.00250373, 0.02894, 0.266073, 0.0435972, 0.000261974, 2.48615e-6, 2.35937e-8, 2.23906e-10, 2.12488e-12, 0.205279, 1.50978, 0.0327468, 157.779, 154.228, 1.45355];
epi =  [0.0226754, 0.0365571, 0.071047, 0.0190197, 0.00838248, 0.0286987, 0.286771, 0.000684817, 0.000157506, 3.62262e-5, 8.33194e-6, 1.91633e-6, 0.0595746, 1.23995, 0.0489362, 463.727, 100.62, 1.54878];

epsilon_m = 0.01*ones(1,18);
epsilon_m(15)=0;


tstart = 0; tend = 5000; dt = 50;
saveat=dt;
time = tstart:dt:tend;
ntstep = length(time);
%initial 200 E cells, max 10000 cells
n_cells = 8000; n_0 = 200;
% track times the cell can divide before terminally differentiated (TD) cell
Div = zeros(n_cells,1);
Div(1:n_0) = randi([2 7],1,n_0); %each cell can divive 2-7 times brfore TD
Cells = zeros(n_cells,18,ntstep); %store gene expression of each cell at each time
Track = zeros(n_cells,ntstep,20); %Track cell,time,tree
tdiv = zeros(n_cells,1); %division time
time_div = zeros(n_cells,1); %time from beginning to the time of last division
time_TD = zeros(n_cells,1);  %keep track of TD cells expression with respect to time
cell_p = zeros(4,ntstep);  %get the population of different cell types at different time pts


%initialize
Cells(1:50,:,1) = epi.*ones(50,18);
Cells(51:100,:,1) = i1.*ones(50,18);
Cells(101:150,:,1) = i2.*ones(50,18);
Cells(151:200,:,1) = mes.*ones(50,18);

Track(1:n_0,1,1) = (1:n_0)'; %track cell lineage

tdiv(1:n_0) =  normrnd(700,200,1,n_0); %division time

%find the time when cell starting to divide
t_div_next = time_div+tdiv;
fastest_div_time = min(t_div_next(t_div_next>0));
i_div = 1; %before the first round of division

% I1 cells do not divide
time_TD(51:100) = 0;
Div(51:100) = 0;



tic
while (fastest_div_time <= tend) && (i_div <= ntstep)
    i_div
    old_i_div = i_div;
    i_div = find(time>=fastest_div_time,1); %first divide time step
    % existing cells before this division
    ind_r = find(sum(Cells(:,:,old_i_div),2)>0);
    if (length(ind_r) == n_cells) || (isempty(ind_r))
        break
    end
    
    %find the division cell within existing cells
    ind_dr = find(tdiv(ind_r)+time_div(ind_r)<=time(i_div));
    ind_d = ind_r(ind_dr);
    
    %find TD cell before this division
    ind_TD = find(time_TD>0);
    
    %existing cells but not TD cells
    ind_c = setdiff(ind_r,ind_TD);
    %what happens before division
    Cells(ind_c,:,old_i_div+1:i_div) = EMT(Cells(ind_c,:,old_i_div),TGF0,k0O,0,time(i_div)-time(old_i_div),epsilon_m,saveat);
    for i = old_i_div+1:i_div %cell expression before division
        Track(ind_c,i,:) = Track(ind_c,i-1,:);
    end
    
    %for division cells
    for i = 1:length(ind_d)
        if ind_r(end)+i > n_cells
            flag =1;
            break
        elseif Div(ind_d(i))>0 %has the ability to divide
            %y_temp= Cells(ind_d(i),:,i_div);
            %Cells(ind_r(end)+i,:,i_div) = y_temp; % daugther cell
            %Track
            a = find(Track(ind_d(i),i_div,:)==0,1);
            Track(ind_d(i),i_div,a) = Track(ind_d(i),i_div,1)+0.1;
            Track(ind_r(end)+i,i_div,:) = Track(ind_d(i),i_div,:);
            Track(ind_r(end)+i,i_div,a) = Track(ind_d(i),i_div,1)+0.2;
            
            %noise at division
            noise_coef = normrnd(0,0.7,1,18);
            noise_coef(15) = 0;
            for ii =1:18
                while abs(noise_coef(ii)) >0.5
                    noise_coef(ii) = normrnd(0,0.7);
                end
            end
            noise_div = Cells(ind_d(i),:,i_div).*noise_coef;
            % daughter cells after noise added
            Cells(ind_d(i),:,i_div) = Cells(ind_d(i),:,i_div) + noise_div;
            Cells(ind_r(end)+i,:,i_div) = Cells(ind_d(i),:,i_div) - noise_div;
            % update how many times left where cell can divide brfore TD
            Div(ind_d(i)) = Div(ind_d(i)) - 1;
            Div(ind_r(end)+i) = Div(ind_d(i));
            %update time from begin to most recent division
            time_div(ind_d(i)) = time(i_div);
            time_div(ind_r(end)+i) = time(i_div);
            
            if Div(ind_d(i))> 0 %still can divide after this division
                
                d1= Cells(ind_d(i),:,i_div);
                d2 =Cells(ind_r(end)+i,:,i_div);
                a = [vecnorm(d1'-epi')',vecnorm(d1'-i1')',vecnorm(d1'-i2')',vecnorm(d1'-mes')'];
                b = [vecnorm(d2'-epi')',vecnorm(d2'-i1')',vecnorm(d2'-i2')',vecnorm(d2'-mes')'];
                [~,a1]=min(a, [], 2);
                [~,a2]=min(b, [], 2);
                if a1==2 %daugher cells are I1 then don't divide any more
                    time_TD(ind_d(i)) = 0;
                    Div(ind_d(i)) = 0;
                    tdiv(ind_d(i)) = 0;
                elseif a2==2
                    time_TD(ind_r(end)+i) = 0;
                    Div(ind_r(end)+i) = 0;
                    tdiv(ind_r(end)+i) = 0;
                else
                    tdiv(ind_d(i)) = normrnd(700,200);
                    tdiv(ind_r(end)+i) = normrnd(700,200);
                end
                    
            else %daughter cells are TD, no division
                tdiv(ind_d(i)) = 0;
                tdiv(ind_r(end)+i) = 0;
                time_TD(ind_d(i)) = time(i_div);
                time_TD(ind_r(end)+i) = time(i_div); %when become TD cells
            end
            t_div_next(ind_d(i)) = time_div(ind_d(i)) + tdiv(ind_d(i));
            t_div_next(ind_r(end)+i) = time_div(ind_r(end)+i) + tdiv(ind_r(end)+i);
            
        elseif Div(ind_d(i))<=0 %no ability to divide
            continue
        end
    end %deal with division cells
    
    index = find(tdiv>0 & Div>0); %cells still can divide
    if isempty(index)
        break
    else
        fastest_div_time = min(t_div_next(index));
    end
end
toc
%% deal with TD cells
%existing cells
%ind_r = find(sum(Cells(:,:,i_div),2)>0);
%ind_div0 = find(Div<=0);
tic
time = tstart:dt:tend;
ind_TD = find(time_TD>0); %TD cells
death_rate = normrnd(700,100,n_cells,1);
% find dead cells
ind_dead = ind_TD(find(time_TD(ind_TD)+death_rate(ind_TD)<time(end)));
for i = 1:length(ind_dead)
    death_ind = find(time>=time_TD(ind_dead(i))+death_rate(ind_dead(i)),1);
    TD_ind = find(time>=time_TD(ind_dead(i)),1);
    Cells(ind_dead(i),:,TD_ind+1:death_ind) = EMT(Cells(ind_dead(i),:,TD_ind),TGF0,k0O,0,time(death_ind)-time(TD_ind),epsilon_m,saveat);
    for ii = TD_ind+1:death_ind
        Track(ind_dead(i),ii,:) = Track(ind_dead(i),TD_ind,:);
    end
    %Cells(ind_dead(i),:,death_ind+1:end) = zeros(1,18); already set 0
end
%live TD cells
ind_live_TD = setdiff(ind_TD,ind_dead);
for i = 1:ind_live_TD
    TD_ind = find(time>=time_TD(ind_live_TD(i)),1);
    Cells(ind_live_TD(i),:,TD_ind+1:end) = EMT(Cells(ind_live_TD(i),:,TD_ind),TGF0,k0O,0,time(end)-time(TD_ind),epsilon_m,saveat);
    for ii = TD_ind+1:ntstep
        Track(ind_live_TD(i),ii,:) = Track(ind_live_TD(i),TD_ind,:);
    end
end

% figure out the population of each cell type during the whole time
for i = 1:ntstep
    ind_r = find(sum(Cells(:,:,i),2)>0);
    a = [vecnorm(Cells(ind_r,:,i)'-epi')',vecnorm(Cells(ind_r,:,i)'-i1')',vecnorm(Cells(ind_r,:,i)'-i2')',vecnorm(Cells(ind_r,:,i)'-mes')'];
    [~,b]=min(a, [], 2); %index of smallest index
    for ii = 1:4
        cell_p(ii,i) = length(find(b==ii));
    end
end
toc                    

save('I1rest.mat')            
    
   