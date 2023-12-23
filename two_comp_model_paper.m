rng(2)
c_m = 10; %nF/mm^2
g_l = 0.003*1e1; % uS/mm^2
E_l_axon = -65; %mV
V_thresh_axon = -63.7;
thresh_to_rest = V_thresh_axon-E_l_axon;
refrac = 4; %ms
dt = .01; %ms

%%% shape of AP
tau_in = 0.3;
tau_out = 0.2;
tran = 0:dt:(7-refrac);
AP = 0.5e3*(exp(-tran/tau_in)-exp(-tran/tau_out));
AP(160:end) = linspace(AP(160),-thresh_to_rest,length(AP(160:end)));
AP_full = [V_thresh_axon+AP E_l_axon*ones(1,round(refrac/dt))];
AP_time = length(AP_full);

Total_time = 100000; %ms
T = Total_time/dt; %total time steps.
E_l_soma = -65; %mV
V_axon  = nan(3,T); %initialize V
V_axon(:,1) = E_l_axon; %mV

V_soma  = nan(3,T); %initialize V
V_soma(:,1) = E_l_soma; %mV

% Add a soma compartment with a high resistance that attenuates the input
% from the axon
radi_by_length = 0.02;
g_comp = radi_by_length/(6*1e-3); % uS/mm^2

I_ext_orig_all = 0.04+0.1*randn(1,T); %nA/mm^2
skip_every = 3/dt;
stim_length = skip_every;
I_ext_orig = zeros(1,T);
nr_curr_input = floor(T/skip_every);
for jj = 0:nr_curr_input
    I_ext_orig(jj*skip_every+1:jj*skip_every+stim_length) = I_ext_orig_all(jj+1);
end

add_for_excit = rand(1,length(I_ext_orig));

E_inh = -65; %mV
E_exc = 0;
%
v_axon_coef = (g_l+g_comp);
tau_mem_axon = c_m/v_axon_coef;

nr_nspk = nan(1,3);
time_nspk = {};
locs_nspk = {};
peaks_nspk = {}; peaks_nspk{1} = nan; peaks_nspk{2} = nan; peaks_nspk{3} = nan;
amp_nspk = {};
mp_nspk = {};
min_prom = 8;

%%
tic
nr_leng = size(V_soma,1);
for ii = 1:nr_leng
    
    if ii == 1
        g_inh = 0; 
    elseif ii == 2  || ii == 3
        g_inh = 0.4*1e0; 
    end
        
    if ii == 1 || ii == 2
        g_exc = 0;    
    elseif ii == 3
      g_exc = 0.0379*1e0; 
    end
    v_soma_coef = (g_l+g_comp+g_inh+g_exc);
    tau_mem_soma = c_m/v_soma_coef;

    sig_noise = 1/sqrt(dt/tau_mem_soma);  
    I_ext_use = sig_noise*I_ext_orig;

    time_nspk{ii} = [];
    cal_exclude = 0;    
    for t=2:T
        
        if ~ismember(t,cal_exclude)
            v_inf_axon = (g_l*E_l_axon+g_comp*V_soma(ii,t-1))/v_axon_coef;
            V_axon(ii,t) = v_inf_axon+exp(-dt/tau_mem_axon)*(V_axon(ii,t-1)-v_inf_axon);
            if V_axon(ii,t)>V_thresh_axon
                time_nspk{ii} = [time_nspk{ii} t];
                if T-t-1 >= AP_time
                    V_axon(ii,t:t+AP_time-1) = AP_full;
                    cal_exclude = t:t+AP_time-1;
                else
                    time_left = T-t+1;
                    V_axon(ii,t:end) = AP_full(1:time_left);
                    cal_exclude = t:T;
                end
            end
        end
        v_inf_soma = (g_l*E_l_soma+g_comp*V_axon(ii,t-1)+g_inh*E_inh+g_exc*E_exc+I_ext_use(t))/v_soma_coef;
        V_soma(ii,t) = v_inf_soma+exp(-dt/tau_mem_soma)*(V_soma(ii,t-1)-v_inf_soma);
        
    end
    [peaks_nspk{ii},locs_nspk{ii},~,~] = findpeaks(V_soma(ii,:),'MinPeakProminence',min_prom);
    nr_nspk(ii) = length(locs_nspk{ii});
    if length(time_nspk{ii})-length(locs_nspk{ii}) == 1
        time_nspk{ii}(end) = [];
    end
    amp_nspk{ii} = peaks_nspk{ii}-V_soma(ii,time_nspk{ii});
    mp_nspk{ii} = V_soma(ii,time_nspk{ii});

end

sprintf('nr nspk bl = %d. nr nspk inh = %d. nr nspk can = %d',nr_nspk(1),nr_nspk(2),nr_nspk(3))
sprintf('mean spk bl = %0.3f. mean spk inh = %0.3f. mean spk can = %0.3f',mean(peaks_nspk{1}),mean(peaks_nspk{2}),mean(peaks_nspk{3}))
toc

%% Figure plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot densities
bspk_thresh = quantile(peaks_nspk{1},0.97);

h = figure('units','inches','position',[6 3 4 2]);
for ii = 1:2
f_all = {};
xi_all = {};
if ii == 1
    data_plot = peaks_nspk;
elseif ii == 2
    data_plot = amp_nspk;
end
[f_all{1},xi_all{1}] = ksdensity(data_plot{1}); 
[f_all{2},xi_all{2}] = ksdensity(data_plot{2}); 
[f_all{3},xi_all{3}] = ksdensity(data_plot{3}); 

subplot(1,2,ii)
p1=plot(xi_all{1},f_all{1},'k','linewidth',1);
hold on
p2=plot(xi_all{2},f_all{2},'r','linewidth',1);
hold on
p3=plot(xi_all{3},f_all{3},'m','linewidth',1);
ax = gca;
y_lim = ax.YLim;
if ii == 1
    hold on
    line([bspk_thresh,bspk_thresh],y_lim,'LineStyle','--')
end
if ii == 1
    xlabel('peak (mV)')
    legend('baseline','inhibition','cancellation','AutoUpdate','off')
    legend('boxoff')

elseif ii == 2   
    xlabel('amp (mV)')
end
ylabel('density')

ax = gca;
ax.XLimSpec = 'tight';
ax.FontSize = 6;
ax.Color = 'none';
ax.XColor = 'k';
ax.YColor = 'k';
ax.Box = 'off';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;

end

%% trace example
h = figure('units','inches','position',[6 3 2 2]);

x_time = [3*1e4,3.6*1e4];
tt = linspace(1,(x_time(2)-x_time(1))*dt,length(V_axon(1,x_time(1):x_time(2))));

plot(tt,V_soma(1,x_time(1):x_time(2)), 'k','linewidth',1.5);
hold on
plot(tt,V_soma(2,x_time(1):x_time(2)), 'r','linewidth',1,'LineStyle','--');
hold on
plot(tt,V_soma(3,x_time(1):x_time(2)), 'c','linewidth',1,'LineStyle','--');
legend('without inh.','with inh.','AutoUpdate','off')
legend('boxoff')
hold on
ax = gca;
x_lim = ax.XLim;
line(x_lim,[bspk_thresh bspk_thresh],'linestyle','--','color',[0.5 0.5 0.5],'linewidth',1.5)
title('V soma')
ylabel('mV')
xlabel('time (ms)')
ax.FontSize = 6;
ax.Color = 'none';
ax.Box = 'off';
ax.XColor = 'k';
ax.YColor = 'k';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;

%% bar graphs
fs = 6;
nr_bspk = nan(1,3);
bspk_thresh = quantile(peaks_nspk{1},0.97);
for ii = 1:nr_leng
    nr_bspk(ii) = sum(peaks_nspk{ii}>bspk_thresh);    
end
%
h = figure('units','inches','position',[3 3 8.1 1.2]);
subplot(1,3,1)
hBar = bar(nr_bspk/(Total_time*1e-3));
hBar.FaceColor = 'flat';
hBar.CData(1,:) = [0 0 0]; hBar.CData(2,:) = [1 0 0]; hBar.CData(3,:) = [1 0 1];
ylabel('broad spike rate (Hz)')
ylim([0 2.2])
ax = gca;
ax.XTickLabel = {'initial','inhibition','cancellation'};
xtickangle(45)
ax.FontSize = 6;
ax.Color = 'none';
ax.Box = 'off';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;
ax.XColor = 'k';
ax.YColor = 'k';

subplot(1,3,2)
hBar = bar(nr_nspk/(Total_time*1e-3));
hBar.FaceColor = 'flat';
hBar.CData(1,:) = [0 0 0]; hBar.CData(2,:) = [1 0 0]; hBar.CData(3,:) = [1 0 1];
ylabel('narrow spike rate (Hz)')
ylim([0 100])
ax = gca;
ax.XTickLabel = {'initial','inhibition','cancellation'};
xtickangle(45)
ax.FontSize = 6;
ax.Color = 'none';
ax.Box = 'off';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;
ax.XColor = 'k';
ax.YColor = 'k';

subplot(1,3,3)
hBar = bar(nr_bspk./nr_nspk);
hBar.FaceColor = 'flat';
hBar.CData(1,:) = [0 0 0]; hBar.CData(2,:) = [1 0 0]; hBar.CData(3,:) = [1 0 1];
ylabel('p')
ax = gca;
ax.XTickLabel = {'initial','inhibition','cancellation'};
xtickangle(45)
ax.FontSize = 6;
ax.Color = 'none';
ax.Box = 'off';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;
ax.XColor = 'k';
ax.YColor = 'k';

%%
fs = 6;

mean_data_for_bar_graph = [nanmean(peaks_nspk{1}) nanmean(peaks_nspk{2}) nanmean(peaks_nspk{3})];
mean_peak_std = mean([nanstd(peaks_nspk{1}) nanstd(peaks_nspk{2}) nanstd(peaks_nspk{3})]);


h = figure('units','inches','position',[3 3 8.1 2]);
subplot(1,3,1)
hBar = bar(mean_data_for_bar_graph);
hBar.FaceColor = 'flat';
hBar.CData(1,:) = [0 0 0]; hBar.CData(2,:) = [1 0 0]; hBar.CData(3,:) = [1 0 1];
ylabel('peak (mV)')
ylim([-51.5, -50.1])
ax = gca;
x_lim = ax.XLim;
hold on
line(x_lim,[bspk_thresh,bspk_thresh],'LineStyle','--')
hold on
line(x_lim,[bspk_thresh-mean_peak_std,bspk_thresh-mean_peak_std],'LineStyle','--')
ax.XTickLabel = {'initial','inhibition','cancellation'};
xtickangle(45)
ax.FontSize = 6;
ax.Color = 'none';
ax.Box = 'off';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;
ax.XColor = 'k';
ax.YColor = 'k';

subplot(1,3,2)
mean_data_for_bar_graph = [nanmean(amp_nspk{1}) nanmean(amp_nspk{2}) nanmean(amp_nspk{3})];
hBar = bar(mean_data_for_bar_graph);
hBar.FaceColor = 'flat';
hBar.CData(1,:) = [0 0 0]; hBar.CData(2,:) = [1 0 0]; hBar.CData(3,:) = [1 0 1];
ylabel('amplitude (mV)')
ylim([12 12.5])
ax = gca;
ax.XTickLabel = {'initial','inhibition','cancellation'};
xtickangle(45)
ax.FontSize = 6;
ax.Color = 'none';
ax.Box = 'off';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;
ax.XColor = 'k';
ax.YColor = 'k';

subplot(1,3,3)
mean_data_for_bar_graph = [nanmean(mp_nspk{1}) nanmean(mp_nspk{2}) nanmean(mp_nspk{3})];
hBar = bar(mean_data_for_bar_graph);
hBar.FaceColor = 'flat';
hBar.CData(1,:) = [0 0 0]; hBar.CData(2,:) = [1 0 0]; hBar.CData(3,:) = [1 0 1];
ylim([-63.8, -62.8])
ylabel('baseline membrane potential')
ax = gca;
ax.XTickLabel = {'initial','inhibition','cancellation'};
xtickangle(45)
ax.FontSize = 6;
ax.Color = 'none';
ax.Box = 'off';
ax.FontName = 'Arial';
ax.LineWidth = 0.8;
ax.XColor = 'k';
ax.YColor = 'k';
