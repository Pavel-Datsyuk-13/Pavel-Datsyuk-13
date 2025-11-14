%RMS simulation script, was used to make figure 2.
folderPath = 'path2\Gly_mI_TE_series\';
cd(folderPath)

clearvars
clc

%% Unpacking the spectra - faster to run one at a time
tic
numpts = 8192;
numloops =length(60:2:88);
store = zeros(8192, numloops);
lb = 3; % Linebroadening in Hz (simulations have inherent 3hz )
time = (1:1:8192)/2000;  % Create a acquisition time vector divided by the spectral width
exp_lb=permute(exp(-(lb*pi*time)),[2 1]);  % Add the linebroadening onto the time vector
indx = 1;

% Pick the metabolite you wanna analyze; rather do this instead of
% making another loop

%met_name = 'on_off_Glc_alpha_te';
%met_name = 'on_off_Glc_beta_te';
%met_name = 'on_off_Gln_ch2ch2ch_te';
%met_name = 'on_off_Gln_nh2_te';
%met_name = 'on_off_Glu_te';
%met_name = 'on_off_Gly_te';
%met_name = 'on_off_Thr_te';
met_name = 'on_off_Ins_te';

for n = 60:2:88
    metab = [met_name, num2str(n), '.mat'];
    test = load(metab);
    %test = load(met_name);
    fields = fieldnames(test);
    offspec= test.(fields{2});
    Frank = zeros(8192, 1);
    x_store = zeros(numpts, 19);
    for x = 1:19
        y_store = zeros(numpts, 19);
        for y = 1:19
            in_off_mI = offspec{x, 1}{1, y}.fids.*exp_lb;
            this_mI_spec = fftshift(fft(in_off_mI));
            this_mI_spec = this_mI_spec(end:-1:1);
            spec_spatial_mI = permute(this_mI_spec,[2 1]);
            spec_spatial_mI = spec_spatial_mI';
            spec_spatial_mI = circshift(spec_spatial_mI, 880);
            spec_spatial_mI = spec_spatial_mI(end:-1:1); 
            y_store(:, y) = spec_spatial_mI(:,1);
        end
        y_store = y_store';
        x_store(:,x) = (sum(y_store));
        
    end
    x_store = x_store';
    store(:,indx) = sum(x_store);
    indx = indx + 1;

end
%%
%temporary storage for summation of same metabolites
if strcmp(met_name, 'on_off_Glc_alpha_te')
    Glc_a_store = store;

elseif strcmp(met_name,'on_off_Gln_ch2ch2ch_te')
    Gln_ch2ch2ch_store = store;
else
end

%save new output to a different folder
if ~exist('processed metabs', 'dir')
    mkdir('processed metabs');
end
cd(fullfile(pwd, 'processed metabs'))

% Scale the simulations to the concentration and save
met2scale = met_name(8:10);
store = scalemetab(store, met2scale);
%sum_store_glc= sum(sum(abs(store)))
save(met2scale, "store")
cd('..')
toc
%% Plot the spectrum from one TE

figure()
hold on
plot(ppm_axis, real(store(:, 3)), 'k') %the number is arbitrary (anything 1-15)
xlabel('ppm')
xlim([3.2 3.8])
set(gca, 'XDir', 'reverse');
grid on

%% compute RMS values for each metabolite
clc
ppm_axis = linspace(-3.1484, 12.5804, 8192);
metab_loop = {'Gly','Ins','Glc','Thr'};

%establish ppm range and echo times used
ppm_start = find(ppm_axis >= 3.5, 1);
ppm_stop = find(ppm_axis > 3.6, 1);
TEs= 60:2:88;
rms_value = zeros(4,15);


% loop for the metabs and compute the rms
for met_idx= 1:length(metab_loop)
    load([metab_loop{met_idx}, '.mat'])
    figure(met_idx)
    for jj = 1:15
        plot(ppm_axis(ppm_start:ppm_stop),real(store(ppm_start:ppm_stop, jj)))
        hold on
    end
    legend(string(num2cell(60:2:88)))
    title(metab_loop(met_idx))

    if met_idx ==1
        load('Ins.mat')
        plot(ppm_axis(ppm_start:ppm_stop),real(store(ppm_start:ppm_stop, end)))
        hold off
    end

    load([metab_loop{met_idx}, '.mat'])
    rms_value(met_idx,:) = rms(real(store(ppm_start:ppm_stop,:)));
    min_te = min(rms_value(met_idx));
    best_te = TEs(rms_value(met_idx,:) == min(rms_value(met_idx,:)));
    [~, best_indx] = min(rms_value(met_idx,:));
    fprintf('%d ms has lowest RMS %s\n', best_te, metab_loop{met_idx});
    check = sum(abs(store(ppm_start:ppm_stop,end)));
end

%scale to the global max and plot values
rms_plot=rms_value./(max(max(rms_value(1,:))));
figure()
for ii= 1:length(metab_loop)
    plot(TEs, rms_plot(ii,:), '-o', "MarkerFaceColor","auto")
    hold on
end
grid on
title('Scaled RMS values at 3.5-3.6 ppm')
xlabel('Echo time (ms)')
ylabel('Scaled RMS value')
hold on
ylim([0 1.7])
legend(metab_loop)
box off
grid off
