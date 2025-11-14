clearvars
cd('path2\GABA_TE_series');
spec_files = dir('*.mat');
spec_files = {spec_files.name};

%% Innitialize axis and loop params
numpts = 8192; % Number of spectral points
spectral_width = 2000; % in Hz. This is how much of the frequency axis you're sampling in your acquisition
bw = linspace(-spectral_width/2, spectral_width/2, numpts);
gyromagnetic_ratio = 42.58;
field_strength = 3; % in Tesla.
ppm_axis = bw/(gyromagnetic_ratio * field_strength);
ppm_axis = ppm_axis + 4.68;
%Loop params
numfiles = length(spec_files);
on_store_allspec = zeros(8192, numfiles);
off_store_allspec = on_store_allspec;
diffspecs = on_store_allspec;
lb = 3; % Linebroadening in Hz CHANGED FOR NET LB = 6HZ
time = (1:1:8192)/2000;  % Create a acquisition time vector divided by the spectral width
exp_lb=permute(exp(-(lb*pi*time)),[2 1]);  % Add the linebroadening onto the time vector
te=60:2:80;
indx = 1; %To move from one file to the next

%% Loop is supposed to go through the files and store the output and sum/average to a 1x8192 array
tic

while indx <= numfiles % for variable number of spectra/metabolites
    this_metabolite = spec_files{indx}; %pulls out file name
    this_metabolite = load(this_metabolite); %loads it in workspace
    fields = fieldnames(this_metabolite);
    this_metabolite_off = this_metabolite.(fields{2});
    this_metabolite_on = this_metabolite.(fields{3});
    %Those two lines above were bc I couldn't figure out how to only inlude
    %The struct w spectra ans leave out the 'more params' field

    %doing on spec first
    x_store = zeros(numpts, 19);
    for x = 1:19
        y_store = zeros(numpts, 19);
        for y = 1:19

            in_off_metabolite = this_metabolite_on{x,1}{1,y}.fids.*exp_lb;
            this_metabolite_spec = fftshift(fft(in_off_metabolite));
            this_metabolite_spec = this_metabolite_spec(end:-1:1);
            spec_spatial_metabolite = permute(this_metabolite_spec,[2 1]);
            spec_spatial_metabolite = spec_spatial_metabolite';
            spec_spatial_metabolite = circshift(spec_spatial_metabolite, 880); 
            spec_spatial_metabolite = spec_spatial_metabolite(end:-1:1);
            y_store(:, y) = spec_spatial_metabolite(:,1);
        end
        y_store = y_store';
        on_spat_store(:,:,x)=y_store();
        x_store(:,x) = (sum(y_store));
    end

    on_Metabolite_spec = x_store';
    on_Metabolite_spec = sum(x_store, 2);
    on_store_allspec(:, indx) = on_Metabolite_spec;

    %same for opn spec
    x_store = zeros(numpts, 19);
    for x = 1:19
        y_store = zeros(numpts, 19);
        for y = 1:19

            in_off_metabolite = this_metabolite_off{x,1}{1,y}.fids.*exp_lb;
            this_metabolite_spec = fftshift(fft(in_off_metabolite));
            this_metabolite_spec = this_metabolite_spec(end:-1:1);
            spec_spatial_metabolite = permute(this_metabolite_spec,[2 1]);
            spec_spatial_metabolite = spec_spatial_metabolite';
            spec_spatial_metabolite = circshift(spec_spatial_metabolite, 880);
            spec_spatial_metabolite = spec_spatial_metabolite(end:-1:1); 
            y_store(:, y) = spec_spatial_metabolite(:,1);
        end
        y_store = y_store';
        off_spat_store(:,:,x)=y_store();
        x_store(:,x) = (sum(y_store));
    end

    off_Metabolite_spec = x_store';
    off_Metabolite_spec = sum(x_store, 2);
    off_store_allspec(:, indx) = off_Metabolite_spec;
    diffspecs(:,indx)= on_Metabolite_spec - off_Metabolite_spec;
    indx = indx + 1;
end

%getting diff spec
GABA_diffspecs = on_store_allspec-off_store_allspec;
toc

% AREA / TE ANALYSIS

GABA_area = zeros(1,size(diffspecs,2));
for aa = 1:size(diffspecs,2)
    GABA_area(1,aa) = sum(diffspecs(3102:3355,aa));
end

GABA_area = GABA_area./max(GABA_area);

%factoring T2-decay
S0 = 1; 
T2 = 88;
t2fun = @(te) S0 .* exp(te ./ T2);
GABA_area_corrected = GABA_area ./ t2fun(te);
GABA_area_corrected = GABA_area_corrected./max(GABA_area_corrected);

% Looking at area plots
figure()
plot(te, GABA_area, 'o-k')
hold on
plot(te, GABA_area_corrected, 'o-b')
ylim([0.3 1.3])
legend('No T2 correction','T2-corrected', 'Location', 'southeast')
title('GABA area TE-evolution')

