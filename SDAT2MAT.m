%%Runs through the SDAT files for all the subjects and creates a .MAT
%%spectra for each, used for analysis in MATLAB
%%NOTE: you will need to comment DIFF save out (or OFF and SUM) and change 
%%The exp LB when processing the data (change in GannetPreInitialise.m)

clearvars
cd ('path2\Volunteer_Data_Anon')
data_path = 'path2\Volunteer_Data_Anon\';
%format: each row is a patient and each column
%is for the scan (in the order of MP64, MP68)

SNR_values = zeros(6, 4); %store the SNR for the means



%extract all folders
folder_path = dir(fullfile(data_path));
all_folders_og = folder_path([folder_path.isdir]);
all_folders_og = all_folders_og(~ismember({all_folders_og.name}, {'.','..','PC_FC_testing'}));
all_folders_og = {all_folders_og.name};
all_folders = repelem(all_folders_og,1, 2); %can change depending on # of scans / subject. 2 in this case
all_folders = convertStringsToChars(all_folders);


%outerloop for per patient
 for pat_idx = 1:length(all_folders_og) %change start to 1 to 2,3,4 etc if you already processed them
    current_folder = all_folders{pat_idx*2};

    %this should automatically pull out the files you need
    current_path = fullfile(data_path, current_folder);
    all_files = dir(fullfile(current_path,'**', '*.SDAT')); % Get all .SDAT files in current/sub folder(s)
    file_names = {all_files.name}; % Extract file names

    %extracting files in this specific order
    %try to maintain this order: Edit 64, edit 68 
    data_files = {file_names(contains(file_names, 'Edit') & contains(file_names, 'TE_64') & contains(file_names, 'act')),...
        file_names(contains(file_names, 'Edit') & contains(file_names, 'TE_68') & contains(file_names, 'act'))};

    wat_ref_files = {file_names(contains(file_names, 'Edit') & contains(file_names, 'TE_64') & contains(file_names, 'ref')),...
        file_names(contains(file_names, 'Edit') & contains(file_names, 'TE_68') & contains(file_names, 'ref'))};

    %change directory for each patient 
    cd([data_path,'P',num2str(pat_idx)]);
    %changes folder name per patient
    folder_prefix = char(data_files{2});
    folder_prefix = [folder_prefix(1:8),'\'];

    for subj_idx = 1:length(data_files)
        all_folders_og=all_folders{subj_idx};
        specinput =  strcat(data_path, current_folder,'\', data_files{subj_idx});
        waterinput = strcat(data_path, current_folder,'\', wat_ref_files{subj_idx});

        MRS_struct = GannetLoad(specinput, waterinput);
        spec = MRS_struct.spec.AllFramesFTrealign;


       if subj_idx == 1
            editon64 = mean(spec(:, 2:2:end),2);
            editoff64 = mean(spec(:,1:2:end),2);
            sum64= editon64+editoff64;

            %edit off
            spec=editoff64;%makes the input smoother

            x =MRS_struct.spec.freq;
            %close gannet figure and the CW
            close all
            clc
            %apply the rest of the adjustments
            spec = spec(end:-1:1);
            mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            SNR_values(pat_idx, subj_idx)=crSNR(spec);

             namesave = 'OFF_64.mat';
            % save(namesave, 'spec')
            

            %Sum
            spec = sum64; %smooth input

            %establish gui inputs
            testspec=spec;
            x =MRS_struct.spec.freq;
           
            %apply the rest of the adjustments
            spec = spec(end:-1:1);
            mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            spec= spec./2;
            SNR_values(pat_idx, subj_idx+1)=crSNR(spec);
            namesave = 'SUM_64.mat';
            % save(namesave, 'spec')

            %DIFF spec
            editon64 = mean(spec(:, 2:2:end),2);
            editoff64 = mean(spec(:,1:2:end),2);
            DIFF64= editon64-editoff64;
            spec=DIFF64;%makes the input smoother
            testspec=spec;
            x =MRS_struct.spec.freq;
            clc
            %apply the rest of the adjustments:
            spec = spec(end:-1:1);
            mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            namesave = 'DIFF_64.mat';
            % save(namesave, 'spec')

            %water data
            spec = MRS_struct.spec.vox1.water;
            spec = spec(end:-1:1);
            water_mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(water_mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            namesave = 'MP_64_water.mat';
            % save(namesave, 'spec')

        else %68 ms

            %edit off
            editon68 = mean(spec(:, 2:2:end),2);
            editoff68 = mean(spec(:,1:2:end),2);
            sum68= editon68+editoff68;
            spec=editoff68;%makes the input smoother

            %establish gui inputs
            testspec=spec;
            x =MRS_struct.spec.freq;
           
            %apply the rest of the adjustments:
            spec = spec(end:-1:1);
            mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            SNR_values(pat_idx, subj_idx+1)=crSNR(spec);

            namesave = 'OFF_68.mat';
            % save(namesave, 'spec')


            %Sum
            spec = sum68; %smooth input

            %establish gui inputs
            testspec=spec;
            x =MRS_struct.spec.freq;
           
            %apply the rest of the adjustments
            spec = spec(end:-1:1);
            mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            spec= spec./2;
            namesave = 'SUM_68.mat';
            SNR_values(pat_idx, subj_idx+2)=crSNR(spec);
            % save(namesave, 'spec')

             %DIFF spec
            editon68 = mean(spec(:, 2:2:end),2);
            editoff68 = mean(spec(:,1:2:end),2);
            DIFF68= editon68-editoff68;
            spec=DIFF68;%makes the input smoother
            testspec=spec;
            x =MRS_struct.spec.freq;
            clc
            %apply the rest of the adjustments:
            spec = spec(end:-1:1);
            mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            namesave = 'DIFF_68.mat';
            % save(namesave, 'spec')

            % water data
            spec = MRS_struct.spec.vox1.water;
            spec = spec(end:-1:1);
            water_mcfid = permute(ifft(ifftshift(spec')),[2,1]);
            spec = fftshift(fft(water_mcfid(end:-1:1)));
            spec = permute(spec,[2,1]);
            namesave = 'MP_68_water.mat';
            % save(namesave, 'spec')
        end

    end

end


