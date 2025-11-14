%script to graph the desired ppm ranges of glycine and mI, as a function 
%of LB and Concentraion ratios
clearvars
clc
basedir = 'path2\Gly_mI_sims\';
scan_count = 4; %1=OFF 64, 2=OFF 68, 3 = SUM 64, 4=SUM 68


numpts = 8192; % Number of spectral points
spectral_width = 2000; % in Hz. This is how much of the frequency axis you're sampling in your acquisition
bw = linspace(-spectral_width/2, spectral_width/2, numpts);
gyromagnetic_ratio = 42.58;
field_strength = 3; % in Tesla.
ppm_axis = bw/(gyromagnetic_ratio * field_strength);
ppm_axis = ppm_axis + 4.68;
ppmaxis = ppm_axis';

%sorting the files down to scan type, LB, and high vs low mI
sim_files = dir(fullfile(basedir,'*.mat'));
sim_files = {sim_files.name};
off64files = sim_files(contains(sim_files,'off') & contains(sim_files, '64'));
off68files = sim_files(contains(sim_files,'off') & contains(sim_files, '68'));
sum64files = sim_files(contains(sim_files,'sum') & contains(sim_files, '64'));
sum68files = sim_files(contains(sim_files,'sum') & contains(sim_files, '68'));


highmI5hzstore = zeros(8192,10);
highmI10hzstore = highmI5hzstore;
lowmI5hzstore = highmI5hzstore;
lowmI10hzstore = highmI5hzstore;


scan = [1 21 41 61];%scan skips by 20
if scan_count ==1
    a1=0;
    a2=10;
elseif scan_count==2
    a1=20;
    a2=30;
elseif scan_count==3
    a1=40;
    a2=50;
else
    a1=60;
    a2=70;
end
for HimI = scan(scan_count):1:scan(scan_count)+9 %high mI
    idx= HimI-a1;
    if rem(HimI,2)==0 %5hz
        load(sim_files{HimI})
        highmI5hzstore(:,idx)=spec;

    else  %10hz

        load(sim_files{HimI})
        highmI10hzstore(:,idx)=spec;
    end

end

for LowmI = scan(scan_count)+10:scan(scan_count)+19 %low mI
    idx= LowmI-a2;
    if rem(LowmI,2)==0 %5hz
        load(sim_files{LowmI})
        lowmI5hzstore(:,(idx))=spec;
    else %10hz
        load(sim_files{LowmI})
        lowmI10hzstore(:,(idx))=spec;
    end
end

lowmI5hzstore(:,1:2:end-1)=[];
highmI5hzstore(:,1:2:end-1)=[];
highmI10hzstore(:,2:2:end)=[];
lowmI10hzstore(:,2:2:end)=[];
%swap back to order of increasing gly conc
lowmI5hzstore=lowmI5hzstore(:, [2 1 3 4 5]);
highmI5hzstore=highmI5hzstore(:, [2 1 3 4 5]);
highmI10hzstore=highmI10hzstore(:, [2 1 3 4 5]);
lowmI10hzstore=lowmI10hzstore(:, [2 1 3 4 5]);
%end

%
figname= {'OFF 64','OFF 68','SUM 64','SUM 68'};
ppmrange = 3217:3740;
testing = lowmI5hzstore;
frank= [-1000000/1.2 2.75e6];
scale_factor = 9e5;
test = figure();
for ll = 1:5
subplot(2,2,1)
plot(ppm_axis(ppmrange)', real(lowmI5hzstore(ppmrange,ll))+scale_factor*ll-1,'o-',MarkerSize=1);
set(gca, 'Xdir','reverse')
hold on
end 

title('3mM myo-Inositol')
xlabel('ppm');
ylabel('5 Hz','FontName','Arial',FontSize=15,FontWeight='bold')
yticks([])
set(gca,'box','off')
ax = gca;
ax.FontSize = 12; 


for ll = 1:5
subplot(2,2,2)
plot(ppm_axis(ppmrange)', real(highmI5hzstore(ppmrange,ll))+scale_factor*ll-1,'o-',MarkerSize=1);
set(gca, 'Xdir','reverse')
hold on
end 

title('10mM myo-Inositol')

xlabel('ppm');
yticks([])
set(gca,'box','off')
ax = gca;
ax.FontSize = 12; 


for ll = 1:5
subplot(2,2,3)
plot(ppm_axis(ppmrange)', real(lowmI10hzstore(ppmrange,ll))+scale_factor*ll-1,'o-',MarkerSize=1);
set(gca, 'Xdir','reverse')
hold on
end 

legend({'0mM','0.5mM','1mM','2mM','3mM'},"Location","north")
xlabel('ppm','FontName','Arial',FontSize=12);
ylabel('10 Hz','FontName','Arial',FontSize=15,FontWeight='bold')
yticks([])
set(gca,'box','off')
ax = gca;
ax.FontSize = 12; 


for ll = 1:5
subplot(2,2,4)
plot(ppm_axis(ppmrange)', real(highmI10hzstore(ppmrange,ll))+scale_factor*ll-1,'o-',MarkerSize=1);
set(gca, 'Xdir','reverse')
hold on
end 

xlabel('ppm');
yticks([])
set(gca,'box','off')
ax = gca;
ax.FontSize = 12; 

sgtitle(figname(scan_count),'FontName','Arial',FontSize=18,FontWeight='bold')
set(test, 'Position', [250 0 1000 1500])

