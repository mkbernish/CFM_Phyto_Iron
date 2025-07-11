
%% just plotting temperature and light fields for supplementary figure
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('\Users\mathp\OneDrive\Desktop\gansett_data\sst\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\m_map1.4\m_map\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\m_mhw1.0-master\'));
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\"))
addpath(genpath("C:\Users\mathp\OneDrive\Documents\paper1\Codes"));
% gridded iron concentration
load('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\new_maggie\Data\dfe\dfe_n_OCIM.mat')
load('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\new_maggie\Data\CTL.mat')
load("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\NO3_Nh4\data\surf_no3_clim.mat")
load('woaptemp.mat')
sst = mean(ptemp(:,:,1,:),4,'omitnan');
fs = 12;
load('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\new_maggie\uploadtodrive\light_new.mat', 'I_OCIM')
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\light_new_NEW.mat")
grid1 = output.grid;
%load transport_v4.mat  %load OCIM1 to get its grid info
lon1 = grid1.XT3d; % XT3d holds lon values for each point in OCIM2 grid
lat1 = grid1.YT3d; % YT3d holds lat values for each point in OCIM2 grid
z = grid1.zt';
lon = lon1(1,:,1);lat = lat1(:,1,1);
fn2 = 'Monthly_dFe_V2.nc';
dfe_n = ncread(fn2,'dFe_RF');lat_n =  ncread(fn2,'Latitude');
lon_n =  ncread(fn2,'Longitude'); z_n = ncread(fn2,'Depth');

depth = 3;
% plotting both
figure('Position',[1 1 900 700])
[ha,pos]=tight_subplot(2,1,[0.02 0.02],[.05 .08],[.03 .12]);
lvl = linspace(0,28,28);
axes(ha(1))
m_proj('robinson','lat',[-88 80],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(mean(ptemp(:,:,1:depth,:),[3 4],'omitnan'),90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''});
cb1 = contourcbar('fontname','times new roman','fontsize',fs);
cb1.Position = [0.8    0.56   0.02    0.3];
cb1.Label.String = ['T (' char(176) 'C)'];
cmocean('thermal',length(lvl))
clim([lvl(1) lvl(end)])

axes(ha(2))
lvl = linspace(0,200,20);
m_proj('robinson','lat',[-88 80],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(mean(I_OCIM(:,:,1:depth),3,'omitnan'),90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12);
cb2 = contourcbar('fontname','times new roman','fontsize',fs);
cb2.Position = [0.8    0.12   0.02    0.3];
cmocean('solar',length(lvl))
clim([lvl(1) lvl(end)])
cb2.Label.String = ['I (\mumol m^{-1} s^{-1})'];

%"C:\Users\mathp\OneDrive\Documents\paper1\Figures"
cd("C:\Users\mathp\OneDrive\Documents\paper1\Figures")
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'temperature_and_light' '.jpeg']);
disp(['Saving figure for:' fig_name]);
export_fig(sprintf(fig_name));