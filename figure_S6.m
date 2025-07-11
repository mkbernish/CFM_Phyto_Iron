
% for supplementary figure 6: comparing two iron fields 
% top: Huang et al., 2022
% bottom: Steady-state (OCIM)
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('C:\Users\mathp\OneDrive\Desktop\FE_opt_Spyd\Optimization\'))
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes"))
addpath(genpath("C:\Users\mathp\OneDrive\Documents\paper1\Codes"));
%%
fs = 12;
load('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\new_maggie\Data\dfe\dfe_n_OCIM.mat')
load('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\new_maggie\Data\CTL.mat')
grid1 = output.grid;
%load transport_v4.mat  %load OCIM1 to get its grid info
lon1 = grid1.XT3d; % XT3d holds lon values for each point in OCIM2 grid
lat1 = grid1.YT3d; % YT3d holds lat values for each point in OCIM2 grid
z = grid1.zt';
lon = lon1(1,:,1);lat = lat1(:,1,1);
fnbp = 'C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\DFe_Pasquier_Holzer_2018.nc';
dfe_bp = ncread(fnbp,'DFe_TYP');
latbp = ncread(fnbp,'lat');lonbp = ncread(fnbp,'lon');
depth = 3;
surf_huang = circshift(mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan'),90,2);
surf_benoit = circshift(mean(dfe_bp(:,:,1:depth),3,'omitnan'),90,2);
% plotting both
figure('Position',[1 1 900 700])
[ha,pos]=tight_subplot(2,1,[0.02 0.02],[.05 .08],[.03 .12]);
lvl = linspace(0,1.2,12);
axes(ha(1))
m_proj('robinson','lat',[-88 80],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan'),90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''});
colormap(gca,cmocean('deep',length(lvl)))
cb1 = contourcbar('horizontal','fontname','times new roman','fontsize',fs);
cb1.Label.String=['[dFe] (nM)'];
clim([0 lvl(end)])
cb1.Position = [0.275    0.925    0.36    0.01];
m_text(lon(end)-180-150,-78,'Huang et al., 2022','fontname','times new roman','fontsize',fs,'color','w')

axes(ha(2))
m_proj('robinson','lat',[-88 80],'lon',[-180 180]);
m_contourf(lonbp-180,latbp,circshift(mean(dfe_bp(:,:,1:depth),3,'omitnan'),90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12);
colormap(gca,cmocean('deep',length(lvl)))
clim([0 lvl(end)])
m_text(lon(end)-180-150,-78,'AWESOME OCIM','fontname','times new roman','fontsize',fs,'color','w')
% saving fig
cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'comparing_iron_fields' '.jpeg']);
disp(['Saving figure for:' fig_name]);
export_fig(sprintf(fig_name));
