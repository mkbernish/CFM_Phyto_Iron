% global plot for all one oceanic/one neritic species
%% global t.oce allocation/stoichiometry
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
%%
% plotting mu of t. oce
t1 = readtable("C:\Users\mathp\OneDrive\Documents\paper1\Codes\optimized_params.xlsx");
afe = t1.a_Fe(3); Afepho = t1.AFe_pho(3); % 3 corresponds to T oceanica (look at excel file)
fe_frac = 1;
Qc=1.00*10^(-12)/12 ; 
Feunit = 1/Qc;
depth = 3;
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,Afepho,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
    .*fe_frac./1000,mean(ptemp(:,:,1:depth,:),[3 4],'omitnan'),mean(I_OCIM(:,:,1:depth),3,'omitnan'));

figure('Position',[1 1 880 740])
[ha,pos]=tight_subplot(3,2,[0.08 0.01],[.03 .08],[.04 .01]);
axes(ha(1))
lvl = linspace(0,2.4,24);
m_proj('robinson','lat',[-82 70],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(mu*86400,90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''});
colormap(gca,flipud(crameri('roma',length(lvl))))
clim([lvl(1) lvl(end)])
cb1=contourcbar('horizontal','fontname','times new roman','fontsize',fs);
cb1.Ruler.MinorTick = 'on';
cb1.Position = [0.145    0.923    0.26    0.012];
cb1.Label.String=['\mu (d^{-1})'];
m_text(30,-74,['a. ' '\it{T. oceanica}'],'fontname','times new roman','fontsize',12,'color','w')

axes(ha(3))
n_c_ratio = Ndna_const_plot+Nconst_protein_plot+Nphoto_plot+Nchl_plot+Nbiosynth_plot+Nrna_plot;
lvl = linspace(0.2,0.24,20);
m_proj('robinson','lat',[-82 70],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(n_c_ratio,90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''});
colormap(gca,cmocean('rain',length(lvl)))
clim([lvl(1) lvl(end)])
cb3=contourcbar('horizontal','fontname','times new roman','fontsize',fs);
cb3.Ruler.MinorTick = 'on';
cb3.Position = [0.145    0.6    0.26    0.012];
cb3.Label.String=['N : C (mol mol^{-1})'];
%
m_text(30,-74,['c. ' '\it{T. oceanica}'],'fontname','times new roman','fontsize',12,'color','w')
%
%
axes(ha(5))
lvl = linspace(0,2e-3,20);
m_proj('robinson','lat',[-82 70],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(qfe.*Feunit.*1e9,90,2),lvl,'linecolor','none')
hold on;
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xtick',[-120 -60 0 60 120]);
%colormap(gca,flipud(crameri('roma',length(lvl))))
colormap(gca,cmocean('rain',length(lvl)))
clim([lvl(1) lvl(end)])
%cmocean('haline')
clim([0 1e-3])
cb5=contourcbar('horizontal','fontname','times new roman','fontsize',fs);
cb5.Ruler.MinorTick = 'on';
cb5.Position = [0.145    0.277    0.26    0.012];
cb5.Label.String=['Fe : C (mol mol^{-1})'];
m_text(30,-74,['e. ' '\it{T. oceanica}'],'fontname','times new roman','fontsize',12,'color','w')
afe = t1.a_Fe(4); Afepho = t1.AFe_pho(4); % 3 corresponds to T weiss(look at excel file)
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,Afepho,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
    .*fe_frac./1000,mean(ptemp(:,:,1:depth,:),[3 4],'omitnan'),mean(I_OCIM(:,:,1:depth),3,'omitnan'));

axes(ha(2))
lvl = linspace(0,1.8,18);
m_proj('robinson','lat',[-82 70],'lon',[-180 180]);
iron_lim = mu./mumax;
m_contourf(lon-180,lat,circshift(mu*86400,90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''},'yticklabel',{''});
colormap(gca,flipud(crameri('roma',length(lvl))))
clim([lvl(1) lvl(end)])
cb2=contourcbar('horizontal','fontname','times new roman','fontsize',fs);
cb2.Ruler.MinorTick = 'on';
cb2.Position = [0.625   0.923    0.26    0.012];
cb2.Label.String=['\mu (d^{-1})'];
m_text(30,-74,['b. ' '\it{T. weissflogii}'],'fontname','times new roman','fontsize',12,'color','w')

axes(ha(4))
n_c_ratio = Ndna_const_plot+Nconst_protein_plot+Nphoto_plot+Nchl_plot+Nbiosynth_plot+Nrna_plot;
lvl = linspace(0.2,0.24,20);
m_proj('robinson','lat',[-82 70],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(n_c_ratio,90,2),lvl,'linecolor','none')
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''},'yticklabel',{''});
colormap(gca,cmocean('rain',length(lvl)))
clim([lvl(1) lvl(end)])
cb4=contourcbar('horizontal','fontname','times new roman','fontsize',fs);
cb4.Ruler.MinorTick = 'on';
cb4.Position = [0.625   0.6    0.26    0.012];
cb4.Label.String=['N : C (mol mol^{-1})'];
m_text(30,-74,['d. ' '\it{T. weissflogii}'],'fontname','times new roman','fontsize',12,'color','w')

axes(ha(6))
lvl = linspace(0,1e-3,20);
m_proj('robinson','lat',[-82 70],'lon',[-180 180]);
m_contourf(lon-180,lat,circshift(qfe.*Feunit.*1e9,90,2),lvl,'linecolor','none')
hold on;
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'yticklabel',{''},'xtick',[-120 -60 0 60 120]);
%colormap(gca,flipud(crameri('roma',length(lvl))))
%colormap(gca,flipud(crameri('roma',length(lvl))))
colormap(gca,cmocean('rain',length(lvl)))
cb6=contourcbar('horizontal','fontname','times new roman','fontsize',fs);
cb6.Ruler.MinorTick = 'on';
cb6.Position = [0.625    0.277    0.26    0.012];
cb6.Label.String=['Fe : C (mol mol^{-1})'];
clim([0 2e-3])
m_text(30,-74,['f. ' '\it{T. weissflogii}'],'fontname','times new roman','fontsize',12,'color','w')

% saving fig
% cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
% set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
% fig_name=([ 'fig6_NEW' '.jpeg']);
% disp(['Saving figure for:' fig_name]);
% export_fig(sprintf(fig_name));
% 
