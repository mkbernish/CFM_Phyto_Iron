% sensitivity analysis - only fe, fe and light, fe and temp, fe-light-temp
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

load('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\new_maggie\uploadtodrive\light_new.mat', 'I_OCIM')
grid1 = output.grid;
%load transport_v4.mat  %load OCIM1 to get its grid info
lon1 = grid1.XT3d; % XT3d holds lon values for each point in OCIM2 grid
lat1 = grid1.YT3d; % YT3d holds lat values for each point in OCIM2 grid
z = grid1.zt';
lon = lon1(1,:,1);lat = lat1(:,1,1);
fn2 = 'Monthly_dFe_V2.nc';
dfe_n = ncread(fn2,'dFe_RF');lat_n =  ncread(fn2,'Latitude');
lon_n =  ncread(fn2,'Longitude'); z_n = ncread(fn2,'Depth');

%
% plotting mu of t. weissflogii
% constant T, I (20 deg C and 500 umol quanta s)
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\ynbest_light_dark_2.csv")
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\ynfbest_light_dark_3.csv")
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\afebest_light_dark_2.csv")
    afe=afebest_light_dark_2(3); %yn = ynbest_light_dark_2(3);
    ynp = ynbest_light_dark_2(3);ynf =  ynfbest_light_dark_3(3);
fe_frac = 1;
Qc=1.00*10^(-12)/12 ; 
Feunit = 1/Qc;
depth = 3;
I_const = ones(size(I_OCIM)).*500;
T_const = ones(size(ptemp)).*20;
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,ynp,ynf,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
    .*fe_frac./1000,mean(T_const(:,:,1:depth,:),[3 4],'omitnan'),mean(I_const(:,:,1:depth),3,'omitnan'));
%
figure('Position',[1 1 1100 500])
[ha,pos]=tight_subplot(2,2,[0.02 0.02],[.05 .08],[.03 .12]);
axes(ha(1))
m_proj('robinson','lat',[-76 70],'lon',[-180 180]);
m_pcolor(lon-180,lat,circshift(mu*86400,90,2))
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''});
colormap(gca,flipud(crameri('roma')))
clim([0 1.8])
m_text(-210,80,['a. ' 'Constant T, I'],'fontname','times new roman','fontsize',12)

axes(ha(2)) %constant T, variable I
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,ynp,ynf,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
    .*fe_frac./1000,mean(T_const(:,:,1:depth,:),[3 4],'omitnan'),mean(I_OCIM(:,:,1:depth),3,'omitnan'));
m_proj('robinson','lat',[-76 70],'lon',[-180 180]);
m_pcolor(lon-180,lat,circshift(mu.*86400,90,2))
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''},'yticklabel',{''});
colormap(gca,flipud(crameri('roma')))
clim([0 1.8])
m_text(-210,80,['b. ' 'Constant T, Variable I'],'fontname','times new roman','fontsize',12)



axes(ha(3)) %constant I, variable T
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,ynp,ynf,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
    .*fe_frac./1000,mean(ptemp(:,:,1:depth,:),[3 4],'omitnan'),mean(I_const(:,:,1:depth),3,'omitnan'));
m_proj('robinson','lat',[-76 70],'lon',[-180 180]);
m_pcolor(lon-180,lat,circshift(mu.*86400,90,2))
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'xtick',[-120 -60 0 60 120]);
colormap(gca,flipud(crameri('roma')))
clim([0 1.8])
m_text(-210,80,['c. ' 'Variable T, Constant I'],'fontname','times new roman','fontsize',12)


axes(ha(4)) %variable T and I
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,ynp,ynf,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
    .*fe_frac./1000,mean(ptemp(:,:,1:depth,:),[3 4],'omitnan'),mean(I_OCIM(:,:,1:depth),3,'omitnan'));
m_proj('robinson','lat',[-76 70],'lon',[-180 180]);
m_pcolor(lon-180,lat,circshift(mu.*86400,90,2))
hold on;
shading interp;
m_coast('patch',[0 0 0],'edgecolor','k');
m_grid('fontname','times new roman','fontsize',12,'yticklabel',{''},'xtick',[-120 -60 0 60 120]);
colormap(gca,flipud(crameri('roma')))
clim([0 1.8])
cb3=contourcbar('fontname','times new roman','fontsize',12);
cb3.Ruler.MinorTick = 'on';
cb3.Position = [0.89    0.24    0.017    0.5];
cb3.Label.String=['\mu (d^{-1})'];
m_text(-210,80,['d. ' 'Variable T, I'],'fontname','times new roman','fontsize',12)

% cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
% set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
% fig_name=([ 'pcal_compound_effects' '.jpeg']);
% disp(['Saving figure for:' fig_name]);
% export_fig(sprintf(fig_name));

