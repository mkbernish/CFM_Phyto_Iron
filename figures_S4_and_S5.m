

%% all 10 modeled species - global carbon allocation

addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('\Users\mathp\OneDrive\Desktop\gansett_data\sst\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\m_map1.4\m_map\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\m_mhw1.0-master\'));
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\"));
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

load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\ynpbest_light_dark_3.csv")
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\ynfbest_light_dark_3.csv")
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\afebest_light_dark_3.csv")
speciesforloop = {'\it E. huxleyi','\it P. calceolata',...
    '\it T. oceanica', '\it T. weissflogii','\it T. pseudonana','\it P. minimum','\it P. micans',...
    ['{\it F. cylindrus} (6 ' char(176) 'C)'],['{\it F. cylindrus} (3 ' char(176) 'C)'],['{\it F. cylindrus} (1 ' char(176) 'C)']};


% C allocation to photosynthesis 
figure('Position',[1 1 780 820])
[ha,pos]=tight_subplot(5,2,[0.01 0.02],[.03 .03],[.05 .11]);
fs = 12;
for i = 1 : 10
    afe=afebest_light_dark_3(i); %yn = ynbest_light_dark_2(3);
    ynp = ynpbest_light_dark_3(i);ynf =  ynfbest_light_dark_3(i);
        species_name = speciesforloop(i);
    fe_frac = 1;
    Qc=1.00*10^(-12)/12 ;
    Feunit = 1/Qc;
    depth = 4;
    [mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
        Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
        Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
        Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
        Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,ynp,ynf,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
        .*fe_frac./1000,mean(ptemp(:,:,1:depth,:),[3 4],'omitnan'),mean(I_OCIM(:,:,1:depth),3,'omitnan'));

    axes(ha(i))
    lvl = linspace(0,100,20);
    m_proj('robinson','lat',[-72 70],'lon',[-180 180]);
    m_contourf(lon-180,lat,circshift(Cpho,90,2),lvl,'linecolor','none')
    hold on;
    shading interp;
    m_coast('patch',[0 0 0],'edgecolor','k');
    colormap(flipud(crameri('roma',length(lvl))))
    clim([lvl(1) lvl(end)])
    j = i;
     if j ==4
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''},'yticklabel',{''})
    elseif j == 3
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''})
    elseif j == 2
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''},'yticklabel',{''})
    elseif j == 1
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''})
    elseif j == 5
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
          m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''})
    elseif j == 6
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
       m_grid('fontname','times new roman','fontsize',12,'xtick',[-120 -60 0 60 120],'xticklabel',{''},'yticklabel',{''})
    elseif j == 7
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
          m_grid('fontname','times new roman','fontsize',12,'xtick',[-120 -60 0 60 120],'xticklabel',{''})
    elseif j == 8
        m_text(-180,79,species_name,'fontname','times new roman','fontsize',fs)
         m_grid('fontname','times new roman','fontsize',12,'yticklabel',{''},'xtick',[-120 0 120],'xticklabel',{''})
    elseif j == 9
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xtick',[-120  0 120])
    else
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
         m_grid('fontname','times new roman','fontsize',12,'yticklabel',{''},'xtick',[-120 0 120])
    end

end
cb3=contourcbar('fontname','times new roman','fontsize',14);
cb3.Ruler.MinorTick = 'on';
cb3.Position = [ 0.9    0.3    0.015    0.4];
cb3.Label.String=['C allocation (%)'];
cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'supp_allCPho_glob_vert' '.jpeg']);
disp(['Saving figure for:' fig_name]);
export_fig(sprintf(fig_name));

%% C allocation to biosynthesis
figure('Position',[1 1 780 820])
[ha,pos]=tight_subplot(5,2,[0.01 0.02],[.03 .03],[.05 .11]);
fs = 12;
for i = 1 : 10
    afe=afebest_light_dark_3(i); %yn = ynbest_light_dark_2(3);
    ynp = ynpbest_light_dark_3(i);ynf =  ynfbest_light_dark_3(i);
        species_name = speciesforloop(i);
    fe_frac = 1;
    Qc=1.00*10^(-12)/12 ;
    Feunit = 1/Qc;
    depth = 4;
    [mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
        Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
        Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
        Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
        Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,ynp,ynf,mean(dfe_n_OCIM(:,:,1:depth),3,'omitnan')...
        .*fe_frac./1000,mean(ptemp(:,:,1:depth,:),[3 4],'omitnan'),mean(I_OCIM(:,:,1:depth),3,'omitnan'));

    axes(ha(i))
    lvl = linspace(0,100,20);
    m_proj('robinson','lat',[-72 70],'lon',[-180 180]);
    m_contourf(lon-180,lat,circshift(Cbio,90,2),lvl,'linecolor','none')
    hold on;
    shading interp;
    m_coast('patch',[0 0 0],'edgecolor','k');
    colormap(flipud(crameri('roma',length(lvl))))
    clim([lvl(1) lvl(end)])
    j = i;
     if j ==4
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''},'yticklabel',{''})
    elseif j == 3
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''})
    elseif j == 2
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''},'yticklabel',{''})
    elseif j == 1
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''})
    elseif j == 5
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
          m_grid('fontname','times new roman','fontsize',12,'xticklabel',{''})
    elseif j == 6
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
       m_grid('fontname','times new roman','fontsize',12,'xtick',[-120 -60 0 60 120],'xticklabel',{''},'yticklabel',{''})
    elseif j == 7
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
          m_grid('fontname','times new roman','fontsize',12,'xtick',[-120 -60 0 60 120],'xticklabel',{''})
    elseif j == 8
        m_text(-180,79,species_name,'fontname','times new roman','fontsize',fs)
         m_grid('fontname','times new roman','fontsize',12,'yticklabel',{''},'xtick',[-120 0 120],'xticklabel',{''})
    elseif j == 9
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
        m_grid('fontname','times new roman','fontsize',12,'xtick',[-120  0 120])
    else
        m_text(-180,80,species_name,'fontname','times new roman','fontsize',fs)
         m_grid('fontname','times new roman','fontsize',12,'yticklabel',{''},'xtick',[-120 0 120])
    end

end
cb3=contourcbar('fontname','times new roman','fontsize',14);
cb3.Ruler.MinorTick = 'on';
cb3.Position = [ 0.9    0.3    0.015    0.4];
cb3.Label.String=['C allocation (%)'];
cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'supp_allCBio_glob_vert' '.jpeg']);
disp(['Saving figure for:' fig_name]);
export_fig(sprintf(fig_name));
