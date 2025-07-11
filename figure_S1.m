%%

% C allocation for all 10 species
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('C:\Users\mathp\OneDrive\Desktop\FE_opt_Spyd\Optimization\'))
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes"))
addpath(genpath("C:\Users\mathp\OneDrive\Documents\paper1\Codes"));

load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\ynbest_light_dark_2.csv")
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\ynfbest_light_dark_3.csv")
load("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes\afebest_light_dark_2.csv")
df = readtable("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\fe_all_new.csv"); % empirical data points for scatter
namesforloop = ["Sunda (1995) hux","Sunda (1995) cal","Sunda (1995) oce",...
    "Sunda (1997) wei","Sunda (1997) pse",...
    "Sunda (1997) min","Sunda (1997) mic",...
    "Jabre (2020) 6","Jabre (2020) 3","Jabre (2020) 1"];
speciesforloop = {'\it E. huxleyi','\it P. calceolata',...
    '\it T. oceanica', '\it T. weissflogii','\it T. pseudonana','\it P. minimum','\it P. micans',...
    ['{\it F. cylindrus} (6 ' char(176) 'C)'],['{\it F. cylindrus} (3 ' char(176) 'C)'],['{\it F. cylindrus} (1 ' char(176) 'C)']};
%% temperature
figure('Position',[1 1 1100 500])
[ha,pos]=tight_subplot(2,5,[0.08 0.04],[.2 .06],[.06 .02]);
fs = 12;
fe_range = [1e-7:1e-6:10/1000]';
for i = 1 :10
    j = i;
    author_name = namesforloop(j);
    rows = find(df.author==author_name); %grabbing the Fe and Mu from each species
    Fe =df.Fe_uM_(rows);% note: Fe is in uM
    Mu_d = df.mu(rows);

    %Removing negative values
    Fe(Mu_d<0) = nan;
    Mu_d(Mu_d<0) = nan;
    if author_name == "Sunda (1995) hux" || author_name == "Sunda (1995) oce" || author_name == "Sunda (1995) cal" || author_name == "Sunda (1997) pse" || ...
            author_name == "Sunda (1997) wei"
        I = 500;T = 20;
    elseif author_name == "Sunda (1997) mic" || author_name == "Sunda (1997) min" || author_name == "Hudson (1990) br"
        I = 50;T = 20;
    elseif author_name == "Hudson (1990)"
        I = 100;T = 20;
    elseif author_name == "Sunda (1995) hux 175" || author_name == "Sunda (1995) oce 175"
        I = 175;T = 20;
    elseif author_name == "Jabre (2020) 6"
        T = 6;I = 50;
    elseif author_name == "Jabre (2020) 3"
        T = 3;I = 50;
    elseif author_name == "Jabre (2020) 1"
        T = 1;I = 50;
    end
    species_name = speciesforloop(j);
    afe=afebest_light_dark_2(i); %yn = ynbest_light_dark_2(3);
    ynp = ynbest_light_dark_2(i);ynf =  ynfbest_light_dark_3(i);
    [mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
        Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
        Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
        Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
        Cess,Cpho,Cbio,Csto,mumax] = CFM_Fe(afe,ynp,ynf,fe_range,T,I);
    % C-based allocation under changing [Fe]
    axes(ha(i))
    Y = [Cess,Cpho,Cbio,Csto];
    ba =area(fe_range.*1000,Y);
    colorscheme = crameri('tokyo');
    cmap2 = crameri('lapaz');
    color_indices = length(colorscheme(:,1))/4;
    lighten = [1.1 1.1 1.1];
    ba(1).FaceColor =  [0.5839    0.3760    0.5355].*lighten; % nconst protein, purple
    ba(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
    ba(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
    ba(4).FaceColor = colorscheme(120,:); %ndna
    %set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''})
    grid on; box on;
    set(gca,'Layer','top')
    xlim([0 max(Fe)*1000+(max(Fe)*100)])
    ylim([0 100])
    %ylabel(['C allocation (%)'])
    % l=legend({'Ess.','Pho.','Bio.','DNA','RNA','Chl.'},'NumColumns',3);
    % l.Position=[0.3576    0.4666    0.3021    0.0473];
    j = i;
    if j ==4
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},...
            'xtick',linspace(0,0.8,5))
        text(0.34,6,species_name,'fontname','times new roman','fontsize',fs)
        text(0.72,95,'d.','fontname','times new roman','fontsize',fs)
        xlim([0 0.8])
    elseif j == 3
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''})

        text(0.03,6,species_name,'fontname','times new roman','fontsize',fs);
        text(0.054,95,'c.','fontname','times new roman','fontsize',fs)
        xlim([0 0.06])
    elseif j == 2
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''})
        text(0.027,6,species_name,'fontname','times new roman','fontsize',fs)
        text(0.054,95,'b.','fontname','times new roman','fontsize',fs)
        xlim([0 0.06])
    elseif j == 1
        set(gca,'fontname','times new roman','fontsize',fs)
        text(0.034,6,species_name,'fontname','times new roman','fontsize',fs)
        text(0.054,95,'a.','fontname','times new roman','fontsize',fs)
        xlim([0 0.06])
    elseif j == 5
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},...
            'xtick',linspace(0,0.8,5))
        text(0.31,6,species_name,'fontname','times new roman','fontsize',fs)
        text(0.72,95,'e.','fontname','times new roman','fontsize',fs)
        xlim([0 0.8])
    elseif j == 6
        set(gca,'fontname','times new roman','fontsize',fs,...
            'xtick',linspace(0,0.8,5))
        text(0.38,6,species_name,'fontname','times new roman','fontsize',fs)
        y = text(-0.21,80,['C allocation (%)'],'rotation',90,'fontname','times new roman','fontsize',fs);
        text(0.72,94.5,'f.','fontname','times new roman','fontsize',fs)
        xlim([0 0.8])
    elseif j == 7
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},...
            'xtick',linspace(0,0.8,5))
        text(0.45,6,species_name,'fontname','times new roman','fontsize',fs)
        text(0.72,95,'g.','fontname','times new roman','fontsize',fs)
        xlim([0 0.8])
    elseif j == 8
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},'xtick',linspace(0,0.8,5))

        text(0.16,6,species_name,'fontname','times new roman','fontsize',fs)
        x = text(0.3,-22,['Fe (nM)'],'fontname','times new roman','fontsize',fs);
        text(0.72,95,'h.','fontname','times new roman','fontsize',fs)
        ylim([0 100])
        xlim([0 0.8])
    elseif j == 9
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},'xtick',linspace(0,0.8,5))
        text(0.16,6,species_name,'fontname','times new roman','fontsize',fs)
        text(0.72,95,'i.','fontname','times new roman','fontsize',fs)
        ylim([0 100])
        xlim([0 0.8])
    else
        set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},'xtick',linspace(0,0.8,5))
        text(0.16,6,species_name,'fontname','times new roman','fontsize',fs)
        text(0.72,95,'j.','fontname','times new roman','fontsize',fs)
        ylim([0 100])
        xlim([0 0.8])

    end

end
l = legend({'Ess.' ,'Pho.','Bio.','Sto'},'NumColumns',4,'Position',[0.0591    0.1078    0.2715    0.0413
],'box','off');

% saving figure
cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'fig_s4' '.jpeg']);
disp(['Saving figure for:' fig_name]);
export_fig(sprintf(fig_name));