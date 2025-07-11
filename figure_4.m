%%
% N:C ratio for 10 modeled species under changing [Fe]
% C allocation for 10 modeled species under iron-deplete and iron-replete
% conditions
clear;
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('C:\Users\mathp\OneDrive\Desktop\FE_opt_Spyd\Optimization\'))
addpath(genpath("C:\Users\mathp\OneDrive\Documents\paper1\Codes"));

t1 = readtable("C:\Users\mathp\OneDrive\Documents\paper1\Codes\optimized_params.xlsx");
afefull = t1.a_Fe; Afepho = t1.AFe_pho;
namesforloop = ["Sunda (1995) hux","Sunda (1995) cal","Sunda (1995) oce",...
    "Sunda (1997) wei","Sunda (1997) pse",...
    "Sunda (1997) min","Sunda (1997) mic",...
    "Jabre (2020) 6","Jabre (2020) 3","Jabre (2020) 1"];
fe_range = [1e-7:1e-6:1/1000]';
figure('Position',[1 1 960 560])
ha=tight_subplot(2,2, [.08 .02],[ .09 .07], [.16 .03]);
% LOW iron (1e-6 uM)
%pre-allocating space
y1 = nan(10,1);y2 = nan(10,1);y3 = nan(10,1);
c1 = nan(10,1);c2 = nan(10,1);c3 = nan(10,1);c4 = nan(10,1);
count = 0;
for i = 1 : 10
    author_name = namesforloop(i);
    count = count + 1;
    afe=afefull(i); 
    Afe_pho = Afepho(i);
    if author_name == "Sunda (1995) hux" || author_name == "Sunda (1995) oce" || author_name == "Sunda (1995) cal" || author_name == "Sunda (1997) pse" || ...
            author_name == "Sunda (1997) wei"
        I = 500;
    elseif author_name == "Sunda (1997) mic" || author_name == "Sunda (1997) min" || author_name == "Jabre (2020) 6" || author_name == "Jabre (2020) 3" || ...
            author_name == "Jabre (2020) 1" || author_name == "Hudson (1990) br"
        I = 50;
    elseif author_name == "Hudson (1990)"
        I = 100;
    elseif author_name == "Sunda (1995) hux 175" || author_name == "Sunda (1995) oce 175"
        I = 175;
    end
    if i == 8
        T = 6;
    elseif i == 9
        T = 3;
    elseif i == 10
        T = 1;
    else
        T = 20;
    end
    [mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
        Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
        Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
        Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
        Cess,Cpho,Cbio,Csto,Ness,Npho,Nbio]= CFM_Fe(afe,Afe_pho,fe_range(5),T,I);
    y1(count) = Ndna_const_plot+Nconst_protein_plot; y2(count) = Nphoto_plot+Nchl_plot;
    y3(count) = Nbiosynth_plot+Nrna_plot;
    c1(count) = Cess; c2(count) = Cpho; c3(count) = Cbio; c4(count) = Csto;
end
Y = [y1,y2,y3]; C = [c1,c2,c3,c4];
axes(ha(1))
ba = barh(1:10,C,'stacked', 'FaceColor','flat','barwidth',1);
colorscheme = crameri('tokyo');
lighten = [1.1 1.1 1.1];
ba(1).FaceColor =  [0.5839    0.3760    0.5355].*lighten; % nconst protein, purple
ba(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba(4).FaceColor = colorscheme(120,:); %sto
ba(1).FaceAlpha = 0.8;
ba(4).FaceAlpha = 0.8;
xlim([0 100])
ylim([0.5 10.5])
box on;
set(gca,'ydir','reverse','fontname','times new roman','fontsize',12,'yticklabel',{'\it E. huxleyi','\it P. calceolata',...
    '\it T. oceanica', '\it T. weissflogii','\it T. pseudonana','\it P. minimum','\it P. micans',...
    ['{\it F. cylindrus} (6 ' char(176) 'C)'],['{\it F. cylindrus} (3 ' char(176) 'C)'],['{\it F. cylindrus} (1 ' char(176) 'C)']})
text(90,11.9,'C allocation (%)','fontname','times new roman','fontsize',12)
text(3,1,'a.','fontname','times new roman','fontsize',12)
title({'Iron-limited'},'FontWeight','Normal')
l=legend({'Ess.','Photo.','Bio.','Sto.'},'NumColumns',4);
l.Position=[  0.0037    0.009    0.3643    0.04];
legend('boxoff')

axes(ha(3))
ba = barh(1:10,Y,'stacked', 'FaceColor','flat','barwidth',1);
lighten = [1.1 1.1 1.1];
ba(1).FaceColor =  [0.5839    0.3760    0.5355].*lighten; % nconst protein, purple
ba(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba(1).FaceAlpha = 0.8;
xlim([0 0.25])
ylim([0.5 10.5])
box on;
set(gca,'ydir','reverse','fontname','times new roman','fontsize',12,'yticklabel',{'\it E. huxleyi','\it P. calceolata',...
    '\it T. oceanica', '\it T. weissflogii','\it T. pseudonana','\it P. minimum','\it P. micans',...
    ['{\it F. cylindrus} (6 ' char(176) 'C)'],['{\it F. cylindrus} (3 ' char(176) 'C)'],['{\it F. cylindrus} (1 ' char(176) 'C)']})
text(0.22,11.9,'N : C (mol mol^{-1})','fontname','times new roman','fontsize',12)
text(0.007,1,'c.','fontname','times new roman','fontsize',12)
y1 = nan(10,1);y2 = nan(10,1);y3 = nan(10,1);
c1 = nan(10,1);c2 = nan(10,1);c3 = nan(10,1);c4 = nan(10,1);
count = 0;
for i = 1 : 10
    author_name = namesforloop(i);

    count = count + 1;
    afe=afefull(i); 
    Afe_pho = Afepho(i);
    if author_name == "Sunda (1995) hux" || author_name == "Sunda (1995) oce" || author_name == "Sunda (1995) cal" || author_name == "Sunda (1997) pse" || ...
            author_name == "Sunda (1997) wei"
        I = 500;
    elseif author_name == "Sunda (1997) mic" || author_name == "Sunda (1997) min" || author_name == "Jabre (2020) 6" || author_name == "Jabre (2020) 3" || ...
            author_name == "Jabre (2020) 1" || author_name == "Hudson (1990) br"
        I = 50;
    elseif author_name == "Hudson (1990)"
        I = 100;
    elseif author_name == "Sunda (1995) hux 175" || author_name == "Sunda (1995) oce 175"
        I = 175;
    end
    if i == 8
        T = 6;
    elseif i == 9
        T = 3;
    elseif i == 10
        T = 1;
    else
        T = 20;
    end
    [mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
        Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
        Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
        Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
        Cess,Cpho,Cbio,Csto,Ness,Npho,Nbio]= CFM_Fe(afe,Afe_pho,fe_range(500),T,I);
    y1(count) = Ndna_const_plot+Nconst_protein_plot; y2(count) = Nphoto_plot+Nchl_plot;
    y3(count) = Nbiosynth_plot+Nrna_plot;
    c1(count) = Cess; c2(count) = Cpho; c3(count) = Cbio; c4(count) = Csto;
end
Y = [y1,y2,y3]; C = [c1,c2,c3,c4];
axes(ha(2))
ba = barh(1:10,C,'stacked', 'FaceColor','flat','barwidth',1);
colorscheme = crameri('tokyo');
lighten = [1.1 1.1 1.1];
ba(1).FaceColor =  [0.5839    0.3760    0.5355].*lighten; % nconst protein, purple
ba(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba(4).FaceColor = colorscheme(120,:); %ndna
ba(1).FaceAlpha = 0.8;
xlim([0 100])
ylim([0.5 10.5])
 box on;
set(gca,'ydir','reverse','fontname','times new roman','fontsize',12,'yticklabel',{})
title({'Iron-replete'},'FontWeight','Normal')
text(3,1,'b.','fontname','times new roman','fontsize',12)
axes(ha(4))
ba = barh(1:10,Y,'stacked', 'FaceColor','flat','barwidth',1);
colorscheme = crameri('tokyo');
lighten = [1.1 1.1 1.1];
ba(1).FaceColor =  [0.5839    0.3760    0.5355].*lighten; % nconst protein, purple
ba(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba(1).FaceAlpha = 0.8;
xlim([0 0.25])
ylim([0.5 10.5])
box on;
set(gca,'ydir','reverse','fontname','times new roman','fontsize',12,'yticklabel',{})
text(0.007,1,'d.','fontname','times new roman','fontsize',12)

% saving figure
% cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
% set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
% fig_name=([ 'c_and_n_c_Functional_alloc' '.jpeg']);
% disp(['Saving figure for:' fig_name]);
% export_fig(sprintf(fig_name));

