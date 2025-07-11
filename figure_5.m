%%

% N:C and P: C under changing [dFe]
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('C:\Users\mathp\OneDrive\Desktop\FE_opt_Spyd\Optimization\'))
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes"))
addpath(genpath("C:\Users\mathp\OneDrive\Documents\paper1\Codes"));
%%
t1 = readtable("C:\Users\mathp\OneDrive\Documents\paper1\Codes\optimized_params.xlsx");
afe = t1.a_Fe(3); Afepho = t1.AFe_pho(3); % 3 corresponds to T oceanica (look at excel file)
fs = 12;
fe_range = [1e-7:1e-7:1/1000]';


%
figure('Position',[1 1 800 600])
[ha,pos]=tight_subplot(2,2,[0.06 0.05],[.1 .1],[.18 .18]);

axes(ha(1))
% high light (500 umol quanta s), low temp (10 deg C)
colororder({'k','k'})
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,Afepho,fe_range,10,500);
Y = [Cess,Cpho,Cbio,Csto];
yyaxis left
ba1 =area(fe_range.*1000,Y,'FaceAlpha',0.8);
colorscheme = crameri('tokyo');
lighten = [1.2 1.2 1.2];
ba1(1).FaceColor =  [0.4839    0.2760    0.4355].*lighten;% nconst protein, purple
ba1(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba1(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba1(4).FaceColor = colorscheme(120,:); %ndna
ba1(4).FaceAlpha = 0.8;
set(gca,'fontname','times new roman','fontsize',fs,'xticklabel',{''})
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.04e-3*1000])
ylim([0 100])
text(0.036,94,'a.','fontname','times new roman','fontsize',12)
text(-0.007,-26,'C allocation (%)','fontname','times new roman','fontsize',12,'Rotation',90)
yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
set(gca,'fontname','times new roman','fontsize',fs,'ytick',[0,0.4,0.8,1.2,1.6,2],...
    'yticklabel',{''})
xlim([0 0.04e-3*1000])
axes(ha(2))
% high light (500), high temp (20)
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,Afepho,fe_range,20,500);
Y = [Cess,Cpho,Cbio,Csto];
yyaxis left
ba2 =area(fe_range.*1000,Y,'FaceAlpha',0.8);
colorscheme = crameri('tokyo');
lighten = [1.2 1.2 1.2];
ba2(1).FaceColor =  [0.4839    0.2760    0.4355].*lighten; % nconst protein, purple
ba2(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba2(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba2(4).FaceColor = colorscheme(120,:); %ndna
ba2(4).FaceAlpha = 0.8;
set(gca,'fontname','times new roman','fontsize',fs,'xticklabel',{''},'yticklabel',{''})
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.04e-3*1000])
ylim([0 100])
text(0.036,94,'b.','fontname','times new roman','fontsize',12)
yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
xlim([0 0.04e-3*1000])
set(gca,'ytick',[0,0.4,0.8,1.2,1.6,2])
axes(ha(3))
% low light (50), low temp
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,Afepho,fe_range,10,50);
Y = [Cess,Cpho,Cbio,Csto];
yyaxis left
ba3 =area(fe_range.*1000,Y,'FaceAlpha',0.8);
colorscheme = crameri('tokyo');
cmap2 = crameri('lapaz');
color_indices = length(colorscheme(:,1))/4;
lighten = [1.2 1.2 1.2];
ba3(1).FaceColor =  [0.4839    0.2760    0.4355].*lighten;
ba3(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba3(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba3(4).FaceColor = colorscheme(120,:); %ndna
ba3(4).FaceAlpha = 0.8;
set(gca,'fontname','times new roman','fontsize',fs)
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.04e-3*1000])
ylim([0 100])
text(0.036,94,'c.','fontname','times new roman','fontsize',12)
yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
xlim([0 0.04e-3*1000])
set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},'xtick',linspace(0,0.04e-3*1000,5),...
    'XTickLabel',{'0','0.01','0.02','0.03','0.04'},'ytick',[0,0.4,0.8,1.2,1.6,2])

axes(ha(4))
% low light (10 umol quanta s), high temp (20 deg C)
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,Afepho,fe_range,20,50);
Y = [Cess,Cpho,Cbio,Csto];
yyaxis left
ba4 =area(fe_range.*1000,Y,'FaceAlpha',0.8);
colorscheme = crameri('tokyo');
cmap2 = crameri('lapaz');
color_indices = length(colorscheme(:,1))/4;
lighten = [1.2 1.2 1.2];
ba4(1).FaceColor =  [0.4839    0.2760    0.4355].*lighten; % nconst protein, purple
ba4(2).FaceColor = [0.6277    0.8116    0.6023]; % nphoto, light green
ba4(3).FaceColor = [0.9997    0.9998    0.8505];  %nbio, light yellow
ba4(4).FaceColor = colorscheme(120,:); %ndna
ba4(4).FaceAlpha = 0.8;
set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''})
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.04e-3*1000])
ylim([0 100])
text(0.036,94,'d.','fontname','times new roman','fontsize',12)
text(-0.008,-16,'[Fe] (nM)','fontname','times new roman','fontsize',12)
text(0.047,120,['\mu ' '(d^{-1})'],'fontname','times new roman','fontsize',12,'rotation',-90)
l=legend([ba1,ba2,ba3,ba4],{'Ess.','Pho.','Bio.','Sto.'},'NumColumns',2);
l.Position=[0.2798    0.5349    0.1917    0.0631];
yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
xlim([0 0.04e-3*1000])
set(gca,'fontname','times new roman','fontsize',fs,'xtick',linspace(0,0.04e-3*1000,5),'XTickLabel',{'0','0.01','0.02','0.03','0.04'},...
    'ytick',[0,0.4,0.8,1.2,1.6,2])

%insert path for saving figure
% cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
% set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
% fig_name=([ 'C_allocation_TandI_toce_fin' '.jpeg']);
% disp(['Saving figure for:' fig_name]);
% export_fig(sprintf(fig_name));