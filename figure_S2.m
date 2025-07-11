
%% qFe as a function of temp and light
% for figure S2
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('C:\Users\mathp\OneDrive\Desktop\FE_opt_Spyd\Optimization\'))
addpath(genpath("C:\Users\mathp\OneDrive\Desktop\final_paper1_codes"))
addpath(genpath("C:\Users\mathp\OneDrive\Documents\paper1\Codes"));
%% loading optimized params, selecting t oceanica (3rd entry)
load("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\ynfbest_light_dark_3.csv");
load("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\ynpbest_light_dark_3.csv");
load("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\afebest_light_dark_3.csv")
afe=afebest_light_dark_3(3); %yn = ynbest_light_dark_2(3);
ynp = ynpbest_light_dark_3(3);
ynf =  ynfbest_light_dark_3(3);
fs = 12;
fe_range = [1e-7:1e-6:1/1000]';
Qc=1.00*10^(-12)/12 ; 
Feunit = 1/Qc; %conversion to mol mol-1


figure('Position',[1 1 800 600])
[ha,pos]=tight_subplot(2,2,[0.06 0.05],[.1 .1],[.18 .18]);
axes(ha(1))
% high light (500 umol quanta s), low temp (10 deg C)
colororder({'k','k'})
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,ynp,ynf,fe_range,10,500);
Y = qfe.*Feunit.*1e9;
yyaxis left
ba1 =area(fe_range.*1000,Y,'FaceAlpha',0.3);

ba1(1).FaceColor =  [0.4839    0.2760    0.4355]; % nconst protein, purple
set(gca,'fontname','times new roman','fontsize',fs,'xticklabel',{''})
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.2e-3*1000])
ylim([0 0.8e-3])
text(0.18,7.4e-4,'a.','fontname','times new roman','fontsize',12)
text(-0.027,-3e-4,'Fe : C (mol mol^{-1})','fontname','times new roman','fontsize',12,'Rotation',90)
yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''})
xlim([0 0.2e-3*1000])
axes(ha(2))
% high light (500), high temp (20)
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,ynp,ynf,fe_range,20,500);
Y = qfe.*Feunit.*1e9;
yyaxis left
ba2 =area(fe_range.*1000,Y,'FaceAlpha',0.3);
ba2(1).FaceColor =  [0.4839    0.2760    0.4355]; % nconst protein, purple
set(gca,'fontname','times new roman','fontsize',fs,'xticklabel',{''},'yticklabel',{''})
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.2e-3*1000])
ylim([0 0.8e-3])
text(0.18,7.4e-4,'b.','fontname','times new roman','fontsize',12)
yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
xlim([0 0.2e-3*1000])
axes(ha(3))
% low light (50), low temp
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,ynp,ynf,fe_range,10,100);
Y = qfe.*Feunit.*1e9;
yyaxis left
ba3 =area(fe_range.*1000,Y,'FaceAlpha',0.3);
ba3(1).FaceColor =  [0.4839    0.2760    0.4355]; % nconst protein, purple
set(gca,'fontname','times new roman','fontsize',fs)
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.2e-3*1000])
ylim([0 0.8e-3])
text(0.18,7.4e-4,'c.','fontname','times new roman','fontsize',12)
yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
xlim([0 0.2e-3*1000])
set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''},'xtick',linspace(0,0.2e-3*1000,5),'XTickLabel',{'0','0.05','0.1','0.15','0.2'})

axes(ha(4))
% low light (10 umol quanta s), high temp (20 deg C)
[mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax]= CFM_Fe(afe,ynp,ynf,fe_range,20,100);
Y = qfe.*Feunit.*1e9;
yyaxis left
ba4 =area(fe_range.*1000,Y,'FaceAlpha',0.3);
ba4(1).FaceColor =  [0.4839    0.2760    0.4355]; % nconst protein, purple
set(gca,'fontname','times new roman','fontsize',fs,'yticklabel',{''})
grid on; box on;
cmocean('balance')
set(gca,'Layer','top')
xlim([0 0.2e-3*1000])
ylim([0 0.8e-3])
text(0.18,7.4e-4,'d.','fontname','times new roman','fontsize',12)

yyaxis right
plot(fe_range.*1000,mu.*86400,'k--','linewidth',1.4)
ylim([0 2])
xlim([0 0.2e-3*1000])
set(gca,'fontname','times new roman','fontsize',fs,'xtick',...
    linspace(0,0.2e-3*1000,5),'XTickLabel',{'0','0.05','0.1','0.15','0.2'})
text(0.232,2.4,['\mu ' '(d^{-1})'],'fontname','times new roman','fontsize',12,'rotation',-90)
text(-0.05,-0.32,'[Fe] (nM)','fontname','times new roman','fontsize',12)
%insert path for saving figure
cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
fig_name=([ 'q_fe_TI_Toce_new' '.jpeg']);
disp(['Saving figure for:' fig_name]);
export_fig(sprintf(fig_name));


