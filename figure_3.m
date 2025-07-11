
%% Showing model output (using optimized afe, A_Pho_Fe values)
% alongside datapoints
% Figure 1 in paper
% author: Maggie Bernish, URI-GSO (2024)
addpath(genpath('C:\Users\mathp\OneDrive\Documents\toolbox\'));
addpath(genpath('C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\'))
addpath(genpath('C:\Users\mathp\OneDrive\Desktop\FE_opt_Spyd\Optimization\'))
addpath(genpath("C:\Users\mathp\OneDrive\Documents\paper1\Codes"));

% load optimized params here
t1 = readtable("C:\Users\mathp\OneDrive\Documents\paper1\Codes\optimized_params.xlsx");
afe = t1.a_Fe; Afepho = t1.AFe_pho;
df = readtable("C:\\Users\\mathp\\OneDrive\\Desktop\\FE_opt_Spyd\\Optimization\\fe_all_new.csv"); % empirical data points for scatter
namesforloop = ["Sunda (1995) hux","Sunda (1995) cal","Sunda (1995) oce",...
    "Sunda (1997) wei","Sunda (1997) pse",...
    "Sunda (1997) min","Sunda (1997) mic",...
    "Jabre (2020) 6","Jabre (2020) 3","Jabre (2020) 1"];
speciesforloop = {'\it E. huxleyi','\it P. calceolata',...
    '\it T. oceanica', '\it T. weissflogii','\it T. pseudonana','\it P. minimum','\it P. micans',...
    ['{\it F. cylindrus} (6 ' char(176) 'C)'],['{\it F. cylindrus} (3 ' char(176) 'C)'],['{\it F. cylindrus} (1 ' char(176) 'C)']};
%%
figure('Position',[1 1 1100 500])
[ha,pos]=tight_subplot(2,5,[0.12 0.04],[.12 .14],[.06 .03]);
for j = 1:10
    author_name = namesforloop(j);
    rows = find(df.author==author_name); %grabbing the Fe and Mu from each species
    Fe =df.Fe_uM_(rows);% note: Fe is in uM
    Mu_d = df.mu(rows);
    %Removing negative values
    Fe(Mu_d<0) = nan;
    Mu_d(Mu_d<0) = nan;
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
    afebest = afe(j);
    Afephobest=  Afepho(j);
    % plotting
    Fe_high_res = [1e-6:1e-6:max(Fe)+(max(Fe)*0.5)];
    species_name = speciesforloop(j);
    Numbertoarray=ones(size(Fe_high_res));
    Mu = calc_growth(afebest,Afephobest,Fe_high_res,I,author_name);
    D = Mu;
    axes(ha(j))
    plot(Fe_high_res.*1000,D*86400,'k-','linewidth',1.2);
    hold on;
    scatter(Fe.*1000,Mu_d,50,'k','filled','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerFaceAlpha',0.8)
    box on; grid on;
    fs = 12;
    set(gca,'fontname','times new roman','fontsize',fs)
    ylim([0 max(Mu_d)+0.4]) 
    xlim([0 max(Fe)*1000+(max(Fe)*100)])

    if j ==4
        text(0.34,0.07,species_name,'fontname','times new roman','fontsize',fs)
        text(0.04,1.18,'d.','fontname','times new roman','fontsize',fs)
    elseif j == 3
        text(0.41,0.12,species_name,'fontname','times new roman','fontsize',fs);
        text(0.04,1.84,'c.','fontname','times new roman','fontsize',fs)
    elseif j == 2
        text(0.28,0.086,species_name,'fontname','times new roman','fontsize',fs)
        text(0.038,1.32,'b.','fontname','times new roman','fontsize',fs)
    elseif j == 1
        text(0.46,0.09,species_name,'fontname','times new roman','fontsize',fs)
        text(0.04,1.54,'a.','fontname','times new roman','fontsize',fs)
    elseif j == 5
        text(11.5,0.13,species_name,'fontname','times new roman','fontsize',fs)
          text(1,2.04,'e.','fontname','times new roman','fontsize',fs)
    elseif j == 6
        text(4.4,0.076,species_name,'fontname','times new roman','fontsize',fs)
        y = text(-2.2,0.8,['Specific growth rate (d^{-1})'],'rotation',90,'fontname','times new roman','fontsize',fs);
        text(0.4,1.16,'f.','fontname','times new roman','fontsize',fs)
    elseif j == 7
        text(2.6,0.046,species_name,'fontname','times new roman','fontsize',fs)
          text(0.2,0.74,'g.','fontname','times new roman','fontsize',fs)
    elseif j == 8
        text(0.32,0.024,species_name,'fontname','times new roman','fontsize',fs)
        x = text(0.38,-0.1,['Fe (nM)'],'fontname','times new roman','fontsize',fs);
          text(0.06,0.368,'h.','fontname','times new roman','fontsize',fs)
          ylim([0 0.4])
          xlim([0 1.5])
    elseif j == 9
        text(0.32,0.024,species_name,'fontname','times new roman','fontsize',fs)
          text(0.062,0.368,'i.','fontname','times new roman','fontsize',fs)
        ylim([0 0.4])
        xlim([0 1.5])
    else
        text(0.34,0.024,species_name,'fontname','times new roman','fontsize',fs)
          text(0.066,0.368,'j.','fontname','times new roman','fontsize',fs)
         ylim([0 0.4])
         xlim([0 1.5])
    end
end
% uncomment the following lines and change the path in 103 to export figure as a jpeg
% cd("C:\Users\mathp\OneDrive\Documents\Fe_opt_MATLAB\figs_iron")
% set(gcf, 'Color', 'w');set(gcf,'PaperPositionMode','auto');%maintains size of fig when saved
% fig_name=([ 'fig2_newc' '.jpeg']);
% disp(['Saving figure for:' fig_name]);
% export_fig(sprintf(fig_name));


%%
% calculates D or mu given optimized parameters, light intensity, iron conc
function D = calc_growth(aFe,Afe_pho,Fe,I_s,author_name)
% parameter sets
Pmax=0.00320513285659728;
OT=0.00863364097132997;
m=3.79146798299876E-19   ;      %(mol C s-1 cell-1) maintenance carbonhydrate consumption (idea from 172-7)
Nstore_max=2.91679384515998E-15  ;       %(molN cell-1) Constant protein pool in nitrogen (193-25)
Ypthylakoid_chl=0.0281633095303638 ;       %((molP cell-1)/(molC chl cell-1)) the shoichiometric ratio for cell phosphorus in thylakoid membrane to chlorophyll (193-26)
Pconst_other=5.44534485638617E-17  ;             %(molP cell-1) Constant part of phosphorus (193-26) * This includes ATP ADP, Phospholipid and DNA RNA
Qp_max=25.26/(3.097e16)             ;                                                                 %(molP cell-1) total phosphorus content in the cell (193-26)
Cessential=1.51786753491048E-15      ;    %(molC cell-1) essential carbon (lipid membrane, etc.) *8.33e-14/10 is 10%
Nconst_protein = 4.45336898828389E-15 ;  %(molN cell-1) Constant protein pool in nitrogen (193-25)


Molar_mass_DNA_AT_average=308.47  ;      %(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_DNA_CG_average=308.97   ;     %(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_AT_average=317.47    ;    %(g mol-1) Molar mass AT average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")
Molar_mass_RNA_CG_average=324.97     ;   %(g mol-1) Molar mass CG average (from "10 Amino acids and nucleic acids stoichiometry.xlsx")

%E coli
CG_Ecoli=0.506  ;        %(dimensionless) from (http://www.ncbi.nlm.nih.gov/genome/167 (accessed 06/18/2016))
AT_Ecoli=1-CG_Ecoli;     %(dimensionless)

Molar_mass_DNA_Ecoli=Molar_mass_DNA_AT_average*CG_Ecoli+Molar_mass_DNA_CG_average*AT_Ecoli ;    %(g mol-1) Molar mass of DNA unit
Molar_mass_RNA_Ecoli=Molar_mass_RNA_AT_average*CG_Ecoli+Molar_mass_RNA_CG_average*AT_Ecoli  ;   %(g mol-1) Molar mass of RNA unit

RNA_DNA_mass_ratio=17.844/6.5239;  %(ug/ug) from values ad D=0 "07 Bremer and Dennis 1996 data plot.xlsx"
RNA_DNA_molar_ratio=RNA_DNA_mass_ratio/Molar_mass_RNA_Ecoli*Molar_mass_DNA_Ecoli  ;  %(mol mol-1)


%Stoichiometric parameters for DNA and RNA
CG=0.563      ;
AT = 1-CG;
AU = 1-CG;
C_CG_RNA = 19/2;
N_CG_RNA = 8/2;
C_AU_RNA = 19/2;
N_AU_RNA = 7/2;

%DNA
C_CG_DNA = 19/2;
N_CG_DNA = 8/2;

C_AT_DNA = 20/2;
N_AT_DNA = 7/2;

YnucacidN_P = N_CG_RNA*CG + N_AU_RNA*AU; %same for RNA and DNA
YnucacidP_N = 1/YnucacidN_P;

YdnaC_N = (C_CG_DNA*CG + C_AT_DNA*AT)/(N_CG_DNA*CG + N_AT_DNA*AT);
YrnaC_N = (C_CG_RNA*CG + C_AU_RNA*AU)/(N_CG_RNA*CG + N_AU_RNA*AU);
DNAmb=2.1269     ;              %(Mb) Megabase pair of synechococcus DNA in mega (million) base pairs (http://www.ncbi.nlm.nih.gov/genome/13522 (accessed 06/18/2016))
Avogadro=6.022*10^23  ;         %(molecules mol-1) Avogadro constant
Pdna_const=DNAmb*2*10^6/Avogadro ;               %(molP cell-1) Constant part of DNA in phosphorus
Prna_const=Pdna_const*RNA_DNA_molar_ratio   ;    %(molP cell-1) Constant part of RNA in phosphorus
%* Make sure to multiply by 2 as they are base PAIRs"
Ndna_const=Pdna_const/YnucacidP_N  ;    %(molN cell-1) Constant part of DNA in nitrogen
Nrna_const=Ndna_const*RNA_DNA_molar_ratio ;  %(molN cell-1) Constatn part of RNA in phosphorus
Ndna=Ndna_const  ;  %(molN cell-1) DNA in nitrogen (here assuming constant)
%Ynphoto_chl=3.56099164557551
%YphotoFe_N=0.001636364 ; %(molFe molN-1) Fe/N ratio in photosystem iron (0.001636364 in Inomura et al., 2022)

E=0.7742;
Qc=1.00*10^(-12)/12 ;     %(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
YchlN_C=4/55;

%Conversion parameters================
CNprotein=4.49 ;  %(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
YcyanoC_N=2  ;                           %(molC molN) C/N molar ratio of cyanophycin
YpgC_P=40   ;                        %(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
Nunit=1/Qc;%*14*10**6/(12*10**3)          %((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
Punit=1/Qc;%*30.97*10**6/(12*10**3)       %((ug P / mgC)/(molP cell-1) unit conversion term (164-20)
Ea = 70000 ;%activation energy
R = 8.3;
A =Ea/R;
%Tt - ambient temperature that you are reading in
if author_name == "Jabre (2020) 1"
    Tt = 1 + 273.15;
elseif author_name == "Jabre (2020) 3"
    Tt = 3 + 273.15;
elseif author_name == "Jabre (2020) 6"
    Tt = 6 + 273.15;
else
    Tt = 20 + 273.15;
end

Tref = 293;
Arr = exp(-A*((1/Tt)-(1/Tref)));

Cnbiosynth=4.34728279914354E-10/Arr  ;      %(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Cnrna_variable=6212.59249917364/Arr   ;    %(s) Constant for Variable part of RNA (193-26)

Pchl=Pmax*(1-exp(-OT*I_s)); %(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
if author_name == "Jabre (2020) 1" || author_name == "Jabre (2020) 3" || author_name == "Jabre (2020) 6"
    ld_cycle = 1;
else
    ld_cycle =     0.71428571428;
end
YphotoFe_N=0.001636364  ;%(molFe molN-1) Fe/N ratio in photosystem iron (0.001636364 in Inomura et al., 2022)clc
Ynphoto_chl = Afe_pho./YphotoFe_N;
Pchl = Pchl*ld_cycle;
Vfe=aFe*Fe    ;%(mol fe cell-1 s-1) iron uptake per cell
A=((1+E)*Qc*Ynphoto_chl)/Pchl+Cnbiosynth;
B=Nconst_protein+(m*Ynphoto_chl)/Pchl;
L = ((1 + E)*Qc*Ypthylakoid_chl)/Pchl;
M = (m*Ypthylakoid_chl)/Pchl;
%========================
% Fe limitation related
%========================
R=((1+E)*Qc*Ynphoto_chl*YphotoFe_N)/Pchl;
S=(m*Qc*Ynphoto_chl*YphotoFe_N)/Pchl;

aFe = R;
bFe = S;
cFe = -Vfe;

DFe=DSolver(aFe,bFe,cFe);
D=DFe;
% DFe=solver_2D(aFe,bFe,cFe)

Chl_const = m/Pchl           ;                   % (molC chl cell-1) cN(i)hlrophyll concentration (193-25)
Chl_D = (1 + E)*Qc/Pchl     ;                      % (molC chl cell-1) cN(i)hlrophyll concentration (193-25)

%Nitrogen related
Nchl_D = Chl_D*YchlN_C                  ;        % (molN chl cell-1) Chlorophyll N concentration
Nphoto_D = Chl_D*Ynphoto_chl            ;        % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nchl_const = Chl_const*YchlN_C          ;        % (molN chl cell-1) Chlorophyll N concentration
Nphoto_const = Chl_const*Ynphoto_chl        ;    % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nbiosynth_D = Cnbiosynth                   ;     % (molN cell-1) various part of biosynthesis related protein in N (193-37)
Nprotein_D = Nphoto_D + Nbiosynth_D            ; % (molN cell-1) All the proteins in N (193-26)
Nprotein_const = Nphoto_const + Nconst_protein;  % (molN cell-1) All the proteins in N (193-26)
Nrna_D = Nprotein_const*Cnrna_variable;          % (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
Nrna_D2 = Nprotein_D*Cnrna_variable  ;           % (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)

%Constant carbon parameters------------------
Cconst_protein = Nconst_protein*CNprotein;  %(molC cell-1) carbon in other protein assumed constant (195-16)
Cdna_const = Ndna_const*YdnaC_N  ;    %(molC cell-1) carbon in constant part of DNA (195-16)

%Calculating factors and Qc_essential-------------

Qc_D2 = A*Cnrna_variable*YrnaC_N;
Qc_D = (1+E)*Qc/Pchl + A*CNprotein + L*YpgC_P + B*Cnrna_variable*YrnaC_N;
Qc_const = m/Pchl + B*CNprotein + M*YpgC_P + Nrna_const*YrnaC_N + Ndna_const*YdnaC_N + Cessential;

Qc_essential = Qc_essential_computation(D,Chl_const, Chl_D, Nchl_const, Nchl_D, Nphoto_const, Nphoto_D, Nbiosynth_D,...
    Nprotein_const, Nprotein_D, Nrna_const, Nrna_D, Nrna_D2, Ypthylakoid_chl, YrnaC_N,...
    CNprotein, YpgC_P, Cdna_const, Cconst_protein, Cessential);

i0 = Qc_essential >= Qc   ;%conditions where Mu max applies
i1 = isnan(Qc_essential) ; %Including this avoide when Qc_essential is nan with high Vn
i22 = i0 + i1;
D(i22==1) = 2*(Qc - Qc_const)/(Qc_D + sqrt(Qc_D*Qc_D + 4*Qc_D2*(Qc - Qc_const)));
end



