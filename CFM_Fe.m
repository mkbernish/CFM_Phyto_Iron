function [mu,qfe,n2p,Qn,Qp,Nconst_protein_plot,Nphoto_plot,...
    Nbiosynth_plot,Ndna_const_plot,Nrna_plot,Nchl_plot,...
    Prna_variable_plot,Pdna_variable_plot,Pthylakoid_plot,...
    Pconst_other_plot,Pdna_const_plot,Prna_const_plot,...
    Cess,Cpho,Cbio,Csto,Chl,mumax,Ness,Npho,Nbio] = CFM_Fe(aFe,A_pho_Fe,Fe,temp,I)
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

% change here
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
Ndna_const=Pdna_const/YnucacidP_N  ;    %(molN cell-1) Constant part of DNA in nitrogen
Nrna_const=Ndna_const*RNA_DNA_molar_ratio ;  %(molN cell-1) Constatn part of RNA in phosphorus
Ndna=Ndna_const  ;  %(molN cell-1) DNA in nitrogen (here assuming constant)

E=0.7742;
Qc=1.00*10^(-12)/12 ;     %(molC/cell) biomass C per cell (196-18)(average of N and P limited cases from Healey 1985)
YchlN_C=4/55;

%Conversion parameters================
CNprotein=4.49 ;  %(molC molN) the ratio of C to N in protein (derived from Brown 1991) calculation in "13 Amino acid composition of different phytoplankton.xlsx"
YcyanoC_N=2  ;                           %(molC molN) C/N molar ratio of cyanophycin
YpgC_P=40   ;                        %(molC molP) C/P molar ratio of PG: Phosphatidyl glycerol (assuming C16 fatty acids (since actually mostly C16 (Huflejt et al., 1990)
Nunit=1/Qc;%*14*10^6/(12*10^3);          %((ug N / mgC)/(molN cell-1) unit conversion term (164-20)
Punit=1/Qc;%*30.97*10^6/(12*10^3)  ;     %((ug P / mgC)/(molP cell-1) unit conversion term (164-20)
Feunit = 1/Qc;
Ea = 70000 ;%activation energy
R = 8.3;
A =Ea/R;
%Tt - ambient temperature that you are reading in
Tt=temp+273.15;
Tref = 293.15;
Arr = exp(-A*((1./Tt)-(1/Tref)));
YphotoFe_N=0.001636364  ;%(molFe molN-1) Fe/N ratio in photosystem iron (0.001636364 in Inomura et al., 2022)
Ynphoto_chl = A_pho_Fe./YphotoFe_N;
Cnbiosynth=4.34728279914354E-10./Arr  ;      %(molN cell-1 s) Constant for varible part of biosynthesis protein nitrogen (193-37)
Cnrna_variable=6212.59249917364./Arr   ;    %(s) Constant for Variable part of RNA (193-26)
Pchl=Pmax*(1-exp(-OT*I)); %(C mol s-1 Chl mol-1) Carbohydrate fixation rate per chlorophyll (167-1)(193-25)
Vfe=aFe*Fe    ;%(mol fe cell-1 s-1) iron uptake per cell
A=((1+E)*Qc*Ynphoto_chl)./Pchl+Cnbiosynth;
B=Nconst_protein+(m*Ynphoto_chl)./Pchl;
L = ((1 + E)*Qc*Ypthylakoid_chl)./Pchl;
M = (m*Ypthylakoid_chl)./Pchl;
%========================
% Fe limitation related
%========================
R=((1+E)*Qc*Ynphoto_chl*YphotoFe_N)./Pchl;
S=(m*Qc*Ynphoto_chl*YphotoFe_N)./Pchl;

aFe = R;
bFe = S;
cFe = -Vfe;

DFe=DSolver(aFe,bFe,cFe);
D=DFe;
% DFe=solver_2D(aFe,bFe,cFe)

Chl_const = m./Pchl           ;                   % (molC chl cell-1) cN(i)hlrophyll concentration (193-25)
Chl_D = (1 + E)*Qc./Pchl     ;                      % (molC chl cell-1) cN(i)hlrophyll concentration (193-25)

%Nitrogen related
Nchl_D = Chl_D*YchlN_C                  ;        % (molN chl cell-1) Chlorophyll N concentration
Nphoto_D = Chl_D*Ynphoto_chl            ;        % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nchl_const = Chl_const*YchlN_C          ;        % (molN chl cell-1) Chlorophyll N concentration
Nphoto_const = Chl_const*Ynphoto_chl        ;    % (molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nbiosynth_D = Cnbiosynth                   ;     % (molN cell-1) various part of biosynthesis related protein in N (193-37)
Nprotein_D = Nphoto_D + Nbiosynth_D            ; % (molN cell-1) All the proteins in N (193-26)
Nprotein_const = Nphoto_const + Nconst_protein;  % (molN cell-1) All the proteins in N (193-26)
Nrna_D = Nprotein_const.*Cnrna_variable;          % (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
Nrna_D2 = Nprotein_D.*Cnrna_variable  ;           % (molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)

%Constant carbon parameters------------------
Cconst_protein = Nconst_protein*CNprotein;  %(molC cell-1) carbon in other protein assumed constant (195-16)
Cdna_const = Ndna_const*YdnaC_N  ;    %(molC cell-1) carbon in constant part of DNA (195-16)

%Calculating factors and Qc_essential-------------

Qc_D2 = A.*Cnrna_variable.*YrnaC_N;
Qc_D = (1+E)*Qc./Pchl + A.*CNprotein + L*YpgC_P + B.*Cnrna_variable*YrnaC_N;
Qc_const = m./Pchl + B*CNprotein + M*YpgC_P + Nrna_const*YrnaC_N + Ndna_const*YdnaC_N + Cessential;

Qc_essential = Qc_essential_computation(D,Chl_const, Chl_D, Nchl_const, Nchl_D, Nphoto_const, Nphoto_D, Nbiosynth_D,...
    Nprotein_const, Nprotein_D, Nrna_const, Nrna_D, Nrna_D2, Ypthylakoid_chl, YrnaC_N,...
    CNprotein, YpgC_P, Cdna_const, Cconst_protein, Cessential);
if ~isvector(D)
    msk = (Qc_essential >= Qc & ~isnan(Qc_essential));
    [idx,idy] = find(msk==1);
    for ii = 1:length(idx)
        DFe(idx(ii),idy(ii)) = 2.*(Qc - Qc_const(idx(ii),idy(ii)))./...
            (Qc_D(idx(ii),idy(ii)) + sqrt(Qc_D(idx(ii),idy(ii)).*Qc_D(idx(ii),idy(ii))...
            + 4.*Qc_D2(idx(ii),idy(ii)).*(Qc - Qc_const(idx(ii),idy(ii)))));
    end
    qfe=(A_pho_Fe.*Chl_D.*(86400.*DFe))+(A_pho_Fe.*Chl_const);
    mu=DFe;
    D = DFe;
else
i0 = Qc_essential >= Qc   ;%conditions where Mu max applies
i1 = isnan(Qc_essential) ; %Including this avoide when Qc_essential is nan with high Vn
i22 = i0 + i1;
D(i22==1) = 2*(Qc - Qc_const)/(Qc_D + sqrt(Qc_D*Qc_D + 4*Qc_D2*(Qc - Qc_const)));
DFe = D;
qfe=(A_pho_Fe.*Chl_D.*(86400.*DFe))+(A_pho_Fe.*Chl_const);
mu=DFe;

end

mumax=2.*(Qc - Qc_const)./...
            (Qc_D + sqrt(Qc_D.*Qc_D...
            + 4.*Qc_D2.*(Qc - Qc_const)));
Numtoarray=ones(size(DFe)) ;
Chl       = Chl_const      + Chl_D.*DFe;% masked aray, ~4.4 e-11
Nchl      = Nchl_const     + Nchl_D.*DFe; %masked array, ~3.2 e-12             %(molN chl cell-1) Chlorophyll N concentration
Nphoto    = Nphoto_const   + Nphoto_D.*DFe; %masked array, ~5.6 e-10    % Nphoto_const should be like e-16  %(molN cell-1) Photosynthesis related protein nitrogen (193-25)6
Nbiosynth =                  Nbiosynth_D.*DFe;% masked array, ~1.7 e-10
Nprotein  = Nprotein_const + Nprotein_D.*DFe ; %masked array, ~7.4 e-10     %(mol N cell-1) all the protein
Nrna      = Nrna_const     + Nrna_D.*DFe + Nrna_D2.*DFe.*DFe;%masked array, 1.8 ~e-6
Nessential= Nchl + Nphoto + Nbiosynth + Nconst_protein + Nrna + Ndna; %masked array, ~1.8 e -6


Nchl=Chl*YchlN_C ;       % %(molN chl cell-1) Chlorophyll N concentration
Nphoto=Chl*Ynphoto_chl;  %%(molN cell-1) Photosynthesis related protein nitrogen (193-25)
Nbiosynth=D.*Cnbiosynth;        %     %(molN cell-1) various part of biosynthesis related protein in N (193-37)
%Nprotein=Nphoto+Nconst_protein+Nbiosynth ;  % %(molN cell-1) All the proteins in N (193-26)
Nrna_variable=Nprotein .* D .* Cnrna_variable;  %      %(molN cell-1) variable part of nitrogen in RNA (193-26)(193-37)
Nrna=Nrna_const + Nrna_variable ;%
Ndna_variable=Ndna_const*D/1.2*(18.3-7.6)/7.6     ;%   %(molN cell-1) variable part of nitrogen in DNA (193-26) Increasing ratio based on Bremmer 1996

Pthylakoid=Chl*Ypthylakoid_chl         ;% %(molP cel-1) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)
Prna_variable=Nrna_variable*YnucacidP_N  ;%   %(molP cell-1) variable part of phosphorus in RNA (193-26)
Pdna_variable=Ndna_variable*YnucacidP_N    ;% %(molP cell-1) variable part of phosphorus in DNA (193-26)

Cchl = Chl  ; %masked array, 4.4 e-11                 %(molC cell-1) carbon in chlorophyll (195-16)
Cphoto = Nphoto*CNprotein; %masked array, 2.5 e-9 %Nphoto should be like e-15     %(molC cell-1) carbon in photosystem protein (195-16)
Cbiosynth = Nbiosynth*CNprotein ;%masked array, 7.99 e-10  %(molC cell-1) carbon in biosynthesis protein (195-16)
Crna = Nrna*YrnaC_N; %masked array, 5.3 e-6     %(molC cell-1) carbon RNA (195-16)
CthylakoidPG = Pthylakoid*YpgC_P ; %masked array, 5.0 e-11         %(molC cell-1) carbon in PG (phosphatidyl glycerol) in thylakoid membranes
Cnstore = nan(size(DFe));%Nstore*YcyanoC_N

Cother = Qc - Cphoto - Cbiosynth - Cconst_protein - Cchl...
    - Crna - Cdna_const...
    - Cessential - CthylakoidPG ;%- Cnstore ;%masked array, -5.2 e-6
percentorratio=100  ;     %100: percent, 1:ratio;
Cphoto_plot=Cphoto/Qc*percentorratio ;  %masked array, 3 e+6  %Cphoto should be like e -14/e -15
Cbiosynth_plot=Cbiosynth/Qc*percentorratio; %masked array, 9e+5
Cconst_protein_plot=Cconst_protein/Qc*percentorratio*Numtoarray; %not masked array, 23 (all same value)
Cchl_plot=Cchl/Qc*percentorratio ;%masked array, 5 e+4
Crna_plot=Crna/Qc.*percentorratio.*Numtoarray ;%masked array, 6.4 e+9
Cdna_const_plot=Cdna_const/Qc*percentorratio*Numtoarray ;%not masked array, 0.09 (all same)
Cother_plot=Cother/Qc*percentorratio ;%masked array, -6.4 e+9
Cessential_plot=Cessential/Qc*percentorratio*Numtoarray;  %not masked array, 1.82 (all same)
CthylakoidPG_plot=CthylakoidPG/Qc*percentorratio ;%masked array, 6 e+4
Nconst_protein_plot=Nconst_protein*Nunit*Numtoarray; %array, .05    %(ug N/ mgC) constant protein pool in nitrogen (193-25)(193-33)
Nphoto_plot=Nphoto*Nunit; %array, 0.01    %(ug N/ mgC) Photosynthesis related protein nitrogen (193-25)(193-33)
Nbiosynth_plot=Nbiosynth*Nunit; %array, 0.0003     %(ug N/ mgC) biosynthesis related protein in N (193-37)
Ndna_const_plot=Ndna_const*Nunit*Numtoarray; %array, 0.0003 (values are same)   %(ug N/ mgC) Nitrogen in constant part of DNA
Nrna_plot=Nrna*Nunit; %array, 0.0008   %(ug N/ mgC) Nitrogen in constant part of RNA
Nchl_plot=Nchl*Nunit ;  % array, 0.0002    %(ug N/ mgC) Chlorophyll nitrogen (actually almost negligiable) (193-33)
Prna_variable_plot=Prna_variable*Punit   ;   %(ug P/mgC) Phosphorus in variable part of RNA (193-37)
Pdna_variable_plot=Pdna_variable*Punit    ;  %(ug P/mgC) Phosphorus in variable part of DNA (193-37)
Pthylakoid_plot=Pthylakoid*Punit  ;   %(ug P/ mgC) Phosphorus in thylakoid membranes: phospholipid, etc. (193-26)(193-33)
Pconst_other_plot=Pconst_other*Punit*Numtoarray    ;   %(ug P/ mgC) PHosphorus in other parts (ex. phospholipid in outer membrane, ATP, ADP, etc. (assuming constant) (193-33)
Pdna_const_plot=Pdna_const*Punit*Numtoarray ;   %(ug P/ mgC) Phosphorus in constant part of DNA
Prna_const_plot=Prna_const*Punit*Numtoarray  ;  %(ug P/ mgC) Phosphorus in constant part of RNA
% Pstore_plot=Pstore*Punit;    %(ug P/ mgC) Phosphorus in phosphorus storage
%  Nstore[i]=Qn_test[i]-Nprotein[i]-Nrna_variable[i]-Nrna_const-Ndna_variable[i]-Ndna_const-Nchl[i]  #(molN cell-1) Nitrogen storage in the cell
Qn=Nprotein+Nrna_variable+Nrna_const+Ndna_variable+Ndna_const+Nchl   ;%         %(molN cell-1) total nitrogen in the cell without storage
Qp=Pconst_other+Pthylakoid+Prna_variable+Prna_const+Pdna_variable+Pdna_const   ;%   %(molP cell-1) total phosphorus in the cell without storage
%Nconst_protein_plot+Ndna_const_plot,Nphoto_plot+Nchl_plot,...
   % Nbiosynth_plot+Nrna_plot
Ness = Nconst_protein_plot+Ndna_const_plot;
Npho = Nphoto_plot+Nchl_plot;
Nbio = Nbiosynth_plot+Nrna_plot;

n2p=Qn./Qp.*14.*10^6./(30.97.*10^6)   ;%     %(ug N /ug P) biomass N to P ratio (164-20)
Cess = Cessential_plot+Cconst_protein_plot+Cdna_const_plot; %ess
Cpho = Cphoto_plot+CthylakoidPG_plot+Cchl_plot; % photo
Cbio = Crna_plot+Cbiosynth_plot; % bio
Csto = Cother_plot; %storage?
end
