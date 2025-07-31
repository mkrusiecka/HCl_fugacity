%C-H-O-S-Cl speciation written by M.K. Rusiecka Department of Earth
%Sciences University of Oxford 
%This program calculates water fugacity based on H2O activity in the melt
%(Burnham 1979); carbon species - Eguchi and Dasgupta 2018; S species -
%Boulliung and Wood 2022, 2023; HCl and Cl2 fugacity - Rusiecka and Wood (in
%preparation); Holloway-type MRK EOS is used to calculate mole fractions of
%fluid species and adjust fugacities accordingly;
%This program can be used to calculate composition of a gas that would be
%in equilibrium with a melt (for example a set of melt inlcusions; as input
%you have to set P, T, composition anf fO2 (fO2 is set starting line 87 -
%you choose a buffer (FMQ, NNO, CCO) and set offset, buffer choice is line
%93 and offset is line 94
function HCl_fugacity_with_S()
    clear; clc;

    %% Load data
    data = load('MI_Gas_Speciation.csv');
    pressure = data(:,18); % bar
    temperature = data(:,1); % °C
    data = data';

    % Reorder bulk composition
    bulkComp = zeros(19, size(data, 2));
    bulkComp([1:3 5:8 11:14, 15, 16, 17, 18, 19], :) = data([2:4 10 5 11 6:9 12 13 14, 15, 16, 17], :);

    T = temperature + 273.15; % K
    P = pressure; % bar
    PG = P / 10000; % GPa
    R = 8.314; % J/mol·K
    %%
    % Calculate oxide mole numbers 
    SiO2_mol=(bulkComp(1,:)./60.08)';
    TiO2_mol=(bulkComp(2,:)./79.866)';
    Al2O3_mol=(bulkComp(3,:)./(101.96/2))';
    FeO_mol=(bulkComp(6,:)./71.844)';
    MgO_mol=(bulkComp(8,:)./40.3044)';
    CaO_mol=(bulkComp(11,:)./56.0774)';
    Na2O_mol=(bulkComp(12,:)./(61.97/2))';
    K2O_mol=(bulkComp(13,:)./(94.2/2))';
    Cr2O3_mol=(bulkComp(9,:)./(75.995/2))';
    MnO_mol=(bulkComp(10,:)./70.9374)';
    P2O5_mol=(bulkComp(14,:)./(283.889/2))';
    H2O_mol=(bulkComp(15,:)./(18.01528/2))';
    Cl_mol=(bulkComp(18,:)./35.45)';
    CatSum = SiO2_mol + TiO2_mol + Al2O3_mol + FeO_mol + MgO_mol + CaO_mol + Na2O_mol + K2O_mol + Cr2O3_mol + MnO_mol + P2O5_mol + H2O_mol;
    ntot = SiO2_mol + TiO2_mol + (bulkComp(3,:)./(101.96))' + FeO_mol + CaO_mol + (bulkComp(12,:)./(61.97))' + (bulkComp(13,:)./(94.2))' + (bulkComp(9,:)./(75.995))' + (bulkComp(14,:)./(283.889))' + (bulkComp(15,:)./(18.01528))' + MnO_mol; %sum of moles for Kress and Carmichael 1991
    % Mole fractions
    XCaO = CaO_mol ./ CatSum;
    XSiO2 = SiO2_mol ./ CatSum;
    XFeO = FeO_mol ./ CatSum;
    XK2O = K2O_mol ./ CatSum;
    XMgO = MgO_mol ./ CatSum;
    XNa2O = Na2O_mol ./ CatSum;
    XAl2O3 = Al2O3_mol ./ CatSum;
    XMnO = MnO_mol./CatSum;
    XCl=Cl_mol./CatSum;
    Cl_mes = bulkComp(18,:)'; % Cl measured in wt%
    totalCO2_wt=bulkComp(16,:)'; %CO2 in wt. %
    totalCO2_wt(totalCO2_wt <= 0 | isnan(totalCO2_wt)) = 1e-20; %if there is no CO2 data or CO2=0 it is getting replaced by CO2(wt.%)=1e-20 so the code can still work
    % --- New: Calculate oxygen moles ---
    O_moles = 2*SiO2_mol + 2*TiO2_mol + 3*(bulkComp(3,:)./(101.96))' + FeO_mol + MgO_mol + CaO_mol + ...
          (bulkComp(12,:)./(61.97))' + (bulkComp(13,:)./(94.2))' + 3*Cr2O3_mol + MnO_mol + 5*(bulkComp(14,:)./(283.889))'+ (bulkComp(15,:)./(18.01528))';

    %% Calculate H2O mole fraction on 8 oxygen basis for Burnham (1979) model
    X_H2O_O = H2O_mol .*(8./O_moles);
    %% Chloride capacity (CCl)
    logC = 1.127 + (4784.*XCaO - 2975.*XSiO2 + 2705.*XFeO + 541.*XMgO - 2661.*XK2O - 605.*PG) ./ T;
    CCl_calc = 10.^logC;

    %% fH2O from Burnham 1979 model
    Po = 1.01325; % bar
    lnkw = 5.83 + log(P) .* (4.481e-8 .* T.^2 - 1.51e-4 .* T - 1.137) + ...
           (log(P)).^2 .* (1.831e-8 .* T.^2 - 4.882e-5 .* T + 0.04656) + ...
           7.80e-3 .* (log(P)).^3 - 5.012e-4 .* (log(P)).^4 + ...
           T .* (4.754e-3 - 1.621e-6 .* T);
    kw = exp(lnkw);
   % --- Modify water activity (aH2O) calculation to use oxygen basis mole fraction X_H2O_O instead of H2O_mol./CatSum ---
aH2O = zeros(size(T));
for idx = 1:length(T)
    if X_H2O_O(idx) >= 0.5
        aH2O(idx) = 0.25 * kw(idx) * exp((6.52 - 2667 / T(idx)) * (X_H2O_O(idx) - 0.5));
    else
        aH2O(idx) = (X_H2O_O(idx))^2 * kw(idx);
    end
end
    fH2O = aH2O .* (P ./ Po);
%% Oxygen fugacity buffer (3 options + user offset)
% Choose one of: 'FMQ', 'NNO', 'CCO'
% You can either:
%   (a) set buffer_choice='FMQ' and buffer_offset=+1;  or
%   (b) use a combined string like 'FMQ+1', 'NNO-0.5', 'CCO+0.3'
% Assumptions: T in Kelvin, P in bar; (P-1)/T term uses 1 bar reference.
buffer_choice = 'NNO';         % 'FMQ' | 'NNO' | 'CCO'  (or 'FMQ+1' style)
buffer_offset = 0;             % numeric Δlog10 fO2 relative to the buffer

% --- input checks (optional) ---
assert(isnumeric(T) && isnumeric(P), 'T and P must be numeric.');
assert(all(T(:) > 0), 'T must be in Kelvin and > 0.');
if ~isscalar(T) && ~isscalar(P)
    assert(isequal(size(T), size(P)), 'T and P must be the same size or one must be scalar.');
end

% --- allow combined syntax like 'FMQ+1' or 'NNO-0.5' (robust to spaces/case/unicode minus) ---
s = regexprep(strtrim(buffer_choice), char(8722), '-');  % U+2212 → '-'
tok = regexp(lower(s), '^\s*(fmq|nno|cco)\s*([+-]\s*\d*\.?\d+)?\s*$', 'tokens', 'once');
if isempty(tok)
    error('buffer_choice must be FMQ, NNO, or CCO, optionally with an offset like FMQ+1.');
end
base = upper(tok{1});
if numel(tok) > 1 && ~isempty(tok{2})
    buffer_offset = buffer_offset + str2double(regexprep(tok{2}, '\s+', ''));
end

% --- base log10(fO2) for each buffer ---
switch base
    case 'FMQ'
        base_logfO2 = (-25096.3 ./ T + 8.735 + 0.110 .* (P - 1) ./ T);
    case 'NNO'
        base_logfO2 = (-24930   ./ T + 9.36  + 0.046 .* (P - 1) ./ T);
    case 'CCO'
        base_logfO2 = (4.325 - (21803 ./ T) + 0.171 .* ((P - 1) ./ T) - 0.2);
end

% --- apply user offset and exponentiate ---
logfO2 = base_logfO2 + buffer_offset;  % log10 fO2 relative to chosen buffer
fO2    = 10.^logfO2;                   % fO2 in bar
    %% fCl2 from Cl capacity
    fCl2 = (Cl_mes ./ CCl_calc).^2 .* (fO2.^0.5);

    %% fHCl from equilibrium
    logKf_H2O = 12850 ./ T - 2.8675;
    KH2O = 10.^logKf_H2O;

    logKf_HCl = 4894.6 ./ T + 0.3436;
    KHCl = 10.^logKf_HCl;

    fHCl = KHCl .* ((fH2O.^0.5) ./ (KH2O.^0.5)) .* (Cl_mes ./ CCl_calc);

    %% fH2 from H2O equilibrium
    fH2 = fH2O ./ (KH2O .* sqrt(fO2));
    %% Kress and Carmichael 1991
    X_SiO2=SiO2_mol/ntot; %mole fraction of SiO2
    X_TiO2=TiO2_mol/ntot; %mole fraction of TiO2
    X_Al2O3=(bulkComp(3,:)./(101.96))'/ntot; %mole fraction of Al2O3
    X_FeOt=FeO_mol/ntot; %mole fraction of TOTAL FeO
    X_MnO=MnO_mol/ntot; %mole fraction of MnO
    X_MgO=MgO_mol/ntot; %mole fraction of MgO
    X_CaO=CaO_mol/ntot; %mole fraction of CaO
    X_Na2O=(bulkComp(12,:)./(61.97))'/ntot; %mole fraction of Na2O
    X_K2O=(bulkComp(13,:)./(94.2))'/ntot; %mole fraction of K2O
    X_P2O5=(bulkComp(14,:)./(283.889))'/ntot; %mole fraction of P2O5
    X_H2O=(bulkComp(15,:)./(18.01528))'./ntot; %mole fraction of H2O
    a_KC=0.196;
    b_KC=11492;
    c_KC=-6.675;
    d_Al=-2.243;
    d_Fe=-1.828;
    d_Na=5.854;
    d_K=6.215;
    d_Ca=3.201;
    lnXFe3XFe2=a_KC.*log(fO2)+b_KC./T+c_KC+(d_Al.*(bulkComp(3,:)./(101.96))'+d_Na.*(bulkComp(12,:)./(61.97))'+d_Fe.*XFeO+d_K.*(bulkComp(13,:)./(94.2))'+d_Ca.*XCaO);
    XFe3XFe2=exp(lnXFe3XFe2);
    n_Fe2O3=(ntot.*XFe3XFe2.*FeO_mol)./(ntot).*(1./(1+(ntot.*XFe3XFe2.*2)./(ntot)));
    n_FeO=FeO_mol-2.*(n_Fe2O3);
    X_Fe2O3=n_Fe2O3./ntot;
    X_FeO=n_FeO./ntot;
    XFeO4=2.*((XFe3XFe2.*X_FeO)./(1.0776-XFe3XFe2.*0.0776));
    XFeO5=(XFeO4.*1.0776)./(XFe3XFe2);
    wtFeO_KC=XFeO5.*(55.845+15.999)./((bulkComp(1,:))'+(bulkComp(2,:))'+(bulkComp(3,:))'+(bulkComp(10,:))'+(bulkComp(8,:))'+(bulkComp(11,:))'+(bulkComp(12,:))'+(bulkComp(13,:))'+(bulkComp(14,:))'+(bulkComp(15,:))'+n_FeO.*(55.845+15.999)+n_Fe2O3.*(55.845*2+15.999*3)).*100;
    wtFe2O3_KC=XFeO4.*(79.846)./((bulkComp(1,:))'+(bulkComp(2,:))'+(bulkComp(3,:))'+(bulkComp(10,:))'+(bulkComp(8,:))'+(bulkComp(11,:))'+(bulkComp(12,:))'+(bulkComp(13,:))'+(bulkComp(14,:))'+(bulkComp(15,:))'+n_FeO.*(55.845+15.999)+n_Fe2O3.*(55.845*2+15.999*3)).*100;
    FeO_mol_KC=wtFeO_KC./(55.845+15.999);
    Fe2O3_mol_KC=wtFe2O3_KC./(159.66/2);
    CatSum_KC = SiO2_mol + TiO2_mol + Al2O3_mol + FeO_mol_KC + Fe2O3_mol_KC + MgO_mol + CaO_mol + Na2O_mol + K2O_mol + Cr2O3_mol + MnO_mol + P2O5_mol + H2O_mol;
    % Mole fractions
    XCaO_KC = CaO_mol ./ CatSum_KC;
    XSiO2_KC = SiO2_mol ./ CatSum_KC;
    XFeO_KC = FeO_mol_KC ./ CatSum_KC;
    XK2O_KC = K2O_mol ./ CatSum_KC;
    XMnO_KC = MnO_mol./CatSum_KC;
    %% Calculate S from Boulliung and Wood 2022 & 2023
%% Sulfur: use BOTH sulfide and sulfate capacities to solve fS2 (base-10 logs)
% Measured total S in the melt (wt%)
S_mes = bulkComp(17,:)'; 
S_mes(S_mes <= 0) = 1e-10; % avoid log(0)

% ---------- (A) Sulfide capacity (as you had, log10 form) ----------
% NOTE: these X_*_KC are on your KC-normalized basis (CatSum_KC). Keep consistent with the calibration you used.
logCS2 = 0.225 + (25237.*XFeO_KC + 5214.*XCaO_KC + 12705.*XMnO_KC + 19828.*XK2O_KC - 1109.*XSiO2_KC - 8879) ./ T ...
                  + ((P-1).*6.2e-3)./(2.303.*R.*T);   % log10 C_S2-

% ---------- (B) Sulfate capacity (your new relation, log10 base) ----------
% Use anhydrous oxide mole-fraction basis (commonly used in capacity calibrations).
CatSum_anh = SiO2_mol + TiO2_mol + Al2O3_mol + FeO_mol + MgO_mol + CaO_mol + ...
             Na2O_mol + K2O_mol + Cr2O3_mol + MnO_mol + P2O5_mol;   % exclude H2O

XNa2O_anh  = Na2O_mol  ./ CatSum_anh;
XCaO_anh   = CaO_mol   ./ CatSum_anh;
XMgO_anh   = MgO_mol   ./ CatSum_anh;
XMnO_anh   = MnO_mol   ./ CatSum_anh;
XAl2O3_anh = Al2O3_mol ./ CatSum_anh;

% Your supplied formula (log10 base):
% logCS6 = -213.65 + (25696 XNa2O + 15076 XCaO + 9543 XMgO + 16158 XMnO + 4316 XAl2O3 + 68254)/T + 55.03 log10(T)
logCS6 = -213.65 + ...
         (25696.*XNa2O_anh + 15076.*XCaO_anh + 9543.*XMgO_anh + 16158.*XMnO_anh + 4316.*XAl2O3_anh + 68254)./T + ...
         55.03.*log10(T);

% ---------- (C) Solve mass balance S_tot = S(2-) + SO4(2-) for fS2 ----------
% Equilibrium constant for SO2 formation (log10)
logK_SO2 = 18880./T - 3.8018;
lt  = @(x) log10(max(x, realmin));   % safe log10
ten = @(x) 10.^x;

n = numel(T);
fS2  = zeros(n,1);
fSO2 = zeros(n,1);

for i = 1:n
    logfO2_i  = lt(fO2(i));
    logCS2_i  = logCS2(i);
    logCS6_i  = logCS6(i);
    S_tot_i   = S_mes(i);          % wt% (ensure same units as capacities)

    % Mass-balance residual in wt% (unknown is log10 fS2)
    % S2- (wt%)  = 10^logCS2 * fS2^(1/2) * fO2^(-1/2)
    % SO4 (wt%)  = 10^logCS6 * fS2^(1/2) * fO2^(+3/2)   [from your definition]
    massbal = @(logfS2) ...
        ( ten(logCS2_i) .* ten( 0.5*logfS2 - 0.5*logfO2_i ) ) + ...
        ( ten(logCS6_i) .* ten( 0.5*logfS2 + 1.5*logfO2_i ) ) - S_tot_i;

    % Bracket and solve
    lo = -30; hi = 5;
    if massbal(lo)*massbal(hi) > 0
        lo = lo - 10; hi = hi + 5;
    end
    logfS2_i = fzero(massbal, [lo, hi]);
    fS2(i)   = ten(logfS2_i);

    % Gas equilibrium linking SO2 (log10 fSO2 = logK_SO2 + 0.5 log fS2 + log fO2)
    fSO2(i)  = ten(logK_SO2(i)) * sqrt(max(fS2(i), realmin)) * max(fO2(i), realmin);
end

% ---------- (D) Now H2S from your equilibrium (unchanged) ----------
logK_H2S = 27103./T - 4.1973;  % (you already defined above, keep it consistent)
logfH2S = lt(fH2O) + lt(fSO2) - 1.5*lt(fO2) - logK_H2S;
fH2S    = 10.^logfH2S;

% ---------- (E) Optional diagnostics you might find handy ----------
S2m_wt   = 10.^logCS2 .* sqrt(fS2) ./ sqrt(fO2);      % predicted melt S2- (wt%)
SO4m_wt  = 10.^logCS6 .* sqrt(fS2) .* (fO2.^(1.5));   % predicted melt SO4(2-) (wt%)
frac_SO4 = SO4m_wt ./ max(S2m_wt + SO4m_wt, eps);     % fraction of melt S as sulfate
    % fSO2 and fH2S from equilibrium (source of K Boulliung and Wood 2023
    % Sulfur oxidation state and solubility in silicate melts
     %% NBO for CO2 model
    ox_total=(bulkComp(1,:))'+(bulkComp(2,:))'+(bulkComp(3,:))'+(bulkComp(10,:))'+(bulkComp(8,:))'+(bulkComp(11,:))'+(bulkComp(12,:))'+(bulkComp(13,:))'+(bulkComp(14,:))'+(bulkComp(15,:))'+wtFe2O3_KC+wtFeO_KC;
    wtSiO2=(bulkComp(1,:))'./ox_total.*100;
    wtTiO2=(bulkComp(2,:))'./ox_total.*100;
    wtAl2O3=(bulkComp(3,:))'./ox_total.*100;
    wtMgO=(bulkComp(8,:))'./ox_total.*100;
    wtMnO=(bulkComp(10,:))'./ox_total.*100;
    wtCaO=(bulkComp(11,:))'./ox_total.*100;
    wtNa2O=(bulkComp(12,:))'./ox_total.*100;
    wtK2O=(bulkComp(13,:))'./ox_total.*100;
    wtP2O5=(bulkComp(14,:))'./ox_total.*100;
    wtH2O=(bulkComp(15,:))'./ox_total.*100;
    %cations
    XSi=wtSiO2./60.08;
    XTi=wtTiO2./79.866;
    XAl=wtAl2O3./101.96.*2;
    XMg=wtMgO./40.4044;
    XMn=wtMnO./70.9374;
    XCa=wtCaO./56.0774;
    XNa=wtNa2O./61.97.*2;
    XK=wtK2O./94.2.*2;
    XP=wtP2O5./283.889.*2;
    XH=wtH2O./18.01528.*2;
    XFe2=FeO_mol_KC;
    XFe3=Fe2O3_mol_KC.*2;
    O_tot=XSi.*2+XTi.*2+wtAl2O3./101.96.*3+XNa./2+XK./2+wtP2O5./283.889.*5+XMg+XMn+XCa+XH./2+XFe2+Fe2O3_mol_KC.*3;
    %% carbon species Eguchi and Dasgupta 2018
    % --- Carbon formation constants (log10) ---
% C(s) + O2(g)   = CO2(g)
logKf_CO2 = 20599./T + 0.0531;
% C(s) + 1/2 O2(g) = CO(g)
logKf_CO  = 5839.7./T + 4.5515;
   % --- Precompute fCO2 from melt model, and capture carbon-saturation phase ---
n = length(T);
fCO2 = zeros(n,1);
phase = ones(n,1);   % 1 = fluid-sat (no C), 2 = graphite, 3 = diamond

for i = 1:n
    [fCO2_i, phase_i, ~] = fco2_from_totalC( ...
        P(i), T(i), totalCO2_wt(i), ...
        wtSiO2(i), wtTiO2(i), wtAl2O3(i), wtFe2O3_KC(i), wtFeO_KC(i), wtMnO(i), wtMgO(i), ...
        wtCaO(i), wtNa2O(i), wtK2O(i), wtP2O5(i), fO2(i));

    % If carbon-saturated, override fCO2 using the formation constant
    if phase_i ~= 1
        fCO2_i = 10.^logKf_CO2(i) * max(fO2(i), realmin);   % fCO2 = Kf_CO2 * fO2
        fCO2_i = min(fCO2_i, 0.95*P(i));                    % cap
    end
    fCO2(i) = fCO2_i;
    phase(i) = phase_i;
end

    %% Calculate mole fractions from MRK EOS
fprintf('\nCalculating mole fractions with MRK + self-consistent H2/H2O coupling...\n');

species = {'HCl','H2O','Cl2','H2','SO2','H2S','CO2','CO'};
crit.Tc = [324.7; 647.3; 416.9;  33.2; 430.8; 373.2; 304.13; 132.9];
crit.Pc = [ 83.1; 220.64; 76.1; 12.97;  78.8;  89.6;  73.77; 34.94];

Nsp  = numel(species);
kIJ  = zeros(Nsp);  % or your actual BIP matrix
opts = struct('tol',1e-8,'maxIter',200,'verbose',false);

iHCl = 1; iH2O = 2; iCl2 = 3; iH2 = 4; iSO2 = 5; iH2S = 6; iCO2 = 7; iCO = 8;

mole_fractions = zeros(length(T), Nsp);
% --- EOS-consistent fugacities to report (computed as y_i * phi_i * P) ---
fHCl_eos = zeros(length(T),1);
fCl2_eos = zeros(length(T),1);
fH2_eos  = zeros(length(T),1);
fSO2_eos = zeros(length(T),1);
fH2S_eos = zeros(length(T),1);
fH2O_eos = zeros(length(T),1);
fCO2_eos = zeros(length(T),1);
fCO_eos  = zeros(length(T),1);
for i = 1:length(T)

    % Equilibrium fugacities not depending on y
    fHCl_i = fHCl(i);
    fCl2_i = fCl2(i);
    fSO2_i = fSO2(i);
    fH2S_i = fH2S(i);
    fCO2_i = fCO2(i);

    % Initial guesses for fH2, fH2O gas, and fCO
    fH2_i  = max(1e-3*P(i), 1.0);
    fH2O_g = max(1e-2*P(i), 10.0);
    fCO_i  = max(1e-3*P(i), 1.0);

    fvec = zeros(Nsp,1);
    fvec([iHCl iH2O iCl2 iH2 iSO2 iH2S iCO2 iCO]) = [fHCl_i; fH2O_g; fCl2_i; fH2_i; fSO2_i; fH2S_i; fCO2_i; fCO_i];

    outer_max = 50; outer_tol = 5e-3; damp = 0.5;
    converged_outer = false;

% keep previous fugacities to monitor convergence
f_prev = fvec;  

for it = 1:outer_max

    % --- EOS step: invert MRK for y, phi at current target fugacities ---
    [y_i, phi_i, Z_i, info_i] = calc_mole_fractions_MRK(fvec, P(i), T(i), crit, kIJ, opts);

    % --- EOS-consistent fugacities from current composition ---
    fH2O_g_new = y_i(iH2O) * phi_i(iH2O) * P(i);
    if ~isfinite(fH2O_g_new) || fH2O_g_new <= 0, fH2O_g_new = max(1e-12, realmin); end

    % --- H2 from H2O–H2–fO2 equilibrium (uses KH2O at this T) ---
    KH2O_i = KH2O(i);               % already precomputed as 10.^(12702/T - 2.8502)
    fH2_new = fH2O_g_new / (KH2O_i * sqrt(fO2(i)));

    % cap and sanitize
    fH2_cap = 0.95 * P(i);
    if ~isfinite(fH2_new) || fH2_new <= 0, fH2_new = max(1e-6*P(i), 1.0); end
    if fH2_new > fH2_cap, fH2_new = fH2_cap; end

    % ---------- CO from formation constants ----------
    fCO2_i = fCO2(i);   % use the precomputed CO2 target for this row
    % Use the ratio when unsaturated; direct Kf when carbon-saturated
    Kratio = 10.^(logKf_CO2(i) - logKf_CO(i));   % = Kf_CO2 / Kf_CO  (log10 base)
    if phase(i) ~= 1
    % carbon present: fCO = Kf_CO * fO2^0.5
        fCO_new = 10.^logKf_CO(i) * sqrt(max(fO2(i), realmin));
    else
     % unsaturated: fCO = fCO2 / [(Kf_CO2/Kf_CO)*fO2^0.5]
        fCO_new = fCO2_i / (Kratio * sqrt(max(fO2(i), realmin)));
    end
    fCO_new = min(max(fCO_new, 0), 0.95*P(i));

    % SO2 from: log10(fSO2) = logK_SO2 + 0.5*log10(fS2) + log10(fO2)
    logK_SO2_i = logK_SO2(i);
    fSO2_new = 10^(logK_SO2_i) * sqrt(max(fS2(i), realmin)) * max(fO2(i), realmin);

    % H2S from: log10(fH2S) = log10(fH2O) + log10(fSO2) - 1.5*log10(fO2) - logK_H2S
    logK_H2S_i = logK_H2S(i);
    fH2S_new = ( fH2O_g_new * max(fSO2_new, realmin) ) / ( (max(fO2(i), realmin))^1.5 * 10^(logK_H2S_i) );

    % HCl from: fHCl = KHCl * sqrt(fH2O/KH2O) * (Cl_mes/CCl)
    KHCl_i = KHCl(i);
    fHCl_new = KHCl_i * sqrt( fH2O_g_new / KH2O_i ) * ( Cl_mes(i) / CCl_calc(i) );

    % Basic sanitization to avoid negatives/NaN
    if ~isfinite(fSO2_new) || fSO2_new < 0, fSO2_new = max(1e-12, realmin); end
    if ~isfinite(fH2S_new) || fH2S_new < 0, fH2S_new = max(1e-12, realmin); end
    if ~isfinite(fHCl_new) || fHCl_new < 0, fHCl_new = max(1e-12, realmin); end

    % (Optional) damping to improve stability, like you do for H2/H2O/CO
    fH2_i  = (1-damp)*fH2_i  + damp*fH2_new;
    fH2O_g = (1-damp)*fH2O_g + damp*fH2O_g_new;
    fCO_i  = (1-damp)*fCO_i  + damp*fCO_new;
    fSO2_i = (1-damp)*fSO2_i + damp*fSO2_new;   % NEW
    fH2S_i = (1-damp)*fH2S_i + damp*fH2S_new;   % NEW
    fHCl_i = (1-damp)*fHCl_i + damp*fHCl_new;   % NEW
    % ========================= END NEW/UPDATED PART ========================

    % --- Refresh target fugacity vector for next MRK inversion ---
    fvec([iHCl iH2O iCl2 iH2 iSO2 iH2S iCO2 iCO]) = ...
        [fHCl_i; fH2O_g; fCl2_i; fH2_i; fSO2_i; fH2S_i; fCO2_i; fCO_i];

    % --- Convergence checks: EOS consistency and chem-coupled f's ---
    % (1) Σ f/(φP) ≈ 1
    err_fp = abs(info_i.sum_f_over_phiP - 1);

    % (2) relative change in key coupled fugacities
    rel = @(x,y) abs(x - y) / max(y, 1e-12);
    err_chem = max([ ...
        rel(fH2_i,  f_prev(iH2)); ...
        rel(fH2O_g, f_prev(iH2O)); ...
        rel(fCO_i,  f_prev(iCO)); ...
        rel(fSO2_i, f_prev(iSO2)); ...
        rel(fH2S_i, f_prev(iH2S)); ...
        rel(fHCl_i, f_prev(iHCl)) ]);

    if (err_fp < 1e-3) && (err_chem < outer_tol)
        converged_outer = true;
        break
    end

    if it == 1
        fprintf('Row %d seed: Σ f/(φP)=%.2f, yH2O=%.3g, yH2=%.3g\n', ...
            i, info_i.sum_f_over_phiP, y_i(iH2O), y_i(iH2));
    end

    f_prev = fvec; % update for next iteration’s error
end

if ~converged_outer
    warning('Row %d: MRK+chem outer loop did not fully converge (Σ f/(φP)=%.2f after %d iters).', ...
        i, info_i.sum_f_over_phiP, it);
end
% Store EOS-consistent fugacities (y * phi * P) for reporting
fHCl_eos(i) = y_i(iHCl) * phi_i(iHCl) * P(i);
fCl2_eos(i) = y_i(iCl2) * phi_i(iCl2) * P(i);
fH2_eos(i)  = y_i(iH2)  * phi_i(iH2)  * P(i);
fSO2_eos(i) = y_i(iSO2) * phi_i(iSO2) * P(i);
fH2S_eos(i) = y_i(iH2S) * phi_i(iH2S) * P(i);
fH2O_eos(i) = y_i(iH2O) * phi_i(iH2O) * P(i);
fCO2_eos(i) = y_i(iCO2) * phi_i(iCO2) * P(i);
fCO_eos(i)  = y_i(iCO)  * phi_i(iCO)  * P(i);
    mole_fractions(i,:) = y_i.';
end
    % Optional: internal consistency check
    if abs(info_i.sum_f_over_phiP - 1) > 0.05
        warning('Row %d: Σ f/(φP) = %.3f at P=%.2f bar, T=%.1f K. Check equilibrium inputs/units.', ...
            i, info_i.sum_f_over_phiP, P(i), T(i));
    end
fprintf('   Temp (K) |   X(HCl) |  X(H2O) | X(Cl2) |  X(H2) |  X(SO2) |  X(H2S) |  X(CO2) | X(CO)\n');
for i = 1:length(T)
   fprintf('%10.1f  | %8.5f | %7.5f | %6.5f | %7.5f | %6.5f | %6.5f | %6.5f | %6.5f\n ', ...
    T(i), mole_fractions(i,:));
end
%% Display fugacity output (using EOS-consistent values)
safeLog10 = @(x) log10(max(x, realmin)); % avoid log10(0)

fprintf(['Temp (K)   |  log fO2 | log fCl2 | log fHCl |  log fH2 |  log fH2O | ', ...
         'log fSO2 | log fH2S | log fCO2 |   log fCO\n']);
for i = 1:length(T)
    fprintf('%10.3f | %8.3f | %8.3f | %8.3f | %9.3f | %9.3f | %8.3f | %8.3f | %9.3f | %9.3f\n', ...
        T(i), ...
        safeLog10(fO2(i)), ...
        safeLog10(fCl2_eos(i)), ...
        safeLog10(fHCl_eos(i)), ...
        safeLog10(fH2_eos(i)), ...
        safeLog10(fH2O_eos(i)), ...
        safeLog10(fSO2_eos(i)), ...
        safeLog10(fH2S_eos(i)), ...
        safeLog10(fCO2_eos(i)), ...
        safeLog10(fCO_eos(i)));
end
%% Plot fugacities (EOS-consistent)
figure; hold on;
plot(P, safeLog10(fO2),     'ko', 'DisplayName', 'log fO2');
plot(P, safeLog10(fCl2_eos),'bo', 'DisplayName', 'log fCl2');
plot(P, safeLog10(fHCl_eos),'ro', 'DisplayName', 'log fHCl');
plot(P, safeLog10(fH2O_eos),'mo', 'DisplayName', 'log fH2O');
plot(P, safeLog10(fH2_eos), 'go', 'DisplayName', 'log fH2');
plot(P, safeLog10(fSO2_eos),'c^', 'DisplayName', 'log fSO2');
plot(P, safeLog10(fH2S_eos),'ms', 'DisplayName', 'log fH2S');
plot(P, safeLog10(fCO2_eos),'bs', 'DisplayName', 'log fCO2');
plot(P, safeLog10(fCO_eos), 'gs', 'DisplayName', 'log fCO');
xlabel('Pressure (bar)'); ylabel('log Fugacity (bar)');
legend('Location', 'best'); title('Fugacity vs Pressure'); grid on;
   %% Plot mole fractions
figure; hold on;
plot(P, mole_fractions(:,1), 'ro', 'DisplayName', 'X HCl');
plot(P, mole_fractions(:,2), 'mo', 'DisplayName', 'X H2O');
plot(P, mole_fractions(:,3), 'bo', 'DisplayName', 'X Cl2');
plot(P, mole_fractions(:,4), 'go', 'DisplayName', 'X H2');
plot(P, mole_fractions(:,5), 'c^', 'DisplayName', 'X SO2'); % SO2
plot(P, mole_fractions(:,6), 'ms', 'DisplayName', 'X H2S'); % H2S
plot(P, mole_fractions(:,7), 'bs', 'DisplayName', 'X CO2'); % CO2
plot(P, mole_fractions(:,8), 'gs', 'DisplayName', 'X CO'); % CO
xlabel('Pressure (bar)');
ylabel('Mole Fraction');
legend('Location', 'best');
title('Gas Mole Fractions vs Pressure');
grid on;
% --- Build EOS fugacity struct from your variables ---
fEOS = struct( ...
    'fO2',  fO2(:), ...
    'fCl2', fCl2_eos(:), ...
    'fHCl', fHCl_eos(:), ...
    'fH2',  fH2_eos(:), ...
    'fH2O', fH2O_eos(:), ...
    'fSO2', fSO2_eos(:), ...
    'fH2S', fH2S_eos(:), ...
    'fCO2', fCO2_eos(:), ...
    'fCO',  fCO_eos(:) ...
    );

% --- Aux with your K(T) and capacity (use yours so bases/conventions match) ---
aux = struct();
aux.KH2O       = KH2O(:);
aux.KHCl       = KHCl(:);
aux.logK_SO2   = logK_SO2(:);
aux.logK_H2S   = logK_H2S(:);
aux.logKf_CO2 = logKf_CO2(:);
aux.logKf_CO  = logKf_CO(:);
aux.CCl_calc   = CCl_calc(:);
aux.Cl_mes     = Cl_mes(:);
aux.phaseFlag = phase;  % 1 = fluid-sat (no C), 2 = graphite, 3 = diamond
aux.opts       = struct('tol',1e-10,'maxIter',300,'verbose',false); % MRK solver opts

% --- Critical constants and kIJ you used for MRK ---
crit = struct();
crit.Tc = [324.7; 647.3; 416.9;  33.2; 430.8; 373.2; 304.13; 132.9];
crit.Pc = [ 83.1; 220.64; 76.1; 12.97;  78.8;  89.6;  73.77; 34.94];
kIJ  = zeros(8);  % or your calibrated matrix

% --- Run the validator ---
report = validate_MRK_Holloway(T(:), P(:), crit, kIJ, fEOS, aux);

% (Optional) act on failures
if ~report.pass
   warning('Validator failed: check report.table for details.');
end
end
function [y, phi, Z, info] = calc_mole_fractions_MRK(fug, P, T, crit, kIJ, opts)
% calc_mole_fractions_MRK
% Invert the Modified Redlich–Kwong (MRK) EOS for a vapor mixture at (P,T)
% to match target fugacities f_i (bar):  f_i = phi_i * y_i * P.
%
% MRK here = Redlich–Kwong cubic with a(T) ~ T^{-1/2}, linear b-mixing,
% quadratic a-mixing with optional binary interaction parameters k_ij.
%
% Inputs
%   fug  : N×1 fugacities (bar), order must match 'crit'
%   P    : pressure (bar)
%   T    : temperature (K)
%   crit : struct with fields (all N×1 vectors) .Tc (K), .Pc (bar), .omega (unused here), .names
%   kIJ  : N×N binary interaction parameters (default zeros)
%   opts : struct with fields: .tol (1e-8), .maxIter (200), .verbose (false)
%
% Outputs
%   y    : N×1 vapor mole fractions
%   phi  : N×1 fugacity coefficients
%   Z    : compressibility factor (vapor root)
%   info : diagnostics (iterations, convergence, sum_f_over_phiP, etc.)
 % Critical constants: [Tc (K), Pc (bar), ω]
    N = numel(fug);
    if nargin < 5 || isempty(kIJ), kIJ = zeros(N); end
    if nargin < 6, opts = struct; end
    if ~isfield(opts,'tol'),     opts.tol     = 1e-8; end
    if ~isfield(opts,'maxIter'), opts.maxIter = 200; end
    if ~isfield(opts,'verbose'), opts.verbose = false; end

    R = 0.08314472; % L·bar/(mol·K)

    Tc = crit.Tc(:); Pc = crit.Pc(:);
    if numel(Tc) ~= N || numel(Pc) ~= N
        error('crit.Tc and crit.Pc must be length N = numel(fug).');
    end
    if any(fug < 0)
        error('Fugacities must be non-negative.');
    end
    if any(diag(kIJ) ~= 0)
        error('Diagonal of kIJ must be zero.');
    end

    % --- Pure-component MRK parameters ---
    % a_i(T) = a0_i / sqrt(T),   a0_i = 0.42748 * R^2 * Tc^(2.5) / Pc
    % b_i    = 0.08664 * R * Tc / Pc
    a0_i = 0.42748 .* (R.^2) .* (Tc.^2.5) ./ Pc;
    a_i  = a0_i ./ sqrt(T);
    b_i  = 0.08664 .* R .* Tc ./ Pc;

    % Initial guess: ideal-gas y ~ fug / sum(fug)
    if all(fug == 0), error('All fugacities are zero.'); end
    y = fug(:) / sum(fug);

    iter = 0; converged = false; max_dy = NaN;

    while iter < opts.maxIter
        iter = iter + 1;
        y_old = y;

        % --- Mixing with BIPs ---
        % a_mix = sum_i sum_j y_i y_j (1 - k_ij) sqrt(a_i a_j)
        sqrt_a = sqrt(a_i);
        a_mix = 0.0;
        a_i_mix = zeros(N,1); % sum_j y_j (1 - k_ij) sqrt(a_i a_j)
        for i = 1:N
            a_i_mix(i) = sum( y(:)'.*(1 - kIJ(i,:)) .* (sqrt_a(i).*sqrt_a(:))' );
            a_mix = a_mix + y(i) * a_i_mix(i);
        end
        b_mix = sum(y .* b_i);

        % Dimensionless MRK parameters
        % NOTE: RK/MRK uses T^(2.5) in A
        A = a_mix * P / (R^2 * T^2.5);
        B = b_mix * P / (R * T);

        % Cubic in Z for RK/MRK:
        % Z^3 - Z^2 + (A - B - B^2) Z - A B = 0
        coeffs = [1, -1, A - B - B^2, -A*B];
        Zr = roots(coeffs);
        Zr = Zr(abs(imag(Zr)) < 1e-12);
        if isempty(Zr)
            error('No real Z root found.');
        end
        Z = max(real(Zr)); % vapor root

        if Z <= B
            % Avoid log singularities/numerical issues
            Z = max(B + 1e-12, Z);
        end

        % Fugacity coefficients (RK/MRK form)
        % ln(phi_i) = (b_i/b)(Z - 1) - ln(Z - B)
        %             - (A/B) * ( 2*(sum_j y_j a_ij)/a - b_i/b ) * ln(1 + B/Z)
        lnTerm = log(1 + B/Z);
        phi = zeros(N,1);
        for i = 1:N
            bi_b  = b_i(i)/b_mix;
            sum_a = a_i_mix(i); % already sum_j y_j (1 - k_ij) sqrt(a_i a_j)
            lnphi = bi_b*(Z - 1) - log(Z - B) ...
                    - (A/B) * ( 2*sum_a/a_mix - bi_b ) * lnTerm;
            phi(i) = exp(lnphi);
        end

        % Update y from f = phi*y*P
        denom = phi * P;
        if any(~isfinite(denom) | denom <= 0)
            error('Non-physical φ_i * P encountered.');
        end
        y = fug(:) ./ denom;
        y(y < 0) = 0;
        s = sum(y);
        if s == 0, error('All updated y are zero.'); end
        y = y / s;

        max_dy = max(abs(y - y_old));
        if max_dy < opts.tol
            converged = true;
            break;
        end
    end

    % Diagnostics
    info.iterations = iter;
    info.converged  = converged;
    info.max_dy     = max_dy;
    info.sum_f_over_phiP = sum(fug(:)./(phi(:)*P));
    info.note = 'Σ f/(φP) ≈ 1 indicates internal consistency of the fugacity set with MRK at (P,T).';

    if ~converged
        warning('MRK mixture did not converge: max|Δy| = %g after %d iters', max_dy, iter);
    end
    if isfield(opts,'verbose') && opts.verbose
        fprintf('MRK mix: iter=%d, conv=%d, max|dy|=%.3e, Σ(f/(φP))=%.4f, Z=%.5f\n', ...
            iter, converged, max_dy, info.sum_f_over_phiP, Z);
    end
end
function [fCO2_bar, phase, details] = fco2_from_totalC( ...
    P_bar, T_K, totalCO2_wt, ...
    wtSiO2, wtTiO2, wtAl2O3, wtFe2O3, wtFeO, wtMnO, wtMgO, wtCaO, wtNa2O, wtK2O, wtP2O5, fO2)

% fco2_from_totalC
% Inverts the original model: given total dissolved CO2 (wt%) in the melt,
% compute the CO2 fugacity (bar) that yields that solubility at (P, T) and composition.
%
% Inputs
%   P_bar  : Pressure in bar (for EOS terms and pressure-dependent coefficients)
%   T_K    : Temperature in Kelvin
%   totalCO2_wt : target total dissolved CO2 in wt%
%   oxide wt%   : anhydrous base OR include water via wtH2O (optional, default = 0)
%   fO2    : oxygen fugacity relative to 1 bar (i.e., fO2 in bar)
%   wtH2O  : optional, default = 0
%
% Outputs
%   fCO2_bar : CO2 fugacity (bar)
%   phase    : 1 = fluid-saturated (melt+fluid); 2 = graphite-saturated; 3 = diamond-saturated
%   details  : struct with intermediate values (NBO, Cmax, etc.)

if nargin < 16 || isempty(wtH2O), wtH2O = 0; end

R = 8.3144; % J/mol/K (used by the original regression)
% ---- Molar masses (g/mol)
M = struct( ...
  'SiO2',60.08, 'TiO2',79.866, 'Al2O3',101.96, 'Fe2O3',159.688, 'FeO',71.844, ...
  'MnO',70.9374, 'MgO',40.3044, 'CaO',56.0774, 'Na2O',61.979, 'K2O',94.2, 'P2O5',141.944);

% ------------------ Normalize composition to 100 wt% (including H2O if provided)
ox_vec = [wtSiO2, wtTiO2, wtAl2O3, wtFe2O3, wtFeO, wtMnO, wtMgO, wtCaO, wtNa2O, wtK2O, wtP2O5];
ox_sum = sum(ox_vec);
if ox_sum <= 0, error('Oxide totals sum to zero.'); end
scale = 100/ox_sum;
[wtSiO2, wtTiO2, wtAl2O3, wtFe2O3, wtFeO, wtMnO, wtMgO, wtCaO, wtNa2O, wtK2O, wtP2O5] = ...
    deal(wtSiO2*scale, wtTiO2*scale, wtAl2O3*scale, wtFe2O3*scale, wtFeO*scale, wtMnO*scale, ...
         wtMgO*scale, wtCaO*scale, wtNa2O*scale, wtK2O*scale, wtP2O5*scale);

% ------------------ Cation "moles" following your original approach
XSi = wtSiO2/M.SiO2;
XTi = wtTiO2/M.TiO2;
XAl = (wtAl2O3/M.Al2O3)*2;
XFe3= (wtFe2O3/M.Fe2O3)*2;
XFe2= wtFeO/M.FeO;
XMn = wtMnO/M.MnO;
XMg = wtMgO/M.MgO;
XCa = wtCaO/M.CaO;
XNa = (wtNa2O/M.Na2O)*2;
XK  = (wtK2O/M.K2O)*2;
XP  = (wtP2O5/M.P2O5)*2;

% Total oxygens (your formula mirrored)
O_tot = XSi*2 + XTi*2 + (wtAl2O3/M.Al2O3)*3 + XNa/2 + XK/2 + (wtP2O5/M.P2O5)*5 + XMg + XMn + XCa + XFe2 + (wtFe2O3/M.Fe2O3)*3;

% Network modifiers and NBO (as in your code)
NM = XFe2 + XMg + XCa + XNa + XK + XMn;
leftover_nms = XAl - NM;
if leftover_nms <= 0
    XAlt = XAl;
else
    XAlt = NM;
end
leftover_nms = -leftover_nms;
Fe3_nms = XFe3 - leftover_nms;
if Fe3_nms <= 0, XFe3t = XFe3; else, XFe3t = leftover_nms; end
tetrahedrals = XAlt + XFe3t + XSi + XTi + XP;
NBO = 2*O_tot - 4*tetrahedrals;

% ------------------ Per-oxygen atomic fractions for fwone (your method)
si1o  = XSi/O_tot; ti1o  = XTi/O_tot; al1o = XAl/O_tot; fe21o = XFe2/O_tot; fe31o = XFe3/O_tot;
mn1o  = XMn/O_tot; mg1o  = XMg/O_tot; ca1o = XCa/O_tot; na1o  = XNa/O_tot; k1o   = XK/O_tot; p1o = XP/O_tot;

fwone = (si1o*28.05 + ti1o*47.867 + al1o*26.9815386 + fe21o*55.845 + fe31o*55.845 + ...
         mn1o*54.938045 + mg1o*24.305 + ca1o*40.078 + na1o*22.98976928 + k1o*39.0983 + ...
         p1o*30.973762 + 15.999);

% ------------------ Oxide mole fractions for gCaO, gNa2O, gK2O terms
n_CaO = wtCaO/M.CaO; n_Na2O = wtNa2O/M.Na2O; n_K2O = wtK2O/M.K2O;
n_all = (wtSiO2/M.SiO2)+(wtTiO2/M.TiO2)+(wtAl2O3/M.Al2O3)+(wtFe2O3/M.Fe2O3)+(wtFeO/M.FeO)+ ...
        (wtMnO/M.MnO)+(wtMgO/M.MgO)+n_CaO+n_Na2O+n_K2O+(wtP2O5/M.P2O5);
X_CaO = n_CaO/max(n_all,eps);
X_Na2O= n_Na2O/max(n_all,eps);
X_K2O = n_K2O/max(n_all,eps);

% ------------------ Graphite/diamond boundary & fCO2 at C-saturation
Pkbar = P_bar/1000;           % kbar
PPa   = P_bar*1e5;            % Pa
trans_P_kbar = 0.025*T_K + 12.592; % Day (2012) crossover (graphite/diamond)

if Pkbar < trans_P_kbar
    % Holloway et al. (1992)
    logK1 = 40.07639 - 2.53932e-2*T_K + 5.27096e-6*T_K.^2 + 0.0267*(P_bar/1000 - 1)./T_K;
else
    % Duncan & Dasgupta (2017)  (note: uses P in Pa)
    logK1 = 2.450e-1 + (2.027e4)/T_K + (-4.730e4)/(T_K.^2) + (1.949e-7*(PPa - 0.001e9))/T_K;
end
fCO2_graph_bar = 10.^logK1 .* fO2;  % bar
phase = (Pkbar < trans_P_kbar)*2 + (Pkbar >= trans_P_kbar)*3; % default if saturated

% ------------------ Regression coefficients (as in your code)
Fco3 = [2.384e-5, -1.6448e5, 1.4732e3, -43.6385, 3.2910, 1.6800e5, 1.7590e5, 2.1085e5];
Fco2 = [1.9244e-5, -9.0212e4, 1.1149e3, -43.0815, -7.0937];

    function [wtCO3, wtCO2, total_wt] = dissolvedCO2_at_f(fco2_bar)
        % Build X array expected by original eqns
        X1 = PPa;     % Pa
        X2 = T_K;     % K
        X3 = log(fco2_bar); % ln fCO2
        X4 = NBO;     % NBO
        X5 = X_CaO; X6 = X_Na2O; X7 = X_K2O;

        % ln mole fractions from original equations
        ln_x_co3 = -(Fco3(1)*X1)/(R*X2) + Fco3(2)/(R*X2) + (X3*Fco3(3))/X2 + Fco3(4)/R + ...
                   (Fco3(5)*X4) + (Fco3(6)*X5 + Fco3(7)*X6 + Fco3(8)*X7)/(R*X2);
        ln_x_co2 = -(Fco2(1)*X1)/(R*X2) + Fco2(2)/(R*X2) + (X3*Fco2(3))/X2 + Fco2(4)/R + ...
                   (Fco2(5)*X4);

        x_co3   = exp(ln_x_co3);
        x_mco2  = exp(ln_x_co2);

        % Convert mole fractions to wt% using "fwone" (as in your code)
        wtCO3 = (44.01*x_co3) ./ (44.01*x_co3 + (1 - (x_co3 + x_mco2)).*fwone) * 100;
        wtCO2 = (44.01*x_mco2)./ (44.01*x_mco2 + (1 - (x_co3 + x_mco2)).*fwone) * 100;

        if NBO > 0.5
            wtCO2 = 0; % your original constraint
        end
        total_wt = wtCO3 + wtCO2;
    end

% ------------------ Check maximum dissolved CO2 at C-saturation
[~, ~, Cmax_wt] = dissolvedCO2_at_f(fCO2_graph_bar);

if totalCO2_wt >= Cmax_wt*(1 - 1e-6)
    % Carbon-saturated: fCO2 is fixed at graph/diamond buffer
    fCO2_bar = fCO2_graph_bar;
    details = struct('NBO',NBO,'fwone',fwone,'Cmax_wt',Cmax_wt,'cap_fCO2_bar',fCO2_graph_bar);
    return
end

% ------------------ Otherwise: solve C_model(fCO2) = C_target in [fmin, f_graph]
fmin = 1e-20; % bar (very reducing / low fCO2)
fun  = @(lnf) (dissolvedCO2_at_f(exp(lnf)) - totalCO2_wt);
% Use log-bracketing for better conditioning
ln_low  = log(max(fmin, realmin));
ln_high = log(fCO2_graph_bar) - 1e-6;

% Ensure the bracket actually brackets the root
C_low = dissolvedCO2_at_f(exp(ln_low));
C_high= dissolvedCO2_at_f(exp(ln_high));
if ~( (C_low - totalCO2_wt) < 0 && (C_high - totalCO2_wt) > 0 )
    % Try expanding the upper bound a bit below the cap
    ln_high = log(fCO2_graph_bar) - 1e-3;
    C_high  = dissolvedCO2_at_f(exp(ln_high));
end

lnfCO2 = fzero(fun, [ln_low, ln_high]);
fCO2_bar = exp(lnfCO2);
phase = 1; % fluid-saturated (no C(s) present)

% Return some diagnostics
[wtCO3, wtCO2, Ccheck] = dissolvedCO2_at_f(fCO2_bar);
details = struct('NBO',NBO,'fwone',fwone,'wtCO3',wtCO3,'wtCO2',wtCO2,'Ccheck',Ccheck, ...
                 'cap_fCO2_bar',fCO2_graph_bar,'Cmax_wt',Cmax_wt);

end
function report = validate_MRK_Holloway(T, P, crit, kIJ, fEOS, aux, varargin)
% validate_MRK_Holloway
% Cross-checks a Holloway-type MRK gas speciation result.
%
% Inputs
%   T, P        : vectors (N×1) of temperature [K] and pressure [bar]
%   crit, kIJ   : same structs/matrices you used in calc_mole_fractions_MRK
%   fEOS        : struct with fields (N×1 each), EOS-consistent fugacities [bar]
%                 fO2, fCl2, fHCl, fH2, fH2O, fSO2, fH2S, fCO2, fCO
%   aux         : struct with optional fields:
%                 - KH2O, KHCl           (N×1): equilibrium constants you used
%                 - logK_SO2, logK_H2S   (N×1)
%                 - logK_COCO2           (N×1)
%                 - CCl_calc, Cl_mes     (N×1): chloride capacity and measured Cl (wt% or matching units)
%                 - opts: struct for MRK (tol, maxIter, verbose)
%
% Optional:
%   'CheckCl2Capacity' (default true): verify fCl2 capacity relation if data provided
%
% Output
%   report: struct with fields
%     .sum_f_over_phiP  (N×1)  : Σ f/(φP)
%     .residuals        (struct of arrays): H2O, COCO2, H2S, SO2 (if available), HCl, Cl2cap
%     .maxAbsResidual   (scalar)
%     .pass             (logical)
%     .table            (table)  : per-row summary
%
% Notes
% - Residuals are computed in log10 space by default. Replace with ln if your K(T) are ln-based.
% - If aux.* K arrays are provided, they are used verbatim (assumed to match your chosen log base).
%
% Example call (after your main routine):
%   fEOS = struct('fO2',fO2,'fCl2',fCl2_eos,'fHCl',fHCl_eos,'fH2',fH2_eos,...
%                 'fH2O',fH2O_eos,'fSO2',fSO2_eos,'fH2S',fH2S_eos,'fCO2',fCO2_eos,'fCO',fCO_eos);
%   aux  = struct('KH2O',KH2O,'KHCl',KHCl,'logK_SO2',logK_SO2,'logK_H2S',logK_H2S,...
%                 'logK_COCO2',logK_CO_CO2,'CCl_calc',CCl_calc,'Cl_mes',Cl_mes);
%   report = validate_MRK_Holloway(T,P,crit,kIJ,fEOS,aux);

% ---------------- Options ----------------
p = inputParser;
addParameter(p,'CheckCl2Capacity',true,@islogical);
parse(p,varargin{:});
doCl2Cap = p.Results.CheckCl2Capacity;

N = numel(T);
if any([numel(P), numel(fEOS.fO2), numel(fEOS.fH2O)] ~= N)
    error('Input arrays must be the same length N.');
end
if nargin < 4 || isempty(kIJ), kIJ = zeros(8); end
if ~isfield(aux,'opts') || isempty(aux.opts)
    aux.opts = struct('tol',1e-10,'maxIter',300,'verbose',false);
end

% ------------- Build default K(T) if not supplied -------------
% IMPORTANT: These assume **log10 K** forms. If your implementation uses ln, pass arrays in aux.
useLog10 = true;
if ~isfield(aux,'KH2O') || isempty(aux.KH2O)
    if useLog10
        aux.KH2O = 10.^(12702./T - 2.8502);  
    else
        aux.KH2O = exp(12702./T - 2.8502);
    end
end
if ~isfield(aux,'KHCl') || isempty(aux.KHCl)
    if useLog10
        aux.KHCl = 10.^(4856.8./T + 0.3405);
    else
        aux.KHCl = exp(4856.8./T + 0.3405);
    end
end
if ~isfield(aux,'logK_SO2') || isempty(aux.logK_SO2)
    aux.logK_SO2 = 18880./T - 3.8018;          % log10 K for: S2 + O2 = 2 SO2  (as used in your code)
end
if ~isfield(aux,'logK_H2S') || isempty(aux.logK_H2S)
    aux.logK_H2S = 27103./T - 4.1973;         % log10 K for: H2S + 3/2 O2 = H2O + SO2  (rearranged in code)
end
if ~isfield(aux,'logK_COCO2') || isempty(aux.logK_COCO2)
    aux.logK_COCO2 = -3395./T + 2.93;         % log10 K for: CO2 + H2 = CO + H2O
end

% ------------- EOS check: Σ f/(φP) ≈ 1 -------------
sum_f_over_phiP = nan(N,1);
Z_out  = nan(N,1);
y_out  = nan(N,8);
phi_out= nan(N,8);

speciesOrder = {'HCl','H2O','Cl2','H2','SO2','H2S','CO2','CO'}; % must match crit and your main code
for i = 1:N
    fvec = [fEOS.fHCl(i); fEOS.fH2O(i); fEOS.fCl2(i); fEOS.fH2(i); ...
            fEOS.fSO2(i); fEOS.fH2S(i); fEOS.fCO2(i); fEOS.fCO(i)];
    [y_i, phi_i, Z_i, info_i] = calc_mole_fractions_MRK(fvec, P(i), T(i), crit, kIJ, aux.opts);
    sum_f_over_phiP(i) = info_i.sum_f_over_phiP;
    Z_out(i)   = Z_i;
    y_out(i,:) = y_i(:).';
    phi_out(i,:)= phi_i(:).';
end

% ------------- Reaction residuals (log10 space by default) -------------
% Helper for safe logs
lt = @(x) log10(max(x, realmin));

res = struct();

% (1) H2 + 1/2 O2 = H2O
%    logK = log fH2O - log fH2 - 0.5 log fO2
lhs = lt(aux.KH2O);   % log10 K
rhs = lt(fEOS.fH2O) - lt(fEOS.fH2) - 0.5*lt(fEOS.fO2);
res.H2O = lhs - rhs;

% Helper
lt = @(x) log10(max(x, realmin));

% Masks
if isfield(aux,'phaseFlag') && numel(aux.phaseFlag)==numel(T)
    mask_sat   = (aux.phaseFlag ~= 1);  % carbon present (graphite/diamond)
    mask_unsat = (aux.phaseFlag == 1);  % fluid-saturated, no C(s)
else
    mask_sat   = true(size(T));   % fallback if no flag provided
    mask_unsat = false(size(T));
end

% Formation checks (apply ONLY where carbon is present)
if isfield(aux,'logKf_CO2')
    res.CO2_form = zeros(size(T));
    m = mask_sat;
    res.CO2_form(m) = lt(fEOS.fCO2(m)) - (aux.logKf_CO2(m) + lt(fEOS.fO2(m)));
end
if isfield(aux,'logKf_CO')
    res.CO_form = zeros(size(T));
    m = mask_sat;
    res.CO_form(m) = lt(fEOS.fCO(m)) - (aux.logKf_CO(m) + 0.5*lt(fEOS.fO2(m)));
end

% Water–gas check for UNSATURATED rows:
% Prefer the thermodynamically-consistent K from your inputs:
if any(mask_unsat)
    res.COCO2 = zeros(size(T));
    m = mask_unsat;

    use_composite = isfield(aux,'logKf_CO') && isfield(aux,'logKf_CO2') && isfield(aux,'KH2O') && ...
                    ~isempty(aux.logKf_CO) && ~isempty(aux.logKf_CO2) && ~isempty(aux.KH2O);

    if use_composite
        log10_Kwg = aux.logKf_CO(m) + lt(aux.KH2O(m)) - aux.logKf_CO2(m);  % log10 K_wg
    else
        % fallback to legacy constant if composite pieces absent
        if ~isfield(aux,'logK_COCO2') || isempty(aux.logK_COCO2)
            aux.logK_COCO2 = -3395./T + 2.93;   % log10
        end
        log10_Kwg = aux.logK_COCO2(m);
    end

    rhs = lt(fEOS.fCO(m)) + lt(fEOS.fH2O(m)) - lt(fEOS.fCO2(m)) - lt(fEOS.fH2(m));
    res.COCO2(m) = log10_Kwg - rhs;
end
% (3) H2S + 3/2 O2 = H2O + SO2   (equivalent to your rearranged form)
lhs = aux.logK_H2S;
rhs = lt(fEOS.fH2O) + lt(fEOS.fSO2) - lt(fEOS.fH2S) - 1.5*lt(fEOS.fO2);
res.H2S = lhs - rhs;

% (4) If you also want to check SO2 from S2+O2=2SO2 and you have fS2, add it to fEOS and uncomment:
if isfield(fEOS,'fS2') && ~isempty(fEOS.fS2)
    lhs = aux.logK_SO2;                  % for S2 + O2 = 2 SO2
    rhs = 2*lt(fEOS.fSO2) - lt(fEOS.fS2) - lt(fEOS.fO2);
    res.SO2 = lhs - rhs;
end

% (5) HCl relation as used in your code:
%     fHCl = KHCl * sqrt(fH2O/KH2O) * (Cl_mes/CCl_calc)
if isfield(aux,'KHCl') && isfield(aux,'KH2O') && isfield(aux,'Cl_mes') && isfield(aux,'CCl_calc')
    lhs = lt(fEOS.fHCl);
    rhs = lt(aux.KHCl) + 0.5*(lt(fEOS.fH2O) - lt(aux.KH2O)) + (lt(aux.Cl_mes) - lt(aux.CCl_calc));
    res.HCl = lhs - rhs;
end

% (6) Cl2 capacity relation (optional): fCl2 = (Cl_mes/CCl)^2 * fO2^(0.5)
if doCl2Cap && isfield(aux,'CCl_calc') && isfield(aux,'Cl_mes')
    lhs = lt(fEOS.fCl2);
    rhs = 2*(lt(aux.Cl_mes) - lt(aux.CCl_calc)) + 0.5*lt(fEOS.fO2);
    res.Cl2cap = lhs - rhs;
end

% ------------- Build summary table -------------
absres = structfun(@(v) abs(v(:)), res, 'UniformOutput', false);
fields = fieldnames(absres);
maxAbsByRow = zeros(N,1);
for i=1:N
    cur = zeros(1,numel(fields));
    for j=1:numel(fields), cur(j) = absres.(fields{j})(i); end
    maxAbsByRow(i) = max(cur);
end

% Pass/fail thresholds (tune as desired)
thr_sum = 1e-3;   % Σ f/(φP) deviation
thr_rxn = 0.02;   % 0.02 in log10 ~ 4.7% multiplicative error

pass_sum = abs(sum_f_over_phiP - 1) < thr_sum;
pass_rxn = maxAbsByRow < thr_rxn;
pass_all = pass_sum & pass_rxn;

% Table
varnames = [{'T_K','P_bar','sum_f_over_phiP','max_rxn_residual_log10','pass'} , fields'];
T_tbl = table();
T_tbl.T_K = T(:);
T_tbl.P_bar = P(:);
T_tbl.sum_f_over_phiP = sum_f_over_phiP;
T_tbl.max_rxn_residual_log10 = maxAbsByRow;
T_tbl.pass = pass_all;
for j=1:numel(fields)
    T_tbl.(sprintf('res_%s',fields{j})) = res.(fields{j})(:);
end

% ------------- Output struct -------------
report = struct();
report.sum_f_over_phiP = sum_f_over_phiP;
report.residuals = res;
report.maxAbsResidual = max(maxAbsByRow);
report.pass = all(pass_all);
report.table = T_tbl;
report.y = y_out;
report.phi = phi_out;
report.Z = Z_out;

% ------------- Print concise summary -------------
fprintf('\n=== MRK/Holloway Validator Summary ===\n');
fprintf('Rows: %d | max |Σ f/(φP) - 1| = %.3e | max reaction residual (log10) = %.3e\n', ...
        N, max(abs(sum_f_over_phiP-1)), report.maxAbsResidual);
fprintf('Overall PASS: %s\n', string(report.pass));
fprintf('Thresholds: Σ f/(φP) within %.1e, reaction residuals within %.3f log10\n\n', thr_sum, thr_rxn);

% Head of table
disp(head(report.table, min(10,N)));

end
