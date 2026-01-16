function Out = HCl_fugacity_with_S_C(~, opts)
% HCl_fugacity_with_S_C
% Strict melt–fluid equilibrium at given (P,T): compute melt-side fugacities f_i,
% then map to gas y_i and φ_i with MRK
%
% H2O_model:  'Moore' | 'Burnham'
%
% CO2_model: 'EguchiDasgupta' (only)
%
% CH4 equilibrium (Ohmoto & Kerrick, 1977):
%   CH4 + 2 O2 = CO2 + 2 H2O
%   log10 K = 41997/T + 0.719*log10(T) - 2.404   (T in K)
%
% Input file: 'MI_Gas_SpeciationGEOROC.csv'
% This version normalizes each row (oxides + H2O + CO2 + S + Cl) to 100 wt%.

    %% ---------------- Options ----------------
    if nargin < 2, opts = struct; end
    if ~isfield(opts,'H2O_model'),     opts.H2O_model     = 'Moore'; end
    if ~isfield(opts,'buffer_choice'), opts.buffer_choice = 'NNO+0.77'; end
    if ~isfield(opts,'doPlots'),       opts.doPlots       = true; end
    if ~isfield(opts,'outCSV'),        opts.outCSV        = 'Fugacities_STRICT.csv'; end

    %% ---------------- Load & normalize ----------------
    data = load('MI_Gas_Speciation2.csv');
    pressure    = data(:,18);    % bar
    temperature = data(:,1);     % °C
    data = data';
    % rows: 1..14 oxides (SiO2..P2O5), 15=H2O, 16=CO2, 17=S, 18=Cl, 19=F (wt%)
    bulkComp = zeros(19, size(data, 2));
    bulkComp([1:3 5:8 11:14, 15, 16, 17, 18, 19], :) = data([2:4 10 5 11 6:9 12 13 14, 15, 16, 17], :);

    % Normalize per column to 100 wt% for rows 1..18 (F ignored)
    comp_rows = 1:18;
    bulkComp(comp_rows, :) = max(bulkComp(comp_rows, :), 0);
    bulkComp(comp_rows, :) = fillmissing(bulkComp(comp_rows, :), 'constant', 0);
    col_sum = sum(bulkComp(comp_rows, :), 1);
    mask = col_sum > 0;
    bulkComp(comp_rows, mask) = bulkComp(comp_rows, mask) .* (100 ./ col_sum(mask));

    %% ---------------- Scalars/vectors ----------------
    T  = temperature + 273.15; % K
    P  = pressure(:);          % bar
    PG = P/1e4;                % GPa

    % Oxide -> moles (note: these are formula-unit moles for oxides,
    % cation-based for some (e.g., Al2O3_mol uses 101.96/2), consistent with
    % original code where CatSum/ntot are defined.)
    SiO2_mol=(bulkComp(1,:)./60.08)';                    TiO2_mol=(bulkComp(2,:)./79.866)';
    Al2O3_mol=(bulkComp(3,:)./(101.96/2))';              FeO_mol=(bulkComp(6,:)./71.844)';
    MgO_mol=(bulkComp(8,:)./40.3044)';                   CaO_mol=(bulkComp(11,:)./56.0774)';
    Na2O_mol=(bulkComp(12,:)./(61.97/2))';               K2O_mol=(bulkComp(13,:)./(94.2/2))';
    Cr2O3_mol=(bulkComp(9,:)./(75.995/2))';              MnO_mol=(bulkComp(10,:)./70.9374)';
    P2O5_mol=(bulkComp(14,:)./(283.889/2))';             H2O_mol=(bulkComp(15,:)./(18.01528/2))';
    H2O_mol_ox=(bulkComp(15,:)./18.01528)';              Cl_mol=(bulkComp(18,:)./35.45)';

    CatSum = SiO2_mol + TiO2_mol + Al2O3_mol + FeO_mol + MgO_mol + CaO_mol + Na2O_mol + K2O_mol + Cr2O3_mol + MnO_mol + P2O5_mol + H2O_mol;
    ntot   = SiO2_mol + TiO2_mol + (bulkComp(3,:)./(101.96))' + FeO_mol + CaO_mol + (bulkComp(12,:)./(61.97))' + (bulkComp(13,:)./(94.2))' + (bulkComp(9,:)./(75.995))' + (bulkComp(14,:)./(283.889))' + (bulkComp(15,:)./(18.01528))' + MnO_mol;

    % Cation-sum mole fractions (used e.g. in Cl capacity)
    XCaO  = CaO_mol ./ CatSum;  XSiO2 = SiO2_mol ./ CatSum;  XFeO = FeO_mol ./ CatSum;
    XK2O  = K2O_mol ./ CatSum;  XNa2O = Na2O_mol ./ CatSum; XMgO = MgO_mol./CatSum;

    % Volatiles / measured
    Cl_mes      = bulkComp(18,:)';                % wt% Cl
    totalCO2_wt = bulkComp(16,:)'; totalCO2_wt(totalCO2_wt<=0 | isnan(totalCO2_wt)) = 1e-40;

    % Oxygen moles for Burnham H2O basis
    O_moles = 2*SiO2_mol + 2*TiO2_mol + 3*(bulkComp(3,:)./(101.96))' + FeO_mol + MgO_mol + CaO_mol + ...
              (bulkComp(12,:)./(61.97))' + (bulkComp(13,:)./(94.2))' + 3*Cr2O3_mol + MnO_mol + 5*(bulkComp(14,:)./(283.889))'+ (bulkComp(15,:)./(18.01528))';
    X_H2O_O = max(0, min(1, H2O_mol_ox .*(8./O_moles)));

   %% Cl capacity
    logC      = 1.15+(4359.*XCaO-3055.*XSiO2+2059.*XFeO-3875.*XK2O+163.*XMgO-514.*(PG))./T;
    CCl_calc  = 10.^logC;
    logKf_H2O = 12850 ./ T - 2.8675;     KH2O = 10.^logKf_H2O;
    logKf_HCl = 4894.6 ./ T + 0.3436;    KHCl = 10.^logKf_HCl;

    %% ---------------- MRK species ------------------
    species = {'HCl','H2O','Cl2','H2','SO2','H2S','CO2','CO','CH4'};
    crit.Tc = [324.7; 647.3; 416.9;  33.2; 430.8; 373.2; 304.13; 132.9; 190.6];   % K
    crit.Pc = [ 83.1; 220.64; 76.1; 12.97;  78.8;  89.6;  73.77; 34.94;  45.99];  % bar
    Nsp  = numel(species);
    kIJ  = zeros(Nsp);
    iHCl = 1; iH2O = 2; iCl2 = 3; iH2 = 4; iSO2 = 5; iH2S = 6; iCO2 = 7; iCO = 8; iCH4 = 9;

    %% ---------------- Utilities ----------------
    function fO2_bar = fO2_from_buffer(Pbar, Tkel, buffer_choice)
        s = regexprep(strtrim(buffer_choice), char(8722), '-');
        tok = regexp(lower(s), '^\s*(fmq|nno|cco)\s*([+-]\s*\d*\.?\d+)?\s*$', 'tokens', 'once');
        if isempty(tok), error('buffer_choice must be FMQ, NNO, or CCO, optionally like NNO+0.3'); end
        base = upper(tok{1}); buffer_offset = 0;
        if numel(tok)>1 && ~isempty(tok{2}), buffer_offset = str2double(regexprep(tok{2},'\s+','')); end
        switch base
            case 'FMQ', base_logfO2 = (-25096.3 ./ Tkel + 8.735 + 0.110 .* (Pbar - 1) ./ Tkel);
            case 'NNO', base_logfO2 = (-24930   ./ Tkel + 9.36  + 0.046 .* (Pbar - 1) ./ Tkel);
            case 'CCO', base_logfO2 = (4.325 - (21803 ./ Tkel) + 0.171 .* ((Pbar - 1) ./ Tkel));
        end
        fO2_bar = 10.^(base_logfO2 + buffer_offset);
        fO2_bar = max(fO2_bar, realmin);
    end

    function phi = phi_pure_MRK_any(Pbar, Tkel, ic)
        % Pure-component MRK fugacity coefficient 
        if ~isfinite(Pbar) || Pbar <= 0, Pbar = 1; end
        if ~isfinite(Tkel) || Tkel <= 0, Tkel = 1000; end
        Rg = 0.08314472; % L·bar/mol·K
        Tc = crit.Tc(ic); Pc = crit.Pc(ic);
        a0 = 0.42748*(Rg^2)*(Tc^2.5)/Pc;   a  = a0 / sqrt(Tkel);
        b  = 0.08664*Rg*Tc/Pc;
        A  = a*Pbar / (Rg^2 * Tkel^2.5);
        B  = b*Pbar / (Rg   * Tkel);
        coeff = [1, -1, A - B - B^2, -A*B];
        if any(~isfinite(coeff)), phi = 1.0; return; end
        Zr = roots(coeff); Zr = Zr(abs(imag(Zr)) < 1e-12);
        if isempty(Zr)
            Z = max(B + 1e-12, 1);
        else
            Z = max(real(Zr)); if Z <= B, Z = B + 1e-12; end
        end
        lnTerm = log(1 + B/Z);
        if B > 1e-12
            lnphi = (Z - 1) - log(Z - B) - (A/B) * lnTerm;
        else
            lnphi = (Z - 1) - log(max(Z - B,1e-12)) - (A/max(Z,1e-12));
        end
        phi = exp(lnphi); if ~isfinite(phi) || phi <= 0, phi = 1.0; end
    end

    %% ---------------- fO2 @ input P ----------------
    fO2 = fO2_from_buffer(P, T, opts.buffer_choice);
    fO2 = max(fO2, realmin);

    %% ---------------- H2O model (Moore / Burnham only) ----------------
    H2O_model = validatestring(opts.H2O_model, {'Moore','Burnham'});
    fH2O_target = zeros(size(T));

    switch H2O_model
        case 'Moore'
            CatSum_an  = SiO2_mol+TiO2_mol+(bulkComp(3,:)./(101.96))'+FeO_mol+ ...
                         MgO_mol+CaO_mol+(bulkComp(12,:)./(61.97))'+(bulkComp(13,:)./(94.2))'+MnO_mol;
            XAl2O3_an = (bulkComp(3,:)./(101.96))'./CatSum_an;
            XFeO_an   = FeO_mol./CatSum_an;
            XH2O_an   = (bulkComp(15,:)./(18.01528))'./CatSum_an;
            XNa2O_an  = (bulkComp(12,:)./(61.97))'./CatSum_an;
            a_M=2565; bAl=-1.997; bFe=-0.9275; bNa=2.736; c=1.171; d=-14.21;
            B_M    = (bAl*XAl2O3_an + bFe*XFeO_an + bNa*XNa2O_an) .* ((P)./T);
            lnfH2O = (2*log(XH2O_an) - a_M./T - B_M - d)./c;
            fH2O_target = exp(lnfH2O);   % bar

        case 'Burnham'
            Po = 1.01325; % bar
            lnkw = 5.00 + log(P) .* (4.481e-8 .* T.^2 - 1.51e-4 .* T - 1.137) + ...
                   (log(P)).^2 .* (1.831e-8 .* T.^2 - 4.882e-5 .* T + 0.04656) + ...
                   7.80e-3 .* (log(P)).^3 - 5.012e-4 .* (log(P)).^4 + ...
                   T .* (4.754e-3 - 1.621e-6 .* T);
            kw = exp(lnkw);
            aH2O = zeros(size(T));
            for idx = 1:length(T)
                if X_H2O_O(idx) >= 0.5
                    aH2O(idx) = 0.25 * kw(idx) * exp((6.52 - 2667 / T(idx)) * (X_H2O_O(idx) - 0.5));
                else
                    aH2O(idx) = (X_H2O_O(idx))^2 * kw(idx);
                end
            end
            % pure-component φ for H2O at (P,T)
            phiH2O_pure = arrayfun(@(Pi,Ti) phi_pure_MRK_any(Pi,Ti,iH2O), P, T);
            fH2O_target = aH2O .* (phiH2O_pure .* P)./Po; % bar
    end

    fH2O_target = max(fH2O_target, realmin);
    Cl_mes   = max(Cl_mes,   realmin);
    CCl_calc = max(CCl_calc, realmin);

    %% ---------------- Kress & Carmichael Fe3+/Fe2+ (anhydrous oxide basis) ----------------
    % Build anhydrous oxides (SiO2..P2O5) from bulkComp, renormalized to 100 wt%
    ox_anh     = bulkComp(1:14,:);                        % major oxides only
    ox_anh_sum = sum(ox_anh,1);
    ox_anh_sum(ox_anh_sum <= 0) = 1;                      % avoid division by zero
    scale_ox_anh = 100 ./ ox_anh_sum;
    ox_anh     = ox_anh .* scale_ox_anh;                  % each column sums to 100 wt% anhydrous

    % Extract anhydrous wt% for K&C-relevant oxides
    wtSiO2_anh   = ox_anh(1,:)';
    wtTiO2_anh   = ox_anh(2,:)';
    wtAl2O3_anh  = ox_anh(3,:)';
    wtFeOtot_anh = ox_anh(6,:)';   % all Fe as FeO* here
    wtMgO_anh    = ox_anh(8,:)';
    wtMnO_anh    = ox_anh(10,:)';
    wtCaO_anh    = ox_anh(11,:)';
    wtNa2O_anh   = ox_anh(12,:)';
    wtK2O_anh    = ox_anh(13,:)';
    wtP2O5_anh   = ox_anh(14,:)';

    % Convert to moles of oxide formula units
    nSiO2  = wtSiO2_anh   ./ 60.08;
    nTiO2  = wtTiO2_anh   ./ 79.866;
    nAl2O3 = wtAl2O3_anh  ./ 101.96;
    nFeOtot= wtFeOtot_anh ./ 71.844;
    nMgO   = wtMgO_anh    ./ 40.3044;
    nMnO   = wtMnO_anh    ./ 70.9374;
    nCaO   = wtCaO_anh    ./ 56.0774;
    nNa2O  = wtNa2O_anh   ./ 61.979;
    nK2O   = wtK2O_anh    ./ 94.2;
    nP2O5  = wtP2O5_anh   ./ 141.944;

    n_tot_ox = nSiO2 + nTiO2 + nAl2O3 + nFeOtot + nMgO + nMnO + nCaO + nNa2O + nK2O + nP2O5;
    n_tot_ox = max(n_tot_ox, realmin);

    % Oxide mole fractions on anhydrous basis (for K&C compositional terms)
    XAl2O3_ox = nAl2O3  ./ n_tot_ox;
    XFeO_ox   = nFeOtot ./ n_tot_ox;    % total Fe as FeO* here
    XCaO_ox   = nCaO    ./ n_tot_ox;
    XNa2O_ox  = nNa2O   ./ n_tot_ox;
    XK2O_ox   = nK2O    ./ n_tot_ox;

    % ln( XFe2O3 / XFeO ) from Kress & Carmichael (1991) with oxide mole fractions
    lnXFe3XFe2 = 0.196*log(fO2) + 11492./T - 6.675 + ...
                 (-2.243*XAl2O3_ox + 5.854*XNa2O_ox - 1.828*XFeO_ox + ...
                  6.215*XK2O_ox + 3.201*XCaO_ox);
    R_Fe3_Fe2_ox = exp(lnXFe3XFe2);          % ≈ XFe2O3 / XFeO (oxide-mole basis)
    R_Fe3_Fe2_ox = max(R_Fe3_Fe2_ox, realmin);

    % Split total Fe (as FeO*) between FeO and Fe2O3 on oxide basis:
    % Let nFeOtot be moles of FeO*; then:
    %   n_Fe2O3 = R * nFeOtot / (2R + 1)
    %   n_FeO   = nFeOtot - 2*n_Fe2O3
    nFeO_tot = nFeOtot;
    n_Fe2O3  = (R_Fe3_Fe2_ox .* nFeO_tot) ./ (2.*R_Fe3_Fe2_ox + 1);
    n_FeO    = nFeO_tot - 2.*n_Fe2O3;

    % Convert FeO and Fe2O3 moles back to wt% in the full bulk (including volatiles)
    M_FeO   = 71.844;
    M_Fe2O3 = 159.688;

    mass_FeO   = n_FeO   .* M_FeO;
    mass_Fe2O3 = n_Fe2O3 .* M_Fe2O3;

    % Mass of all other components taken from current bulkComp
    % (1..18 are normalized to 100; remove original FeO row 6, then add new FeO & Fe2O3)
    mass_other = sum(bulkComp([1:5 7:14 15:18],:),1)';  % everything except original FeO
    total_mass = mass_other + mass_FeO + mass_Fe2O3;
    total_mass = max(total_mass, realmin);

    wtFeO_KC   = mass_FeO   ./ total_mass * 100;
    wtFe2O3_KC = mass_Fe2O3 ./ total_mass * 100;

    FeO_mol_KC   = wtFeO_KC   ./ M_FeO;
    Fe2O3_mol_KC = wtFe2O3_KC ./ M_Fe2O3;

    %% ---------------- Sulfur capacities → fS2 → fSO2 (needs Fe redox) ----------------
    % 1-O basis for capacity regressions
    n_SiO2  = (bulkComp(1,:)/(60.08/2))';   n_TiO2  = (bulkComp(2,:)/(79.866/2))';
    n_Al2O3 = (bulkComp(3,:)/(101.96/3))'; n_FeO2p = FeO_mol_KC; n_Fe2O3 = Fe2O3_mol_KC;
    n_MnO   = MnO_mol; n_MgO = MgO_mol;  n_CaO = CaO_mol;
    n_Na2O  = (bulkComp(12,:)/(61.97))';  n_K2O = (bulkComp(13,:)/(94.2))';
    O_tot1  = n_SiO2+n_TiO2+n_Al2O3+n_FeO2p+n_Fe2O3+n_MnO+n_MgO+n_CaO+n_Na2O+n_K2O;

    X_Si05O    = n_SiO2  ./ O_tot1; 
    X_Ti05O    = n_TiO2  ./ O_tot1;
    X_Al067O   = n_Al2O3 ./ O_tot1;
    X_FeO_2p   = n_FeO2p ./ O_tot1;
    X_Fe067O3p = n_Fe2O3 ./ O_tot1;
    X_MnO_1O   = n_MnO   ./ O_tot1;
    X_MgO_1O   = n_MgO   ./ O_tot1;
    X_CaO_1O   = n_CaO   ./ O_tot1;
    X_Na2O_1O  = n_Na2O  ./ O_tot1;
    X_K2O_1O   = n_K2O   ./ O_tot1;

    logCS2_base = 0.65 + (-3368.*X_Si05O -1233.*X_Al067O +1295.*X_CaO_1O +44885.*X_K2O_1O + ...
                          10914.*(X_FeO_2p .* X_Si05O) -871864.*(X_FeO_2p .* X_K2O_1O) - ...
                          225569.*(X_FeO_2p .* X_Na2O_1O) + 54392.*(X_FeO_2p .* X_Al067O) - 7585)./T + 3.9*erf(X_FeO_2p);
    logCS6_base = 196.0 + (-7633.*X_Si05O + -10670.*X_Ti05O + -6901.*X_Al067O + -4625.*X_FeO_2p + ...
                           19114.*X_MnO_1O + 11740.*X_CaO_1O + 33267.*X_Na2O_1O - 4095)./T - 57.3.*log10(T);

    % ΔV terms: add (P-1)*ΔV/(2.303 R T)
    logCS2 = logCS2_base + ((P - 1) .* 0.62) ./ (T.*8.314.*2.303);
    logCS6 = logCS6_base + ((P - 1) .* 2.92) ./ (T.*8.314.*2.303);

    S_mes = bulkComp(17,:)'; S_mes(S_mes <= 0 | isnan(S_mes)) = 1e-12;
    % S equilibrium constants (Ohmoto-style parameterization; T in K; base-10 logs)
    logK_SO2 = 18880./T - 3.8018;  % SO2 = 0.5 S2 + O2
    logK_H2S = 27103./T - 4.1973;  % H2S + 1.5 O2 = H2O + 0.5 S2
    fS2  = zeros(size(P)); fSO2 = zeros(size(P)); fH2S = zeros(size(P));
    lt = @(x) log10(max(x, realmin)); ten = @(x) 10.^x;
    for i = 1:numel(P)
        logfO2_i = lt(fO2(i));
        massbal = @(logfS2) ( ten(logCS2(i)) .* ten( 0.5*logfS2 - 0.5*logfO2_i ) ) + ...
                            ( ten(logCS6(i)) .* ten( 0.5*logfS2 + 1.5*logfO2_i ) ) - S_mes(i);
        lo = -40; hi = 10; tries=0;
        while massbal(lo)*massbal(hi) > 0 && tries < 8
            lo = lo - 5; hi = hi + 5; tries=tries+1;
        end
        if massbal(lo)*massbal(hi) > 0
            logfS2_i = -30; % fallback
        else
            logfS2_i = fzero(massbal, [lo, hi]);
        end
        fS2(i)  = ten(logfS2_i);
        fSO2(i) = 10^(logK_SO2(i)) * sqrt(max(fS2(i),realmin)) * max(fO2(i),realmin);
        % fH2S computed after fH2O (below)
    end

    %% ---------------- CO2 (Eguchi & Dasgupta 2018) ----------------
    ox_total = sum(bulkComp([1 2 3 10 8 11 12 13 14 15],:),1)' + (wtFe2O3_KC+wtFeO_KC);
    wSi  = (bulkComp(1,:)'./ox_total)*100; wTi=(bulkComp(2,:)'./ox_total)*100; wAl=(bulkComp(3,:)'./ox_total)*100;
    wMg  = (bulkComp(8,:)'./ox_total)*100; wMn=(bulkComp(10,:)'./ox_total)*100; wCa=(bulkComp(11,:)'./ox_total)*100;
    wNa  = (bulkComp(12,:)'./ox_total)*100; wK =(bulkComp(13,:)'./ox_total)*100; wP =(bulkComp(14,:)'./ox_total)*100;
    wFeO = wtFeO_KC; wFe2O3 = wtFe2O3_KC;

    fCO2 = zeros(size(P)); phaseC = ones(size(P));
    for i = 1:numel(P)
        [fCO2_i, phase_i, ~] = fco2_from_totalC( ...
            P(i), T(i), totalCO2_wt(i), ...
            wSi(i), wTi(i), wAl(i), wFe2O3(i), wFeO(i), wMn(i), wMg(i), ...
            wCa(i), wNa(i), wK(i), wP(i), fO2(i), []); % last arg (wtH2O) unused internally
        if phase_i ~= 1
            % carbon-saturated: cap via CO2 formation constant
            fCO2_i = min( 10.^(20599./T(i) + 0.0531) * fO2(i), 0.95*P(i) );
        end
        fCO2(i)  = max(fCO2_i, realmin);
        phaseC(i)= phase_i;
    end

    %% ---------------- Other volatiles at (P,T) ----------------
    % Chlorine system (Rusiecka & Wood 2026)
    fHCl = KHCl .* sqrt( fH2O_target ./ KH2O ) .* ( Cl_mes ./ CCl_calc );
    fCl2 = (Cl_mes ./ CCl_calc).^2 .* sqrt(fO2);

    % H2 from H2O–H2 equilibrium
    fH2 = fH2O_target ./ (KH2O .* sqrt(fO2));

    % H2S now that fH2O known
    logK_SO2 = 18880./T - 3.8018;  % SO2 = 0.5 S2 + O2
    logK_H2S = 27103./T - 4.1973;  % H2S + 1.5 O2 = H2O + 0.5 S2
    fH2S = ( fH2O_target .* max(fSO2,realmin) ) ./ ( (max(fO2,realmin)).^1.5 .* 10.^(logK_H2S) );

    % CO via CO2/CO equilibrium
    logKf_CO2 = (20599./T + 0.0531);
    logKf_CO  = (5839.7./T + 4.5515);
    Kratio    = 10.^(logKf_CO2 - logKf_CO);
    fCO       = fCO2 ./ (Kratio .* sqrt(fO2));
    fCO       = max(fCO, realmin);

    % CH4 via Ohmoto & Kerrick (1977): CH4 + 2 O2 = CO2 + 2 H2O
    log10K_OK77 = 41997./T + 0.719*log10(T) - 2.404;
    K_OK        = 10.^log10K_OK77;
    fCH4        = ( fCO2 .* (fH2O_target.^2) ) ./ ( K_OK .* (fO2.^2) );
    fCH4        = max(fCH4, realmin);

    %% ---------------- MRK inversion (strict mapper) ----------------
    optsMRK = struct('tol',1e-8,'maxIter',10000,'verbose',false);

    % fugacity vector (columns by species order)
    F = [fHCl(:), fH2O_target(:), fCl2(:), fH2(:), fSO2(:), fH2S(:), fCO2(:), fCO(:), fCH4(:)];
    mole_fr   = nan(numel(P), Nsp);
    f_eos_out = nan(numel(P), Nsp);
    phi_out   = nan(numel(P), Nsp);
    sum_cl    = nan(numel(P), 1);
    maxdiff   = nan(numel(P), 1);

    for i = 1:numel(P)
        fvec = max(F(i,:).', realmin);  % column vector
        [y_i, phi_i, ~, info_i] = calc_mole_fractions_MRK(fvec, P(i), T(i), crit, kIJ, optsMRK);
        mole_fr(i,:)  = y_i.';
        phi_out(i,:)  = phi_i.';
        f_check       = (y_i(:).*phi_i(:)*P(i)).';
        f_eos_out(i,:)= f_check;
        sum_cl(i)     = info_i.sum_f_over_phiP;
        maxdiff(i)    = max(abs(log10( max(f_check,realmin) ./ max(fvec.',realmin) )));
    end

    %% ---------------- Output ----------------
    Out = table(temperature(:), P(:), fO2(:), ...
        fH2S(:), fSO2(:), fCl2(:), fHCl(:), fH2O_target(:), ...
        fCO2(:), fCO(:), fCH4(:), Cl_mes(:), ...
        mole_fr(:,1), mole_fr(:,2), mole_fr(:,3), mole_fr(:,4), mole_fr(:,5), mole_fr(:,6), ...
        mole_fr(:,7), mole_fr(:,8), mole_fr(:,9), ...
        totalCO2_wt(:), logCS6(:), logC(:), sum_cl(:), maxdiff(:), ...
        'VariableNames', {'T_C','P_bar','fO2_bar', ...
                          'fH2S_bar','fSO2_bar','fCl2_bar','fHCl_bar','fH2O_bar', ...
                          'fCO2_bar','fCO_bar','fCH4_bar','Cl_wt_pct', ...
                          'X_HCl','X_H2O','X_Cl2','X_H2','X_SO2','X_H2S','X_CO2','X_CO','X_CH4', ...
                          'TotalCO2_wt','logCS6', 'logCCl', 'sum_f_over_phiP','maxLog10Diff'});

    writetable(Out, opts.outCSV);
    fprintf('\nWrote %s (STRICT mapping; H2O=%s; CO2=EguchiDasgupta; CH4=Ohmoto&Kerrick77).\n', ...
        opts.outCSV, opts.H2O_model);

    %% ---------------- Plots (optional) ----------------
    if opts.doPlots
        safeLog10 = @(x) log10(max(x, realmin));
        figure; hold on;
        plot(P, safeLog10(fO2),     'ko', 'DisplayName', 'log fO_2');
        plot(P, safeLog10(fCl2),    'bo', 'DisplayName', 'log fCl_2');
        plot(P, safeLog10(fHCl),    'ro', 'DisplayName', 'log fHCl');
        plot(P, safeLog10(fH2O_target),'mo','DisplayName','log fH_2O');
        plot(P, safeLog10(fH2),     'go', 'DisplayName', 'log fH_2');
        plot(P, safeLog10(fSO2),    'c^', 'DisplayName', 'log fSO_2');
        plot(P, safeLog10(fH2S),    'ms', 'DisplayName', 'log fH_2S');
        plot(P, safeLog10(fCO2),    'bs', 'DisplayName', 'log fCO_2');
        plot(P, safeLog10(fCO),     'gs', 'DisplayName', 'log fCO');
        plot(P, safeLog10(fCH4),    'k^', 'DisplayName', 'log fCH_4');
        xlabel('Pressure (bar)'); ylabel('log Fugacity (bar)');
        legend('Location','best'); grid on;

        figure; hold on;
        plot(P, mole_fr(:,1), 'ro', 'DisplayName', 'X HCl');
        plot(P, mole_fr(:,2), 'mo', 'DisplayName', 'X H_2O');
        plot(P, mole_fr(:,3), 'bo', 'DisplayName', 'X Cl_2');
        plot(P, mole_fr(:,4), 'go', 'DisplayName', 'X H_2');
        plot(P, mole_fr(:,5), 'c^', 'DisplayName', 'X SO_2');
        plot(P, mole_fr(:,6), 'ms', 'DisplayName', 'X H_2S');
        plot(P, mole_fr(:,7), 'bs', 'DisplayName', 'X CO_2');
        plot(P, mole_fr(:,8), 'gs', 'DisplayName', 'X CO');
        plot(P, mole_fr(:,9), 'k^', 'DisplayName', 'X CH_4');
        xlabel('Pressure (bar)'); ylabel('Mole Fraction');
        legend('Location','best'); grid on;
    end
end

% ======================= MRK mixture inversion =======================
function [y, phi, Z, info] = calc_mole_fractions_MRK(fug, P, T, crit, kIJ, opts)
    N = numel(fug);
    if nargin < 5 || isempty(kIJ), kIJ = zeros(N); end
    if nargin < 6, opts = struct; end
    if ~isfield(opts,'tol'),     opts.tol     = 1e-8; end
    if ~isfield(opts,'maxIter'), opts.maxIter = 10000; end
    if ~isfield(opts,'verbose'), opts.verbose = false; end

    R = 0.08314472; % L·bar/(mol·K)
    Tc = crit.Tc(:); Pc = crit.Pc(:);
    if numel(Tc) ~= N || numel(Pc) ~= N, error('crit.Tc and crit.Pc must be length N'); end
    if any(fug < 0), error('Fugacities must be non-negative.'); end
    if any(diag(kIJ) ~= 0), error('Diagonal of kIJ must be zero.'); end

    a0_i = 0.42748 .* (R.^2) .* (Tc.^2.5) ./ Pc;
    a_i  = a0_i ./ sqrt(T);
    b_i  = 0.08664 .* R .* Tc ./ Pc;

    if all(fug == 0), error('All fugacities are zero.'); end
    y = fug(:) / sum(fug);

    iter = 0; converged = false; max_dy = NaN;
    while iter < opts.maxIter
        iter = iter + 1;
        y_old = y;

        sqrt_a = sqrt(a_i);
        a_mix = 0.0;
        a_i_mix = zeros(N,1);
        for ii = 1:N
            a_i_mix(ii) = sum( y(:)'.*(1 - kIJ(ii,:)) .* (sqrt_a(ii).*sqrt_a(:))' );
            a_mix = a_mix + y(ii) * a_i_mix(ii);
        end
        b_mix = sum(y .* b_i);

        A = a_mix * P / (R^2 * T^2.5);
        B = b_mix * P / (R * T);
        coeffs = [1, -1, A - B - B^2, -A*B];
        if any(~isfinite(coeffs))
            Z = max(B + 1e-12, 1);
            phi = ones(N,1);
            break
        end
        Zr = roots(coeffs); Zr = Zr(abs(imag(Zr)) < 1e-9);
        if isempty(Zr)
            Z = max(B + 1e-12, 1);
        else
            Z = max(real(Zr)); if Z <= B, Z = B + 1e-12; end
        end

        lnTerm = log(1 + B/Z);
        phi = zeros(N,1);
        for ii = 1:N
            bi_b  = b_i(ii)/max(b_mix,realmin);
            sum_a = a_i_mix(ii);
            lnphi = bi_b*(Z - 1) - log(Z - B) - (A/max(B,realmin)) * ( 2*sum_a/max(a_mix,realmin) - bi_b ) * lnTerm;
            phi(ii) = exp(lnphi);
        end

        denom = phi * P;
        if any(~isfinite(denom) | denom <= 0), error('Non-physical φ_i * P encountered.'); end
        y = fug(:) ./ denom; y(y < 0) = 0; y = y / sum(y);

        max_dy = max(abs(y - y_old));
        if max_dy < opts.tol, converged = true; break; end
    end

    info.iterations = iter;
    info.converged  = converged;
    info.max_dy     = max_dy;
    info.sum_f_over_phiP = sum(fug(:)./(phi(:)*P));

    if ~converged
        warning('MRK mixture did not converge: max|Δy| = %g after %d iters', max_dy, iter);
    end
end

% ============ CO2 solubility inversion (Eguchi & Dasgupta 2018) =============
function [fCO2_bar, phase, details] = fco2_from_totalC( ...
    P_bar, T_K, totalCO2_wt, ...
    wtSiO2, wtTiO2, wtAl2O3, wtFe2O3, wtFeO, wtMnO, wtMgO, wtCaO, wtNa2O, wtK2O, wtP2O5, fO2, wtH2O)

    if nargin < 16 || isempty(wtH2O), wtH2O = 0; end 
    R = 8.3144; % J/mol/K
    M = struct('SiO2',60.08, 'TiO2',79.866, 'Al2O3',101.96, 'Fe2O3',159.688, 'FeO',71.844, ...
               'MnO',70.9374, 'MgO',40.3044, 'CaO',56.0774, 'Na2O',61.979, 'K2O',94.2, 'P2O5',141.944);

    % Normalize oxide vector to 100 wt%
    ox_vec = [wtSiO2, wtTiO2, wtAl2O3, wtFe2O3, wtFeO, wtMnO, wtMgO, wtCaO, wtNa2O, wtK2O, wtP2O5];
    scale = 100/sum(ox_vec);
    [wtSiO2, wtTiO2, wtAl2O3, wtFe2O3, wtFeO, wtMnO, wtMgO, wtCaO, wtNa2O, wtK2O, wtP2O5] = ...
        deal(wtSiO2*scale, wtTiO2*scale, wtAl2O3*scale, wtFe2O3*scale, wtFeO*scale, wtMnO*scale, ...
             wtMgO*scale, wtCaO*scale, wtNa2O*scale, wtK2O*scale, wtP2O5*scale);

    XSi = wtSiO2/M.SiO2; XTi = wtTiO2/M.TiO2; XAl = (wtAl2O3/M.Al2O3)*2;
    XFe3= (wtFe2O3/M.Fe2O3)*2; XFe2= wtFeO/M.FeO; XMn = wtMnO/M.MnO; XMg = wtMgO/M.MgO;
    XCa = wtCaO/M.CaO; XNa = (wtNa2O/M.Na2O)*2; XK = (wtK2O/M.K2O)*2; XP = (wtP2O5/M.P2O5)*2;

    O_tot = XSi*2 + XTi*2 + (wtAl2O3/M.Al2O3)*3 + XNa/2 + XK/2 + (wtP2O5/M.P2O5)*5 + ...
            XMg + XMn + XCa + XFe2 + (wtFe2O3/M.Fe2O3)*3;

    NM = XFe2 + XMg + XCa + XNa + XK + XMn;
    leftover_nms = XAl - NM;  XAlt = (leftover_nms <= 0)*XAl + (leftover_nms > 0)*NM;
    leftover_nms = -leftover_nms;
    Fe3_nms = XFe3 - leftover_nms;  XFe3t = (Fe3_nms <= 0)*XFe3 + (Fe3_nms > 0)*leftover_nms;
    tetrahedrals = XAlt + XFe3t + XSi + XTi + XP;
    NBO = 2*O_tot - 4*tetrahedrals;

    si1o=XSi/O_tot; ti1o=XTi/O_tot; al1o=XAl/O_tot; fe21o=XFe2/O_tot; fe31o=XFe3/O_tot;
    mn1o=XMn/O_tot; mg1o=XMg/O_tot; ca1o=XCa/O_tot; na1o=XNa/O_tot; k1o=XK/O_tot; p1o=XP/O_tot;
    fwone = (si1o*28.05 + ti1o*47.867 + al1o*26.9815386 + fe21o*55.845 + fe31o*55.845 + ...
             mn1o*54.938045 + mg1o*24.305 + ca1o*40.078 + na1o*22.98976928 + k1o*39.0983 + ...
             p1o*30.973762 + 15.999);

    n_CaO = wtCaO/M.CaO; n_Na2O = wtNa2O/M.Na2O; n_K2O = wtK2O/M.K2O;
    n_all = (wtSiO2/M.SiO2)+(wtTiO2/M.TiO2)+(wtAl2O3/M.Al2O3)+(wtFe2O3/M.Fe2O3)+(wtFeO/M.FeO)+ ...
            (wtMnO/M.MnO)+(wtMgO/M.MgO)+n_CaO+n_Na2O+n_K2O+(wtP2O5/M.P2O5);
    X_CaO = n_CaO/max(n_all,eps); X_Na2O = n_Na2O/max(n_all,eps); X_K2O = n_K2O/max(n_all,eps);

    Pkbar = P_bar/1000; PPa = P_bar*1e5;
    trans_P_kbar = 0.025*T_K + 12.592;

    if Pkbar < trans_P_kbar
        logK1 = 40.07639 - 2.53932e-2*T_K + 5.27096e-6*T_K.^2 + 0.0267*(P_bar/1000 - 1)./T_K;
    else
        logK1 = 2.450e-1 + (2.027e4)/T_K + (-4.730e4)/(T_K.^2) + (1.949e-7*(PPa - 0.001e9))/T_K;
    end
    fCO2_graph_bar = 10.^logK1 .* fO2;
    phase = 2 + double(Pkbar >= trans_P_kbar);   % 2 graphite, 3 diamond

    Fco3 = [2.384e-5, -1.6448e5, 1.4732e3, -43.6385, 3.2910, 1.6800e5, 1.7590e5, 2.1085e5];
    Fco2 = [1.9244e-5, -9.0212e4, 1.1149e3, -43.0815, -7.0937];

    function [wtCO3, wtCO2, total_wt] = dissolvedCO2_at_f(fco2_bar)
        X1 = PPa; X2 = T_K; X3 = log(fco2_bar); X4 = NBO; X5 = X_CaO; X6 = X_Na2O; X7 = X_K2O;
        ln_x_co3 = -(Fco3(1)*X1)/(R*X2) + Fco3(2)/(R*X2) + (X3*Fco3(3))/X2 + Fco3(4)/R + ...
                   (Fco3(5)*X4) + (Fco3(6)*X5 + Fco3(7)*X6 + Fco3(8)*X7)/(R*X2);
        ln_x_co2 = -(Fco2(1)*X1)/(R*X2) + Fco2(2)/(R*X2) + (X3*Fco2(3))/X2 + Fco2(4)/R + ...
                   (Fco2(5)*X4);
        x_co3  = exp(ln_x_co3); x_mco2 = exp(ln_x_co2);
        wtCO3  = (44.01*x_co3) ./ (44.01*x_co3 + (1 - (x_co3 + x_mco2)).*fwone) * 100;
        wtCO2  = (44.01*x_mco2)./ (44.01*x_mco2 + (1 - (x_co3 + x_mco2)).*fwone) * 100;
        if NBO > 0.5, wtCO2 = 0; end
        total_wt = wtCO3 + wtCO2;
    end

    [~, ~, Cmax_wt] = dissolvedCO2_at_f(fCO2_graph_bar);
    if totalCO2_wt >= Cmax_wt*(1 - 1e-6)
        fCO2_bar = fCO2_graph_bar;
        details = struct('NBO',NBO,'fwone',fwone,'Cmax_wt',Cmax_wt,'cap_fCO2_bar',fCO2_graph_bar);
        return
    end

    % Robust bracket & solve
    fmin = 1e-100;
    fun  = @(lnf) (dissolvedCO2_at_f(exp(lnf)) - totalCO2_wt);

    ln_low  = log(max(fmin, realmin));
    ln_high = log(max(fCO2_graph_bar, 1e-50)) - 1e-9;

    if ~isfinite(ln_high) || ln_high <= ln_low
        fCO2_bar = exp(ln_low); phase = 1;
        [wtCO3, wtCO2, Ccheck] = dissolvedCO2_at_f(fCO2_bar); 
        details = struct('NBO',NBO,'fwone',fwone,'wtCO3',wtCO2,'wtCO2',wtCO2,'Ccheck',Ccheck, ...
                         'cap_fCO2_bar',fCO2_graph_bar,'Cmax_wt',Cmax_wt);
        return
    end

    L = linspace(ln_low, ln_high, 80);
    V = arrayfun(fun, L); V(~isfinite(V)) = sign(V(~isfinite(V))) .* 1e3;

    if all(V < 0)
        fCO2_bar = fCO2_graph_bar; phase = 2 + double((P_bar/1000) >= trans_P_kbar);
        details  = struct('NBO',NBO,'fwone',fwone,'Cmax_wt',Cmax_wt,'cap_fCO2_bar',fCO2_graph_bar);
        return
    end
    if all(V > 0)
        fCO2_bar = exp(ln_low); phase = 1;
        [wtCO3, wtCO2, Ccheck] = dissolvedCO2_at_f(fCO2_bar);
        details  = struct('NBO',NBO,'fwone',fwone,'wtCO3',wtCO3,'wtCO2',wtCO2,'Ccheck',Ccheck, ...
                          'cap_fCO2_bar',fCO2_graph_bar,'Cmax_wt',Cmax_wt);
        return
    end

    idx = find(V(1:end-1).*V(2:end) <= 0, 1, 'first');
    if isempty(idx)
        obj = @(lnf) (fun(lnf)).^2;
        lnfCO2 = fminbnd(obj, ln_low, ln_high);
    else
        lnfCO2 = fzero(fun, [L(idx), L(idx+1)]);
    end

    fCO2_bar = exp(lnfCO2);
    phase    = 1;

    [wtCO3, wtCO2, Ccheck] = dissolvedCO2_at_f(fCO2_bar); 
    details = struct('NBO',NBO,'fwone',fwone,'wtCO3',wtCO3,'wtCO2',wtCO2,'Ccheck',Ccheck, ...
                     'cap_fCO2_bar',fCO2_graph_bar,'Cmax_wt',Cmax_wt);
end
