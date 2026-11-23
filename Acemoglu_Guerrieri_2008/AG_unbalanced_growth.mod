//----------------------------------------------------------------------------
/* Discrete-time version of Acemoglu & Guerrieri (2008),
 * "Capital Deepening and Non-Balanced Economic Growth", Section 3 model,
 * calibrated as in Section 4.
 * Published version: https://www.journals.uchicago.edu/doi/10.1086/589523
 *
 * Note on solution method:
 *  AG discretize the continuous-time ODEs with the Euler method and
 *  approximate the 3000-year path by a global polynomial in t, fitted by
 *  nonlinear least squares (see footnote 19). Here the discrete-time system is 
 *  solved exactly (up to Newton tolerance) for all T periods by Dynare's 
 *  stacked-Newton perfect-foresight solver; the only approximation is the
 *  terminal condition (affects roughly the last 100 periods). lambda, kappa and
 *  sigma_K match Figure 3 closely, but r declines monotonically from r(0)=0.136
 *  to r*~0.104 instead of the non-monotonic path in the paper. Since
 *  r = rho + theta*g_C and g_C > g* along the transition, r should stay above
 *  r*; the paper's r path is likely an artefact of the polynomial approximation.
 *
 * Model (notation: a_i = capital share in sector i = 1 - alpha_i of the paper)
 *   Y   = [gam*Y1^((eps-1)/eps) + (1-gam)*Y2^((eps-1)/eps)]^(eps/(eps-1))     (15)
 *   Yi  = A_i * Ki^a_i * Li^(1-a_i),  A_i = M_i^(1/(nu-1)) grows at m_i_hat  (27)
 *   L grows at n, K1+K2 = K, L1+L2 = L
 *   max sum_t beta^t (C_t^(1-theta)-1)/(1-theta)              (13)
 *   K_t = (1-delta)K_{t-1} + Y_t - C_t,  delta = 0            (16)
 *
 * Normalisation (eps < 1, sector 1 asymptotically dominant: m1hat/(1-a1) <= m2hat/(1-a2))
 *   Aggregates are divided by X_t = A1_t^(1/(1-a1)) * L_t, which grows at
 *   G = (1+m1hat)^(1/(1-a1))*(1+n), i.e. g* = n + m1hat/alpha1 = 2% (see p.30).
 *   The remaining non-stationary object is the exogenous index
 *       Omega_t = A2_t / A1_t^((1-a2)/(1-a1)),
 *   which drives non-balanced growth: z = (Y2/Y1)^((eps-1)/eps) -> 0,
 *   kappa = K1/K -> 1, lambda = L1/L -> 1 (Theorem 2).
 *   The allocation conditions (29)-(32) reduce to closed forms in z,
 *   equivalent to eqs. (33)-(34).
 *
 * Initial conditions: 1948 NIPA data (paper p.30): L, K, nominal Y1, Y2. As in the
 * paper, only the composite technology terms matter, so gam=0.5 is a normalisation
 * and A1(0), A2(0) are backed out from the data.
 * Terminal condition: conditional steady state with Omega fixed at Omega_{T+1}.
 * Omega moves very slowly, so this only distorts the last ~100 periods: use T large
 * and discard the end of the simulation.
 *
 * (C) 2026 Max Gödl, University of the Bundeswehr Munich. Claude Opus 5.5 
 * was used in developing this code.
 */
//----------------------------------------------------------------------------


@#ifndef T
  @#define T = 3000
@#endif
@#ifndef EPS_LIST
  @#define EPS_LIST = [0.5, 0.4, 0.6]   // first entry = benchmark
@#endif


var k        ${k}$          (long_name='capital K/X (end of period)')
    c        ${c}$          (long_name='consumption C/X')
    inv      ${i}$          (long_name='net investment I/X')
    y        ${y}$          (long_name='output Y/X')
    y1       ${y_1}$        (long_name='sector-1 output Y1/X')
    R        ${r}$          (long_name='rental rate of capital (interest rate)')
    z        ${z}$          (long_name='(Y2/Y1)^((eps-1)/eps)')
    Phi      ${\Phi}$       (long_name='gam+(1-gam)z')
    kappa    ${\kappa=K_1/K}$     (long_name='capital share sector 1')
    lambda   ${\lambda=L_1/L}$    (long_name='employment share sector 1')
    nom_share1 ${\sigma_1}$ (long_name='nominal output share sector 1, p1Y1/Y')
    cap_share  ${\sigma_K}$ (long_name='capital share in national income')
    Y2_Y1    ${Y_2/Y_1}$    (long_name='relative real output')
    L2_L1    ${L_2/L_1}$    (long_name='relative employment')
    log_p2_p1  ${\log(p_2/p_1)}$ (long_name='log relative price sector 2')
    g_Y      ${g_Y}$        (long_name='growth aggregate output')
    g_C      ${g_C}$        (long_name='growth aggregate consumption')
    g_Y1     ${g_{Y1}}$     (long_name='growth sector-1 output')
    g_Y2     ${g_{Y2}}$     (long_name='growth sector-2 output')
    ;

varexo log_Omega ${\log\Omega}$ (long_name='log A2/A1^((1-a2)/(1-a1))');

parameters a1 a2 eps gam theta beta delta m1hat m2hat n G gO k_init logOmega0;

//----------------------------------------------------------------------------
// Calibration, Section 4 (annual)
//----------------------------------------------------------------------------
a1     = 1-0.72;      // capital share sector 1 (paper: labour share alpha1 = 0.72)
a2     = 1-0.52;      // capital share sector 2 (paper: alpha2 = 0.52)
gam    = 0.5;         // normalisation (only B_i composites matter, see p.29)
n      = 0.01;        // population growth
m1hat  = 0.0072;      // Hicks-neutral TFP growth sector 1, = m1/(nu-1)
m2hat  = 0.0072;      // benchmark: equal to m1hat
theta  = 4;           // from g* = (r*-rho)/theta with r*=0.10, rho=0.02, g*=0.02
beta   = 1/1.02;      // rho = 0.02
delta  = 0;           // paper has no depreciation

G  = (1+m1hat)^(1/(1-a1))*(1+n);                   // gross growth of X (= 1+g*)
gO = log(1+m2hat) - (1-a2)/(1-a1)*log(1+m1hat);    // growth of log Omega

// 1948 NIPA data (paper p.30), millions of current $ / thousands FTE
L0  = 37169;
K0 = 560000;
Y1N0 = 85885;
Y2N0 = 108473;
YN0 = Y1N0 + Y2N0;
sig1_0   = Y1N0/YN0;  // p1Y1/Y with P=1

model;
    [name='Law of motion capital']
    k = (1-delta)*k(-1)/G + inv;
    [name='Resource constraint']
    c + inv = y;
    [name='Euler equation, aggregate C']
    c^(-theta) = beta*G^(-theta)*c(+1)^(-theta)*(1 + R(+1) - delta);
    [name='Sector-1 production']
    y1 = (kappa*k(-1)/G)^a1 * lambda^(1-a1);
    [name='Relative output index']
    z = ( exp(log_Omega)*(k(-1)/G)^(a2-a1)
          *(1-kappa)^a2*(1-lambda)^(1-a2)
          /(kappa^a1*lambda^(1-a1)) )^((eps-1)/eps);
    [name='Aggregator helper']
    Phi = gam + (1-gam)*z;
    [name='Final output']
    y = y1*Phi^(eps/(eps-1));
    [name='Interest rate, eq. (31)']
    R = gam*a1*Phi^(1/(eps-1))*y1/(kappa*k(-1)/G);
    [name='Capital allocation, eqs. (31)-(32)']
    kappa = 1/(1 + (1-gam)*a2*z/(gam*a1));
    [name='Labour allocation, eqs. (29)-(30)']
    lambda = 1/(1 + (1-gam)*(1-a2)*z/(gam*(1-a1)));

    // --- reporting variables ---
    nom_share1 = gam/Phi;
    [name='Capital share in national income, eq. (38)']
    cap_share  = a1*nom_share1 + a2*(1-nom_share1);
    Y2_Y1 = z^(eps/(eps-1));
    L2_L1 = (1-lambda)/lambda;
    log_p2_p1 = log((1-gam)/gam) - log(z)/(eps-1);
    g_Y  = y/y(-1)*G - 1;
    g_C  = c/c(-1)*G - 1;
    g_Y1 = y1/y1(-1)*G - 1;
    g_Y2 = y1/y1(-1)*Y2_Y1/Y2_Y1(-1)*G - 1;
end;

// eps for the initial guess: first entry of EPS_LIST
eps = @{EPS_LIST[1]};

//----------------------------------------------------------------------------
// Initial guess, set once before the loop (Dynare does not allow initval after
// endval). Only a starting value: k(0) is overwritten with the 1948 capital stock
// in each run, and the other period-0 values only enter growth rates in period 1.
//----------------------------------------------------------------------------
initval;
    log_Omega = 0;
    R      = G^theta/beta - 1 + delta;
    k      = G*(a1*gam^(eps/(eps-1))/R)^(1/(1-a1));
    kappa  = 0.5;
    lambda = 0.6;
    z      = 0.5;
    Phi    = gam + (1-gam)*z;
    y1     = (kappa*k/G)^a1*lambda^(1-a1);
    y      = y1*Phi^(eps/(eps-1));
    inv    = (1-(1-delta)/G)*k;
    c      = y - inv;
    nom_share1 = gam/Phi;
    cap_share  = a1*nom_share1 + a2*(1-nom_share1);
    Y2_Y1 = z^(eps/(eps-1));
    L2_L1 = (1-lambda)/lambda;
    log_p2_p1 = log((1-gam)/gam) - log(z)/(eps-1);
    g_Y = G-1; g_C = G-1; g_Y1 = G-1; g_Y2 = G-1;
end;

steady;


iv   = @(s) find(strcmp(s, M_.endo_names));
sims = {};
@#for eps_val in EPS_LIST
//----------------------------------------------------------------------------
// Run for eps = @{eps_val}
//----------------------------------------------------------------------------
eps = @{eps_val};

// eps-dependent initial values: technology levels backed out from 1948 data
kappa_0  = 1/(1 + (a2/a1)*(1-sig1_0)/sig1_0);            // eq. (60): 0.32
lambda_0 = 1/(1 + ((1-a2)/(1-a1))*(1-sig1_0)/sig1_0);    // eq. (34): 0.52
Y1_0 = YN0*(sig1_0/gam)^(eps/(eps-1));                   // from p1 = gam (Y1/Y)^(-1/eps)
Y2_0 = YN0*((1-sig1_0)/(1-gam))^(eps/(eps-1));
A1_0 = Y1_0/((kappa_0*K0)^a1*(lambda_0*L0)^(1-a1));
A2_0 = Y2_0/(((1-kappa_0)*K0)^a2*((1-lambda_0)*L0)^(1-a2));
X_0  = A1_0^(1/(1-a1))*L0;
k_init    = G*K0/X_0;                                    // so that k(-1)/G = K_1948/X_1948
logOmega0 = log(A2_0) - (1-a2)/(1-a1)*log(A1_0);

//----------------------------------------------------------------------------
// Exogenous path: period 1 = 1948, Omega grows at gO
//----------------------------------------------------------------------------
Omega_path = logOmega0 + gO*(0:@{T}-1);

shocks(overwrite);
    var log_Omega;
    periods 1:@{T};
    values (Omega_path);
end;

//----------------------------------------------------------------------------
// Terminal condition: conditional steady state at Omega_{T+1}
//----------------------------------------------------------------------------
endval;
    log_Omega = logOmega0 + gO*@{T};
end;
steady;

perfect_foresight_setup(periods=@{T});

// some hand-holding for initial values
if ~isempty(sims)
    oo_.endo_simul(:, 2:end-1) = sims{1}(:, 2:end-1);
end
//----------------------------------------------------------------------------
// Initial capital stock from 1948 data (histval cannot be combined with endval,
// so overwrite the initial condition directly)
//----------------------------------------------------------------------------
oo_.endo_simul(strcmp('k', M_.endo_names), 1) = k_init;

perfect_foresight_solver;

// Period 1 = 1948, period 57 = 2004 (compare Table 1 of the paper)
fprintf('\neps = %4.2f      1948     2004\n', eps);
fprintf('Y2/Y1     %6.3f   %6.3f\n', oo_.endo_simul(iv('Y2_Y1'),[2 58]));
fprintf('L2/L1     %6.3f   %6.3f\n', oo_.endo_simul(iv('L2_L1'),[2 58]));
fprintf('sigma_K   %6.3f   %6.3f\n', oo_.endo_simul(iv('cap_share'),[2 58]));
fprintf('r         %6.3f   %6.3f\n', oo_.endo_simul(iv('R'),[2 58]));
fprintf('Y2/Y1 change: %+.1f%%,  L2/L1 change: %+.1f%%\n', ...
    100*(oo_.endo_simul(iv('Y2_Y1'),58)/oo_.endo_simul(iv('Y2_Y1'),2)-1), ...
    100*(oo_.endo_simul(iv('L2_L1'),58)/oo_.endo_simul(iv('L2_L1'),2)-1));

sims{end+1} = oo_.endo_simul;
@#endfor


//----------------------------------------------------------------------------
// Figure 3: lambda, kappa, r, sigma_K for all eps in EPS_LIST
// (start at column 2 = 1948)
//----------------------------------------------------------------------------
eps_vals  = [@{EPS_LIST}];
plot_vars = {'lambda','kappa','R','cap_share'};
plot_ttl  = {'Sector 1 labor share, \lambda = L_1/L','Sector 1 capital share, \kappa = K_1/K','Interest rate, r','Capital income share, \sigma_K'};
styles    = {'-','--',':','-.'};
Tplot = @{T} - 100;
t     = 0:Tplot-1;
figure('Name','AG (2008), Figure 3');
for p = 1:numel(plot_vars)
    subplot(2,2,p); hold on;
    for j = 1:numel(sims)
        plot(t, sims{j}(iv(plot_vars{p}), 2:Tplot+1), styles{mod(j-1,4)+1}, 'LineWidth', 1.2);
    end
    hold off; box off; title(plot_ttl{p}); xlabel('years after 1948'); xlim([0 Tplot-1]);
end
subplot(2,2,1);
legend(arrayfun(@(e) sprintf('\\epsilon = %.1f', e), eps_vals, 'UniformOutput', false), 'Location', 'southeast');
print(gcf, 'AG_figure3_replication.png', '-dpng', '-r300');

figure('Name','AG (2008), Aggregate growth rate')
hold on;
for j = 1:numel(sims)
    plot(1:Tplot-1, sims{j}(iv('g_Y'), 3:Tplot+1), styles{mod(j-1,4)+1}, 'LineWidth', 1.2); 
end
hold off; box off; title('Aggregate output growth, g_Y'); xlabel('years after 1948'); xlim([1 Tplot]);
legend(arrayfun(@(e) sprintf('\\epsilon = %.1f', e), eps_vals, 'UniformOutput', false), 'Location', 'northeast');
print(gcf, 'AG_agg_growth.png', '-dpng', '-r300');