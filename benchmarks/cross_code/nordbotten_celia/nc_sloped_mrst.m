%% nc_sloped_mrst.m
% MRST co2lab mirror of wombat test/tests/nordbotten_celia/nc_sloped_fv.i
% Verification-ladder step 2: 1D sloped aquifer, constant props, no trap, S_wr=0.
%
% Case A: CO2 forced INCOMPRESSIBLE (constant rho = 700) to match wombat's
% ConstantFluidProperties. makeVEFluid otherwise builds a real compressible CO2
% EOS (rho = 700 only near ~93 bar; ~1.8 kg/m3 at 0.3 bar), which floods the
% domain when run at the low p=0-referenced pressures used here.
%
% Analytical anchors at t = 2e6 s:
%   injected mass  = Q*rho*t      = 4e-6*700*2e6 = 5600 kg   (per m width)
%   mass integral  = integral sat_n dx = Q*t/(phi*H) = 2.0 m
%   gravity nose   = v_N*t,  v_N = K*drho*g*sin(theta)/(phi*mu_n) = 1.096e-5 m/s
%                  = 21.9 m
%
% Run from MRST after: mrstModule add co2lab ad-core ad-blackoil

mrstModule add co2lab ad-core ad-blackoil
gravity reset on

%% ── GRID: 1D strip, 100 m, 200 cells, unit width, H = 20 m ───────────────────
nx = 200; Lx = 100; Wy = 1.0; H = 20.0;
G3D = cartGrid([nx, 1, 1], [Lx, Wy, H]);

% wombat z_top = 0.034899*x is an ELEVATION (updip = +x). MRST uses depth
% (positive down), so depth decreases with x. The +200 datum keeps depth
% positive; with incompressible fluid the absolute value is irrelevant -- only
% grad(z_top) drives buoyancy.
sgn = 0.034899;                       % sin(2 deg)
xn  = G3D.nodes.coords(:,1);
G3D.nodes.coords(:,3) = (200 - sgn*xn) + G3D.nodes.coords(:,3);
G3D = computeGeometry(G3D);
Gt  = topSurfaceGrid(G3D);

%% ── ROCK ─────────────────────────────────────────────────────────────────────
rock2D = makeRock(Gt, 1e-12, 0.20);    % K = 1e-12 m^2, phi = 0.20

%% ── FLUID (forced incompressible -- Case A) ──────────────────────────────────
fluid = makeVEFluid(Gt, rock2D, 'sharp_interface_simple', ...
    'co2_mu_ref', 5e-5, 'wat_mu_ref', 8e-4, ...
    'co2_rho_ref', 700, 'wat_rho_ref', 1020, ...
    'residual', [0.0, 0.0], 'dissolution', false);   % S_wr = 0

% Override formation-volume factors to constants => rho is pressure-independent
% (rhoG = rhoGS = 700, rhoW = rhoWS = 1020), matching wombat ConstantFluidProperties.
% '1 + 0*p' keeps AD type with a zero derivative (true zero compressibility).
fluid.bG = @(p, varargin) 1 + 0*p;
fluid.bW = @(p, varargin) 1 + 0*p;

model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);

%% ── IC / BC ──────────────────────────────────────────────────────────────────
% Incompressible => absolute pressure level is irrelevant; keep wombat's p=0
% reference at the updip (x=Lx) open boundary.
g   = norm(gravity);
xcc = Gt.cells.centroids(:,1);
state0       = initResSol(Gt, 1020*g*sgn*(Lx - xcc), [1, 0]);   % [water, CO2]
state0.sGmax = state0.s(:,2);

bf = boundaryFaces(Gt);
rightFaces = bf(abs(Gt.faces.centroids(bf,1) - Lx) < 1e-3);
bc = addBC([], rightFaces, 'pressure', 0.0, 'sat', [1, 0]);     % brine inflow only

%% ── INJECTION: CO2 flux across the x=0 boundary, Q = 4e-6 m^3/s ──────────────
% Mirrors wombat's FVNeumannBC (flux INTO the domain across the downdip face),
% NOT a volumetric cell source -- a cell source raises local pressure and adds a
% viscous updip drive on top of buoyancy, which a boundary flux does not. With
% ny=1 the left boundary is a single face, so 'flux' = 4e-6 is the total rate.
leftFaces = bf(abs(Gt.faces.centroids(bf,1) - 0) < 1e-3);
bc = addBC(bc, leftFaces, 'flux', 4e-6, 'sat', [0, 1]);         % pure CO2 in

%% ── STAGE 1: short closed-system check (t = 2e5 s, nose stays << 100 m) ──────
sched_short = simpleSchedule(repmat(1e4, 20, 1), 'bc', bc);
[wellSols_s, sts_s] = simulateScheduleAD(state0, model, sched_short);

t_short = 2e5;
co2_vol = sum(sts_s{end}.s(:,2) .* poreVolume(Gt,rock2D) .* Gt.cells.H);  % H put back
nose_s  = max([0; xcc(sts_s{end}.s(:,2) > 1e-3)]);

fprintf('\n=== STAGE 1: short closed check (t = %.0e s) ===\n', t_short);
fprintf('  domain CO2 vol : %.4e m3   (injected R*t = %.4e)\n', co2_vol, 4e-6*t_short);
fprintf('  ratio          : %.3f       (want 1.0)\n', co2_vol/(4e-6*t_short));
fprintf('  nose           : %.2f m     (analytic 2.2; must be << 100)\n', nose_s);

%% ── STAGE 2: full run (t = 2e6 s) -- compare to the three analytical anchors ─
sched = simpleSchedule(repmat(1e4, 200, 1), 'bc', bc);
[~, states] = simulateScheduleAD(state0, model, sched, 'Verbose', true);

sfin      = states{end}.s(:,2);
dx        = Lx/nx;
mass_int  = sum(sfin) * dx;                                  % expect 2.0 m
mass_kg   = sum(sfin .* poreVolume(Gt,rock2D) .* Gt.cells.H) * 700;  % expect 5600 kg
nose      = max([0; xcc(sfin > 1e-3)]);                     % expect 21.9 m
[~, ord]  = sort(xcc);                                       % inlet = smallest x
inlet_sat = sfin(ord(1));                                    % wombat FV: 0.0913

fprintf('\n=== STAGE 2: full run (t = 2e6 s) vs wombat FV ===\n');
fprintf('  mass integral : %.3f m    (analytic 2.0)\n',  mass_int);
fprintf('  injected mass : %.1f kg   (analytic 5600)\n', mass_kg);
fprintf('  inlet sat     : %.4f      (wombat FV 0.0913)  <-- key check\n', inlet_sat);
fprintf('  nose position : %.1f m    (analytic 21.9; diffuse in FV)\n', nose);

%% ── CENTRELINE OUTPUT (mirror wombat VectorPostprocessors/centreline_sat) ────
% Write CO2 saturation vs x at every step, for overlay against wombat's
% anticline_a1_csv_centreline_sat_*.csv-style profiles.
for i = 1:numel(states)
    T = table(xcc, states{i}.s(:,2), 'VariableNames', {'x_m','sat_co2'});
    writetable(T, sprintf('nc_mrst_sat_%04d.csv', i));
end
fprintf('\nWrote %d centreline CSVs (nc_mrst_sat_*.csv).\n', numel(states));
