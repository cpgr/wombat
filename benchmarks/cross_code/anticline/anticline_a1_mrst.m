%% anticline_a1_mrst.m
% MRST co2lab mirror of benchmarks/cross_code/anticline/anticline_a1.i
% Cross-code benchmark Tier A1: conceptual anticline, constant density, sharp
% interface, S_wr = 0.2, no trapping/dissolution. 2D, 10x10 km, dome at (7,5) km.
%
% This carries EVERY fix learned on the nordbotten_celia rung -- see
% benchmarks/cross_code/nordbotten_celia/ and doc/content/benchmarks/cross_code.md:
%   1. CO2 forced INCOMPRESSIBLE (bG=bW=1) to match wombat ConstantFluidProperties
%      -- makeVEFluid's CO2 is a real EOS; at the ~25 bar here it is ~1.7x off.
%   2. Phase order black-oil: s(:,1)=water, s(:,2)=CO2; init [1,0]; set sGmax.
%   3. CO2 volume / mass uses poreVolume(Gt,rock).*Gt.cells.H (poreVolume omits H).
%   4. Capillary spreading term: MRST's VE model includes Pc=drho*g*h intrinsically.
%      The wombat side MUST run with capillary=true (now wired in anticline_a1.i)
%      for this to be like-for-like -- it governs dome fill/spill.
%
% Compare against wombat's anticline_a1_csv.csv: mobile_co2_mass and footprint vs
% time are the robust cross-code metrics (NOT crest saturation).
%
% Run in MATLAB after: mrstModule add co2lab ad-core ad-blackoil

mrstModule add co2lab ad-core ad-blackoil
gravity reset on

%% ── GRID: 10x10 km, 50x50, H=50 m, regional dip + Gaussian dome ──────────────
nx = 50; ny = 50; Lx = 1e4; Ly = 1e4; H = 50.0;
G3D = cartGrid([nx, ny, 1], [Lx, Ly, H]);

% wombat z_top is an ELEVATION: -800 + 0.0261859*x + 40*exp(-r^2/(2*1200^2)).
% MRST uses depth (positive down): depth_top = -elevation.
slope = 0.0261859;          % tan(1.5 deg) regional dip, updip = +x
xc = 7000; yc = 5000; amp = 40; sig = 1200;   % dome
xn = G3D.nodes.coords(:,1); yn = G3D.nodes.coords(:,2);
depth_top = 800 - slope*xn - amp*exp(-((xn-xc).^2 + (yn-yc).^2)/(2*sig^2));
G3D.nodes.coords(:,3) = depth_top + G3D.nodes.coords(:,3);   % top->depth_top, bottom->+H
G3D = computeGeometry(G3D);
Gt = topSurfaceGrid(G3D);

%% ── ROCK: 200 mD isotropic, phi = 0.20 ──────────────────────────────────────
rock2D = makeRock(Gt, 200*9.869e-16, 0.20);     % 1.9738e-13 m^2

%% ── FLUID: constant rho/mu, sharp interface, S_wr = 0.2; forced incompressible
fluid = makeVEFluid(Gt, rock2D, 'sharp_interface_simple', ...
    'co2_mu_ref', 6e-5, 'wat_mu_ref', 8e-4, ...
    'co2_rho_ref', 700, 'wat_rho_ref', 1000, ...
    'residual', [0.0, 0.2], 'dissolution', false);
fluid.bG = @(p, varargin) 1 + 0*p;     % rho_CO2  = 700  everywhere (incompressible)
fluid.bW = @(p, varargin) 1 + 0*p;     % rho_brine= 1000 everywhere

% RELPERM: use MRST's native 'sharp_interface_simple' linear curves -- NO override.
% These are krG(sG)=sG, krW(sw)=sw, which never reach krW=0 (at a filled column
% sG=0.8 -> krW(0.2)=0.2), so the brine equation stays well-conditioned and the
% solve completes. The wombat side is matched to THIS by setting its
% VERelPermSharpUO S_wr=0 (giving kr_n=sat_n, kr_w=1-sat_n) while keeping the
% capillary column height at S_wr=0.2 (residual sets the height, not the relperm).
%   Trade-off: the renormalised form kr_n=sat_n/(1-S_wr) is the more accurate VE
%   upscaling, but it drives krW->0 at a filled column and stalls MRST's stock
%   direct solver (wombat survives that via PETSc '-pc_factor_shift_type NONZERO').
%   Using the simple form in BOTH codes gives a consistent, solvable comparison.

model = CO2VEBlackOilTypeModel(Gt, rock2D, fluid);

% Limit the Newton saturation step for robustness through the dome-filling phase.
model.dsMaxAbs = 0.1;

% Native simple relperm: krG(sG)=sG. wombat is matched by setting S_wr=0 there,
% which gives kr_n(sat_n)=sat_n -- so these should read 0.8 and 0.4.
fprintf('krG(0.8) = %.4f   (wombat S_wr=0: kr_n(0.8)=0.8)\n', fluid.krG(0.8));
fprintf('krG(0.4) = %.4f   (wombat S_wr=0: kr_n(0.4)=0.4)\n', fluid.krG(0.4));

%% ── IC: hydrostatic brine, zero gauge at updip (x=Lx) open boundary ──────────
g   = norm(gravity);
xcc = Gt.cells.centroids(:,1);
ycc = Gt.cells.centroids(:,2);
state0       = initResSol(Gt, 1000*g*slope*(Lx - xcc), [1, 0]);   % [water, CO2]
state0.sGmax = state0.s(:,2);

%% ── BOUNDARY: open updip (right) at p=0; CO2 leaves naturally ────────────────
bf = boundaryFaces(Gt);
rightFaces = bf(abs(Gt.faces.centroids(bf,1) - Lx) < 1e-3);
bc_open = addBC([], rightFaces, 'pressure', 0.0, 'sat', [1, 0]);

%% ── INJECTION: 1 Mt/yr CO2 flux across the x=0 patch, y in [4800,5200] ───────
% Mirrors wombat's FVNeumannBC on the injection_zone. Total mass 1 Mt/yr ->
% volume rate (incompressible) = (1e9 kg/yr / 3.1557e7 s/yr) / 700 = 0.04527 m3/s.
yr = 3.1557e7;
Q_vol = (1e9 / yr) / 700;                                   % 0.04527 m^3/s total
leftInj = bf(abs(Gt.faces.centroids(bf,1) - 0) < 1e-3 & ...
             Gt.faces.centroids(bf,2) >= 4800 & Gt.faces.centroids(bf,2) <= 5200);
% addBC 'flux' assigns the value PER FACE, so divide the total by the face count
% (otherwise N injection faces deliver N x too much -- this gave 20 Mt at 2 faces).
bc_inj = addBC(bc_open, leftInj, 'flux', Q_vol/numel(leftInj), 'sat', [0, 1]);
fprintf('Injection: %d faces, %.4e m3/s each, %.4e m3/s total\n', ...
        numel(leftInj), Q_vol/numel(leftInj), Q_vol);

%% ── SCHEDULE: 10 yr injection, then migration to 200 yr ──────────────────────
% Two controls: ctrl1 carries the injection BC, ctrl2 drops it (open updip only).
% BC sets (unlike wells) may change between controls, so no facility issue here.
dt_inj  = repmat(0.5*yr, 20, 1);     % 0-10 yr
dt_post = repmat(2.0*yr, 95, 1);     % 10-200 yr
schedule.control(1) = struct('W',[],'bc',bc_inj, 'src',[]);
schedule.control(2) = struct('W',[],'bc',bc_open,'src',[]);
schedule.step.val     = [dt_inj; dt_post];
schedule.step.control = [ones(numel(dt_inj),1); 2*ones(numel(dt_post),1)];

%% ── RUN (with crash-recovery output handler) ─────────────────────────────────
nls = NonLinearSolver('useRelaxation', true, 'maxTimestepCuts', 12, 'maxIterations', 30);
oh  = ResultHandler('storeInMemory', true, 'clearAfterRead', false);
try
    simulateScheduleAD(state0, model, schedule, 'NonLinearSolver', nls, ...
                       'OutputHandler', oh, 'Verbose', true);
catch err
    fprintf('\nStopped at step %d: %s\n', numel(oh.getValidIds()), err.message);
end
% Defensive collection: getValidIds can report an index past what is stored when
% the run stops mid-step, so stop at the first gap rather than indexing past it.
ids = oh.getValidIds();
states = {};
for k = ids(:)'
    try, states{end+1} = oh{k}; catch, break; end %#ok<AGROW>
end
ids = ids(1:numel(states));          % keep ids aligned with what was collected
fprintf('Completed %d / %d steps.\n', numel(states), numel(schedule.step.val));

%% ── INJECTION GATE: CO2 in place at end of injection (step 20) should be ~10 Mt
pv = poreVolume(Gt, rock2D) .* Gt.cells.H;     % full 3D pore volume (H put back)
if numel(states) >= 20
    m10 = sum(states{20}.s(:,2) .* pv) * 700;
    fprintf('CO2 in place at t=10yr: %.2f Mt (expect ~10)\n', m10/1e9);
end

%% ── OBSERVABLES vs wombat anticline_a1_csv.csv (time, mobile mass, footprint) ─
tall = cumsum(schedule.step.val);
t = tall(ids);  t = t(:);                                      % column
mobile = cellfun(@(s) sum(s.s(:,2).*pv)*700, states);          % kg
foot   = cellfun(@(s) sum(Gt.cells.volumes(s.s(:,2)>1e-3)), states);  % m^2
mobile = mobile(:);  foot = foot(:);                           % columns to match t
T = table(t, mobile, foot, 'VariableNames', {'time_s','mobile_co2_kg','footprint_m2'});
writetable(T, 'anticline_mrst_observables.csv');
fprintf('Wrote anticline_mrst_observables.csv (%d rows).\n', height(T));

% Centreline CO2 saturation (y=4900 row) per step, to overlay on wombat's
% anticline_a1_csv_centreline_sat_*.csv if a profile comparison is wanted.
row = abs(ycc - 4900) < 1e-3;
[~, ord] = sort(xcc(row)); rc = find(row); rc = rc(ord);
for i = 1:numel(states)
    writetable(table(xcc(rc), states{i}.s(rc,2), 'VariableNames', {'x_m','sat_co2'}), ...
               sprintf('anticline_mrst_centreline_%04d.csv', i));
end

fprintf('\nCompare anticline_mrst_observables.csv against wombat anticline_a1_csv.csv\n');
fprintf('(mobile_co2_mass and footprint vs time -- the robust cross-code metrics).\n');
