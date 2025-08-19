%============NUMERICAL VORTEX PANEL METHOD CODE============%
clc
clear

%========== Define Knowns ==========%
U = 1; %free stream velocity
alphaD = 5; %angle of attach in degrees
c = 1; %chord length of foil
rho = 1.035; %density of air 

%========== Define Airfoil Profile ==========%
def_foil = 'Use .dat File'; %variable to control how profile is generated

[XB, YB, XC, YC, phiR, betaR, S, numPan] = LoadPanels(def_foil, c, alphaD)

%======= Determine Normal & Tangentential Infuence Coefficients =======%
[K, L] = VPM_InfluenceCoeff(XC, YC, XB, YB, phiR, S, numPan); %function solving for influence coefficients

%========== Solve Linear System of Equations ==========%
[gamma, Vt, Cp, Circulation, Lift, Cl_tot] = SolveVortexPanels(K, L, U, betaR, numPan, S, rho)


%========== Plot Streamlines ==========%
[Nxx , Nyy, Vxy, rp, psi, THETA, Cpxy_mask] = PM_streamlines(XC, YC, XB, YB, phiR, S, gamma, U, alphaD, Cp, numPan)

%========== Plot Coefficients =========%
PlotC_(XC, Cp, alphaD, numPan)