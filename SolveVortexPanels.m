%Function to solve for strength of the vortex panels (gamma), the surface
%velocity at each panel, the pressure distribution across the airfoil, and
%the circulation + resultant lift

function [gamma, Vt, Cp, Circulation, Lift, Cl_tot] = SolveVortexPanels(K, L, U, betaR, numPan, S, rho)
A = zeros(numPan, numPan); %Influence coefficient storage matrix

for i = 1:numPan
    for j = 1:numPan
        if (i == 1 && j == 1) || (i == 1 && j == numPan) %enforce the Kutta condition at lower TE 
            A(i,j) = 1;
        elseif (i == j && i ~= 1)
            A(i,j) = 0; %no self influence 
        elseif (i~= j && i~= 1)
            A(i,j) = -K(i,j); %assign coefficient value
        end
    end
end

b = zeros(numPan,1); %normal free-stream terms vector (right hand side of equation)

for n = 1:numPan
    if n ~= 1
        b(n) = -U*2*pi*cos(betaR(n)); %normal compute free-stream terms
    else
        b(n) = 0; %enforce Kutta Condition at lower TE
    end
end

gamma = A\b; %SOLVE LINEAR SYSTEM OF EQUATIONS

gamma_ds = gamma(:).*S(:); %sum of all the contributions to circulation from each panel
Circulation = sum(gamma_ds); %CALCULATE NET CIRCULATION
Lift = rho*U*Circulation; %CALCULATE NET LIFT

%========== Surface Velocity and Aerodynamic Loads ==========%
Vt = zeros(numPan,1);  %Panel tangential velocity vector
Cp = zeros(numPan,1); %Panel pressure coefficient vector
Cl = zeros(numPan,1); %Panel lift coefficient vector

for i = 1:numPan
    gamma_vals = 0;
    for j = 1:numPan
        if i ~= j
            gamma_vals = gamma_vals + (gamma(j)/(2*pi))*(L(i,j)); %cotntribution of all other jth panels to tangent velocity of the ith panel
    end
    Vt(i) = U*sin(betaR(i))-gamma_vals; %tangential velocity of the ith panel
    Cp(i) = 1-(Vt(i)/U)^2; %pressure coefficient evaluated at ith control point
    Cl(i) = (2*Lift)/(rho*Vt(i)^2*S(i)); %lift coefficient evaluated at the ith panel
    end
end

S_tot = sum(S); %total length of curve bounding airfoil
F_lift = sum(Cl(:).*S(:)); %comparison lift formula
dp = 0.5*rho*U; %dynamic pressure term
Cl_tot = F_lift/(dp*S_tot); %NET LIFT COEFFICIENT OF AIRFOIL