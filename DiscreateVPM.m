%============NUMERICAL VORTEX PANEL METHOD CODE============%
clc
clear

%========== Define Knowns ==========%
U = 1; %free stream velocity
alphaD = 10; %angle of attach in degrees
c = 1; %chord length of foil

%========== Define Airfoil Profile ==========%
def_foil = 'Use .dat File'; %variable to control how profile is generated

switch def_foil
    case 'Cosine'
        N = 71; %select number of boundary points

        theta_U = linspace(pi,0,N)'; %analytical angle for upper surface cosine spacing
        theta_L = linspace(0,pi,N)'; %analytical angle for lower surface cosine spacing

        xU = (c/2)*(1 - cos(theta_U)); %upper surface x-coords
        xL = (c/2)*(1 - cos(theta_L)); %lower surface y-coords
        XB = vertcat(xL, xU(2:end)); %x-coords starting at TE and traveling clockwise

        for i = 1:N
            yU(i,1) = 0.594689181*(0.298222773*sqrt(xU(i)) - 0.127125232.*xU(i) - 0.357907906.*xU(i).^2 + 0.291984971.*xU(i).^3 - 0.105174606.*xU(i).^4);
        end
     
        yL = -yU; %lower surface y-coords
        YB = vertcat(yL, yU(2:end)); %y-coords starting at trailing edge and traveling clockwise

    case 'Use .dat File'
        load('foilData.dat');
        N = length(foilData(:,1));
        XB = foilData(:,1);
        YB = -foilData(:,2);

end

%========== Define Panel Orientation and Geometry =========%
numPan = N-1; %N-1 panels for N coordinate points
XC = zeros(N-1, 1); %X control point coordinate vector
YC = zeros(N-1, 1); %Y control point coordinate vector
S = zeros(N-1, 1); %panel lengths vector
phiD = zeros(N-1, 1); %panel orientation vector

dX = zeros(N-1, 1); %XB(i+1)-XB(i) --- for use in panel orientation and length
dY = zeros(N-1, 1); %YB(i+1)-YB(i) --- for use in panel orientation and length

for i = 1:numPan
    XC(i) = 0.5*(XB(i+1)+XB(i)); %X-coordinate of control point of ith panel
    YC(i) = 0.5*(YB(i+1)+YB(i)); %Y-coordinate of control point of ith panel
    dX(i) = (XB(i+1)-XB(i)); %change in x along ith panel
    dY(i) = (YB(i+1)-YB(i)); %change in y along ith panel
    S(i) = sqrt(dX(i)^2 + dY(i)^2); %length of ith panel
    phiD(i) = atan2d(dY(i), dX(i)); %angle of ith panel with respect to positive x-axis
    if (phiD(i) < 360)
        phiD(i) = phiD(i) + 360; %makes all angles positive angles
    end
end

deltaD = phiD + 90; %angle from positive x-axis to outward facing unit normal
betaD = deltaD-alphaD; %angle between freestream velocity vector and outward facing unit normal
betaD(betaD>360) = betaD(betaD>360)-360; %make angles less than 360

phiR = phiD.*(pi/180); %angle phi converted to radians
betaR = betaD.*(pi/180); %angle beta converted to radians

%======= Determine Normal & Tangentential Infuence Coefficients =======%
[K, L] = VPM_InfluenceCoeff(XC, YC, XB, YB, phiR, S); %function solving for influence coefficients
for i = 1:length(L(1,:))
    L(i,i) = 0.5
end
%========== Solve Linear System of Equations ==========%
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

RHS = zeros(numPan,1); %normal free-stream terms vector (right hand side of equation)

for n = 1:numPan
    if n ~= 1
        RHS(n) = -U*2*pi*cos(betaR(n)); %normal compute free-stream terms
    else
        RHS(n) = 0; %enforce Kutta Condition at lower TE
    end
end

gamma = A\RHS; %SOLVE LINEAR SYSTEM OF EQUATIONS
Circulation = sum(gamma(:).*S(:)); %sum of each panel strength (per unit length) multiplied by panel length

%========== Surface Velocity and Aerodynamic Loads ==========%

Vt = zeros(numPan,1);  %Panel tangential velocity vector
Cp = zeros(numPan,1); %Panel pressure coefficient vector

for i = 1:numPan
    gamma_vals = 0;
    for j = 1:numPan
        if i ~= j
            gamma_vals = gamma_vals + (gamma(j)/(2*pi))*(L(i,j)); %cotntribution of all other jth panels to tangent velocity of the ith panel
    end
    Vt(i) = U*sin(betaR(i))-gamma_vals; %tangential velocity of the ith panel
    Cp(i) = 1-(Vt(i)/U)^2; %pressure coefficient evaluated at ith control point
    end
end


%========== Plot Streamlines ==========%
[Nx , Ny] = VPM_xyCoeff(XC, YC, XB, YB, phiR, S); %Compute X and Y velocity influence coefficients for Meshgrid

[Nxx , Nyy, Vxy, rp, psi, THETA, Cpxy_mask] = PM_streamlines(XC, YC, XB, YB, phiR, S, gamma, U, alphaD, Cp, Nx, Ny, numPan)

half_x = floor(numPan/2);
figure; hold on;
axis([0 1 -1.1 1.1]);
plot(XC(1:half_x), Cp(1:half_x), 'ro');
plot(XC(half_x+1:end), Cp(half_x+1:end), 'bo'); 
plot(XC, Cp, 'k')
set(gca, 'YDir','reverse')
title(['Pressure Distribution on Airfoil Surface ($\alpha = ', num2str(alphaD), ')$'], 'Interpreter','latex');
xlabel('X-Coordinate of Airfoil');
ylabel('Coefficient of Pressure (Cp)');
legend('Bottom Cp', 'Top Cp');