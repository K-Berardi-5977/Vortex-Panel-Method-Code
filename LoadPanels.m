%========== FUNCTION TO DEFINE PANELS AND SURFACE GEOMETRY ==========%

function [XB, YB, XC, YC, phiR, betaR, S, numPan] = LoadPanels(user_input, chord, alpha)
%user_input: feeds the switch case variable to determine whether a
%pre-defined profile is used or a custom one

switch user_input
    case 'Cosine'
        N = 71; %select number of boundary points

        theta_U = linspace(pi,0,N)'; %analytical angle for upper surface cosine spacing
        theta_L = linspace(0,pi,N)'; %analytical angle for lower surface cosine spacing

        xU = (chord/2)*(1 - cos(theta_U)); %upper surface x-coords
        xL = (chord/2)*(1 - cos(theta_L)); %lower surface y-coords
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
betaD = deltaD-alpha; %angle between freestream velocity vector and outward facing unit normal
betaD(betaD>360) = betaD(betaD>360)-360; %make angles less than 360

phiR = phiD.*(pi/180); %angle phi converted to radians
betaR = betaD.*(pi/180); %angle beta converted to radians



end