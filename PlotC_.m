function PlotC_(XC, Cp, alphaD, numPan)
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
end