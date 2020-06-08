function kc = getKcFromHalfFilledCavity(a, b, c, h, epr, sigmaCond, m, n, axisLimits, resolution)
%getKcFromHalfFilledCavity - Plots (and solve) the trascendental equation you get
% from solving the half-filled cavity shown in TAP_prism, Chatterjee, etc.
% Geometry is a x b x c, where the cavity is filled with dielectric from 0 to h in the
% z direction.
% The examples in Chatterjee's paper are a = 0.01, b = 0.001 (to discard solutions
% x1x), c = 0.01 and h = 0.005.
% If guess is not provided, a plot of the two trascendental equations is shown.
% It is the equation from (8-115) in Balanis' Advanced Engineering Electromagnetics.
% resolution is the argument to the digits function which is used to obtain an arbitrary accuracy.
% Author: aamorm. Please write to aamor89@gmail.com for concerns.

    % Initialize variable.
    kc = 0;

    syms kcsquared kcsym

    % Lossy material
    if (sigmaCond > eps)
        lhsEquation = sqrt(kcsym^2-(m*pi/a)^2-(n*pi/b)^2)*cot(sqrt(kcsym^2-(m*pi/a)^2-(n*pi/b)^2)*(c-h));
        eprMod = epr - 1i*sigmaCond*4*pi*1e-7*299792458/kcsym;
        rhsEquation = -sqrt(eprMod*kcsym^2-(m*pi/a)^2-(n*pi/b)^2)*cot(sqrt(eprMod*kcsym^2-(m*pi/a)^2-(n*pi/b)^2)*h);
    else
        axisLimits = axisLimits.^2;
        lhsEquation = sqrt(kcsquared-(m*pi/a)^2-(n*pi/b)^2)*cot(sqrt(kcsquared-(m*pi/a)^2-(n*pi/b)^2)*(c-h));
        eprMod = epr;
        rhsEquation = -sqrt(eprMod*kcsquared-(m*pi/a)^2-(n*pi/b)^2)*cot(sqrt(eprMod*kcsquared-(m*pi/a)^2-(n*pi/b)^2)*h);
    end

    if (nargin > 9)

        digitsOld = digits(resolution);
        equation = lhsEquation-rhsEquation == 0;
        if (sigmaCond > eps)
            kc = (vpasolve(equation, kcsym, axisLimits));
        else
            kc = (sqrt(vpasolve(equation, kcsquared, axisLimits)));
        end
        digits(digitsOld);
    else
        close all
        fplot(real(lhsEquation),axisLimits)
        hold on
        fplot(real(rhsEquation),axisLimits);
        fplot(imag(rhsEquation),axisLimits);
        xlabel('k_c')
        ylabel('Equation')
        legend('LHS','RHS','Location','best')
        set(gca, 'fontsize', 22)
    end