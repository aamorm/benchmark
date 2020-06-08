function zeroValue = getZerosBesselFunction(m, n, derivativeFlag, resolution)
    %getZerosBesselFunction - Gets the zeros from bessel function (and his derivative) with arbitrary resolution
    % using symbolic toolbox.
    % m: parameter of the bessel function. Only to m = 11 is supported
    % n: number of the root. Only 5 zeros are supported.
    % Author: aamorm. Please write to aamor89@gmail.com for concerns.

    if (nargin == 2)
        derivativeFlag = false;
    end

    if (nargin <= 3)
        resolution = 16;
    end

    if (m>11) || (n>5)
        error ('Only values up to m = 11 or n = 5 are allowed.')
    end

    syms z
    if (derivativeFlag)
        equation = besselj(m-1,z) - besselj(m+1,z) == 0;
    else
        equation = besselj(m, z) == 0;
    end
    digitsOld = digits(resolution);
    zeroValue = (vpasolve(equation, z, getGuessForBesselFunction (m, n, derivativeFlag)));
    digits(digitsOld);

end

% Obtained from Balanis' book Advanced Engineering Electromagnetics
function guess = getGuessForBesselFunction (m, n, derivativeFlag)

    if (derivativeFlag)

        switch m
        case (0)
            switch n
            case (1)
                guess = 3.8318;
            case (2)
                guess = 7.0156;
            case (3)
                guess = 10.1735;
            case (4)
                guess = 13.3237;
            case (5)
                guess = 16.4706;
            end
        case (1)
            switch n
            case (1)
                guess = 1.8412;
            case (2)
                guess = 5.3315;
            case (3)
                guess = 8.5363;
            case (4)
                guess = 11.7060;
            case (5)
                guess = 14.8636;
            end
        case (2)
            switch n
            case (1)
                guess = 3.0542;
            case (2)
                guess = 6.7062;
            case (3)
                guess = 9.9695;
            case (4)
                guess = 13.1704;
            case (5)
                guess = 16.3475;
            end
        case (3)
            switch n
            case (1)
                guess = 4.2012;
            case (2)
                guess = 8.0153;
            case (3)
                guess = 11.3459;
            case (4)
                guess = 14.5859;
            case (5)
                guess = 17.7888;
            end
        case (4)
            switch n
            case (1)
                guess = 5.3175;
            case (2)
                guess = 9.2824;
            case (3)
                guess = 12.6819;
            case (4)
                guess = 15.9641;
            case (5)
                guess = 19.1960;
            end
        case (5)
            switch n
            case (1)
                guess = 6.4155;
            case (2)
                guess = 10.5199;
            case (3)
                guess = 13.9872;
            case (4)
                guess = 17.3129;
            case (5)
                guess = 20.5755;
            end
        case (6)
            switch n
            case (1)
                guess = 7.5013;
            case (2)
                guess = 11.7349;
            case (3)
                guess = 15.2682;
            case (4)
                guess = 18.6375;
            case (5)
                guess = 21.9317;
            end
        case (7)
            switch n
            case (1)
                guess = 8.5777;
            case (2)
                guess = 12.9324;
            case (3)
                guess = 16.5294;
            case (4)
                guess = 19.9419;
            case (5)
                guess = 23.2681;
            end
        case (8)
            switch n
            case (1)
                guess = 9.6474;
            case (2)
                guess = 14.1155;
            case (3)
                guess = 17.7740;
            case (4)
                guess = 21.2291;
            case (5)
                guess = 24.5872;
            end
        case (9)
            switch n
            case (1)
                guess = 10.7114;
            case (2)
                guess = 15.2867;
            case (3)
                guess = 19.0046;
            case (4)
                guess = 22.5014;
            case (5)
                guess = 25.8913;
            end
        case (10)
            switch n
            case (1)
                guess = 11.7708;
            case (2)
                guess = 16.4479;
            case (3)
                guess = 20.2230;
            case (4)
                guess = 23.7607;
            case (5)
                guess = 27.1820;
            end
        case (11)
            switch n
            case (1)
                guess = 12.8264;
            case (2)
                guess = 17.6003;
            case (3)
                guess = 21.4309;
            case (4)
                guess = 25.0085;
            case (5)
                guess = 28.4609;
            end
        end

    else

        switch m
        case (0)
            switch n
            case (1)
                guess = 2.4049;
            case (2)
                guess = 5.5201;
            case (3)
                guess = 8.6537;
            case (4)
                guess = 11.7915;
            case (5)
                guess = 14.9309;
            end
        case (1)
            switch n
            case (1)
                guess = 3.8318;
            case (2)
                guess = 7.0156;
            case (3)
                guess = 10.1735;
            case (4)
                guess = 13.3237;
            case (5)
                guess = 16.4706;
            end
        case (2)
            switch n
            case (1)
                guess = 5.1357;
            case (2)
                guess = 8.4173;
            case (3)
                guess = 11.6199;
            case (4)
                guess = 14.7960;
            case (5)
                guess = 17.9598;
            end
        case (3)
            switch n
            case (1)
                guess = 6.3802;
            case (2)
                guess = 9.7610;
            case (3)
                guess = 13.0152;
            case (4)
                guess = 16.2235;
            case (5)
                guess = 19.4094;
            end
        case (4)
            switch n
            case (1)
                guess = 7.5884;
            case (2)
                guess = 11.0647;
            case (3)
                guess = 14.3726;
            case (4)
                guess = 17.6160;
            case (5)
                guess = 20.8269;
            end
        case (5)
            switch n
            case (1)
                guess = 8.7715;
            case (2)
                guess = 12.3386;
            case (3)
                guess = 15.7002;
            case (4)
                guess = 18.9801;
            case (5)
                guess = 22.2178;
            end
        case (6)
            switch n
            case (1)
                guess = 9.9361;
            case (2)
                guess = 13.5893;
            case (3)
                guess = 17.0038;
            case (4)
                guess = 20.3208;
            case (5)
                guess = 23.5861;
            end
        case (7)
            switch n
            case (1)
                guess = 11.0864;
            case (2)
                guess = 14.8213;
            case (3)
                guess = 18.2876;
            case (4)
                guess = 21.6415;
            case (5)
                guess = 24.9349;
            end
        case (8)
            switch n
            case (1)
                guess = 12.2251;
            case (2)
                guess = 16.0378;
            case (3)
                guess = 19.5545;
            case (4)
                guess = 22.9452;
            case (5)
                guess = 26.2668;
            end
        case (9)
            switch n
            case (1)
                guess = 13.3543;
            case (2)
                guess = 17.2412;
            case (3)
                guess = 20.8071;
            case (4)
                guess = 24.2339;
            case (5)
                guess = 27.5838;
            end
        case (10)
            switch n
            case (1)
                guess = 14.4755;
            case (2)
                guess = 18.4335;
            case (3)
                guess = 22.0470;
            case (4)
                guess = 25.5095;
            case (5)
                guess = 28.8874;
            end
        case (11)
            switch n
            case (1)
                guess = 15.5898;
            case (2)
                guess = 19.6160;
            case (3)
                guess = 23.2759;
            case (4)
                guess = 26.7733;
            case (5)
                guess = 30.1791;
            end
        end
    end
end
