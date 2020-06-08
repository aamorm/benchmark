function [kc, indices] = getAllKcFromCircCavity(a, h, number_kc, resolution)
%getAllKcFromCircCavity - Returns a sorted list of kc and the respective mode in indices.
% Also, in indices, if 4th index is -1 is a TE, if 4th index is -2, is a TM mode.
% Arbitrary accuracy (specified by resolution) is provided thanks to the symbolic toolbox.
%
% Syntax: kc = getAllKcFromCircCavity(a,b,c,number_kc)
% Author: aamorm. Please write to aamor89@gmail.com for concerns.
%

    if (nargin == 3)
        resolution = 16;
    end

    digitsOld = digits(resolution);

    limit_kc = 5;

    kc = sym(zeros(1, 2*limit_kc^3));
    indices = zeros(4, 2*limit_kc^3);

    counter = 1;

    for m = 0:limit_kc
        for n = 1:limit_kc
            for p = 0:limit_kc
                if (p>0)
                    % TE mode
                    kc(counter) = sqrt((debug.getZerosBesselFunction(m,n,true,resolution)/a)^2+(p*pi/h)^2);
                    indices(:,counter) = [m, n, p, -1];
                    counter = counter + 1;
                end
                % TM mode
                kc(counter) = sqrt((debug.getZerosBesselFunction(m,n,false,resolution)/a)^2+(p*pi/h)^2);
                indices(:,counter) = [m, n, p, -2];
                counter = counter + 1;
            end
        end
    end

    % To remove zeros.
    indices = indices(:,1:counter-1);
    kc = kc(:,1:counter-1);

    [kc, sorted_indices] = sort(kc);
    indices = indices(:,sorted_indices);

    kc = kc(1:number_kc);
    indices = indices(:,1:number_kc);

    digits(digitsOld)
