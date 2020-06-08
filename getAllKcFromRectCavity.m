function [kc, indices] = getAllKcFromRectCavity(a, b, c, numKc, resolution)
%getAllKcFromRectCavity - Returns a sorted list of kc and the respective mode in indices.
% Also, in indices, if 4th index is -1 is a TE, if 4th index is -2, is a TM mode.
% Arbitrary accuracy (specified by resolution) is provided thanks to the symbolic toolbox.
%
% Syntax: kc = getAllKcFromRectCavity(a,b,c,numKc)
% Author: aamorm. Please write to aamor89@gmail.com for concerns.

    % Default resolution
    if (nargin == 4)
        resolution = 16;
    end
    % Change the resolution of the symbolic package.
    digitsOld = digits(resolution);
    kc = sym(zeros(1, 2*(numKc+1)^3));
    indices = zeros(4, 2*(numKc+1)^3);

    counter = 1;

    for m = 0:numKc
        for n = 0:numKc
            if ((m==0) && (n==0))
                continue
            end
            for p = 0:numKc
                if (p>0)
                    % TE mode
                    kc(counter) = sqrt((m*pi/a)^2+(n*pi/b)^2+(p*pi/c)^2);
                    indices(:,counter) = [m, n, p, -1];
                    counter = counter + 1;
                end
                if ((m>0) && (n>0))
                    % TM mode
                    kc(counter) = sqrt((m*pi/a)^2+(n*pi/b)^2+(p*pi/c)^2);
                    indices(:,counter) = [m, n, p, -2];
                    counter = counter + 1;
                end
            end
        end
    end

    % To remove zeros.
    indices = indices(:,1:counter-1);
    kc = vpa(kc(:,1:counter-1),resolution);

    [kc, sorted_indices] = sort(kc);
    indices = indices(:,sorted_indices);

    kc = kc(1:numKc);
    indices = indices(:,1:numKc);

    % Restore previous resolution.
    digits(digitsOld)
