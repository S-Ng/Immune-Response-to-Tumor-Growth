function [a, b, sigma, f, h, w, m, k ,q, r, c, g, d, lambda, s, j] = getParameters(model, variation)
% outputs the correct variables for each mouse variation, nn, nl, ln, and ll
% human variation not input yet

if variation ~= 1 && variation ~= 2 && variation ~= 3 && variation ~= 4 % check input to avoid graphing without correct inputs
    error('variation must be 1-4. nn = 1, nl = 2, ln = 3, ll = 4');
end

if model == 0 % mouse
    a = 5.14*10^-1 ; % tumor growth rate
    b = 1.02*10^-9 ; % 1/b is tumor carrying capacity
    sigma = 1.3*10^4;
    f = 4.12*10^-2;
    h = 2.02*10^7;
    w = 1*10^-7;
    m = 2*10^-2;
    k = 2.02*10^7;
    q = 3.42*10^-10;
    r = 1.1*10^-7;
    if variation == 1 || variation == 2 % nn or nl
        c = 3.23*10^-7;
        g = 2.5*10^-2;
        if variation == 1 % nn
            d = 1.43;
            lambda = 5.8*10^-1;
            s = 2.73;
            j = 3.75*10^-2;
        else % nl
            d = 3.6;
            lambda = 4.6*10^-1;
            s = 1.61;
            j = 3.75*10^-2;
        end
    else % ln or ll
        c = 3.5*10^-6;
        g = 2*10^-1; %4*g(n)
        if variation == 3 % ln
            d = 3.51;
            lambda = 9*10^-1;
            s = 5.07;
            j = 1.13*10^-1; % 3*j(nn)
        else % ll
            d = 7.17;
            lambda = 7.5*10^-1;
            s = 4*10^-1;
            j = 3.0*10^-1; %8*j(nn)
        end
    end
elseif model == 1 % human
    
    a = 5.14*10^-1;
    b = 1.02*10^-9;
    c = 3.23*10^-1;
    if variation ~= 5 % patient 9
    d = 5.80;
    elseif variation ~= 6 % patient 10
    d = 4.23;
    end
    sigma = 1.3*10^4;
    if variation ~= 5 % patient 9
    lambda = 1.36;
    elseif variation ~= 6 % patient 10
    lambda = 1.43;
    end
    f = 4.12*10^-2;
    g = 2.5*10^-2;
    h = 2.02*10^7;
    j = 3.75*10^-2;
    k = 2.0*10^7;
    m = 2.00*10^-2;
    q = 3.42*10^-10;
    p = 1.00*10^-7;
    if variation ~= 5 % patient 9
    s = 2.5*10^-1;
    elseif variation ~= 6 % patient 10
    s = 3.6*10^-1;
    end
    r = 1.1*10^-7;
 end
end
