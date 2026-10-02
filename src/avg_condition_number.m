function avgCond = avg_condition_number(Nt, Nr, rhoRange, numTrials)
%AVG_CONDITION_NUMBER  Mean condition number of the channel matrix vs. correlation.
%
%   avgCond = AVG_CONDITION_NUMBER(Nt, Nr, rhoRange, numTrials) draws
%   numTrials channel realisations for every correlation factor in rhoRange
%   and returns the average COND(H) for each one.
%
%   A larger condition number means a more ill-conditioned channel, i.e.
%   a harder detection problem (especially for the ZF detector).
%
%   See also CORRELATION_FACTORS, GENERATE_CHANNEL.

    avgCond = zeros(size(rhoRange));
    for i = 1:numel(rhoRange)
        C = correlation_factors(Nt, Nr, rhoRange(i));
        c = zeros(1, numTrials);
        for t = 1:numTrials
            c(t) = cond(generate_channel(C));
        end
        avgCond(i) = mean(c);
    end
end
