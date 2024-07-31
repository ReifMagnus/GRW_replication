function r = BetaCoef(x,mosd)
al  = x(1);
bet = x(2);

mode = mosd(1);
sd   = mosd(2);

r(1) = mode-(al-1)/(al+bet-2);
r(2) = sd-sqrt(al*bet/((al+bet)^2*(al+bet+1)));
end