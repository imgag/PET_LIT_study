* Encoding: UTF-8.
data list free /rho n.
begin data
0,72 64
0,54 64
0,69 60
0,67 60
end data.
compute zaehler = rho*(n-2)**0.5.
compute nenner = (1-rho**2)**0.5.
execute.
compute t = zaehler/nenner.

execute.
compute p = (1-cdf.t(t,n-2))*2.
execute.

*https://datatab.de/tutorial/spearman-korrelation.
