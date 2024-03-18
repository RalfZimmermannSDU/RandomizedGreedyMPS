function [ROM] = prepare_PODDEIMROM(FOM, ROM, rDEIM)
%
%ROM struct gets enhanced with DEIM bases and magic point lists
%
[Qa,Ea,Ra] = svd(FOM.a_nonlin(:,1:FOM.n_tsteps));
ROM.Qa     = Qa(:,1:rDEIM);  % DEIM basis for a-terms
[Qb,Eb,Rb] = svd(FOM.b_nonlin(:,1:FOM.n_tsteps));
ROM.Qb     = Qb(:,1:rDEIM);  % DEIM basis for b-terms
%figure;
%semilogy(1:N, diag(Eb), 1:N, diag(Ea));

% perform DEIM, output is point lists Pa, Pb
ROM.Pa0 = deim(ROM.Qa);
ROM.Pb0 = deim(ROM.Qb);

ROM.dimDEIM = rDEIM;
end