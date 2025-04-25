% Generate surrogate data with matching amplitude spectrum and
% signal distribution (Schreiber and Schmitz, 1996).
%
% H0: signal is a realisation of a linear Gaussian stochastic process, followed
%	by a static nonlinearity (observation function)
%
% set random seed using rand('seed',s)
%
% Usage: [Xs E] = generate_iAAFT (X, specflag);
%	specflag	exact amplitude spectrum (1, default), otherise signal distr
%	E		error plot

function [Xs,E] = generate_iAAFT (X, specflag)

%El algoritmo consiste en crear una surrogada con la misma distribución de
%probabilidad (mismos valores, pero mezclados) y mismo espectro. Para ellos
%hace dos pasos en cada iteración: 1- igualar espectro; 2- igualar
%distribución. Cada paso implica deshacer en parte lo anterior, pero en
%cada iteración hace que se parezcan más ambos factores. Cuando llegamos a
%un error aceptable, acabamos las iteraciones y obtendremos nuestra
%surrogada.


if (nargin<1)
	Xs = [];
	E = [];
	return;
end
if (nargin<2)
	specflag = 1;
end
if (isempty(X))
	return;
end

max_it = 500;
pp = length(X);

% Initial Conditions
rn = X(randperm(pp));
Xsorted = sort(X);	% Desired signal distribution
Yamp = abs(fft(X));	% Desired amplitude spectrum

%prev_err = 0;
E = zeros(1,max_it);
c = 1;
prev_err = 1e10;
err = prev_err - 1;
while (prev_err>err) && (c<max_it)

	% Match Amplitude Spec
	Yrn = fft(rn);
	Yang = angle(Yrn);
	sn = real(ifft(Yamp.*exp(sqrt(-1).*Yang)));

	% Scale to Original Signal Distribution
	[~, INDs] = sort(sn);
	rn(INDs) = Xsorted;

	% Eval Convergence
	prev_err = err;
	A2 = abs(Yrn);
	err = mean(abs(A2-Yamp));
	E(c) = err;

	c = c+1;
end
E = E(1:c-1);
if (specflag==1)
	Xs = sn;	% Exact Amplitude Spectrum
else
	Xs = rn;	% Exact Signal Distribution
end


