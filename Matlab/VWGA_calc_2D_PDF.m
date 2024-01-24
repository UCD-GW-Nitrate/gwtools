function PDF = VWGA_calc_2D_PDF(X, Y, s, varargin)
%% Calculates 2D propability density function

in = ~isnan(X) & ~isnan(Y) & ~isinf(X) & ~isinf(Y);
X = X(in);
Y = Y(in);

if nargin == 5
    xlm = varargin{1};
    ylm = varargin{2};
else
    xlm = [min(X) max(X)];
    ylm = [min(Y) max(Y)];
end

dx = linspace(xlm(1), xlm(2), 101);
dy = linspace(ylm(1), ylm(2), 101);
dx_n = linspace(0, 100, 101);
dy_n = linspace(0,100, 101);
% rescale
X_n = interp1(dx, dx_n, X);
Y_n = interp1(dy, dy_n, Y);

[Xg, Yg] = meshgrid(0:100);
%
% Radial basis function
nrm = @(x,y) sqrt((Xg - x).^2 + (Yg - y).^2);
phi = @(r,s) exp(-(r./s).^2);
V = zeros(size(Xg));
for ii = 1:length(X_n)
    xm = X_n(ii);
    ym = Y_n(ii);
    V = V + phi(nrm(xm, ym),s);
end

%
V = V./max(max(V));
PDF.X = dx;
PDF.Y = dy;
PDF.V = V;
[nx, ny] = ndgrid(PDF.X, PDF.Y);
PDF.F = griddedInterpolant(nx, ny, PDF.V','linear','nearest');



% clf
% surf(dx,dy,PDF,'edgecolor','none')
% hold on 
% plot(X,Y,'.')
% alpha(0.5)
% view(0,90)