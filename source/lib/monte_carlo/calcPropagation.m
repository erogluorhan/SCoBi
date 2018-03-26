% Mehmet Kurum - Nov. 2nd, 2007

function calcPropagation

%% GET GLOBAL DIRECTORIES
dir_afsa = SimulationFolders.getInstance.afsa;


%% GET GLOBAL PARAMETERS 
% Satellite Parameters
f_Hz = SatParams.getInstance.f_MHz * Constants.MHz2Hz ;
% Vegetation Parameters
dim_layers_m = VegParams.getInstance.dim_layers_m;
TYPKND = VegParams.getInstance.TYPKND;
dsty = VegParams.getInstance.dsty;
dim1_m = VegParams.getInstance.dim1_m;
dim2_m = VegParams.getInstance.dim2_m; 
dim3_m = VegParams.getInstance.dim3_m;
epsr = VegParams.getInstance.epsr;
parm1_deg = VegParams.getInstance.parm1_deg;
parm2_deg = VegParams.getInstance.parm2_deg;


%% CALCULATIONS
% Layer
[Nlayer, Ntype] = size(TYPKND) ;
NkindMax = max(max(TYPKND)) ;

% Angle
ANGDEG = 0 : 89 ;
% ANGDEG = 36 : 40 ;
Na = length(ANGDEG) ;

fXAmp = zeros(2, Na, NkindMax, Ntype, Nlayer) ;

for ii = 1 : Nlayer
    
    for jj = 1 : Ntype
        
        Nkind = TYPKND(ii, jj) ;
        
        for kk = 1 : Nkind
            
            if jj == VegParams.getInstance.TYPES.L % Elliptic Disk (Leaf)
                
                A = dim1_m(kk, jj, ii) ;
                B = dim2_m(kk, jj, ii) ;
                T = dim3_m(kk, jj, ii) ;
                EPS = epsr(kk, jj, ii);
                TH1 = parm1_deg(kk, jj, ii) ;
                TH2 = parm2_deg(kk, jj, ii) ;
                PROB_deg = [TH1, TH2] ;
                
                Object = Object_mode1(PROB_deg, A, B, T, EPS, 'L') ;
                
            else  % Circular Cylinder, i.e., dim1=dim2
                
                RAD = dim1_m(kk, jj, ii) ;
                LEN = dim3_m(kk, jj, ii) ;
                EPS = epsr(kk, jj, ii);
                TH1 = parm1_deg(kk, jj, ii) ;
                TH2 = parm2_deg(kk, jj, ii) ;
                PROB_deg = [TH1, TH2] ;
                
                Object = Object_mode1(PROB_deg, RAD, RAD, LEN, EPS, 'B') ;
                
            end % Cylinder or disk ?
                
            disp(strcat('Layer:', num2str(ii), '-Type:',...
                num2str(jj), '-Kind:', num2str(kk)))
            disp('calculating...')
            tic ;

            for aa = 1 : Na

                tho = ANGDEG(aa) * Constants.deg2rad ;
                pho = 0 ;
                % Average Forward Scattering Amplitude
                fXAmp(:, aa, kk, jj, ii) ...
                    = compute_avfscatamp(tho, pho, f_Hz, Object) ;
            end

            toc ;
            
        end % Nkind
        
    end % Ntype
    
end % Nlayer


%% SAVE
filename = 'ANGDEG' ;
writeVar(dir_afsa, filename, ANGDEG)

filename = 'fXAmp' ;
writeComplexVar(dir_afsa, filename, fXAmp)


%% PROPAGATION AND ATTENUATION
% Incremental Propagation Constant
dddKz = zeros(2, Na, NkindMax, Ntype, Nlayer) ;
ddKz = zeros(2, Na, Ntype, Nlayer) ;
dKz = zeros(2, Na, Nlayer) ;

% Attenuation
ATPIQS = zeros(2, Na, NkindMax, Ntype, Nlayer) ;
ATTPIQS = zeros(2, Na, Ntype, Nlayer) ;
ATTENPIQS = zeros(2, Na, Nlayer) ;
ATTENH = zeros(Na, 1) ;
ATTENV = zeros(Na, 1) ;

% Propagation Constants
for aa = 1 : Na
    
    tho = ANGDEG(aa) * Constants.deg2rad ;
    
    for ii = 1 : Nlayer
        
        Di = dim_layers_m(ii, 1) ;
        
        for jj = 1 : Ntype
            
            Nkind = TYPKND(ii, jj) ;
            
            for kk = 1 : Nkind
                
                RHO = dsty(kk, jj, ii) ;
                fXAmp0 = fXAmp(:, aa, kk, jj, ii) ;
               
                [~, dddKz0] = compute_dKzn(tho, f_Hz, RHO, fXAmp0) ;
                
                dddKz(:, aa, kk, jj, ii) = dddKz0 ;
                % +++++++++++++++++++++
                % Determine incremental propagation constant for each
                % scatter type
                ddKz(:, aa, jj, ii) = ddKz(:, aa, jj, ii) + dddKz0 ;
                
                % 20 log10(exp(imag(dkz)*d))
                ATPIQS(:, aa, kk, jj, ii) = 20 * 0.4343 * imag(dddKz0) * Di;
           
            end % Nkind
            
            ddKz0 = ddKz(:, aa, jj, ii) ;
            % +++++++++++++++++++++
            % Calculate attenuation due to each scatter type in all layers
            % ATTPIQS(4, Ntype, Nlayer) : Incoming attenuation for P
            % poln and outgoing attenuation for Q poln
            ATTPIQS(:, aa, jj, ii) = 20 * 0.4343 * imag(ddKz0) * Di;
            
            % +++++++++++++++++++++
            % Determine incremental propagation contant for each layer
            % Sum over the scatter types assigned to each layer
            dKz(:, aa, ii) = dKz(:, aa, ii) + ddKz0 ;
            
        end % Ntype
        
        dKz0 = dKz(:, aa, ii) ;
        % ++++++++++++++++++++++++
        % Calculate total incoming and outgoing attenuation for each layer
        ATTENPIQS(:, aa, ii) = 20 * 0.4343 * imag(dKz0) * Di;
        
        % ++++++++++++++++++++++++
        % Calculate overall attenuation
        ATTENH(aa, 1) = ATTENH(aa, 1) + ATTENPIQS(1, aa, ii) ;
        ATTENV(aa, 1) = ATTENV(aa, 1) + ATTENPIQS(2, aa, ii) ;
        
    end % Nlayer
    
end


%% SAVE ALL
% Attenuation
filename = 'ATPIQS' ;
writeVar(dir_afsa, filename, ATPIQS) ;
filename = 'ATTPIQS' ;
writeVar(dir_afsa, filename, ATTPIQS) ;
filename = 'ATTENPIQS' ;
writeVar(dir_afsa, filename, ATTENPIQS) ;
filename = 'ATTENH' ;
writeVar(dir_afsa, filename, ATTENH) ;
filename = 'ATTENV' ;
writeVar(dir_afsa, filename, ATTENV) ;
% Propagation Constant
filename = 'dKz' ;
writeComplexVar(dir_afsa, filename, dKz) ;


end


%
% ----------------- compute_dKzn ------------------
%
function [kz, dKz] = compute_dKzn(tin, f_Hz, rho, fXAmp)
% COMPUTE_DKZN Calculates the contribution of each scatterer type to the attenuation

ko = 2 * pi * f_Hz / 3e+08;
kz = ko * cos(tin);
cfac = 2 * pi / kz;

% reshape into a n x 2 matrix as hh vv
% fXAmp = reshape(fXAmp.',[2, length(Sigma)]).';

dKz = cfac * rho * fXAmp;

end


%% compute_avfscatamp
% compute average forward scattering amplitudes

function avFScatAmp = compute_avfscatamp(tin, pin, f_Hz, Object)

switch lower(Object.category)
    case 'cylinder'
        
        scatFunc = @CYLINDER ;
        pdffun = Object.pdf ;
        
    case 'edisc'
        
        scatFunc = @EDISC ;
        pdffun = Object.pdf ;
        
    otherwise
        
        error('Unknown scatterer category!')
        
end

nth = length(tin) ;
nph = length(pin) ;

% hwaitbar = waitbar(0,'Please wait...','Name','Calculating average forward scattering amplitudes');

avFScatAmp = zeros(nth, nph, 4) ;

for m = 1 : nth % loop over ts
    
    for n = 1 : nph % loop over ps
        
        avFScatAmp(m, n, :) = compute_PhTh_average(pdffun, @compute_fscatamp, Object.theta1, ...
            Object.theta2, scatFunc, tin(m), pin(n), f_Hz, Object) ;
    end
    %         waitbar(m / nth, hwaitbar)
end
% close(hwaitbar)

% Select only hh and vv components, also eliminate singleton dimensions
avFScatAmp = squeeze(avFScatAmp(:, :, [1, 4])) ;

end


%% compute_fscatamp
%

function fScatAmp = compute_fscatamp(th, ph, scatFun, tin, pin, f_Hz, Object)

ts = pi - tin ;
ps = pi + pin ;

cobject = struct2cell(Object) ;
fScatAmp = feval(scatFun, tin, pin, ts, ps, th, ph, f_Hz, cobject{5 : end}) ;
fScatAmp = fScatAmp(:) ; % transform to vector from 2X2 matrix

end


%% compute_PhTh_average
%

function afun = compute_PhTh_average(pdf, fun, th1, th2, varargin)

tol = 1e-02 ;
if isequal(class(pdf), 'function_handle'),  pdf = func2str(pdf) ; end
if isequal(class(fun), 'function_handle'),  fun = func2str(fun) ; end

afun = qsimp(@compute_Th_average, 0, 2 * pi, tol, 'pdfxfun', pdf, fun, ...
    th1, th2, varargin{:}) ;
afun = afun / 2 / pi ;

end


%% compute_Th_average
%

function avf = compute_Th_average(phi, pxf, pdf, fun, th1, th2, varargin)

tol = 1e-02 ;
if th1 == th2 | pdf == 1 %#ok<OR2>
    
    avf = feval(fun, th1, phi, varargin{:}) ;
    
else
    
    avf = qsimp(pxf, th1, th2, tol, pdf, fun, th1, th2, phi, varargin{:}) ;
    
end

end


%% pdfxfun
%

function pxf = pdfxfun(x, pdf, fun, pdfarg1, pdfarg2, varargin)

pdf = feval(pdf, x, pdfarg1, pdfarg2) ;
fun = feval(fun, x, varargin{:}) ;
pxf = pdf .* fun ;


end


%% qsimp
%

function [q, n] = qsimp(fun, a, b, tol, varargin)

% TO-DO: Magic numbers?
nmax = 15 ;
q2old = -1e-30 ;
qold  = -1e-30 ;
q2 = 0 ;
it = 1 ;
for n = 1 : nmax
    
    [q2, it] =  trapzd(fun, a, b, q2, n, it, varargin{:}) ;
    q = (4 * q2 - q2old) / 3 ;
    if (abs(norm(q(:) - qold(:))) < abs(norm(qold(:))) * tol), return, end
    qold = q ;
    q2old = q2 ;
end


end


%% trapzd
%

function [s, it] = trapzd(fun, a, b, s, n, it, varargin)

if n == 1
    
    s = 0.5 * (b - a) * (feval(fun, a, varargin{:}) + feval(fun, b, varargin{:})) ;
    it = 1 ;
    
else
    
    del = (b - a) / it ;
    x = a + 0.5 * del ;
    sum = 0 ;
    
    for k = 1 : it
        
        sum = sum + feval(fun, x, varargin{:}) ;
        x = x + del ;
        
    end
    
    s = 0.5 * (s + sum * del) ;
    it = 2 * it ;
    
end

end


%% unifpdf (from Matlab statistics toolbox)
%UNIFPDF Uniform (continuous) probability density function (pdf).
%   Y = UNIFPDF(X,A,B) returns the continuous uniform pdf on the
%   interval [A,B] at the values in X. By default A = 0 and B = 1.
%
%   The size of Y is the common size of the input arguments. A scalar input
%   functions as a constant matrix of the same size as the other inputs.

%   Reference:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.34.

%   Copyright 1993-2000 The MathWorks, Inc.
%   $Revision: 2.9 $  $Date: 2000/05/26 18:53:54 $

function y = unifpdf(x, a, b)

if nargin < 1
    error('Requires at least one input argument.') ;
end

if nargin == 1
    a = 0 ;
    b = 1 ;
end

[errorcode x a b] = distchck(3, x, a, b) ;

if errorcode > 0
    error('Requires non-scalar arguments to match in size.') ;
end

% Initialize Y to zero.
y = zeros(size(x)) ;

k1 = find(a >= b) ;
if any(k1)
    tmp   = NaN ;
    y(k1) = tmp(ones(size(k1))) ;
end

k = find(x >= a & x <= b & a < b) ;
if any(k)
    y(k) = 1 ./ (b(k) - a(k)) ;
end

end


% ------------ distchck (from Matlab statistics toolbox)-----------
% checks the argument list for the probability functions.

%   B.A. Jones  1-22-93
%   Copyright 1993-2000 The MathWorks, Inc.
%   $Revision: 2.9 $  $Date: 2000/05/26 17:28:46 $

function [errorcode, varargout] = distchck(nparms, varargin)

errorcode = 0 ;
n = nargout - 1 ;
varargout = cell(1, n) ;

if nparms == 1
    varargout{1} = varargin{1} ;
    return;
end

% Get size of each input, check for scalars, copy to output
sz = cell(1, n) ;
isscalar = logical(zeros(1, n)); %#ok<LOGL>
for j = 1 : n
    s = size(varargin{j}) ;
    sz{j} = s ;
    isscalar(j) = (prod(s) == 1) ;
    varargout{j} = varargin{j} ;
end

% Done if all inputs are scalars.  Otherwise fetch their common size.
if (all(isscalar)), return ; end
t = sz(~isscalar) ;
size1 = t{1} ;

% Scalars receive this size.  Other arrays must have the proper size.
for j = 1 : n
    sizej = sz{j} ;
    if (isscalar(j))
        t = zeros(size1) ;
        t(:) = varargin{j} ;
        varargout{j} = t ;
    elseif (~isequal(sizej, size1))
        errorcode = 1 ;
        return;
    end
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%% MODE 1 FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% initialize_mode1
% initialize parameters for mode 1

function Object = Object_mode1(prob_deg, dim1_m, dim2_m, dim3_m, eps_r, filename)

theta1_deg = prob_deg(1, 1) ;
theta2_deg = prob_deg(1, 2) ;

theta1_rad = theta1_deg * Constants.deg2rad ;
theta2_rad = theta2_deg * Constants.deg2rad ;

EPSILON = eps_r ;
var1 = EPSILON ;

if filename(1, 1) == 'B' || filename(1, 1) == 'T' ...
        || filename(1, 1) == 'N'
    
    OBJECT_CATEGORY = 'cylinder' ;
    
else
    
    OBJECT_CATEGORY = 'edisc' ;
    
end

switch lower(OBJECT_CATEGORY)
    
    case 'cylinder'
        
        RADIUS = dim1_m ;
        LENGTH = dim3_m ;
        
        var2 = RADIUS ; % meters
        var3 = LENGTH ; % meters
        
        field2 = {'pdf' @unifpdf};
        extrafields = {'radius', 'length', 'epsilon'} ;
        orderindex = [1 5 2 4 3 6] ; % this will give fields the order: 'radius', 'length', 'epsilon'
        
    case 'edisc'
        
        SEMI_MAJORX = dim1_m ;
        SEMI_MINORX = dim2_m ;
        THICKNESS = dim3_m ;
        
        var2 = THICKNESS ; % meters
        var3 = SEMI_MAJORX ; % meters
        xmnr = SEMI_MINORX ; % meters
        
        extrafields = {'thickness', 'semiMajorx', 'semiMinorx', num2cell(xmnr), 'epsilon'};
        orderindex = [1 7 2 6 3 4 5 8] ; % this will give fields the order: 'thickness', 'semiMajorx', 'semiMinorx', 'epsilon'
        
        field2 = {'pdf' @unifpdf} ;
        
end

var1 = num2cell(var1) ;
var2 = num2cell(var2) ;
var3 = num2cell(var3) ;
extrafields = [extrafields{:}, var3, var2, var1] ;
extrafields = extrafields(orderindex) ; % reverse order
clear var*

% generate the object
Object = struct('category', OBJECT_CATEGORY, field2{:}, ...
    'theta1', theta1_rad, 'theta2', theta2_rad, extrafields{:}) ;

end
