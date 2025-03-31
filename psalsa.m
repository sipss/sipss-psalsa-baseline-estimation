function [Xc,Z]= psalsa(X,varargin)
%  Baseline correction using psalsa with inputparser for flexible arguments

% This function is an alternative baseline correction method inspired by the
% asymetric least squares method (ALS) and the Adaptive Iteratively 
% Reweighted Penalized Least Squares (airPLS). Unlike these methods,  
% psalsa modifies the penalty weighting scheme by introducing an adaptive  
% exponential decay factor, allowing better handling of intense peaks.
% Some parts of the code are inspired by the MATLAB implementation described in reference (2) 
% and also implemented in an airPLS MATLAB function by Zhimin Zhang @ Central South University on Mar 30,2011. 

% Syntax:
% [Xc,Z]=psalsa(X)
% [Xc,Z]=psalsa(X,'Lambda',Lambda_val, 'order', order_val, 'p', p_val, 'k', kval, 'itermax', itermax_val)
% 
% Input 
%   X:row matrix of spectra or chromatogram (size m*n)
% Optional input: name-value pair parameter
%   'Lambda'  - adjustable smoothing parameter: The larger lambda, the smoother Z (default 1e6) 
%   'order'   - an integer indicating the order of the difference of penalties (default 2)
%   'p'       - asymmetry parameter for the start and end (default 0.4)
%   'k'       - parameter that controls the exponential decay of the weights (default 10% of the maximum value of X)
%   'itermax' - maximum iteration times (default 20)
%        
%  Output
%    Xc: the corrected spectra or chromatogram vector (size m*n)
%    Z: the fitted vector
%  Examples:
%         Xc=psalsa(X); %uses default values and returns the corrected spectra only
%         [Xc,Z]=psalsa(X,'Lambda', 1e7,'p', 0.01); %custom lambda and p, returns the corrected spectra and the estimated baseline
%         
%  Main Reference:
%       (1) Oller-Moreno, et al. Adaptive Asymmetric Least Squares baseline estimation for analytical instruments
%  Additional rerferences:
%       (2) Eilers, P. H. C., A perfect smoother. Analytical Chemistry 75 (14), 3631 (2003).
%       (3) Eilers, P. H. C., Baseline Correction with Asymmetric Least Squares Smoothing, http://www.science.uva.nl/~hboelens/publications/draftpub/Eilers_2005.pdf
%       (4) Gan, Feng, Ruan, Guihua, and Mo, Jinyuan, Baseline correction by improved iterative polynomial fitting with automatic threshold. Chemometrics and Intelligent Laboratory Systems 82 (1-2), 59 (2006).
% 
%  Authors: Antonio Pardo, Luis Fernández, Sergio Oller, Santiago Marco @ Universitat de Barcelona 2025

 parser = inputParser;

    % Add parameters with their validation
    addRequired(parser,'X', @(x) ismatrix(x) && size(x, 1) > 0 && size(x, 2) > 0); %This is needed

    %These are optionals
    addParameter(parser, 'Lambda', 1e6, @isscalar);
    addParameter(parser, 'order', 2, @(x) isscalar(x) && x >= 1);
    addParameter(parser, 'p', 0.4, @(x) isscalar(x) && x >= 0 && x <= 1);
    addParameter(parser, 'k', max(max(X)) / 10, @(x) isscalar(x) && x >= 0);
    addParameter(parser, 'itermax', 20, @(x) isscalar(x) && x > 0);

    % Parse the input arguments
    parse(parser, X, varargin{:});
    
    % Extract values
    lambda = parser.Results.Lambda;
    order = parser.Results.order;
    p = parser.Results.p;
    k = parser.Results.k;
    itermax = parser.Results.itermax;

[m,n]=size(X); %Dimensions of the input matrix X (m rows, n columns)
D = diff(speye(n), order); %Create the finite difference matrix of order "order"
DD = lambda*D'*D; %Compute the smoothness penalty matrix by multiplying by the lambda parameter

%Iterate over each row of matrix X (each spectrum or chromatogram)
for i=1:m
    w=ones(n,1); %Initialize weights to 1 (all points initially have equal importance)
    d = 2*ones(n,1)'; %Initialize residuals
    x=X(i,:); %Extract the current row from X (spectrum or chromatogram)

    for j=1:itermax    %Iterative process to adjust the baseline using the psalsa algorithm
        W=spdiags(w, 0, n, n); %Construct a sparse diagonal matrix with the current weight
        %Solve the system of equations to find the adjusted baseline z
        C = chol(W + DD); %Cholesky decomposition
        z = (C\(C'\(w .* x')))'; %Solve the system using Cholesky factorization
        d_geq_old = d >=0; %For the stop criterion
        d = x-z; %Compute the residue between the original signal and the adjusted baseline
        
        %Update weights for the next iteration
        w(d>=0) = p*exp(-d(d>=0)/k); %for positive residuals a pondered exponential is applied to p
        w(d<0)  = 1-p; %for negative residue a fixed value (1-p) is used
         if ( all((d>=0) == d_geq_old)) %Stop criterion. If the residuals do not change, stop
            break;
        end
    end
    if j > itermax %Check if the maximum number of iterations has been reached 
        disp('psalsa: iterations exhausted');
    end
    Z(i,:)=z; %Store the adjusted baseline for this row
end
Xc=X-Z; % Compute the corrected signal by subtracting the baseline from the original signal