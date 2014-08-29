function [sol,info] = solvebilevel(varargin)
%SOLVEBILEVEL Simple global bilevel solver
%
%   min        Objective(Pouter)(x,y)
%   subject to Constraints(Pouter)(x,y)>0
%              y = arg min Objective(Pinner)(x,y)
%              subject to Constraints(Pinner)(x,y)>0
%
%   [DIAGNOSTIC,INFO] = SOLVEBILEVEL(Pouter, Pinner, y, options)
%
%   diagnostic : Struct with standard YALMIP diagnostics
%   info       : Bilevel solver specific information
%
%   Input
%      Pouter   : Outer problem 
%      Pinner   : Inner problem
%      y        : Inner variables
%      options  : solver options from SDPSETTINGS.
%
%   The behaviour of the bilevel solver can be controlled
%   using the field 'bilevel' in SDPSETTINGS
%
%      bilevel.outersolver : Solver for outer problems with inner KKT removed
%      bilevel.innersolver : Solver for inner problem
%      bilevel.rootcut     : Number of cuts (based on complementary
%                            constraints) added in root (experimental)
%      bilevel.relgaptol   : Termination tolerance
%      bilevel.compslacktol: Tolerance for accepting complementary slackness
%      bilevel.feastol     : Tolerance for feasibility in outer problem
%
%
%   See also SDPVAR, SDPSETTINGS, SOLVESDP

if isa(varargin{2},'optimizer')
    Pout = varargin{1};
    s = struct(varargin{2});
    if nargin < 3
        options = sdpsettings;
    else
        options = varargin{3};
    end     
    [sol,info] = solvebilevel(Pout.Constraints,Pout.Objective,s.F,s.h,s.output.expression,options)    
else
    Pout = varargin{1};
    Pinn = varargin{2};
    y = varargin{3};
    
    if nargin < 4
        options = sdpsettings;
    else
        options = varargin{4};
    end
    
    [sol,info] = solvebilevel(Pout.Constraints,Pout.Objective,Pinn.Constraints,Pinn.Objective,y,options)
end
