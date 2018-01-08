%YALMIPDEMO Brief tutorial and examples.
%
% See also YALMIPTEST

disp('The examples here are obsolete.')
disp(sprintf('Please check out the <a href="http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Tutorials.Tutorials">on-line tutorials</a> on the YALMIP Wiki instead'));
return

% Check for paths
if ~(exist('socpex')==2)
    disp('You have to set the path to the demo library (...\yalmip\demos\)')
    return;
end


i = 1;
problems{i}.class = 0;
problems{i}.info = 'Getting started, the basics';
problems{i}.call = 'basicsex';i = i+1;

problems{i}.class = 1;
problems{i}.info = 'Linear and quadratic programming';
problems{i}.call = 'regressex';i = i+1;

problems{i}.class = 1;
problems{i}.info = 'Second order cone programming';
problems{i}.call = 'socpex';i = i+1;

problems{i}.class = 2;
problems{i}.info = 'Lyapunov stability (SDP)';
problems{i}.call = 'stabilityex';i = i+1;

problems{i}.class = 0;
problems{i}.info = 'Model predictive control (LP,QP,SDP)';
problems{i}.call = 'mpcex';i = i+1;

problems{i}.class = 2;
problems{i}.info = 'Determinant maximization (MAXDET)';
problems{i}.call = 'maxdetex';i = i+1;

problems{i}.class = 2;
problems{i}.info = 'Decay-rate estimation (SDP)';
problems{i}.call = 'decayex';i = i+1;

problems{i}.class = 0;
problems{i}.info = 'Mixed integer programming (MILP,MIQP,MICP)';
problems{i}.call = 'milpex';i = i+1;

problems{i}.class = 3;
problems{i}.info = 'Working with polynomial expressions';
problems{i}.call = 'nonlinex';i = i+1;

problems{i}.class = 3;
problems{i}.info = 'Working with nonlinear operators';
problems{i}.call = 'nonlinopex';i = i+1;

problems{i}.class = 3;
problems{i}.info = 'Nonlinear semidefinite programming using PENBMI (BMI)';
problems{i}.call = 'bmiex1';i = i+1;

problems{i}.class = 3;
problems{i}.info = 'Decay-rate estimation revisited with PENBMI (BMI)';
problems{i}.call = 'decaybmiex';i = i+1;

problems{i}.class = 3;
problems{i}.info = 'Simultaneous stabilization with PENBMI (BMI)';
problems{i}.call = 'simstabex';i = i+1;

problems{i}.class = 4;
problems{i}.info = 'Sum-of-squares decompositions';
problems{i}.call = 'sosex';i = i+1;

problems{i}.class = 4;
problems{i}.info = 'Polynomial programming using moment-relaxations';
problems{i}.call = 'momentex';i = i+1;

problems{i}.class = 4;
problems{i}.info = 'Global nonlinear programming';
problems{i}.call = 'globalex';i = i+1;

problems{i}.class = 5;
problems{i}.info = 'Multi-parametric programming';
problems{i}.call = 'mptex';i = i+1;

problems{i}.class = 5;
problems{i}.info = 'KYP problems (SDP)';
problems{i}.call = 'kypdex';i = i+1;

problems{i}.class = 5;
problems{i}.info = 'Posynomial geometric programming';
problems{i}.call = 'geometricex';i = i+1;

problems{i}.class = 5;
problems{i}.info = 'Complex-valued problems';
problems{i}.call = 'complexex';i = i+1;

problems{i}.class = 5;
problems{i}.info = 'Dual variables';
problems{i}.call = 'dualex';i = i+1;


while (1)
    clc
    echo off
    
    disp(' ')
    disp(' ')
    disp('                    YALMIP DEMO')
    disp(' ')
    oldclass = 0;
    for i = 1:length(problems)
    %    if problems{i}.class == oldclass
    %        fprintf('\n');
    %    end
        
        fprintf(['       %1.2d) ' problems{i}.info '\n'],i);
    %    oldclass = problems{i}.class;
    end
    
    disp('         ');
    disp('         0) quit')
    inp = input('Select demo: ');
    try
        if ~isempty(inp)
            switch(inp)
                case 0
                    return
                otherwise
                    if inp<=length(problems)
                        eval(problems{inp}.call);
                    end
            end
        end
    catch
        disp(lasterr)
        pause
    end
end

