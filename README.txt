*****************************************    
6 steps towards a successful installation
*****************************************	


1)
Remove any old version of YALMIP

2)
unzip yalmip.zip. This should create the structure

/yalmip
/yalmip/@sdpvar
/yalmip/extras
/yalmip/solvers
/yalmip/modules
/yalmip/operators 

3)
Put the following paths in your MATLAB path 

/yalmip
/yalmip/extras
/yalmip/solvers
/yalmip/modules
/yalmip/modules/parametric
/yalmip/modules/moment
/yalmip/modules/global
/yalmip/modules/robust
/yalmip/modules/sos
/yalmip/operators 

Most easily done either via the gui or addpath(genpath('yourlocation/yalmip')) 
or as explained here https://yalmip.github.io/tutorial/installation/

4)	      
Make sure to have the desired solvers in your path. 

5)
Restart Matlab if you are installing a new version, or at least type "clear classes".

6)
Run yalmiptest.m and everything should work (as long as you have the necessary solvers).

7) Learn more at yalmip.github.io

Forum
https://github.com/yalmip/YALMIP/discussions
https://groups.google.com/forum/?fromgroups#!forum/yalmip

*****************************************

Comments and bug-reports are higly appreciated.

Johan Löfberg, Linköping University
johan.lofberg@liu.se