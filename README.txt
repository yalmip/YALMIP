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
/yalmip/demos
/yalmip/solvers
/yalmip/modules
/yalmip/operators 

3)
Put the following paths in your MATLAB path 

/yalmip
/yalmip/extras
/yalmip/demos
/yalmip/solvers
/yalmip/modules
/yalmip/modules/parametric
/yalmip/modules/moment
/yalmip/modules/global
/yalmip/modules/robust
/yalmip/modules/sos
/yalmip/operators 

Most easily done either via the gui or addpath(genpath('yourlocation/yalmip'))

4)	      
Make sure to have the desired solvers in your path. 

5)
Restart Matlab, or at least type "clear classes".

6)
Run yalmiptest.m and everything should work (as long as you have the 
necessary solvers).

Learn more at
http://users.isy.liu.se/johanl/yalmip

Forum
https://groups.google.com/forum/?fromgroups#!forum/yalmip

*****************************************

Comments and bug-reports are higly appreciated.

Johan Löfberg, Linköping University
johanl@isy.liu.se