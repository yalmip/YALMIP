% YALMIP
% Version 3 (R14) 23-Apr-2004
%   and                       - AND (overloaded)
%   binary                    - BINARY Constrains a set of SDPVAR-variables to be binary (0/1)
%   blkdiag                   - BLKDIAG (overloaded)
%   clean                     - CLEAN Remove terms with small coefficients
%   clearsdpvar               - CLEARSDPVAR Clear solution
%   cone                      - CONE Defines a second order cone constraint ||z||<x
%   conj                      - CONJ (overloaded)
%   ctranspose                - CTRANSPOSE (overloaded)
%   cut                       - CUT Defines a cut constraint
%   degreduce                 - DEGREDUCE Remove higher order terms
%   degree                    - DEGREE Polynomial degree
%   depends                   - DEPENDS Returns indicies to variables used in an SDPVAR object
%   det                       - DET (overloaded)
%   diag                      - DIAG (overloaded)
%   diff                      - DIFF (overloaded)
%   display                   - DISPLAY (overloaded)
%   double                    - DOUBLE Returns current numerical value 
%   eig                       - EIG Generate eigenvalue object (currently not used)
%   end                       - END (overloaded)
%   eq                        - EQ (overloaded)
%   exp                       - EXP (overloaded, currently not used) 
%   extractkyp                - EXTRACTKYP Returns (A,B,P,M) from KYP object
%   find                      - FIND (overloaded)
%   fliplr                    - FLIPLR (overloaded)
%   flipud                    - FLIPUD (overloaded)
%   ge                        - GE (overloaded)
%   generateAB                - GENERATEAB Internal function to generates linear equation system
%   getbase                   - GETBASE Internal function to extract all base matrices
%   getbasematrix             - GETBASEMATRIX Internal function to extract basematrix for variable IND
%   getbasematrixwithoutcheck - GETBASEMATRIXWITHOUTCHECK Internal function to extract basematrix for variable IND
%   getbasevectorwithoutcheck - GETBASEVECTORWITHOUTCHECK Internal function to extract basematrix for variable ind
%   gethackflag               - GETHACKFLAG Internal function to extract constraint type
%   getvariables              - GETVARIABLES Returns variable indicies to variables used in a SDPVAR object
%   gt                        - GT (overloaded)
%   hankel                    - HANKEL (overloaded)
%   homogenize                - HOMOGENIZE Homogenize polynomial
%   horzcat                   - HORZCAT (overloaded)
%   imag                      - IMAG (overloaded)
%   integer                   - INTEGER Constrains variables to be integer
%   is                        - IS Check property of variable.
%   isequal                   - ISEQUAL (overloaded)
%   ishermitian               - ISHERMITIAN Check if variable is Hermitian
%   isinteger                 - ISINTEGER Check if (part of) a variable is integer
%   islinear                  - ISLINEAR Check if variable is linear       
%   ismember                  - ISMEMBER (overloaded)
%   isreal                    - ISREAL (overloaded)
%   issymmetric               - ISSYMMETRIC Check if variable is symmetric
%   kron                      - KRON (overloaded)
%   kyp                       - KYP Create KYP matrix variable
%   le                        - LE (overloaded)
%   length                    - LENGTH (overloaded)
%   loadobj                   - LOADOBJ (overloaded)
%   lt                        - LT (overloaded)
%   minus                     - MINUS (overloaded)
%   mldivide                  - MLDIVIDE (overloaded)
%   mpower                    - MPOWER (overloaded)
%   mrdivide                  - MRDIVIDE (overloaded)
%   mtimes                    - MTIMES (overloaded)
%   nonlineartocone           - NONLINEARTOCONE Convert nonlinear constraint to second order cone
%   norm                      - NORM Currently not used
%   not                       - NOT (overloaded)
%   numel                     - NUMEL (overloaded)
%   or                        - OR (overloaded)
%   plot                      - PLOT (overloaded)
%   plus                      - PLUS (overloaded)
%   power                     - POWER (overloaded)
%   prod                      - PROD (overloaded)
%   quaddecomp                - QUADDECOMP Internal function to decompose quadratic expression
%   rank                      - RANK (overloaded)
%   rcone                     - RCONE Defines a rotated second order cone constraint ||z||^2<2xy, x+y>0
%   rdivide                   - RDIVIDE (overloaded)
%   real                      - REAL (overloaded)
%   relaxdouble               - RELAXDOUBLE Return numerial value treating nonlinear variables as independant
%   replace                   - REPLACE Substitutes variables
%   repmat                    - REPMAT (overloaded)
%   reshape                   - RESHAPE (overloaded)
%   rot90                     - ROT90 (overloaded)
%   saveobj                   - SAVEOBJ (overloaded)
%   sdp2mat                   - SDP2MAT Obsolete, use double
%   sdpvar                    - SDVPAR Create symbolic decision variable
%   sdpvarfun                 - SDPVARFUN Applies operator on matrix variable
%   see                       - SEE Displays internal structure of variable
%   set                       - SET Defines a constraint (the feasible set)
%   setsdpvar                 - SETSDPVAR Assigns a numerical value to an sdpvar
%   setsos                    - SETSOS Internal function
%   size                      - SIZE (overloaded)
%   sos                       - SOS Declare sum-of-squares structure
%   sosd                      - SOSD Returns sum-of-squares decomposition
%   sparse                    - SPARSE (overloaded)
%   spy                       - SPY (overloaded)
%   sqrt                      - SQRT (overloaded) 
%   subsasgn                  - SUBASGN (overloaded)
%   subsref                   - SUBSREF (overloaded)
%   sum                       - SUM (overloaded)
%   times                     - TIMES (overloaded)
%   toeplitz                  - TOEPLITZ (overloaded)
%   trace                     - TRACE (overloaded)
%   transpose                 - TRANSPOSE (overloaded)
%   tril                      - TRIL (overloaded)
%   triu                      - TRIU (overloaded)
%   uminus                    - UMINUS (overloaded)
%   unique                    - UNIQUE (overloaded)
%   uplus                     - UPLUS (overloaded)
%   vertcat                   - VERTCAT (overloaded)
