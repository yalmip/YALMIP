function ei=eyev(n,i)
%EYEV Internal function to generate unit vector

% Author Johan Löfberg
% $Id: eyev.m,v 1.2 2004-07-02 08:17:30 johanl Exp $

if i>n
  disp('Error in eyev')
  return
  else
ei = zeros(n,1);ei(i)=1;
end;

