function ei=eyev(n,i)
%EYEV Internal function to generate unit vector

if i>n
  disp('Error in eyev')
  return
  else
ei = zeros(n,1);ei(i)=1;
end;

