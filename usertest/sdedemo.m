function sdedemo;
%sdeDemo tests SDE with a simple spiral consisting of 100 points.
%please accustom the code to your SDP solver (CSDP/SEDUMI)
%
%function sdeDemo;
%
%tt=linspace(0,4*pi,100);
%X=[tt.*sin(tt);tt.*cos(tt)];
%figure;
%subplot(2,1,1);
%scatter(X(1,:),X(2,:),'o',tt,'filled'); axis equal;
%title('Original');
%drawnow; 
%try
% [Y,D]=sdeCSDP(X);       % CSDP 
% fprintf('Seems like sdeCSDP is working!\n');
%catch
%  [Y,D]=sdeNT(X);           % SEDUMI 
%  fprintf('Seems like sdeSeDuMi is working!\n');
%end;
%subplot(2,1,2);
%scatter(Y(1,:),Y(2,:),'o',tt,'filled'); axis equal
%title('Reduced Dimensionality');

tt=linspace(0,4*pi,100);
X=[tt.*sin(tt);tt.*cos(tt)];
%X=X-repmat(mean(X,2),1,length(X));
figure;
pars.slack=1;

try
 fprintf('Computing distances...');
 Dis=distance(X);  
 fprintf('done\n');
catch
 error('ERROR! Are you sure distance.m is in the path?');
  
end;

subplot(2,1,1);
scatter(X(1,:),X(2,:),60,tt,'filled'); axis equal;
title('Original');
drawnow; 
try
 pars.solver=0;
 [Y,D]=sde(Dis,3,pars);       % CSDP 
 fprintf('\n\nCSDP is working!\n');
catch
  pars.solver=1;
  [Y,D]=sde(Dis,3,pars);           % SEDUMI 
  fprintf('\n\nCSDP does not seem to be installed correctly.\n');
  fprintf('SeDuMi is working!\n');
end;
subplot(2,1,2);
scatter(Y(1,:),Y(2,:),60,tt,'filled'); axis equal;
title('Reduced Dimensionality');
