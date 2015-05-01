


clear all;close all;clc;

%geometry = readAllOut('NURBSAllOut.txt', 'NURBSInput.txt');
geometry = readAllOut('NURBSIncOut.txt', 'NURBSInput.txt');
nframe = max(size(fieldnames(geometry)));
mov(1:nframe)= struct('cdata',[],'colormap',[]);
set(gca,'nextplot','replacechildren');

for k=1%:nframe
  %------------------------------------------------------------------------  
  % Create the figure
  %------------------------------------------------------------------------
  if k==1
      % Patch 1
      nrb1 = nrbmak(geometry.f0.coefs{1},geometry.f0.knots{1});
      nrbplotx(nrb1,[100 100]); hold on;
      % Patch 2
      nrb2 = nrbmak(geometry.f0.coefs{2},geometry.f0.knots{2});
      nrbplotx(nrb2,[100 100]); hold off;
  else
      % Patch 1
      nrb11 = nrbmak(geometry.(['f' num2str(k-1)]).coefs{1},geometry.(['f' num2str(k-1)]).knots{1});
      nrb12 = nrbmak(geometry.('f0').coefs{1},geometry.('f0').knots{1});
      nrbplotx_disp(nrb12,nrb11,[100 100],4); hold on;
      % Patch 2
      nrb21 = nrbmak(geometry.(['f' num2str(k-1)]).coefs{2},geometry.(['f' num2str(k-1)]).knots{2});
      nrb22 = nrbmak(geometry.('f0').coefs{2},geometry.('f0').knots{2});
      nrbplotx_disp(nrb22,nrb21,[100 100],4); hold off;
      
  end
  %------------------------------------------------------------------------
  colorbar;
  view(2);
  %axis([-0.5 5.5 -0.5 4.5])%PatchTest
  axis([-0.5 12.5 -0.5 6])%Ironing
  title(['frame ' num2str(k)]);
  mov(k)=getframe(gcf);
  pause(0.5)
end


%movie2avi(mov, '1moviename.avi', 'compression', 'None');


% nframe=10;
% x=rand(100,nframe);
% mov(1:nframe)= struct('cdata',[],'colormap',[]);
% set(gca,'nextplot','replacechildren')
% for k=1:nframe
%   plot(x(:,k))
%   mov(k)=getframe(gcf);
% end
% movie2avi(mov, '1moviename.avi', 'compression', 'None');

%% Represent inital conditions of the problem

clear all;close all; clc;


fid2 = fopen('NURBSInput.txt','r');
[ncpi, ordersi, knotsi, coefsi] = readNURBS_MP(fid2);
patch1 = nrbmak(coefsi{1},knotsi{1});
patch2 = nrbmak(coefsi{2},knotsi{2});

figure;
nrbplotx(patch1,[100 100]);hold on;
nrbplotx(patch2,[100 100]);
axis([-0.5 12.5 -0.5 6])%Ironing
view(2);


