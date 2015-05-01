%% Criar novo ficheiro Ironing apartir da analise 2-Ironing
clear all;close all; clc;


fid2 = fopen('NURBSInput.txt','r');
[ncpi, ordersi, knotsi, coefsi] = readNURBS_MP(fid2);
patch1 = nrbmak(coefsi{1},knotsi{1});
patch2 = nrbmak(coefsi{2},knotsi{2});
fclose(fid2);

%
%
% figure;
% nrbplotx(patch1,[100 100]);hold on;
% nrbplotx(patch2,[100 100]);
% plot3(squeeze(patch1.coefs(1,:,:)),squeeze(patch1.coefs(2,:,:)),ones(size(squeeze(patch1.coefs(2,:,:)))),'ob');
% plot3(squeeze(patch2.coefs(1,:,:)),squeeze(patch2.coefs(2,:,:)),ones(size(squeeze(patch2.coefs(2,:,:)))),'ob');
% axis([-0.5 12.5 -0.5 6])%Ironing
% view(2);
%
%

%
% Refinement
%
newnurbs = OrderElevation(patch1,1,1);
patch1 = OrderElevation(newnurbs,1,2);

newnurbs = OrderElevation(patch2,1,1);
patch2 = OrderElevation(newnurbs,1,2);
%
%
MP_conn = BuildMPconn(patch1,patch2);
%patch = patch1;
%xi = xi+1;
%eta = eta +1;
%
%
% figure;
% nrbplotx(patch1,[100 100]);hold on;
% nrbplotx(patch2,[100 100]);
% plot3(squeeze(patch1.coefs(1,:,:)),squeeze(patch1.coefs(2,:,:)),ones(size(squeeze(patch1.coefs(2,:,:)))),'ob');
% plot3(squeeze(patch2.coefs(1,:,:)),squeeze(patch2.coefs(2,:,:)),ones(size(squeeze(patch2.coefs(2,:,:)))),'ob');
% axis([-0.5 12.5 -0.5 6])%Ironing
% view(2);
%
%

%
%clear all;clc
% fid = fopen('Ironing_p2','r');
% data = ICOreadMP(fid);
% fclose(fid);
% 
% %tst=data.conn(1,:);
% coefs_list_total = [];
% pos = 0;
% for j = 1:data.patch1.number(2)
%     for i = 1:data.patch1.number(1)
%         pos = pos + 1;
%         coefs_list_total(pos,1) = data.patch1.coefs(1,i,j);
%         coefs_list_total(pos,2) = data.patch1.coefs(2,i,j);
%         coefs_list_total(pos,3) = data.patch1.coefs(3,i,j);
%     end
% end
% for j = 1:data.patch2.number(2)
%     for i = 1:data.patch2.number(1)
%         pos = pos + 1;
%         coefs_list_total(pos,1) = data.patch2.coefs(1,i,j);
%         coefs_list_total(pos,2) = data.patch2.coefs(2,i,j);
%         coefs_list_total(pos,3) = data.patch2.coefs(3,i,j);
%     end
% end

coefs_list_total = [];
pos = 0;
for j = 1:patch1.number(2)
    for i = 1:patch1.number(1)
        pos = pos + 1;
        coefs_list_total(pos,1) = patch1.coefs(1,i,j);
        coefs_list_total(pos,2) = patch1.coefs(2,i,j);
        coefs_list_total(pos,3) = patch1.coefs(3,i,j);
    end
end
for j = 1:patch2.number(2)
    for i = 1:patch2.number(1)
        pos = pos + 1;
        coefs_list_total(pos,1) = patch2.coefs(1,i,j);
        coefs_list_total(pos,2) = patch2.coefs(2,i,j);
        coefs_list_total(pos,3) = patch2.coefs(3,i,j);
    end
end

tst = MP_conn(:,2:end);
mov(1:size(tst,1))= struct('cdata',[],'colormap',[]);


for ciclo=1:size(tst,1)
    hold off;
    nrbplotx(patch1,[100 100]);hold on;
    nrbplotx(patch2,[100 100]);
    plot3(squeeze(patch1.coefs(1,:,:)),squeeze(patch1.coefs(2,:,:)),ones(size(squeeze(patch1.coefs(2,:,:)))),'ob');
    plot3(squeeze(patch2.coefs(1,:,:)),squeeze(patch2.coefs(2,:,:)),ones(size(squeeze(patch2.coefs(2,:,:)))),'ob');
    axis([-0.5 12.5 -0.5 6])%Ironing
    view(2);
    hold on;
    plot3(coefs_list_total(tst(ciclo,:),1),coefs_list_total(tst(ciclo,:),2),coefs_list_total(tst(ciclo,:),3),'dg','MarkerSize',12,'MarkerFaceColor','r');
    hold off;
    mov(ciclo)=getframe(gcf);
    pause(0.5);
end

% for ciclo=1:size(tst,1)
%     hold off;
%     nrbplotx(patch1,[100 100]);hold on;
%     nrbplotx(patch2,[100 100]);
%     plot3(squeeze(patch1.coefs(1,:,:)),squeeze(patch1.coefs(2,:,:)),ones(size(squeeze(patch1.coefs(2,:,:)))),'ob');
%     plot3(squeeze(patch2.coefs(1,:,:)),squeeze(patch2.coefs(2,:,:)),ones(size(squeeze(patch2.coefs(2,:,:)))),'ob');
%     axis([-0.5 12.5 -0.5 6])%Ironing
%     view(2);
%     hold on;
%     plot3(coefs_list_total(tst(ciclo,2:end),1),coefs_list_total(tst(ciclo,2:end),2),coefs_list_total(tst(ciclo,2:end),3),'dg','MarkerSize',12,'MarkerFaceColor','r');
%     hold off;
%     mov(ciclo)=getframe(gcf);
%     pause(0.5);
% end

%movie2avi(mov, '1moviename.avi', 'compression', 'None');