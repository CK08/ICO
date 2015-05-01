function InterpZ = NaNmatrix(xx,yy,GP,stress,Stype)
%
%
%

% xx e yy do mesmo tamanho
% Se GP e stress s?o celulas vai correr celula a celula caso contrario
%  corre uma vez

[linhas,colunas] = size(xx);

x = xx(:,1);
y = yy(1,:);

if Stype == 1
    type = 1; 
elseif Stype == 2
    type = 2;
elseif Stype == 3
    type = 3;
elseif Stype == 4
    VM = @(t1,t2,t3) (sqrt(t1.^2+t2.^2-t1.*t2+3.*t3.^2));
end

for k=1:max(size(GP))
    % Read GaussPoint and Stress cell
    gp = GP{k};
    tensao = stress{k};
    if Stype==4
        tensao = VM(tensao(:,1),tensao(:,2),tensao(:,3));
        type = 1;
    end
    
    % Initialise InterpZ cell
    InterpZ{k} = zeros(linhas,colunas);
    InterpZ{k}(:,:) = NaN;
    
    % Alloc values
    for i = 1:max(size(gp))
        px = gp(i,1);
        py = gp(i,2);
        value = tensao(i,type);
        
        lin = knnsearch(x  ,px);
        col = knnsearch(y.',py);
        InterpZ{k}(lin,col) = value;
    end
end

end