function MP_conn = BuildMPconn(varargin)
%
% Builds the 2-D connectivity matrix for multiple patches (ICO type file).
%
% Copyright (C) Diogo Cardoso, 2015
%
nargs = nargin;
if nargs<2
    error('Need at least 2 patches to build connectivity matrix');
end

index = 1;
MP_conn = [];

for eachpatch = 1:nargs
    patch = varargin{eachpatch};
    p1 = patch.order(1);p2 = patch.order(2);
    n1 = patch.number(1);n2 = patch.number(2);
    nele1 = sum(ceil(diff(patch.knots{1})));
    nele2 = sum(ceil(diff(patch.knots{2})));
    % Indexing the value of MP CP
    if eachpatch==1
        MPindex=0;
    else
        MPindex=max(max(MP_conn));
    end
    % Ciclo por elemento
    for eta = 1:nele2
        for xi = 1:nele1
            % Construir o 1o elemento
            if xi==1 && eta==1
                line = 1:(p1);
                matrix1 = repmat(line,[p2,1]);
                addi = 0;
                for pontos2 = 1:(p2)
                    matrix1(pontos2,:) = matrix1(pontos2,:)+addi;
                    addi = addi + n1;
                end
                %Correction to CP index in MP
                matrix1 = matrix1 + MPindex;
                %Add to MP-conn
                MP_conn(index,:) = reshape(matrix1.',[1,numel(matrix1)]);
                index = index + 1;
                % Incremental matrix (addM)
                addM = ones(size(matrix1)).*(p1-2);
                continue;
            elseif xi==1 && eta>1
                %inicial = min(min(matrix)) + (p1-1);
                %line = inicial:(inicial+p1-1);
                %matrix1 = repmat(line,[p2,1]);
                %addi = 0;
                %for pontos2 = 1:(p2)
                %    matrix1(pontos2,:) = matrix1(pontos2,:)+addi;
                %    addi = addi + n1;
                %end
                matrix1 = matrix1 + n1*(p2-2);
                
                %Correction to CP index in MP
                matrix1 = matrix1 + MPindex;
                %Add to MP-conn
                MP_conn(index,:) = reshape(matrix1.',[1,numel(matrix1)]);
                index = index + 1;
                % Incremental matrix (addM)
                addM = ones(size(matrix1)).*(p1-2);
                continue;
            end
            
            matrix = matrix1 + addM;
                          % Ter em ATENCAO este valor (p1-2) - verificar!!!  
            addM = addM + (p1-2);
            %Add to MP-conn
            MP_conn(index,:) = reshape(matrix.',[1,numel(matrix)]);
            index = index + 1;
        end
    end
end
% Add connectivity index to Output matrix
array = 1:size(MP_conn,1);
MP_conn = [array;MP_conn.'].';
end