function [b_mat] = BlkGenerate(Phi, varargin)
% generalize 1D block structure; run once and save the matrix;
% Input:    Phi, N*M sensitivity matrix;
%           type == 1, generate one dimensional block
%           type == 2, generate two dimensional block
%           Len, blcok size; Len = 4 for 1D default; window length = Len*2+1 for 2D default;
%           s11, 2D mapping.
% Output:   b_mat, embedding structure
% Author: Jinxi Xiang, 27/12/2019

% Default Parameter Values for Any Cases
type    = 1;
len     = 4;

if mod(length(varargin),2)==1
    error('Error! Pars should go by pairs\n');
else
   for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'type'
                type = varargin{i+1};
            case 'len'
                len = varargin{i+1}; 
            case 's11'
                s11 = varargin{i+1};     
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
    end 
end

if type==1
    M = size(Phi,2);
    C0 = sqrtm(eye(len));
    p = M - len + 1;
    b_mat = zeros(M,p*len);
    for i = 1 : p
        Ci = zeros(M,len);
        Ci( i:i+len-1, : ) = C0;
        b_mat(:,(i-1)*len+1:i*len) = Ci;
    end
end 

if type==2
    M = size(Phi, 2);
    N = size(s11,2);
    k = 1;
    for j = 1:N
        for i = 1:N
            if isnan(s11(i,j))==0
                s11(i,j) = k;
                k = k+1;
            end
        end
    end
   b_mat = [];
   
   [row, col] = find(~isnan(s11));
   
   blk_list = [];   % store the elements of each block
   k = 0;
   %%%%%%%%%% find block indices in the 2D map %%%%%%%%%%
    for i = 1:length(col)
        j = i;
        % the block should lay within [1 N]
        row_ul = row(j)-len; c1 = row_ul>0 & row_ul<N+1; 
        col_ul = col(i)-len; c2 = col_ul>0 & col_ul<N+1; 
        row_ur = row(j)-len; c3 = row_ur>0 & row_ur<N+1; 
        col_ur = col(i)+len; c4 = col_ur>0 & col_ur<N+1; 
        row_dl = row(j)+len; c5 = row_dl>0 & row_dl<N+1;
        col_dl = col(i)-len; c6 = col_dl>0 & col_dl<N+1;
        row_dr = row(j)+len; c7 = row_dr>0 & row_dr<N+1;
        col_dr = col(i)+len; c8 = col_dr>0 & col_dr<N+1;

        if c1*c2*c3*c4*c5*c6*c7*c8>0
            % no NaN should exist in the block
            c9 = isnan(s11(row_ul,col_ul))==0;
            c10 = isnan(s11(row_ur,col_ur))==0;
            c11 = isnan(s11(row_dl,col_dl))==0;
            c12 = isnan(s11(row_dr,col_dr))==0;
            if c9*c10*c11*c12>0
                k = k+1;
                s11(row(j),col(i));
                row_map = s11(row_ul:row_dr,col_ul:col_dr);
                blk_list(:,k) = row_map(:);                
            end
        end
        
    end
    blk_num = size(blk_list,2);     % total number of block
    blk_size = size(blk_list,1);    % block size 
    %%%%%%%%%% extract each block and make row manipulation %%%%%%%%%%
    M = size(Phi,2);
    b_mat = zeros(M,blk_num*blk_size);
    for i = 1 : blk_num
        bem = zeros(M,blk_num*blk_size);
        bem(blk_list(:,i),:) = [zeros(blk_size,(i-1)*blk_size)...
            eye(blk_size) zeros(blk_size,(blk_num-i)*blk_size)]; % extract block
        b_mat = b_mat + bem;
    end
end
