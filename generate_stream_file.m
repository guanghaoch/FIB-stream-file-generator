close all
clear all
clc

tic;

nm = 1e-9;
um = 1e-6;
ns = 1e-9;
us = 1e-6;


unit_dwt = 100*ns; % unit dwell time
ADC_steps = 2^16;

dwell_time = 5*us;
n_pass = 100;

d_circ = 0.9*nm;

Nx = 100; % number of periods
Ny = 100;
px = 100*nm; % period
py = 100*nm;

dirName = './';
fName = ['array_',num2str(Nx),'x',num2str(Ny),'_px',num2str(fix(px/nm)),...
    'nm_py',num2str(fix(py/nm)),'nm_',num2str(dwell_time/us),'usDwt_',num2str(n_pass),'passes'];

w_offset = 2.5*um;
h_offset = 2.5*um;
w_pattern = 15*um;
h_pattern = 15*um;
dx = w_pattern / ADC_steps;
dy = h_pattern / ADC_steps;

n_circ_x = d_circ/dx;
n_circ_y = d_circ/dy;
if n_circ_x<1/2 || n_circ_y<1/2
    n_circ_x = 1/2;
    n_circ_y = 1/2;
end

[nx,ny] = meshgrid(-n_circ_x/2:n_circ_x/2,-n_circ_y/2:n_circ_y/2);
ind = sqrt(nx.^2+ny.^2)<fix(n_circ_x/2);
pos_circ(:,1) = fix(nx(ind)+w_offset/dx);
pos_circ(:,2) = fix(ny(ind)+h_offset/dy);


pos_circ = raster2serp(pos_circ);

pos_circ = fix([w_offset/dx, h_offset/dy]);

Np = size(pos_circ,1);

pos = nan(Nx*Ny*Np*n_pass,3);
cc = 0;
for j = 0 : Ny-1
    j+1
    if mod(j,2)==0
        G0 = 0;
    else
        G0 = Nx-1;
    end
    
    for i = 0 : Nx-1
        offset_x = fix(px/dx) * abs(G0-i);
        offset_y = fix(py/dy) * j;
        new_pos_circ = pos_circ + [offset_x, offset_y];
        
        pos(cc*Np*n_pass+1:cc*Np*n_pass+Np*n_pass,:) =...
            repmat([dwell_time/unit_dwt*ones(Np,1), new_pos_circ],n_pass,1);
        cc = cc + 1;
        
%         pos(end+1,:) = [0, pos(end,2:3)];

    end
end


pos = pos(~all(isnan(pos),2),:);


%*** Prepare matrix format ************************************************
nrLines = size(pos,1); %Nr of coordinate lines
%*** Write header *********************************************************
if ~isfolder([fileparts(mfilename('fullpath')),'/Output_StreamFiles.nosync/'])
    mkdir([fileparts(mfilename('fullpath')),'/Output_StreamFiles.nosync/']);
end
fNameFull = [fileparts(mfilename('fullpath')),'/Output_StreamFiles.nosync/',fName,'.txt'];
fID = fopen(fNameFull,'w');
fprintf(fID, '%s\r\n%s\r\n%s\r\n', 's16', '1', num2str(nrLines));
fclose(fID);
writematrix(pos,fNameFull,"Delimiter","\t","Encoding","UTF-8","LineEnding","\r\n",...
    "WriteMode","append");
[~,fn,ext] = fileparts(fNameFull);
movefile(fNameFull,[fileparts(mfilename('fullpath')),'/Output_StreamFiles.nosync/',fName,'.str']);

toc;


function serp_pos = raster2serp(pos)
    
    [~,I] = sort(sqrt(sum(pos.^2,2)), 'ascend');
    pos = pos(I,:);
    pos = unique(pos, 'rows');
    pos_max = max(pos(:));
    n_pt = size(pos,1);
    serp_pos = nan(n_pt,2);
    k = 1;
    for r = linspace(0,pos_max,30)
        for theta = linspace(0,2*pi,50)
            c0 = r * [cos(theta), sin(theta)];
            distance = sqrt(sum((pos-c0).^2, 2));
            [~,I] = min(distance);
            serp_pos(k,:) = pos(I,:); k=k+1;
            pos(I,:) = [];
            if k-1 == n_pt
                return
            end
        end
    end
    
end

















