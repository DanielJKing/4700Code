%% Daniel King 100921117, Part 2
% Using FD code to calculate E field and potential with bottleneck. Called
% with bottleneck conductivity of 0.01, and conductivity of 1 elsewhere.
% Voltage BC of 0.8V on left.

function [output] = fields(Acond, Bcond)
global im fig fc map
ngeo = 2;
Max = 20;
nx=100;
ny=100;

pgeo = zeros(ngeo, 2);
pgeo(1,1) = 90;
pgeo(1,2) = 50;
pgeo(2,1) = 10;
pgeo(2,2) = 50;

rgeo = ones(ngeo,1)*Max;

dgeo = ones(ngeo,1)*Max;
cMap = zeros(nx,ny);


    
    for i = 1:nx
        for j = 1:ny
            cMap(j,i) = Acond;
            for p = 1:ngeo
                dx = abs(pgeo(p, 1) - i);
                dy = abs(pgeo(p, 2) - j);
                if (dx < dgeo(p)) && (dy < dgeo(p))
                    cMap(j, i) = Bcond;
                end
            end
        end
   end
G = zeros(nx*ny);
B = zeros(1,nx*ny);



for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;

        if i == 1
            G(n, :) = 0;
            G(n, n) = 1;
            B(n) = 0.8;
        elseif i == nx
            G(n, :) = 0;
            G(n, n) = 1;
        elseif j == 1
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nyp = j + 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            ryp = (cMap(i, j) + cMap(i, j + 1)) / 2.0;

            G(n, n) = -(rxm+rxp+ryp);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nyp) = ryp;

        elseif j ==  ny
            nxm = j + (i - 2) * ny;
            nxp = j + (i) * ny;
            nym = j - 1 + (i - 1) * ny;

            rxm = (cMap(i, j) + cMap(i - 1, j)) / 2.0;
            rxp = (cMap(i, j) + cMap(i + 1, j)) / 2.0;
            rym = (cMap(i, j) + cMap(i, j - 1)) / 2.0;

            G(n, n) = -(rxm + rxp + rym);
            G(n, nxm) = rxm;
            G(n, nxp) = rxp;
            G(n, nym) = rym;
        else
            nxm = j + (i-2)*ny;
            nxp = j + (i)*ny;
            nym = j-1 + (i-1)*ny;
            nyp = j+1 + (i-1)*ny;

            rxm = (cMap(i,j) + cMap(i-1,j))/2.0;
            rxp = (cMap(i,j) + cMap(i+1,j))/2.0;
            rym = (cMap(i,j) + cMap(i,j-1))/2.0;
            ryp = (cMap(i,j) + cMap(i,j+1))/2.0;

            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp;
        end

    end
end

V = G\B';

Vmap = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j + (i - 1) * ny;
        Vmap(i, j) = V(n);
    end
end

for i = 1:nx
    for j = 1:ny
        if i == 1
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i, j));
        elseif i == nx
            Ex(i, j) = (Vmap(i, j) - Vmap(i - 1, j));
        else
            Ex(i, j) = (Vmap(i + 1, j) - Vmap(i - 1, j)) * 0.5;
        end
        if j == 1
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j));
        elseif j == ny
            Ey(i, j) = (Vmap(i, j) - Vmap(i, j - 1));
        else
            Ey(i, j) = (Vmap(i, j + 1) - Vmap(i, j - 1)) * 0.5;
        end
    end
end

Ex = -2*10^7*Ex;
Ey = -10^7*Ey;
output = [Ex,Ey];
 
%     figure(1);
%     quiver(Ex', Ey');
%     title('E Field')
%     figure(2);
%     surf(cMap');
%     view(2);
%     title('Conductivity')
%     colorbar;
%     figure(3);
%     surf(Vmap')
%     title('Potential Map')
%     colorbar;

end