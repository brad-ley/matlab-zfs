nspins = 3;

Lx = cell(1, nspins); Ly = cell(1, nspins); Lz = cell(1, nspins);
sx = [0 1; 1 0];
sy = [0 -1i; 1i, 0];
sz = [1, 0; 0, -1];
unit = [1,0;0,1];

for n=1:nspins
    lxc=1;lyc=1;lzc=1;
    for k=1:nspins
        if k==n
            lxc=kron(lxc,sx);
            lyc=kron(lyc,sy);
            lzc=kron(lzc,sz);
        else
            lxc=kron(lxc,unit);
            lyc=kron(lyc,unit);
            lzc=kron(lzc,unit);
        end
    end
    Lx{n}=lxc; Ly{n}=lyc;Lz{n}=lzc;
end

cell2mat(Lx)
cell2mat(Ly)
cell2mat(Lz)