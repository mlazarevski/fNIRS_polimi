function [ dod ] = BPF_NIRSToolbox( dod,hpf,lpf )

f = waitbar(0, sprintf('BPF filter:  %u / %u', 0, length(dod)));

for ind=1:length(dod)

    waitbar(ind/length(dod), f, sprintf('BPF filter:  %u / %u', ind, length(dod)));

    val = dod(ind).data;
    SD = nirs.util.probe2sd( dod(ind).probe );
    t = dod(ind).time;
    
%     dodFilt = hmrHighPassFiltConc(val, t, SD, hpf );
    dodFilt = hmrBandpassFilt(val, t, hpf, lpf);
    
    dod(ind).data = dodFilt;
end
close(f);

end
