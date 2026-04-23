function [ dod ] = HPF_NIRSToolbox( dod,hpf )

f = waitbar(0, sprintf('HPF filter:  %u / %u', 0, length(dod)));

for ind=1:length(dod)
    
    waitbar(ind/length(dod), f, sprintf('HPF filter:  %u / %u', ind, length(dod)));

    val = dod(ind).data;
    SD = nirs.util.probe2sd( dod(ind).probe );
    t = dod(ind).time;
    
    dodFilt = hmrHighPassFiltConc(val, t, SD, hpf );
    
    dod(ind).data = dodFilt;
end
close(f);

end

