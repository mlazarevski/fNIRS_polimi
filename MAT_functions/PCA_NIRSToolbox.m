function [ dod ] = PCA_NIRSToolbox( dod,ncomp )

if length(dod) == 1         % custom PCA: no need for waitbar

    tInc = dod.time;
    tInc = ones( 1,length(tInc) );
    SD = nirs.util.probe2sd( dod.probe );

    job = nirs.modules.Run_HOMER2();
    job.fcn = 'enPCAFilter';
    job.vars.nSV = ncomp;
    job.vars.tInc = tInc;

    dod = job.run( dod );

else                        % standard PCA
    
    f = waitbar(0, sprintf('PCA filter:  %u / %u', 0, length(dod)));

    for ind=1:length(dod)        
        tInc = dod(ind).time;
        tInc = ones( 1,length(tInc) );
        SD = nirs.util.probe2sd( dod(ind).probe );

        job = nirs.modules.Run_HOMER2();
        job.fcn = 'enPCAFilter';
        job.vars.nSV = ncomp;
        job.vars.tInc = tInc;

        waitbar(ind/length(dod), f, sprintf('PCA filter:  %u / %u', ind, length(dod)));

        dod(ind) = job.run( dod(ind) );
    end

    close(f);

end

end

