function dimResults = dimAnalysis(pattern)
%% get list of files to analyse
folderlist = dir(pattern);

saveTo = 'dimAnalysis.csv';

%% create Cell to store extracted Information
n = length(folderlist);
dimResults = cell(n+1,11);
dimResults{1,1} = 'Foldername';
dimResults{1,2} = 'dk';
dimResults{1,3} = 'D';
dimResults{1,4} = 'dopt';
dimResults{1,5} = 'L';
dimResults{1,6} = 'maxTrustsite';
dimResults{1,7} = 'loops';
dimResults{1,8} = 'Time';
dimResults{1,9} = 'D_max';
dimResults{1,10}= 'increasedk';
dimResults{1,11}= 'reDim';

%% Start extracting information from sweeps
for k=2:(n+1)
    filename = strcat(folderlist(k-1).name,'/results.mat');
    if exist(filename)
        load(filename);
        dimResults{k,1} = folderlist(k-1).name;
        [dimResults{k,2},dimResults{k,3},dimResults{k,4},dimResults{k,5}] = strread(dimResults{k,1}, '%*u-%*u-%*[SspinBos]-%*[dk]%u%*[D]%u%*[dopt]%u%*[L]%u');
        dimResults{k,6} = max(para.trustsite);
        dimResults{k,7} = para.loop;
        if max(cellfun(@(x) strcmp(x,'time'), fieldnames(results)))
            dimResults{k,8} = results.time;
        end
        dimResults{k,9} = max([results.D{:}]);
        dimResults{k,10} = para.increasedk;
        dimResults{k,11} = mat2str(find(cellfun(@isempty,results.D)==0));
    end
end
end
%% write into file
% fid = fopen(saveTo, 'w') ;
% fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\r\n', dimResults{1,:}) ;
% fprintf(fid, '%s\t%u\t%u\t%u\t%u\t%u\r\n', dimResults{2:end,1:end-1}) ;
% fclose(fid) ;
