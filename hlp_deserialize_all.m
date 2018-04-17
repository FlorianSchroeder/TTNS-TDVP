function s = hlp_deserialize_all(s)
%% Deserialize all variables in a struct or the base namespace
if nargin>0
    % input must be a struct
    Vars = fieldnames(s);
    for ii = 1:size(Vars)
        if strcmp(class(s.(Vars{ii})),'uint8')
            s.(Vars{ii}) = hlp_deserialize(s.(Vars{ii}));
        end
    end
else
    % with no input, go through everything in base matlab workspace
    Vars = evalin('base','whos');
    for ii = 1:size(Vars)
        if strcmp(Vars(ii).class,'uint8')
            evalin('base',sprintf('%s = hlp_deserialize(%s);',Vars(ii).name,Vars(ii).name));
        end
    end
end
end