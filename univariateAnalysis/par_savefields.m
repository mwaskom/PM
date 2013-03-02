function par_savefields(fnam,p)
% helperfunction for batchscript (writing seg output)

fn = fieldnames(p);
if numel(fn)==0, return; end;
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
save(fnam,fn{:});

% %----------------