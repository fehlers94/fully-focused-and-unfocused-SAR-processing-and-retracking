function [out,RemainingVarargin] = parseVarargin(UserVarargin,DefinedVarargin,varargin)

% [out,RemainingVarargin] = parseVarargin(UserVarargin,DefinedVarargin,varargin)
%
% helper function for varargin
%
% argin
%       UserVarargin    ... cell array containing parameter / value pairs ('help' ... produces help)
%       DefinedVarargin ... cell array containing parameter / value pairs (defaults or keywords '$required' or '$normal')
%                       ... optional check of the validity of the values possible (mechanism of validateattributes())
%                           keywords are '$class' and '$attributes'
%                       ... optional description, keyword is '$description'
%       varargin
%                   'strict'   ... <false>, true ... miswritten or unknown
%                                  parameter check, if true: unused UserVargins will be returned in RemainingVarargin
%                   'SETargin' ... struct containing parameter / value pairs
%
% UserVarargin has highest priority, followed by SETargin.
% DefinedVarargin holds default values AND optional $class and $attributes (compare validateattributes())
%
% example for parsing varargin inside a user function / method
%       argin = parseVarargin(varargin,{'Mode',1,'Offset',10,'Label','Test','$class',{'char'}},'strict',true)
%
% GS 05-2010

%% help mode
HelpMode = false;
if length(UserVarargin)==1
    if strcmp(UserVarargin{:},'help'),
        HelpMode = true;
    end
end

%% '$normal' ... means: no key will be given, keys will be inserted automatically
NormalArgins = DefinedVarargin(find(strcmp('$normal',DefinedVarargin))-1);

% if ~isempty(NormalArgins) & isempty(UserVarargin) % ... fill in empty values 
%     x = cell(size(NormalArgins));
% else
%     x = UserVarargin(1:length(NormalArgins));
% end

tmp          = [NormalArgins;UserVarargin(1:length(NormalArgins))];
UserVarargin = [tmp(:)' UserVarargin(length(NormalArgins)+1:end)];

%% check inputs
if mod(length(UserVarargin),2) && length(UserVarargin)>=2,
    error('... UserVarargin has to follow ''key''+''value'' pairs!');
end
if mod(length(DefinedVarargin),2),
    error('... DefinedVarargin has to follow ''key''+''value'' pairs!');
end

%%
tmp = dbstack;
callingFunctionname = tmp(end).name;

%% ensure valid varargin for internal use
argin.strict   = false;
argin.SETargin = [];

for arginName = fieldnames(argin)'
    if ~isempty(strmatch(arginName{:},varargin(1:2:end))),
        argin.(arginName{:}) = varargin{strmatch(arginName{:},varargin(1:2:end))+1};
    end
end

%% DefinedVarargin -> DefinedStruct (transform DefinedVarargin to struct and add KeyWords if requiered)
ReservedKeyWordsAndDefaults = {'description',[],'class',[],'attributes',[]};

StartIndexes = [(find(~any(cell2mat(cellfun(@(x)cellfun(@(y)~isempty(y),x),cellfun(@(z)regexp(DefinedVarargin(1:2:end)',z),ReservedKeyWordsAndDefaults(1:2:end)','UniformOutput',false)','UniformOutput',false)),2))-1)'*2+1 length(DefinedVarargin)+1];
for Idx =  1:length(StartIndexes)-1
    tmp = DefinedVarargin(StartIndexes(Idx):StartIndexes(Idx+1)-1);
    
    DefinedStruct(Idx).name = tmp{1};
    DefinedStruct(Idx).value = tmp{2};
    
    for KeyWord = ReservedKeyWordsAndDefaults(1:2:end),
        KeyWordIdx = strmatch(['$' KeyWord{:}],tmp(1:2:end),'exact');
        if ~isempty(KeyWordIdx),
            DefinedStruct(Idx).(KeyWord{:}) = tmp{(KeyWordIdx-1)*2+2};
        else
            KeyWordIdx = strmatch(KeyWord{:},ReservedKeyWordsAndDefaults(1:2:end),'exact');
            DefinedStruct(Idx).(KeyWord{:}) = ReservedKeyWordsAndDefaults{(KeyWordIdx-1)*2+2};
        end
    end
end

%%  ...
if HelpMode,
    stack = dbstack;
    if length(stack)>2,
        %PrintHelp = @HTMLHelp;
        PrintHelp = @(x)printHTMLTable(x,'emphasizeFirstColumn',true,...
            'emphasizeFirstRow',true,'FontSize',12,'FontColor','blue','Title',stack(2).name);
    else
        PrintHelp = @ConsoleHelp;
    end
    
    tmp = arrayfun(@(x)struct2cell(DefinedStruct(x))',1:length(DefinedStruct),'UniformOutput',false)';
    tmp = reshape([fieldnames(DefinedStruct)' tmp{:}],length(tmp{1}),length(tmp)+1)';
    % fill in '-' @empty
    tmp(cellfun(@(x)isempty(x),tmp)) = {'-'};
    % convert cells to strings
    tmp(cellfun(@(x)iscell(x),tmp)) = cellfun(@(x)sprintf('%s ',x{:}),tmp(cellfun(@(x)iscell(x),tmp)),'UniformOutput',false);
    % convert anything else than strings to strings
    tmp(cellfun(@(x)~ischar(x),tmp)) = cellfun(@(x)mat2str(x,'class'),tmp(cellfun(@(x)~ischar(x),tmp)),'UniformOutput',false);
    
    PrintHelp(tmp);
    
    out = [];
    RemainingVarargin = {};
    
    return;
end

%% initialize out based on DefinedStruct
out               = cell2struct({DefinedStruct.value}',{DefinedStruct.name}');
RemainingVarargin = {};

% UserVarargin -> UserStruct (transform UserVarargin to struct)
UserStruct = cell2struct(UserVarargin(2:2:end)',UserVarargin(1:2:end)');

% keep already SET argin and overwrite defaults
if ~isempty(argin.SETargin),
    for Key = fieldnames(argin.SETargin)'
        if ~isempty(strmatch(Key{:},fieldnames(out),'exact')),
            out.(Key{:}) = argin.SETargin.(Key{:});
        end
    end
end

% built out and RemainingVarargin dependend on strict state
for Key = fieldnames(UserStruct)',
    IsValidKey = ~isempty(strmatch(Key{:},fieldnames(out),'exact'));
    switch true,
        case argin.strict && ~IsValidKey,
            error('... unknown varargin: ''%s''',Key{:});
        case ~argin.strict && ~IsValidKey,
            % collect unused varargins
            RemainingVarargin = [RemainingVarargin Key {UserStruct.(Key{:})}];
        otherwise,
            out.(Key{:}) = UserStruct.(Key{:});
    end
end

% check for required arguments
for Key = fieldnames(out)',
    if strcmp(out.(Key{:}),'$required'), error('%s ... required argin ''%s'' missing!',callingFunctionname,Key{:}); end
end

% validate types
for Idx = 1:length(DefinedStruct),
    Key = DefinedStruct(Idx).name;
    switch true,
        case ~isempty(DefinedStruct(Idx).class) && ~isempty(DefinedStruct(Idx).attributes),
            validateattributes(out.(Key),DefinedStruct(Idx).class,DefinedStruct(Idx).attributes,callingFunctionname,Key);
        case ~isempty(DefinedStruct(Idx).class) && isempty(DefinedStruct(Idx).attributes),
            validateattributes(out.(Key),DefinedStruct(Idx).class,{},callingFunctionname,Key);
        otherwise,
    end
end

%% helper
    function ConsoleHelp(help_cell)
        
        help_cell = help_cell';
        
        % coutn howmany colums are left
        col_cnt = size(help_cell,1);
        
        % get size of the bigest string in each column
        col_width = arrayfun(@(x)max(cellfun(@length,help_cell(x,:))),1:col_cnt);
        
        % function-handle to create a horizontal line
        horz_line = @()({arrayfun(@(x)fprintf('+-%s-',repmat('-',1,col_width(x))) ,1:col_cnt,'UniformOutput',false ); fprintf('+\n')});
        % function-handle to print a row with correct width
        print_row = @(row)({arrayfun(@(x)fprintf('| %-*s ',col_width(x),row{x}) ,1:col_cnt,'UniformOutput',false ); fprintf('|\n')});
        
        %Print Help
        horz_line();
        
        % Header
        print_row(help_cell(:,1));
        horz_line();
        
        % Rows
        arrayfun(@(y)print_row(help_cell(:,y)),2:size(help_cell,2),'UniformOutput',false);
        horz_line();
    end

%     function HTMLHelp
%         color = [[0.95 0.95 0.95];[0.95 0.95 0.98]] ;
%         M = [];
%         for Idx=1:size(help_cell,2)-1
%             if mod(Idx,2),
%                 M = [M;ones(1,size(help_cell,1)-1)];
%             else
%                 M = [M;zeros(1,size(help_cell,1)-1)];
%             end
%         end
%         
%         % http://tomheller.de/theholycymbal/html-farben.html
%         colors=repmat({'white';'bisque'},size(help_cell,2)+1,1);colors = colors(2:size(help_cell,2));
%         fprintf('<html>\n');
% 
%         fprintf('<table cellpadding="4" cellspacing="0" frame="box" style="font-family:arial; font-size:14px; color:darkred">\n');
%         fprintf('<tr><th colspan=5 style="font-family:comic sans ms; font-size:20px; color:bisque; background-color:darkred">available varargins for ''%s''</th></tr>\n',stack(2).name); % functionname
%         fprintf('<tr>%s</tr>\n',...
%             cell2mat(cellfun(@(th)sprintf('<th>%s</th>',th),help_cell(1:end,1)','UniformOutput',false))); % header
%         cellfun(@(tr,color)fprintf('<tr style="background-color:%s">%s</tr>\n',color,[tr{:}]),...
%             num2cell([...
%             cellfun(@(td)sprintf('<td style="min-width:200px"><b>%s</b></td>',td),help_cell(1,2:end)','UniformOutput',false) ... first column
%             cellfun(@(td)sprintf('<td>%s</td>',td),help_cell(2:end,2:end)','UniformOutput',false) ... content
%             ],2),colors,'UniformOutput',false);
%         fprintf('\n</table>\n');
% 
%         fprintf('</html>\n');
%     end
end
