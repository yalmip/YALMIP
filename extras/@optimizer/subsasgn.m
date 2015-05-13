function self = subsasgn(self,subs,data)

if isequal(subs(1).type,'.')
    if strcmp(subs(1).subs,'options')
        if length(subs)==1
            if isa(data,'struct') && isfield(data,'beeponproblem')
                self.model.options = data;
            else
                error('The P.options can only be assigned an sdpsettings structure')
            end
            varargout{1} = self;
        else
            s = ['self.model.' subs(1).subs '.'];
            for i = 2:length(subs)
                s = [s subs(i).subs];
            end
            try
                eval([s '=data;']);
            catch
                error('Field not found in options');
            end
        end
        
    else
        error('You can only maniplate P.options');
    end
    
else
    error('You can only maniplate P.options');
end
