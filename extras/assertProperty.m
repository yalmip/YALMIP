function properties = assertProperty(properties,checkfor,default)
if ~isfield(properties,checkfor)
    properties = setfield(properties,checkfor,default);
end