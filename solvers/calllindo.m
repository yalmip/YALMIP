function output = calllindo(interfacedata)

switch interfacedata.solver.tag

    case {'lindo-NLP'}
        output = calllindo_nlp(interfacedata);
    case {'lindo-MIQP'}
        output = calllindo_miqp(interfacedata);
    otherwise
        error;
end