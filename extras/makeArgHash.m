function key = makeHash(obj)
    % Convert MATLAB data to hashable string
    bytes = getByteStreamFromArray(obj);
    md = java.security.MessageDigest.getInstance('MD5');
    md.update(bytes);
    hash = typecast(md.digest, 'uint8');
    key = sprintf('%02x', hash);
end