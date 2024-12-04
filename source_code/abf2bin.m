function binFilename = abf2bin(abfFilename, format)

fileProperties = loadABF(abfFilename);
binFilename = writeBinary(fileProperties.sweep, [abfFilename(1:end-3) 'dat'], format);