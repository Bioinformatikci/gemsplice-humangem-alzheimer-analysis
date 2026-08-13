function ensID =  refseq2ensembl(refseqID, ID_conversion)
ensID = '';
ix = find(strcmp(ID_conversion.RefSeqMRNA,refseqID));
if ~isempty(ix)
    ensID = ID_conversion.EnsemblTranscriptID{ix};
end

end