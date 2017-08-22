USE SNVBox_dev;
SELECT x.UID, x.Exon, SUM(x.UniprotSum) as UniprotSum, x.aaEnd-x.aaStart+1 as aaLen, SUM(x.UniprotSum) / (x.aaEnd-x.aaStart+1) as UniprotDensity
FROM ( 
    SELECT te.UID, te.Exon, te.aaStart, te.aaEnd, 
           BINDING+NP_BIND+METAL+DNA_BIND+ACT_SITE+SITE+LIPID+CARBOHYD+CA_BIND+MOD_RES+DISULFID+SE_CYS+PROPEP+SIGNALP+TRANSMEM+MOTIF+ZN_FING+REGIONS+PPI+RNABD+TF+LOC+MMBRBD+Chrom+PostModRec+PostModEnz as UniprotSum
    FROM Transcript_Exon te, Uniprot_Xref ux, Uniprot_features uf
    WHERE te.UID=ux.UID AND te.aaStart<=ux.Pos AND te.aaEnd>=ux.Pos AND uf.Pos=ux.UniprotPos AND uf.Acc=ux.Uniprot
) x
GROUP BY x.UID, x.Exon;
