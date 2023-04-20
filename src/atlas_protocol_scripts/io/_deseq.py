from anndata import AnnData

# def fix_sample_ids(pb):
#     """Homogenize all sample ids."""
#     repl = {}
#     for k, v in dict(zip(pb.obs["condition"].index, "_" + pb.obs["condition"].values)).items():
#         repl[k] = k.replace(v, "")

#     return repl


def write_deseq_tables(pb: AnnData, samplesheet_path: str, counts_path: str):
    """Save AnnData as counts matrix and samplesheet as read by the DESeq2 script."""
    samplesheet = pb.obs.copy()
    samplesheet.reset_index(inplace=True)
    # sample_ids_repl = fix_sample_ids(pb)
    # bulk_df = pb.to_df().T.rename(columns=sample_ids_repl)
    bulk_df = pb.to_df().T
    bulk_df.index.name = "gene_id"
    samplesheet.to_csv(samplesheet_path, index=False)
    bulk_df.to_csv(counts_path)
