def validate_obs(adata, ref_meta_dict, keys_to_ignore=None):
    """
    Validates the .obs slot of an AnnData object against a reference metadata dictionary.

    Args:
        adata (anndata.AnnData): AnnData object to be validated.
        ref_meta_dict (dict): Reference metadata dictionary.
        keys_to_ignore (list, optional): List of keys to ignore during validation. Defaults to None.

    Raises:
        AssertionError: If invalid values are found in the AnnData object.

    Returns:
        None

    Note:
        This function validates the metadata information in the `.obs` attribute of the AnnData object `adata` against a
        reference metadata dictionary `ref_meta_dict`. The `ref_meta_dict` should be a dictionary where the keys represent
        the metadata columns to be validated, and the values are dictionaries with a 'values' key containing a list of
        allowed values for that metadata column. The function raises an `AssertionError` if invalid values are found in the
        metadata columns of the AnnData object. The `keys_to_ignore` parameter can be used to specify a list of keys that
        should be ignored during validation.
    """
    if keys_to_ignore is None:
        keys_to_ignore = []

    # Check if keys are present as columns and verify values if present (except keys_to_ignore)
    for key, value in ref_meta_dict.items():
        if key in keys_to_ignore:
            continue  # Skip keys to ignore

        if key not in adata.obs.columns:
            continue  # Skip keys not present in adata.obs

        # Verify values in corresponding column
        allowed_values = value["values"]
        column_values = adata.obs[key].unique()
        invalid_values = [val for val in column_values if val not in allowed_values]
        assert not invalid_values, f"Invalid values found in column '{key}': {invalid_values}"


def find_missing_obs_columns(datasets, ref_columns):
    missing_columns = {}
    for key, dataset in datasets.items():
        missing_cols = set(ref_columns) - set(dataset.obs.columns)
        if missing_cols:
            missing_columns[key] = missing_cols
    return missing_columns


def search_dict(my_dict, columns, search=None):
    values = {}
    for column in columns:
        if column in my_dict:
            column_values = {}
            for key, value in my_dict[column].items():
                if not search or key in search:
                    if key in ["values", "description", "type"]:
                        column_values[key] = value
            values[column] = column_values
    return values
