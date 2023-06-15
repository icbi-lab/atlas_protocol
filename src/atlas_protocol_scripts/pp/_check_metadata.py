import warnings
from typing import Dict, List, Optional

import pandas as pd


def validate_obs(
    adata_obs: pd.DataFrame,
    ref_meta_dict: Dict[str, Dict[str, List]],
    keys_to_ignore: Optional[List[str]] = None,
) -> None:
    """
    Validates the metadata information in the `.obs` attribute of an AnnData object or a Pandas DataFrame against a
    reference metadata dictionary.

    Parameters
    ----------
    adata_obs : pd.DataFrame
    Pandas DataFrame or `.obs` attribute of an AnnData object to be validated.
    ref_meta_dict : dict
    Reference metadata dictionary.
    keys_to_ignore : list[str], optional
    List of keys to ignore during validation. Defaults to None.

    Raises
    ------
    ValueError: If missing columns or invalid values are found in the metadata.

    Returns
    -------
    None

    Note:
    -----
    This function validates the metadata information in the `.obs` attribute of an AnnData object or a Pandas DataFrame
    against a reference metadata dictionary `ref_meta_dict`. The `ref_meta_dict` should be a dictionary where the keys
    represent the metadata columns to be validated, and the values are dictionaries with a 'values' key containing a list
    of allowed values for that metadata column. The function raises a `ValueError` if any missing columns or invalid values
    are found in the metadata columns of the AnnData object. The `keys_to_ignore` parameter can be used to specify a
    list of keys that should be ignored during validation. If missing columns are found, the error message will include
    the missing columns in the order of the input dictionary. If invalid values are found, the error message will
    include the invalid values for the corresponding column.
    """
    if keys_to_ignore is None:
        keys_to_ignore = []

    # Check if all expected keys are present in adata_obs.columns
    expected_cols = [k for k in ref_meta_dict.keys() if k not in keys_to_ignore]
    missing_cols = [c for c in expected_cols if c not in adata_obs.columns]
    if missing_cols:
        missing_cols_str = ", ".join(missing_col for missing_col in expected_cols if missing_col in missing_cols)
        raise ValueError(f"Missing columns: {missing_cols_str}")

    # Check if keys are present as columns and verify values if present (except keys_to_ignore)
    for key, value in ref_meta_dict.items():
        if key in keys_to_ignore:
            continue  # Skip keys to ignore

        # Raise a warning if the column contains only NaN values
        if adata_obs[key].dropna().empty:
            warnings.warn(f"Column '{key}' contains only missing values")

        if key not in adata_obs.columns:
            raise ValueError(f"Missing columns: {key}")

        # Verify type of corresponding column
        column_type = adata_obs[key].dtype
        expected_type = value.get("type", None)

        if expected_type is not None and column_type != expected_type:
            offending_value = adata_obs[key][adata_obs[key].apply(lambda x: type(x) != expected_type)].iloc[0]
            raise ValueError(
                f"Unexpected data type found in column '{key}'. Expected '{expected_type}', but found '{offending_value}'."
            )

        # Verify values in corresponding column
        allowed_values = value.get("values", None)
        column_values = adata_obs[key].unique()

        if "min" in value and "max" in value:
            min_value = value["min"]
            max_value = value["max"]
            invalid_values = [val for val in column_values if not (min_value <= val <= max_value)]
            invalid_values = [
                val for val in invalid_values if not (pd.isna(val) or val == "nan")
            ]  # Add check for non-NA values
            if invalid_values:
                raise ValueError(f"Invalid values found in column '{key}': {invalid_values}")
        elif allowed_values is not None:
            invalid_values = [val for val in column_values if val not in allowed_values]
            invalid_values = [
                val for val in invalid_values if not (pd.isna(val) or val == "nan")
            ]  # Add check for non-NA values
            if invalid_values:
                raise ValueError(f"Invalid values found in column '{key}': {invalid_values}")

        # Verify values in corresponding column
        allowed_values = value.get("values", None)
        column_values = adata_obs[key].unique()

        if "min" in value and "max" in value:
            min_value = value["min"]
            max_value = value["max"]
            invalid_values = [val for val in column_values if not (min_value <= val <= max_value)]
            if invalid_values:
                raise ValueError(f"Invalid values found in column '{key}': {invalid_values}")
        elif allowed_values is not None:
            invalid_values = [
                val for val in column_values if val not in allowed_values and not (pd.isna(val) or val == "nan")
            ]
            if invalid_values:
                raise ValueError(f"Invalid values found in column '{key}': {invalid_values}")
