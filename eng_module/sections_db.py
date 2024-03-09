import pandas as pd
import pathlib

DB_PATH = pathlib.Path(__file__).parent / "section_data"


def pfi_sections():
    df = pd.read_csv(DB_PATH / "eu_pf_sections.csv")
    scaled_columns = {
        "A": 2,
        "Iy": 4,
        "Wel_y": 3,
        "Wpl_y": 3,
        "iy": 1,
        "Avz": 2,
        "Iz": 4,
        "Wel_z": 3,
        "Wpl_z": 3,
        "iz": 1,
        "Ss": 1,
        "It": 4,
        "Iw": 3 + 6,
    }
    for column_name, power in scaled_columns.items():
        df[column_name] = df[column_name] * 10**power
    return df


def aisc_w_sections(si_units: bool = True,  unscaled_si: bool = True):
    if not si_units:
        df = pd.read_csv(DB_PATH / "aisc_w_sections_us.csv")
        return df
    else:
        df = pd.read_csv(DB_PATH / "aisc_w_sections_si.csv")
        scaled_columns = {
            "Ix": 6,
            "Zx": 3,
            "Sx": 3,
            "Iy": 6,
            "Zy": 3,
            "Sy": 3,
            "Cw": 9,
            "J": 3,
        }
        if unscaled_si:
            for column_name, power in scaled_columns.items():
                df[column_name] = df[column_name] * 10**power
        return df
    

def aisc_hss_sections(si_units: bool = True, unscaled_si: bool = True):
    if not si_units:
        df = pd.read_csv(DB_PATH / "aisc_hss_sections_us.csv")
        return df
    else:
        df = pd.read_csv(DB_PATH / "aisc_hss_sections_si.csv")
        scaled_columns = {
            "Ix": 6,
            "Zx": 3,
            "Sx": 3,
            "Iy": 6,
            "Zy": 3,
            "Sy": 3,
            "J": 3,
        }
        if unscaled_si:
            for column_name, power in scaled_columns.items():
                df[column_name] = df[column_name] * 10**power
        return df


def section_filter(sections_df: pd.DataFrame, operator: str, **kwargs) -> pd.DataFrame:
    """
    Returns a new DataFrame representing 'sections_df' but filtered with 'operator'
    according to the given 'kwargs'
    sections_df: A DataFrame with records of structural sections
    operator: str, either {"ge", "le"}, greater-than-or-equal-to and less-than-or-equal-to
    'kwargs': The kwargs provided should correspond to column names in 'sections_df'
        (which should all be valid Python identifiers)
        e.g. if there is a column in 'sections_df' called Ix then the DataFrame
        could be filtered by calling the function as such:
            section_filter(sections_df, 'ge', Ix=400e6)
    """
    sub_df = sections_df.copy()
    for key, value in kwargs.items():
        if operator.lower() == "le":
            sub_df = sub_df.loc[sub_df[key] <= value]
        elif operator.lower() == "ge":
            sub_df = sub_df.loc[sub_df[key] >= value]
        if sub_df.empty:
            print(f"No records match all of the parameters: {kwargs}")
    return sub_df


def sort_by_weight(sections_df: pd.DataFrame, ascending: bool = True) -> pd.DataFrame:
    """
    Returns 'sections_df' sorted by weight
    """
    return sections_df.sort_values("W", ascending=ascending)
