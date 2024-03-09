from dataclasses import dataclass
from typing import Optional
import math
from math import sqrt
from eng_module.utils import str_to_float, read_csv_file
from eng_module import load_factors


@dataclass
class Column:
    h: float
    E: float
    A: float
    Ix: float
    Iy: float
    kx: float
    ky: float

    def critical_buckling_load(self, axis: str):
        """
        Returns the critical buckling load for the axis under consideration.

        axis: One of {'x', 'y'}
        """
        if axis.lower() not in ("x", "y"):
            raise ValueError("Axis argument must be one of {'x', 'y'}.")
        elif axis.lower() == "x":
            return euler_buckling_load(self.h, self.E, self.Ix, self.kx)
        else:
            return euler_buckling_load(self.h, self.E, self.Iy, self.ky)

    def radius_of_gyration(self, axis: str):
        """
        Returns the radius of gyration for the axis under consideration.

        axis: One of {'x', 'y'}
        """
        if axis.lower() not in ("x", "y"):
            raise ValueError("Axis argument must be one of {'x', 'y'}.")
        elif axis.lower() == "x":
            return radius_of_gyration(self.Ix, self.A)
        else:
            return radius_of_gyration(self.Iy, self.A)


@dataclass
class SteelColumn(Column):
    """
    A class to describe an axially loaded, doubly or singly symmetric steel section
    subject to axial compression.

    'fy': yield stress of the steel material
    'n', a factor required that describes the fabrication method of the section
        Acceptable values are either:
        1.34 for hot-rolled sections (e.g. HSS, W sections), DEFAULT
        2.24 for built-up plate sections (e.g. WT sections)
    'phi': material resistance factor (default = 0.9)
    """
    fy: float
    n: float = 1.34
    phi: float = 0.9
    tag: Optional[str] = None

    def factored_compressive_resistance(self) -> float:
        """
        Returns the factored axial capacity of self calculated using CSA S16:19
        in accordance with Cl. 13.3.1 about the governing axis.
        'n', a factor required that describes the fabrication method of the section
            Acceptable values are either:
                1.34 for hot-rolled sections (e.g. HSS, W sections)
                2.24 for built-up plate sections (e.g. WT sections)
        """
        F_e_x = self.critical_buckling_load("x") / self.A
        F_e_y = self.critical_buckling_load("y") / self.A
        F_e = min(F_e_x, F_e_y)
        buckling_factor = steel_buckling_factor(self.fy, F_e, self.n)
        crushing_load = gross_axial_resistance(self.fy, self.A, self.phi)
        return crushing_load * buckling_factor

    def factored_crushing_load(self) -> float:
        """
        Returns the factored axial capacity of self calculated using CSA S16:19
        in accordance with Cl. 13.3.1 about the governing axis.
        'n', a factor required that describes the fabrication method of the section
            Acceptable values are either:
                1.34 for hot-rolled sections (e.g. HSS, W sections)
                2.24 for built-up plate sections (e.g. WT sections)
        """
        return gross_axial_resistance(self.fy, self.A, self.phi)
    

@dataclass
class USSteelColumn(Column):
    """
    A class to describe an axially loaded, doubly or singly symmetric steel section
    subject to axial compression designed according to AISC

    'fy': yield stress of the steel material
    'phi': material resistance factor (default = 0.85)
    """
    fy: float
    phi: float = 0.85
    tag: Optional[str] = None

    def factored_compressive_resistance(self) -> float:
        """
        Returns the factored axial capacity of self calculated using AISC
        """
        lambda_c_x = calculate_lambda_c(self.h, self.Ix, self.A, self.kx, self.fy, self.E)
        lambda_c_y = calculate_lambda_c(self.h, self.Iy, self.A, self.ky, self.fy, self.E)
        lambda_c = max(lambda_c_x, lambda_c_y)

        if lambda_c <= 1.5:
            f_cr = (0.658 ** lambda_c**2) * self.fy
        else:
            f_cr = (0.877 / lambda_c**2) * self.fy

        print(lambda_c, f_cr)

        buckling_load = self.phi * self.A * f_cr
        crushing_load = self.phi * self.A * self.fy
        
        return min(buckling_load, crushing_load)

    def factored_buckling_load(self) -> float:
        """
        Returns the Euler buckling load multiplied by self.phi
        """
        p_cr = min(self.critical_buckling_load("x"), self.critical_buckling_load('y'))
        return self.phi * p_cr

    def factored_crushing_load(self) -> float:
        """
        Returns the factored axial capacity of self calculated using CSA S16:19
        in accordance with Cl. 13.3.1 about the governing axis.
        'n', a factor required that describes the fabrication method of the section
            Acceptable values are either:
                1.34 for hot-rolled sections (e.g. HSS, W sections)
                2.24 for built-up plate sections (e.g. WT sections)
        """
        return gross_axial_resistance(self.fy, self.A, self.phi)


def steel_buckling_factor(f_y: float, F_e: float, n: float) -> float:
    """
    Returns the buckling factor, lambda, an intermediate value
    calculated during S16-19, Cl. 13.3.1.
    """
    lamb = sqrt(f_y / F_e)
    factor = (1 + lamb ** (2 * n)) ** (-1 / n)
    return factor


def gross_axial_resistance(yield_stress: float, area: float, phi: float):
    """
    Returns the gross axial resistance of a section made of a material
    with a 'yield_stress' based on its 'area' and material reduction
    factor, 'phi'
    """
    return phi * area * yield_stress


def euler_buckling_load(l: float, E: float, I: float, k: float) -> float:
    """
    Returns the critical Euler buckling load for an axially
    loaded member.

    'l': length
    'E': Elastic modulus
    'I': Moment of inertia
    'k': Effective length factor

    This function assumes consistent units are used for each parameter.
    """
    P_cr = math.pi**2 * E * I / (k * l) ** 2
    return P_cr


def radius_of_gyration(
    moment_of_inertia: float,
    area: float,
) -> float:
    """
    Returns the radius of gyration from the 'moment_of_inertia' and 'area'.
    """
    return (moment_of_inertia / area) ** 0.5


def euler_buckling_moi(
    P_cr: float,
    l: float,
    E: float,
    k: float,
) -> float:
    """
    Returns the moment of inertia required to resist buckling against the
    given axial load, 'P_cr' using the Euler critical buckling load equation.

    'P_cr': Applied axial load
    'l': length
    'E': Elastic modulus
    'k': Effective length factor

    This function assumes consistent units are used for each parameter.
    """
    I_req = P_cr * (k * l**2) / (E * math.pi**2)
    return I_req


def calculate_lambda_c(
    l: float, I: float, a: float, k: float, fy: float, E: float
) -> float:
    """
    Returns the intermediate value, lambda_c, as found in the AISC design code.

    No units are assumed so consistent dimensions must be used.

    'l': length of member
    'I': second moment of area about the axis of study
    'a': cross-sectional area
    'k': effective length factor
    'fy': yield stress
    'E': elastic modulus
    """
    r = radius_of_gyration(I, a)
    lambda_c = k * l / (r * math.pi) * math.sqrt(fy / E)
    return lambda_c
 

def csv_record_to_steelcolumn(record: list[str], **kwargs) -> SteelColumn:
    """
    Returns a SteelColumn instance populated with the data from 'record'. It is assumed
    that the fields in 'record' are in the following order:

    - Tag
    - Area
    - Height
    - Moment of inertia x
    - Moment of inertia y
    - fy
    - E
    - Effective length factor x
    - Effective lenght factor y
    - Dead load
    - Live load
    """
    area = str_to_float(record[1])
    height = str_to_float(record[2])
    moix = str_to_float(record[3])
    moiy = str_to_float(record[4])
    fy = str_to_float(record[5])
    E = str_to_float(record[6])
    kx = str_to_float(record[7])
    ky = str_to_float(record[8])
    tag = record[0]

    sc = SteelColumn(
        h=height, E=E, A=area, Ix=moix, Iy=moiy, kx=kx, ky=ky, fy=fy, tag=tag, **kwargs
    )
    return sc


def convert_csv_data_to_steelcolumns(
    csv_data: list[list[str]], **kwargs
) -> list[SteelColumn]:
    """
    Returns a list of SteelColumn populated with the data provided in 'csv_data
    """
    steel_columns = []
    for csv_record in csv_data:
        sc = csv_record_to_steelcolumn(csv_record, **kwargs)
        steel_columns.append(sc)
    return steel_columns


def calculate_factored_csv_load(record: list[str]) -> float:
    """
    Returns a facotred load based on the raw data inside the 'record'.
    It is assumed that the dead load is located at index 9 and the live
    load is at index 10
    """
    load_combos = load_factors.nbcc_2020_uls_combos()
    loads = {"D_load": str_to_float(record[9]), "L_load": str_to_float(record[10])}
    return load_factors.max_factored_load(loads, load_combos)


def run_all_columns(csv_filename: str, **kwargs) -> list[SteelColumn]:
    """
    Returns a list of SteelColumn populated with the data in 'csv_filename' that have
    had two additional attributes added to each instance dynamically:
        - .factored_load
        - .demand_capacity_ratio

    The factored load is a quantity calculated from loading data contained within the file
    at 'csv_filename'.

    The demand/capacity ratio is calculated from the factored load / factored compressive resistance.

    The file at 'csv_filename' is expected to have the following fields (in this order):
    - Tag
    - Area
    - Height
    - MoIx
    - MoIy
    - fy
    - E
    - kx
    - ky
    - Dead load
    - Live load
    """
    csv_data = read_csv_file(csv_filename)
    steel_columns = convert_csv_data_to_steelcolumns(csv_data[1:], **kwargs)
    factored_loads = [calculate_factored_csv_load(record) for record in csv_data[1:]]
    for idx, steel_column in enumerate(steel_columns):
        steel_column.factored_load = factored_loads[idx]
        dcr = (
            steel_column.factored_load / steel_column.factored_compressive_resistance()
        )
        steel_column.demand_capacity_ratio = dcr
    return steel_columns


def export_steelcolumn_results(
    steel_columns: list[SteelColumn], export_filename: str
) -> None:
    """
    Writes a csv file to disk with name, 'export_filename'. Writes all of the data contained in each
    SteelColumn instance of 'steel_columns'
    """
    with open(export_filename, "w") as csv_export:
        header_row = f"Tag,Factored Load,Factored Compressive Resistance,Utilization\n"
        csv_export.write(header_row)
        for sc in steel_columns:
            tag = sc.tag
            factored_load = sc.factored_load
            factored_resistance = sc.factored_compressive_resistance()
            dcr = sc.demand_capacity_ratio
            line_to_write = f"{tag},{factored_load},{factored_resistance},{dcr}\n"
            csv_export.write(line_to_write)
