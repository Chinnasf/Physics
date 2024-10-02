import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os



def get_file(filename):
    """
    returns raw dataframe from any type of file
    """
    with open(filename, "r") as file:
        file_content = file.readlines()
    df = pd.DataFrame(file_content, columns=['data'])
    return df


def get_structure(ionomer,slices):
    data = ionomer.loc[slices[0]:slices[1]].copy().reset_index(drop=True)
    data["data"] = data["data"].str.split('\s').apply(lambda x: list(filter(None, x)))
    data["atom id"] = data["data"].apply(lambda x: x[0]).astype(int)
    data["molecule id"] = data["data"].apply(lambda x: x[1]).astype(int)
    data["atom type"] = data["data"].apply(lambda x: x[2]).astype(int)
    data["charge"] = data["data"].apply(lambda x: x[3]).astype(float)
    data["x"] = data["data"].apply(lambda x: x[4]).astype(float)
    data["y"] = data["data"].apply(lambda x: x[5]).astype(float)
    data["z"] = data["data"].apply(lambda x: x[6]).astype(float)
    # REMARKS: Based on DOI:10.1063/1.4894813 
    data["elements"] = data["atom type"].map({1:"CB",2:"C",3:"F",4:"O",5:"S",6:"H"})
    data = data.sort_values("atom id").reset_index(drop=True)
    return data

def get_bonds(ionomer,slices):
    data = ionomer[slices[0]:slices[1]].copy().reset_index(drop=True)
    data["data"] = data["data"].str.split('\s').apply(lambda x: list(filter(None, x)))
    data["bond id"] = data["data"].apply(lambda x: x[0]).astype(int)
    data["bond type"] = data["data"].apply(lambda x: x[1]).astype(int)
    data["atom1"] = data["data"].apply(lambda x: x[2]).astype(int)
    data["atom2"] = data["data"].apply(lambda x: x[3]).astype(int)
    return data


def get_angles(data, slices):
    angles = data[slices[0]:slices[1]].copy().reset_index(drop=True)
    angles["data"] = angles["data"].str.split('\s').apply(lambda x: list(filter(None, x)))
    angles["angle id"] = angles["data"].apply(lambda x: x[0]).astype(int)
    angles["angle type"] = angles["data"].apply(lambda x: x[1]).astype(int)
    angles["atom1"] = angles["data"].apply(lambda x: x[2]).astype(int)
    angles["atom2"] = angles["data"].apply(lambda x: x[3]).astype(int)
    angles["atom3"] = angles["data"].apply(lambda x: x[4]).astype(int)
    return angles


def get_dihedrals(data,slices):
    dihedrals = data[slices[0]:slices[1]].copy().reset_index(drop=True)
    dihedrals["data"] = dihedrals["data"].str.split('\s').apply(lambda x: list(filter(None, x)))
    dihedrals["dihedral id"] = dihedrals["data"].apply(lambda x: x[0]).astype(int)
    dihedrals["dihedral type"] = dihedrals["data"].apply(lambda x: x[1]).astype(int)
    dihedrals["atom1"] = dihedrals["data"].apply(lambda x: x[2]).astype(int)
    dihedrals["atom2"] = dihedrals["data"].apply(lambda x: x[3]).astype(int)
    dihedrals["atom3"] = dihedrals["data"].apply(lambda x: x[4]).astype(int)
    dihedrals["atom4"] = dihedrals["data"].apply(lambda x: x[5]).astype(int)
    return dihedrals