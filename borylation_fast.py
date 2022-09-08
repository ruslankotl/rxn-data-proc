# easier version

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from collections import deque
from itertools import chain
from wrapt_timeout_decorator import *
from numpy import nan

# determines if molecules are identical
def molecules_identical(mol1, mol2):
    return mol1.HasSubstructMatch(mol2) & mol2.HasSubstructMatch(mol1)

def smiles_identical(mol1, mol2):
    return Chem.MolToSmiles(mol1)==Chem.MolToSmiles(mol2)

#get number of boron atoms in a molecule
getboroncount = lambda mol : sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum()==5)


# timeout is UNIX-specific, returns borylation substrate and product
@timeout(0.1)
def substrate_product(rxn, sanitize:bool=True)->tuple(str):
    products = [mol for mol in rxn.GetProducts()]
    reactants = [mol for mol in rxn.GetReactants()]
    if sanitize:
        try:
            [Chem.SanitizeMol(prod) for prod in products]
            [Chem.SanitizeMol(reac) for reac in reactants]
        except:
            return None
    for product in products:
        for reactant in reactants:
            if product.HasSubstructMatch(reactant) and \
                getboroncount(product)>getboroncount(reactant):
                return reactant, product
    return None

# get substrate, boron source, and product for mapped reaction
def substrate_product_boron(rxn):
    try:
        substrate, product = substrate_product(rxn)
    except (TypeError,TimeoutError):
        return None, None, None
    for reactant in rxn.GetReactants():
        if molecules_identical(substrate, reactant) == True:
            continue
        if sum(1 for atom in reactant.GetAtoms() if ((atom.GetAtomicNum()==5) and (atom.GetAtomMapNum()>0)))>0:
            return substrate, reactant, product
    return substrate, None, product

# does the reaction look like a formal C–H borylation?
def is_borylation(rxn, sanitize=True):
    return substrate_product(rxn, sanitize) is not None
