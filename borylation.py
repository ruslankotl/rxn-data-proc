# determine whether the reaction is a single-step C–H borylation, rdkit style

from rdkit import Chem
from rdkit.Chem import rdChemReactions
from collections import deque
from itertools import chain
from wrapt_timeout_decorator import *
from numpy import nan

# determines if molecules are identical
def molecules_identical(mol1, mol2):
    return mol1.HasSubstructMatch(mol2) & mol2.HasSubstructMatch(mol1)

# a relevant molecular template
deborylation = rdChemReactions.ReactionFromSmarts('[#6:1][B]>>[#6:1]')

#get number of boron atoms in a molecule
getboroncount = lambda mol : sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum()==5)

# returns all possible deborylation products with ONE break of C–B bond
def unique_deborylation_single_step(mol):
    deborylation_products = list(chain.from_iterable(deborylation.RunReactants([mol])))
    deborylation_unique_products = []
    
    for prod in deborylation_products:
        unique = True
        for u_prod in deborylation_unique_products:
            unique ^= molecules_identical(prod, u_prod)
        if unique:
            deborylation_unique_products.append(prod)
    
    return deborylation_unique_products
    
# generate all possible deborylation products
def deque_deborylation(mol):
    
    if getboroncount(mol)==1:
        return unique_deborylation_single_step(mol)
    
    result = [mol]
    d = deque(result)
    
    while (len(d)>0):
        current_deborylated = unique_deborylation_single_step(d.popleft())
        for prod in current_deborylated:
            
            unique = True
            for u_prod in result:
                unique ^= molecules_identical(prod, u_prod)
            if unique:
                result.append(prod)
                d.append(prod)
        
    
    return result[1:]

# timeout is UNIX-specific, returns borylation substrate and product
@timeout(0.1)
def substrate_product(rxn, sanitize=True):
    borylation = False
    products = [prod for prod in rxn.GetProducts() if getboroncount(prod)>0]
    reactants = rxn.GetReactants()
    if sanitize:
        try:
            [Chem.SanitizeMol(prod) for prod in products]
            [Chem.SanitizeMol(reac) for reac in reactants]
        except:
            return None
    for product in products:
        deborylated_products = deque_deborylation(product)
        for deborylated_product in deborylated_products:
            for reactant in reactants:
                borylation = molecules_identical(reactant, deborylated_product)
                if borylation:
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
