from pathlib import Path
from sys import argv
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from borylation import substrate_product

def parent_centres(reactant_smiles, product_smiles)->list:
        reactant = Chem.MolFromSmiles(reactant_smiles)
        product = Chem.MolFromSmiles(product_smiles)
        match = product.GetSubstructMatch(reactant)

        reaction_sites = [atom.GetIdx() for atom in product.GetAtoms() \
         if (atom.GetIdx() in match) and \
         not all((n.GetIdx() in match) for n in atom.GetNeighbors())]

        centres = [reactant_index for reactant_index, prod_index in enumerate(match)
            if prod_index in reaction_sites]

        return centres

def forgot_to_annotate(rxn_smiles:str):

    substrate_smi, _,product_smi = rxn_smiles.split('>')
    rxn_centres = parent_centres(substrate_smi,product_smi)
    ch_bonds = [match[0] for match in Chem.QuickSmartsMatch(substrate_smi,'[cH:1]')]
    substrate = Chem.MolFromSmiles(substrate_smi)
    borylate = rdChemReactions.ReactionFromSmarts('[cH:1]>>[c:1]-B1OC(C)(C)C(C)(C)O1')
    prods = [prod[0] for prod in borylate.RunReactants([substrate])]
    # dict comprehension: if reactive site, classified as 1, else 0, unique prods taken care as dict keys
    outcomes = {Chem.MolToSmiles(prod):int(site in rxn_centres) for site,prod in zip(ch_bonds, prods)}
    # output would return reaction SMILES with product replaced
    output = {rxn_smiles.replace(product_smi,site):outcome for site,outcome in outcomes.items()}
    return output

    
def main(folder):
    '''folder: path to a T5 dataset. To add: 0 yields'''
    dataset_path = Path(folder)
    yields_path = dataset_path/'yield'
    class_path = (dataset_path/'classification')
    new_yields = dataset_path/'new_yield'
    class_path.mkdir(parents=True,exist_ok=True)
    new_yields.mkdir(parents=True,exist_ok=True)

    files_to_read = zip(sorted(yields_path.glob('*.source')),sorted(yields_path.glob('*.target')))
    for file, yields in files_to_read:
        filename = file.parts[-1].split('.')[0]
        cls_source_path = class_path/f'{filename}.source'
        cls_target_path = class_path/f'{filename}.target'
        yld_source_path = new_yields/f'{filename}.source'
        yld_target_path = new_yields/f'{filename}.target'
        with open(file) as input_file, \
        open(yields) as yields_file,\
        open(cls_source_path,'w') as class_source, \
        open(cls_target_path,'w') as class_target,\
        open(yld_source_path,'w') as yield_source,\
        open(yld_target_path,'w') as yield_target:
            for rxn, yields in zip(input_file, yields_file):
                output = forgot_to_annotate(rxn)
                for site, reacts in output.items():
                    class_source.write(site+'\n')
                    yield_source.write(site+'\n')
                    class_target.write(str(reacts)+'\n')
                    yield_target.write(str(float(yields)*reacts)+'\n')

if __name__=='__main__':
    main(argv[1])