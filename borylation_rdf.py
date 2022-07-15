from CGRtools.files import RDFRead, SMILESRead
from CGRtools.reactor import CGRReactor 
from CGRtools.containers import QueryContainer, MoleculeContainer, ReactionContainer
from io import StringIO
from collections import deque

getboroncount = lambda mol : sum(1 for atom in mol.atoms() if atom[1].atomic_number==5)

boron = QueryContainer()
boron.add_atom('C')
boron.add_atom('B')
boron.add_bond(1, 2, 1)
deborylated = QueryContainer()
deborylated.add_atom('C')
template = ReactionContainer([boron], [deborylated])
reactor = CGRReactor(template)

def deborylation(mol):
    result = set([mol])
    d = deque([mol])
    while len(d)>0:
        current_deborylated = set(r for r in reactor(d.popleft()))
        for prod in current_deborylated:
            if prod not in result:
                result.add(prod)
                d.append(prod)
    result.discard(mol)    
    return list(result)

#takes in ReactionContainers from a successfully parsed RDF file
def substrate_product(rxn):
    products = [prod for prod in rxn.products if getboroncount(prod)>0]
    reactants = rxn.reactants
    for product in products:
        deborylated = deborylation(product)
        for deborylated_prod in deborylated:
            if deborylated_prod in reactants:
                    return deborylated_prod, product
    return None


def is_borylation(rxn):
    return substrate_product(rxn) is not None