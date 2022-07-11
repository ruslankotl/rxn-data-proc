from sys import argv
from csv import DictReader, DictWriter
from urllib.error import HTTPError, URLError
from cirpy import resolve
from pubchempy import get_compounds
from tqdm import tqdm
from lxml.etree import XMLSyntaxError

'''A program which parses reaction csv file to generate a complete reaction smiles string.'''
def main(input_file, output_file, smiles, columns_to_resolve):
    
    # read the csv and arguments
    with open(input_file, 'r') as open_input:
        reader = DictReader(open_input)
        field_list = reader.fieldnames
        compound_cache = {'nan':None}
        
        with open(output_file, 'w') as open_output:
            writer = DictWriter(open_output, fieldnames = field_list)
            writer.writeheader()
            for row in tqdm(reader):
                write_flag = True
                agents = list()
                for column in columns_to_resolve:
                    compounds_to_resolve = row[column].split('|')
                    if len(compounds_to_resolve[0])<1:
                        continue
                    for compound in compounds_to_resolve:
                        if compound not in compound_cache:
                            try:
                                compound_cache[compound] = resolve(compound,'smiles')
                                if compound_cache[compound] is None:
                                    compound_cache[compound] = get_compounds(compound, 'name')[0].canonical_smiles
                            except:
                                compound_cache[compound] = None
                            
                        if compound_cache[compound] is None:
                            write_flag = False
                            break
                        else:
                            agents.append(compound_cache[compound])
                if write_flag == True:
                    agents_string = '.'.join(agents)
                    smiles_split = row[smiles].split('>')
                    smiles_split[1] = agents_string
                    smiles_joined = ('>').join(smiles_split)
                    row[smiles] = smiles_joined
                    writer.writerow(row)
                                
    print('did them all')
                    
if __name__=='__main__':
    main(argv[1], argv[2], argv[3], argv[4:-1])