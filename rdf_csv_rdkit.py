from sys import argv
import csv
from tqdm import tqdm
import rdkit_rdf_io

def main(filename, destination):
    with rdkit_rdf_io.RDFileReader(filename, sanitize=True) as rdf:
        data = list()
        for rxn in tqdm(rdf, desc='Sanitising data...'):
            try:
                data.append(rdkit_rdf_io.RxnToDict(rxn)())
            except:
                continue
    fields = ['RX_ID', 'SMILES', 'RX_NVAR', 'RX_RANK', 'RX_RTYP', 'RX_RCT', 'RX_PRO', 'RX_BIN', 'CL', 'LB', 'TI', 'TXT', 'STP', 'SNR', 'MID', 'MTEXT', 'TIM', 'T', 'P', 'PH', 'COND', 'TYP', 'SUB', 'PRT', 'RXDES', 'LCN', 'COM', 'YPRO', 'YD', 'NYD', 'YDO', 'SRCT', 'RGT', 'CAT', 'SOL', 'citation', 'COPYRIGHT']
    keys = dict()
    for rxn in data:
        for entry in rxn:
            for key in entry.keys():
                if key not in fields:
                    keys[key] = None
                    fields.append(key)
    
    with open(destination, 'w+') as csvfile:
        writer = csv.DictWriter(csvfile, restval='',fieldnames=fields, quoting=csv.QUOTE_NONNUMERIC)
    
        writer.writeheader()
    
        for rxn in tqdm(data, desc='Converting...'):
            try:
                for entry in rxn:
                    writer.writerow(entry)
            except KeyError:
                continue

if __name__ == "__main__":
    main(argv[1], argv[2])