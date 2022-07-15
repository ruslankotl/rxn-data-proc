from CGRtools.files import RDFRead, RDFWrite
from tqdm import tqdm
from sys import argv
from wrapt_timeout_decorator import *
from borylation_rdf import is_borylation

def main(rdf_files, destination):
    
    @timeout(0.1)
    def timed_borylation(entry):
        return is_borylation(entry)

    def rdf_borylation_filter(filename):
        if '.rdf' not in filename:
            yield
        else:
            logname = filename.rsplit('.', maxsplit=1)[0]+'.log'
            log = open(logname, 'w')
        
        with RDFRead(filename) as f:
            for entry in tqdm(f, mininterval=0.5, desc=f'{filename}'):
                try:
                    entry.canonicalize()
                except:
                    continue
                try:
                    if timed_borylation(entry):
                        yield entry
                except TypeError:
                    continue
                except TimeoutError:
                    log.write(str(entry)+'\n')
        log.close()
    
    with RDFWrite(destination) as f:
        for num, piece in enumerate(rdf_files):
            print('File {0} out of {1}'.format(num+1, len(rdf_files)))
            for rxn in rdf_borylation_filter(piece):
                try:
                    f.write(rxn)
                except TimeoutError:
                    continue

if __name__ == "__main__":
    main(argv[1:-1], argv[-1])