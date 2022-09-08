import re
from pathlib import Path
from time import strftime
from io import TextIOWrapper,StringIO
from rdkit import Chem, RDLogger
from rdkit.Chem import rdChemReactions

class RDFileReader:
    def __init__(self, file, sanitize=False,**kwargs):
        lg = RDLogger.logger()
        lg.setLevel(RDLogger.CRITICAL)
        if isinstance(file, str):
            self._file = open(file)
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        self.__file = iter(self._file.readline,'')
        self._data=self.__reader(sanitize)

    

    def read(self):
        """
        Parse whole file
        :return: list of parsed molecules
        """
        return list(iter(self))

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def close(self, force=False):
        """
        Close opened file
        :param force: force closing of externally opened file or buffer
        """
        if not self._is_buffer or force:
            self._file.close()

    def __iter__(self):
        return (x for x in self._data)

    def __next__(self):
        return next(iter(self))

    def __reader(self, sanitize=True):
        record = parser = mkey = pos = None
        failed = False
        file = self._file
        buffer = ''
        def rxn_from_buffer(buffer, sanitize=sanitize):
            rxn_extractor = re.compile(r"^\$RXN[\s\S]*M  END")
            data_extractor = re.compile(r"^\$DTYPE (.+)\n\$DATUM ([^$]+)", flags=re.M|re.I)
            rxn_block = rxn_extractor.search(buffer)
            rxn_data = data_extractor.findall(buffer)
            try:
                rxn = rdChemReactions.ReactionFromRxnBlock(rxn_block.string)
                if sanitize:
                    [Chem.SanitizeMol(reac) for reac in rxn.GetReactants()]
                    [Chem.SanitizeMol(prod) for prod in rxn.GetProducts()]
                [rxn.SetProp(datum[0], datum[1].strip()) for datum in rxn_data]
                return rxn
            except:
                return

        #TODO: get the last reaction, or start the first one properly - try another generator?
        for line in self.__file:
            if line.startswith('$RFMT'):
                if record is not None:
                    rxn = rxn_from_buffer(buffer)
                    yield rxn
                buffer = ''
                record = True
            else:
                buffer += line
        rxn = rxn_from_buffer(buffer)           
        yield rxn

#TODO: make FileWriter, have WriteToRDF and WriteToCSV as subclasses
class FileWriter:
    def __init__(self, file, append:bool=False):

        if isinstance(file, str):
            self._file = open(file, 'a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open('a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. '
                            'TextIOWrapper, StringIO, BytesIO, BufferedReader and BufferedIOBase subclasses possible')

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        self.write = self.__write_closed

        if not self._is_buffer or force:
            self._file.close()

    @staticmethod
    def __write_closed(_):
        raise ValueError('I/O operation on closed writer')

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()


class _WriteToRDF:
    def __init__(self, file, *, append=False):
        super().__init__(file)
        if not append or not (self._is_buffer or self._file.tell() != 0):
            self.write = self.__write

    def __write(self,data):
        del self.write
        self._file.write(strftime('$RDFILE 1\n$DATM    %d/%m/%y %H:%M\n'))
        self.write(data)

class WriteToRDF(_WriteToRDF, FileWriter):
    def write(self, data):
        if isinstance(data, rdChemReactions.ChemicalReaction):
            rxn_block = rdChemReactions.ReactionToRxnBlock(data)+'\n'
            self._file.write('$RFMT\n')
            self._file.write(rxn_block)
            
        else:
            raise TypeError('not a reaction')
        
        names = (x for x in data.GetPropNames())
        props = (data.GetProp(x) for x in names)
        for k, v in zip(names, props):
            self._file.write(f'$DTYPE {k}\n$DATUM {v}\n')

class _WriteToCSV:
    '''
    not implemented
    '''
    def __init__(self, file, *, append=False):
        super().__init__(file)
        if not append or not (self._is_buffer or self._file.tell() != 0):
            self.write = self.__write

    def __write(self,data):
        del self.write
        self.write(data)

class WriteToCSV(_WriteToCSV, FileWriter):
    '''
    not implemented
    '''
    def write(self, data):
        if isinstance(data, rdChemReactions.ChemicalReaction):
            # convert into a SMILES String

            names = (x for x in data.GetPropNames())
            props = (data.GetProp(x) for x in names)
            for k, v in zip(names, props):
                #TODO: implement Path-like navigation
                k = k.split(':')
                
            pass
    pass

class RxnToDict(rdChemReactions.ChemicalReaction):
    
    def __call__(self) -> list[dict]:
        prop_names = (prop for prop in self.GetPropNames())
        props = [[*name.split(':')[1:], self.GetProp(name)] for name in prop_names]
        common = {'SMILES': rdChemReactions.ReactionToSmiles(self)}
        locals = list()
        var = None
    
        for prop in props:
            if len(prop)==2:
                common[prop[0]]=prop[1]
            elif len(prop)==3:
                if var != prop[0]:
                    var = prop[0]
                    local = dict()
                    locals.append(local)
                local[prop[1]]=prop[2]
        return [common|local for local in locals][:-1]
      