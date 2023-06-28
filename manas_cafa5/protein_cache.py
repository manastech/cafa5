from .protein import Protein
import os, io

class ProteinCache:
    def __init__(self, cache_dir):
        self.cache_dir = cache_dir

    def load(self, name, **kwargs):
        file = os.path.join(self.cache_dir, name + '.xml')
        if os.path.exists(file):
            return Protein.from_file(name, file)
        protein = Protein(name, autoload=False)
        data = protein._fetch_xml_url(f'https://rest.uniprot.org/uniprotkb/{name}.xml')
        wh = io.open(file,'w')
        wh.write(data)
        wh.close()
        return Protein.from_data(name, data)
