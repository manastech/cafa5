from Bio.PDB import PDBParser
import requests, ftplib, gzip, os, re
from urllib.parse import urlparse
from io import StringIO
from .utils import fetch_url
import numpy as np

STRUCTURE_TERMS = set([ 'pdb', 'alphafolddb' ])

class Structure:
    def __init__(self, term):
        self.type = term.get('type')
        self.id = term.get('id')
        self.properties = term.get('properties')
        self.structure = None

    def __repr__(self):
        return f'<manas_cafa5.Structure id={self.id} type={self.type} properties={self.properties}>'

    def load(self):
        if self.type == 'pdb':
            self.load_url(f'ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb{self.id}.ent.gz')
        elif self.type == 'alphafolddb':
            self.load_url(f'https://alphafold.ebi.ac.uk/files/AF-{self.id}-F1-model_v4.pdb')
        else:
            raise RuntimeError(f'default endpoint for type {self.type} not available')

    def load_file(self, file, id=None):
        data = None
        if file.split('.')[-1] == 'gz':
            data = gzip.open(file, 'r').read()
        else:
            data = open(file, 'r', encoding='UTF-8').read()
        if id is None:
            id = re.sub(r'(\.(ent|pdb))?\.gz$', '', os.path.split(file)[-1])
        if isinstance(data, bytes):
            data = str(data, 'utf-8')
        parser = PDBParser()
        self.structure = parser.get_structure(id, StringIO(data))

    def load_url(self, url, id=None):
        data = fetch_url(url)
        if data[0] == 0x1f and data[1] == 0x8b: # gzip magic number
            data = gzip.decompress(data)
        parser = PDBParser()
        if id is None:
            u = urlparse(url)
            id = re.sub(r'(\.(ent|pdb))?\.gz$', '', os.path.split(u.path)[-1])
        if isinstance(data, bytes):
            data = str(data, 'utf-8')
        self.structure = parser.get_structure(id, StringIO(data))

    # https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
    def contact_map(self, chain_one, chain_two, threshold, index=0):
        if self.structure is None:
            raise RuntimeError(f'structure for type {self.type} not loaded. call load() first')
        if len(self.structure) == 0:
            raise RuntimeError('no structures available')
        model = self.structure[index]
        c1 = model[chain_one]
        c2 = model[chain_two]
        cmap = np.zeros((len(c1), len(c2)), np.float)
        for row, residue_one in enumerate(c1) :
            for col, residue_two in enumerate(c2) :
                diff  = residue_one['CA'].coord - residue_two['CA'].coord
                cmap[row,col] = np.sqrt(np.sum(diff * diff))
        return np.where(cmap < threshold, 1.0, 0.0)

    def get_chain_ids(self, index=0):
        return [ ch.id for ch in self.structure[index].get_list() ]
