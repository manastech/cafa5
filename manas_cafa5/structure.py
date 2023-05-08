from Bio.PDB import PDBParser
import requests, ftplib, gzip, os, re
from urllib.parse import urlparse
from io import StringIO
from .utils import fetch_url

class Structure:
    def __init__(self):
        self.structure = None

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

    def load_pdb(self, name):
        self.load_url(f'ftp://ftp.wwpdb.org//pub/pdb/data/structures/all/pdb/pdb{name}.ent.gz')

