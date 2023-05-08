import xml.parsers.expat as xml_parser
import requests
import re
import numpy as np

AMINO_ACID_LIST = 'ARNDCEQGHILKMFPSTWYV'
AMINO_ACID_INDEX = { a: AMINO_ACID_LIST.find(a) for a in AMINO_ACID_LIST }

class Protein:
    def __init__(self, name):
        self.name = name
        self.entries = None

    def __len__(self):
        return self.entries and len(self.entries) or 0

    def load_uniprot(self):
        self.load_url(f'https://rest.uniprot.org/uniprotkb/{self.name}.xml')

    def load_url(self, url):
        data = self._fetch_xml_url(url)
        self.entries = Protein.parse_xml(data)

    def load_file(self, file):
        data = open(file, encoding='utf-8').read()
        self.entries = Protein.parse_xml(data)

    def _fetch_xml_url(self, url):
        r = requests.get(url)
        if r.status_code != 200:
            raise RuntimeError(f'unexpected status code fetching xml for {self.name}: {r.status_code}')
        ctype = r.headers.get('content-type')
        if ctype == None:
            raise RuntimeError('content-type header not sent')
        if re.match(r'^(application|text)/xml\s*(;|$)', ctype) == None:
            raise RuntimeError(('unexpected content-type: expected application/xml or text/xml, '
                f'received: {r.headers["content-type"]}'))
        return r.text

    def _auto_load(self):
        if self.entries is None:
            self.load_uniprot()

    def go_terms(self, index=0):
        self._auto_load()
        entry = self.entries[index]
        return entry and entry['terms']['go']

    def get_sequence(self, index=0):
        self._auto_load()
        entry = self.entries[index]
        return entry and entry.get('sequence')

    def one_hot_sequence(self, index=0):
        self._auto_load()
        seq = self.get_sequence(index)
        if seq is None:
            return None
        n = len(seq)
        nd_seq = np.ndarray(shape=(n,20), dtype=float, order='C')
        nd_seq.fill(0.0)
        for i in range(0,n):
            j = AMINO_ACID_INDEX.get(seq[i])
            if j is not None:
                nd_seq[i][j] = 1.0
        return nd_seq

    def parse_xml(xml_data):
        cursor = {
            'entries': [],
            'current_name': None,
        }

        def start_element(cursor, name, attrs):
            atype = attrs.get('type')
            lname = name and name.lower()
            latype = atype and atype.lower()
            if lname == 'entry':
                cursor['entries'].append({
                    'terms': { 'go': [] },
                    'sequence': None,
                })
            elif lname == 'dbreference' and latype == 'go':
                cursor['entries'][-1]['terms']['go'].append({
                    'id': attrs.get('id'),
                    'properties': {},
                })
            else:
                cursor['current_name'] = name
                entries = cursor['entries']
                terms = len(entries) > 0 and entries[-1]['terms']
                if terms is not None and lname == 'property' and latype is not None:
                    terms['go'][-1]['properties'][latype] = attrs.get('value')

        def end_element(cursor, name):
            cursor['current_name'] = None

        def char_data(data):
            if cursor.get('current_name') == 'sequence':
                cursor['entries'][-1]['sequence'] = data

        p = xml_parser.ParserCreate()
        p.StartElementHandler = lambda name, attrs: start_element(cursor, name, attrs)
        p.EndElementHandler = lambda name: end_element(cursor, name)
        p.CharacterDataHandler = char_data
        p.Parse(xml_data)

        return cursor['entries']
