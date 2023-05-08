import xml.parsers.expat as xml_parser
import requests
import re
import numpy as np

AMINO_ACID_LIST = 'ARNDCEQGHILKMFPSTWYV'
AMINO_ACID_INDEX = { a: AMINO_ACID_LIST.find(a) for a in AMINO_ACID_LIST }

class Protein:
    def __init__(self, name):
        self.name = name
        self.terms = None
        self.sequence = None

    def load_uniprot(self):
        self.load_url(f'https://rest.uniprot.org/uniprotkb/{self.name}.xml')

    def load_url(self, url):
        data = self._fetch_xml_url(url)
        self._apply_parsed(Protein.parse_xml(data))

    def load_file(self, file):
        data = open(file, encoding='utf-8').read()
        self._apply_parsed(Protein.parse_xml(data))

    def _apply_parsed(self, parsed):
        self.terms = parsed.get('terms')
        self.sequence = parsed.get('sequence')

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

    def go_terms(self):
        if self.terms == None:
            self.load_uniprot()
        return self.terms.get('go')

    def one_hot_sequence(self):
        n = len(self.sequence)
        seq = np.ndarray(shape=(n,20), dtype=float, order='C')
        seq.fill(0.0)
        for i in range(0,n):
            j = AMINO_ACID_INDEX.get(self.sequence[i])
            if j is not None:
                seq[i][j] = 1.0
        return seq

    def parse_xml(xml_data):
        cursor = {
            'terms': { 'go': [] },
            'dbref': None,
            'current_name': None,
            'sequence': None,
        }

        def start_element(cursor, name, attrs):
            atype = attrs.get('type')
            dbref = cursor.get('dbref')
            cursor['current_name'] = name
            if dbref != None and name == 'property' and atype != None:
                dbref['properties'][atype] = attrs.get('value')
            if name == 'dbReference' and atype == 'GO':
                dbref = {
                    'id': attrs.get('id'),
                    'properties': {},
                }
                cursor['dbref'] = dbref
                cursor['terms']['go'].append(dbref)

        def end_element(cursor, name):
            cursor['current_name'] = None

        def char_data(data):
            if cursor.get('current_name') == 'sequence':
                cursor['sequence'] = data

        p = xml_parser.ParserCreate()
        p.StartElementHandler = lambda name, attrs: start_element(cursor, name, attrs)
        p.EndElementHandler = lambda name: end_element(cursor, name)
        p.CharacterDataHandler = char_data
        p.Parse(xml_data)

        return { key: cursor.get(key) for key in ['terms','sequence'] }
