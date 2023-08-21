import xml.parsers.expat as xml_parser
import sys
import requests
import re
import numpy as np
import obonet, networkx
from functools import reduce


# this is a pointer to the module object instance itself.
this = sys.modules[__name__]
this.uniprot_entries_path = '/data/uniprot/entries/'
this.uniprot_entries_url = 'https://rest.uniprot.org/uniprotkb/'

AMINO_ACID_LIST = 'ARNDCEQGHILKMFPSTWYV'
AMINO_ACID_INDEX = { a: AMINO_ACID_LIST.find(a) for a in AMINO_ACID_LIST }

class Protein:
    def __init__(self, name):
        self.name = name
        self.terms = None
        self.sequence = None

    def load_uniprot(self):
        try: 
            self.load_file(this.uniprot_entries_path+self.name[-2:]+"/"+self.name+".xml")
        except FileNotFoundError:
            self.load_url(f'{this.uniprot_entries_url}{self.name}.xml')

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

    def structure_terms(self):
        if self.terms == None:
            self.load_uniprot()
        return self.terms.get('pdb')

    def go_terms_children(self, graph, max_distance):
        term_set = set()
        for dist in range(1,max_distance+1):
            term_set = reduce(
                lambda terms, term: terms.union(
                    networkx.descendants_at_distance(graph, term['id'], dist)
                ),
                self.go_terms(),
                term_set
            )
        return list(term_set)

    def go_terms_parents(self, graph, max_distance):
        term_set = set()
        last_terms = self.go_terms()
        for _ in range(max_distance):
            new_last_terms = set()
            for term in last_terms:
                ancestors = []
                try:
                    ancestors = networkx.ancestors(graph, term['id'])
                except:
                    try:
                        networkx.ancestors(graph, term)
                    except:
                        pass
                for ancestor in ancestors:
                    new_last_terms.add(ancestor)
            term_set = term_set.union(new_last_terms)
            last_terms = new_last_terms
        return list(term_set)

    @staticmethod
    def build_graph(url_or_file):
        # example url to use: https://current.geneontology.org/ontology/go-basic.obo
        return obonet.read_obo(url_or_file)

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
            'terms': {'go': [], 'pdb': []},
            'dbref': None,
            'current_name': None,
            'sequence': None,
        }

        def start_element(cursor, name, attrs):
            name = name.lower()
            atype = attrs.get('type')
            atype = atype and atype.lower()
            dbref = cursor.get('dbref')
            cursor['current_name'] = name
            if dbref is not None and name == 'property' and atype is not None:
                dbref['properties'][atype] = attrs.get('value')
            if name == 'dbreference' and atype is not None:
                dbref = {
                    'type': atype,
                    'id': attrs.get('id'),
                    'properties': {},
                }
                cursor['dbref'] = dbref
                terms = cursor['terms'].get(atype)
                if terms is None:
                    terms = []
                    cursor['terms'][atype] = terms
                terms.append(dbref)

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
