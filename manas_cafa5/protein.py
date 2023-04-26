import xml.parsers.expat as xml_parser
import requests
import re

class Protein:
    def __init__(self, name):
        self.name = name
        self.terms = None

    def load_uniprot(self):
        self.load_url(f'https://rest.uniprot.org/uniprotkb/{self.name}.xml')

    def load_url(self, url):
        data = self._fetch_xml_url(url)
        self.terms = Protein.parse_xml(data)

    def load_file(self, file):
        data = open(file, encoding='utf-8').read()
        self.terms = Protein.parse_xml(data)

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

    def parse_xml(xml_data):
        cursor = {
            'terms': { 'go': [] },
            'dbref': None,
        }

        def start_element(cursor, name, attrs):
            atype = attrs.get('type')
            dbref = cursor.get('dbref')
            if dbref != None and name == 'property' and atype != None:
                dbref['properties'][atype] = attrs.get('value')
            if name == 'dbReference' and atype == 'GO':
                dbref = {
                    'id': attrs.get('id'),
                    'properties': {},
                }
                cursor['dbref'] = dbref
                cursor['terms']['go'].append(dbref)

        def end_element(name):
            pass
        def char_data(data):
            pass

        p = xml_parser.ParserCreate()
        p.StartElementHandler = lambda name, attrs: start_element(cursor, name, attrs)
        p.EndElementHandler = end_element
        p.CharacterDataHandler = char_data
        p.Parse(xml_data)

        return cursor.get('terms')
