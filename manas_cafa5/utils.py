import networkx
import requests, ftplib
from urllib.parse import urlparse
from functools import reduce

def fetch_url(url):
    if url.find('ftp:') == 0:
        u = urlparse(url)
        ftp = ftplib.FTP()
        ftp.connect(u.hostname, u.port or 21)
        ftp.login()
        chunks = []
        ftp.retrbinary(f'RETR {u.path}', lambda chunk: chunks.append(chunk))
        ftp.quit()
        return b''.join(chunks)
    else:
        r = requests.get(url)
        if r.status_code != 200:
            raise RuntimeError(f'unexpected status code while fetching {url}:  {r.status_code}')
        return r.content

def ancestors_within_distance(graph, term, max_distance):
    parents = set([term])
    for dist in range(1,max_distance+1):
        for parent in parents:
            parents = parents.union(networkx.ancestors(graph, parent))
    parents.discard(term)
    return parents
