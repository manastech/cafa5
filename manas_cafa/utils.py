import requests, ftplib
from urllib.parse import urlparse

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
