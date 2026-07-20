import urllib.request, re
url = 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233049'
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
html = urllib.request.urlopen(req).read().decode('utf-8')
matches = re.findall(r'<a href="/geo/query/acc.cgi\?acc=(GSM\d+)">(.*?)</a>', html)
for m in matches:
    print(f'{m[0]}\t{m[1]}')
