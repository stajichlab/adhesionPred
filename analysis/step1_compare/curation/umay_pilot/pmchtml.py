import sys
import urllib.request

from fetch import html2txt

for pmc in sys.argv[1:]:
    req = urllib.request.Request(
        f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/",
        headers={
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 Chrome/120 Safari/537.36"
        },
    )
    try:
        h = urllib.request.urlopen(req, timeout=60).read().decode("utf8", "ignore")
    except Exception as e:
        print(pmc, "ERR", e)
        continue
    t = html2txt(h)
    open(f"texts/{pmc}.html.txt", "w").write(t)
    print(pmc, len(t))
