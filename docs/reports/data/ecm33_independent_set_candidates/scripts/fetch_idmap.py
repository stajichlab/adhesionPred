# ruff: noqa
# Scoping script kept exactly as run on 2026-10-23 (stdlib only).
import urllib.request,urllib.parse,json,time,sys
ids=[l.strip() for l in open(sys.argv[1]) if l.strip()]
H={"User-Agent":"adhesionPred-scoping/1.0"}
d=urllib.parse.urlencode({'from':'UniProtKB_AC-ID','to':'UniRef50','ids':','.join(ids)}).encode()
j=json.load(urllib.request.urlopen(urllib.request.Request('https://rest.uniprot.org/idmapping/run',data=d,headers=H),timeout=60))['jobId']
while True:
    s=json.load(urllib.request.urlopen(urllib.request.Request('https://rest.uniprot.org/idmapping/status/'+j,headers=H),timeout=60))
    if 'results' in s or s.get('jobStatus')=='FINISHED': break
    if s.get('jobStatus') in ('ERROR','FAILED'): sys.exit(str(s))
    time.sleep(3)
u='https://rest.uniprot.org/idmapping/uniref/results/stream/'+j+'?format=tsv&fields=id,accessions'
# stream to uniref gives cluster rows, not mapping; use plain results
u='https://rest.uniprot.org/idmapping/uniref/results/stream/'+j+'?format=json'
r=json.load(urllib.request.urlopen(urllib.request.Request(u,headers=H),timeout=120))
with open(sys.argv[2],'w') as f:
    f.write('accession\tuniref50\n')
    for x in r['results']: f.write(x['from']+'\t'+x['to']['id']+'\n')
print(len(r['results']),r.get('failedIds'))
