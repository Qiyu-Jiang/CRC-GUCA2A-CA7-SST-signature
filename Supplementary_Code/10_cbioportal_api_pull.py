# -*- coding: utf-8 -*-
"""cBioPortal pull v4: per-gene expression (18 genes), discrete-cna-data, mutations; robust."""
import requests, csv, os, time

H = {'User-Agent': 'Mozilla/5.0 (research)', 'Accept': 'application/json'}
base = 'https://www.cbioportal.org/api'
outdir = r'.\outputs/phase1'
os.makedirs(outdir, exist_ok=True)

entrez = {'GUCA2A': 2980, 'CA7': 766, 'SST': 6750, 'PDCD1': 5133, 'CD274': 29126,
          'CTLA4': 1493, 'LAG3': 3902, 'TIGIT': 201633, 'HAVCR2': 84868, 'BTLA': 151888,
          'IDO1': 3620, 'CD276': 80381, 'VSIR': 64115, 'SIGLEC15': 284266,
          'TNFRSF9': 3604, 'ICOS': 29851, 'CD27': 939, 'CD28': 940}

def get(url, params=None):
    r = requests.get(url, params=params, headers=H, timeout=300)
    r.raise_for_status()
    return r.json()

# sample list size
try:
    sl = get(f'{base}/sample-lists/coadread_tcga_pub_all')
    n_all = sl.get('sampleCount')
except Exception:
    n_all = None
print('all sample count:', n_all)

# 1) per-gene expression
prof = 'coadread_tcga_pub_rna_seq_mrna'
expr = {}
for sym, eid in entrez.items():
    for attempt in range(3):
        try:
            rows = get(f'{base}/molecular-profiles/{prof}/molecular-data',
                       params={'sampleListId': 'coadread_tcga_pub_all', 'entrezGeneId': str(eid)})
            d = {}
            for row in rows:
                d[row['sampleId']] = row['value']
            expr[sym] = d
            print(sym, 'expr samples:', len(d))
            break
        except Exception as e:
            print(sym, 'retry', attempt, str(e)[:80]); time.sleep(2)
samples = sorted({s for d in expr.values() for s in d})
with open(os.path.join(outdir, 'tcga_expr_wide.csv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['gene'] + samples)
    for g, d in expr.items():
        w.writerow([g] + [d.get(s, '') for s in samples])
print('expression matrix:', len(expr), 'x', len(samples))

# 2) sample->patient
srows = get(f'{base}/studies/coadread_tcga_pub/samples', params={'sampleListId': 'coadread_tcga_pub_all', 'projection': 'SUMMARY'})
s2p = {s['sampleId']: s['patientId'] for s in srows}
with open(os.path.join(outdir, 'tcga_sample_patient.tsv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['sampleId', 'patientId'])
    for s, p in s2p.items(): w.writerow([s, p])
print('samples:', len(s2p))

# 3) mutations + 4) CNA (discrete-cna-data)
mc = {}
for sym, eid in [('GUCA2A', 2980), ('CA7', 766), ('SST', 6750)]:
    time.sleep(1)
    try:
        rows = get(f'{base}/molecular-profiles/coadread_tcga_pub_mutations/mutations',
                   params={'sampleListId': 'coadread_tcga_pub_all', 'entrezGeneId': str(eid), 'projection': 'SUMMARY'})
        mc.setdefault(sym, {})['mutationCount'] = len(rows)
        mc[sym]['mutatedSamples'] = len({r.get('sampleId') for r in rows}) if isinstance(rows, list) else 0
    except Exception as e:
        mc.setdefault(sym, {})['mutationError'] = str(e)[:100]
    time.sleep(1)
    try:
        rows = get(f'{base}/molecular-profiles/coadread_tcga_pub_gistic/discrete-cna-data',
                   params={'sampleListId': 'coadread_tcga_pub_all', 'entrezGeneId': str(eid), 'projection': 'SUMMARY'})
        amp = sum(1 for r0 in rows if r0.get('alteration') == 2)
        dele = sum(1 for r0 in rows if r0.get('alteration') == -2)
        mc[sym]['CNA_amp'] = amp; mc[sym]['CNA_del'] = dele; mc[sym]['CNA_total'] = len(rows)
    except Exception as e:
        mc[sym]['cnaError'] = str(e)[:100]
    print(sym, mc[sym])

with open(os.path.join(outdir, 'tcga_mut_cna.csv'), 'w', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['gene', 'mutationCount', 'mutatedSamples', 'CNA_amp', 'CNA_del', 'CNA_total'])
    for sym in ('GUCA2A', 'CA7', 'SST'):
        d0 = mc.get(sym, {})
        w.writerow([sym, d0.get('mutationCount', ''), d0.get('mutatedSamples', ''),
                    d0.get('CNA_amp', ''), d0.get('CNA_del', ''), d0.get('CNA_total', '')])
print('done')
