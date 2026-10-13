# -*- coding: utf-8 -*-
"""Charoentong 2017 (PMC5215856) supplementary files via Europe PMC (Range-resumable)."""
import requests, re, os, subprocess, json

OUT = r"F:\232机器学习联合神经网络\CRC\99.P1\09.charoentong"
os.makedirs(OUT, exist_ok=True)
PMCID = "PMC5215856"

H = {"User-Agent": "research"}
# 1) fullTextXML -> supplementary file list
r = requests.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/fullTextXML", headers=H, timeout=120)
xml = r.text
print("fullTextXML chars:", len(xml))
sups = re.findall(r'<supplementary-material[^>]*id="([^"]*)"[^>]*>(.*?)</supplementary-material>', xml, re.S)
if not sups:
    sups = re.findall(r'<supplementary-material[^>]*>(.*?)</supplementary-material>', xml, re.S)
    sups = [("", s) for _, s in sups]
files = []
for sid, body in sups:
    href = re.search(r'xlink:href="([^"]+)"', body) or re.search(r'href="([^"]+)"', body)
    label = re.search(r'<label>([^<]+)</label>', body)
    cap = re.search(r'<caption>(.*?)</caption>', body, re.S)
    captext = re.sub(r"<[^>]+>", " ", cap.group(1)).strip()[:90] if cap else ""
    if href:
        files.append((sid, href.group(1), label.group(1) if label else "", captext))
print("supplementary files:", len(files))
for f in files:
    print("  ", f)

json.dump(files, open(os.path.join(OUT, "supp_files.json"), "w"), ensure_ascii=False, indent=1)

# 2) download each (resumable via Range)
for sid, href, label, cap in files:
    name = os.path.basename(href)
    dest = os.path.join(OUT, name)
    url = f"https://europepmc.org/articles/{PMCID}/bin/{name}"
    print("downloading:", name, "->", dest)
    for attempt in range(6):
        rc = subprocess.run(["curl", "-sL", "-C", "-", "--connect-timeout", "15",
                             "--max-time", "900", "--retry", "3", "-o", dest, url],
                            capture_output=True).returncode
        sz = os.path.getsize(dest) if os.path.exists(dest) else 0
        print("   attempt", attempt + 1, "->", sz, "bytes (rc", rc, ")")
        if sz > 5000:  # xlsx/pdf files are much larger than 5KB
            break
        time.sleep(3)
print("downloads complete")
