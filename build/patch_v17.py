#!/usr/bin/env python3
"""Website v1.7 patch (round 3): version bump, honest reproducibility wording
on the Methods page, footer version-DOI + analysis-code links, audit v1.7
download source, GSE198950/GSE36830 sensitivity downloads.
Run from web/. The token ZENODO_V17_DOI in build.py is replaced with the real
version DOI after the Zenodo v1.7 release is published."""
import pathlib, sys

P = pathlib.Path('build/build.py')
t = P.read_text(encoding='utf-8')
edits = []
def E(old, new, n=1): edits.append((old, new, n))

E('VERSION = "v1.6"', 'VERSION = "v1.7"')

E('"Reproduction audit (release v1.6): 85 rows',
  '"Reproduction audit (release v1.7): 85 rows')
E('"复现审计（v1.6）：85 行',
  '"复现审计（v1.7）：85 行')

# Methods data-source paragraph: honest wording + analysis link + version DOI
E('build scripts are public with the code repository, so every number is reproducible. Code: <a href=\\"https://github.com/polypatlas/polypatlas\\">github.com/polypatlas/polypatlas</a> \\u00b7 Data archive: <a href=\\"https://doi.org/10.5281/zenodo.23448779\\">Zenodo DOI 10.5281/zenodo.23448779</a>.',
  'build scripts are public with the code repository, and every number on this site traces to a frozen result file in the data archive. Analysis code: <a href=\\"https://github.com/polypatlas/polypatlas/tree/main/analysis\\">github.com/polypatlas/polypatlas/tree/main/analysis</a> \\u00b7 Data archive, release v1.7: <a href=\\"https://doi.org/10.5281/zenodo.ZENODO_V17_DOI\\">Zenodo DOI 10.5281/zenodo.ZENODO_V17_DOI</a> (concept DOI 10.5281/zenodo.23448779).')
E('全部数字可复现。代码仓库：<a href=\\"https://github.com/polypatlas/polypatlas\\">github.com/polypatlas/polypatlas</a> · 数据存档：<a href=\\"https://doi.org/10.5281/zenodo.23448779\\">Zenodo DOI 10.5281/zenodo.23448779</a>。',
  '本站每个数字均可溯源至数据存档中的冻结结果文件。分析代码：<a href=\\"https://github.com/polypatlas/polypatlas/tree/main/analysis\\">github.com/polypatlas/polypatlas/tree/main/analysis</a> · 数据存档（v1.7 版）：<a href=\\"https://doi.org/10.5281/zenodo.ZENODO_V17_DOI\\">Zenodo DOI 10.5281/zenodo.ZENODO_V17_DOI</a>（概念 DOI 10.5281/zenodo.23448779）。')

# Footer: version DOI + analysis code links
E("""FOOTER = (f'<div class="footer">{WORKING_NAME} '
          f'· {VERSION} · '""",
  """FOOTER = (f'<div class="footer">{WORKING_NAME} '
          f'· {VERSION} · '
          f'<a href="https://doi.org/10.5281/zenodo.ZENODO_V17_DOI">Zenodo release v1.7</a> · '
          f'<a href="https://github.com/polypatlas/polypatlas/tree/main/analysis">analysis code</a> · '""")

# Download labels for the new sensitivity files
E('_def("dl_strict",',
  '''_def("dl_sens198950",
     "GSE198950-exclusion sensitivity (release v1.7): consensus lists recomputed without the tissue-site-confounded cohort (9 primary cohorts; 6 strict).",
     "GSE198950 剔除敏感性（v1.7）：剔除组织部位混杂队列后重算的共识基因表（9 个主队列；6 个严格队列）。")
_def("dl_sens36830",
     "GSE36830 donor-handling sensitivity (release v1.7): consensus with the cohort excluded, consensus with its polyp-only contrast, and the marker-panel recount without the cohort.",
     "GSE36830 供者处理敏感性（v1.7）：整队列剔除共识、息肉单侧对比共识，以及剔除该队列后的标志物面板复算。")
_def("dl_strict",''')

# Downloads list entries
E('''        ("did_de.csv", "dl_did"),''',
  '''        ("sens198950_main9_consensus.csv", "dl_sens198950"),
        ("sens198950_strict6_consensus.csv", "dl_sens198950"),
        ("sens36830_excluded_consensus.csv", "dl_sens36830"),
        ("sens36830_polypOnly_consensus.csv", "dl_sens36830"),
        ("sens36830_marker_panel.csv", "dl_sens36830"),
        ("did_de.csv", "dl_did"),''')

# Download copy block: audit v1.7 + sensitivity files
E('''    # v1.6 download files (audit v1.6, strict-tissue sensitivity, DiD)''',
  '''    # v1.7 download files (audit v1.7, strict-tissue sensitivity, DiD, round-3 sensitivities)''')
E('''        (os.path.join(s6, "reproduction_audit_v1.6.csv"), "reproduction_audit.csv"),''',
  '''        (os.path.join(s6, "reproduction_audit_v1.7.csv"), "reproduction_audit.csv"),
        (os.path.join(s2, "sensitivity_gse198950_out", "main9_consensus.csv"), "sens198950_main9_consensus.csv"),
        (os.path.join(s2, "sensitivity_gse198950_out", "strict6_consensus.csv"), "sens198950_strict6_consensus.csv"),
        (os.path.join(s2, "sensitivity_gse36830_donor", "consensus_A_main9_no36830.csv"), "sens36830_excluded_consensus.csv"),
        (os.path.join(s2, "sensitivity_gse36830_donor", "consensus_B_main10_polypOnly.csv"), "sens36830_polypOnly_consensus.csv"),
        (os.path.join(s2, "sensitivity_gse36830_donor", "marker_panel_no36830.csv"), "sens36830_marker_panel.csv"),''')

fails = []
for old, new, n in edits:
    c = t.count(old)
    if c != n: fails.append((c, n, old[:70])); continue
    t = t.replace(old, new)
if fails:
    for c, n, o in fails: print(f"FAIL count={c} expected={n}: {o}")
    sys.exit(1)
P.write_text(t, encoding='utf-8')
print(f"patch_v17 OK: {len(edits)} edits")
