#!/usr/bin/env python3
"""v1.6 website patch for build.py. Asserted single-hit replacements."""
from pathlib import Path

P = Path(__file__).with_name("build.py")
t = P.read_text(encoding="utf-8")
EDITS = []
def E(old, new): EDITS.append((old, new))

E('VERSION = "v1.5"', 'VERSION = "v1.6"')

# --- gene page: extended meta fields -------------------------------------
E('_def("g_f1_k", "cohorts", "队列数")',
  '_def("g_f1_k", "cohorts", "队列数")\n'
  '_def("g_f1_tau2", "\\u03c4\\u00b2 (between-cohort variance)", "\\u03c4\\u00b2（队列间方差）")\n'
  '_def("g_f1_pi", "95% prediction interval", "95% 预测区间")\n'
  '_def("g_f1_loo", "Leave-one-cohort-out range", "留一队列范围")')
E('"g_f1_t", "g_f1_combined", "g_f1_ci", "g_f1_p", "g_f1_i2", "g_f1_k",',
  '"g_f1_t", "g_f1_combined", "g_f1_ci", "g_f1_p", "g_f1_i2", "g_f1_k",\n            "g_f1_tau2", "g_f1_pi", "g_f1_loo",')
E("""  html += '<p><b>'+gs('g_f1_combined')+':</b> '+fmt(m.lfc,3)+' &nbsp; <b>'+gs('g_f1_ci')+':</b> ['+fmt(m.lo,3)+', '+fmt(m.hi,3)+'] &nbsp; <b>'+gs('g_f1_p')+':</b> '+fmtP(m.p)+' &nbsp; <b>'+gs('g_f1_i2')+':</b> '+fmt(m.i2,1)+'% &nbsp; <b>'+gs('g_f1_k')+':</b> '+m.k+'</p>';""",
  """  html += '<p><b>'+gs('g_f1_combined')+':</b> '+fmt(m.lfc,3)+' &nbsp; <b>'+gs('g_f1_ci')+':</b> ['+fmt(m.lo,3)+', '+fmt(m.hi,3)+'] &nbsp; <b>'+gs('g_f1_p')+':</b> '+fmtP(m.p)+' &nbsp; <b>'+gs('g_f1_i2')+':</b> '+fmt(m.i2,1)+'% &nbsp; <b>'+gs('g_f1_k')+':</b> '+m.k+'</p>';
  if(m.t2!==undefined){ html += '<p><b>'+gs('g_f1_tau2')+':</b> '+fmt(m.t2,3)+' &nbsp; <b>'+gs('g_f1_pi')+':</b> ['+fmt(m.plo,3)+', '+fmt(m.phi,3)+'] &nbsp; <b>'+gs('g_f1_loo')+':</b> ['+fmt(m.lmin,3)+', '+fmt(m.lmax,3)+']</p>'; }""")
E('"Meta-analysis: per-gene combined logFC, 95% CI, combined p, I\\u00b2, k (DerSimonian\\u2013Laird, 10 primary contrasts; genes with k\\u22653)."',
  '"Meta-analysis: per-gene combined logFC, 95% CI, combined p, I\\u00b2, \\u03c4\\u00b2, 95% prediction interval, leave-one-cohort-out range, k (DerSimonian\\u2013Laird, 10 primary contrasts; genes with k\\u22653)."')
E('"Meta 分析：逐基因合并 logFC、95% CI、合并 p、I\\u00b2、k（DerSimonian\\u2013Laird，10 个主对比；k\\u22653 的基因）。"',
  '"Meta 分析：逐基因合并 logFC、95% CI、合并 p、I\\u00b2、\\u03c4\\u00b2、95% 预测区间、留一队列范围、k（DerSimonian\\u2013Laird，10 个主对比；k\\u22653 的基因）。"')

# --- methods page: correct low-power list --------------------------------
E('''     "Six bulk cohorts have n\\u22647 per group or low power (GSE72713, GSE150846, GSE174681, GSE330761, "
     "GSE36830, GSE194282) \\u2014 their weight in the evidence matrix should be discounted.",
     "6 个 bulk 队列 n≤7/组或功效低（GSE72713、GSE150846、GSE174681、GSE330761、GSE36830、GSE194282），证据矩阵中权重应下调。"),''',
  '''     "Five bulk cohorts have fewer than five samples in at least one primary-contrast arm "
     "(GSE198950, GSE174681, GSE150846, GSE330761, GSE72713) \\u2014 their contrasts fall outside the "
     "\\\\u22655-per-arm denominator of the marker-recovery validation and are flagged wherever they appear.",
     "5 个 bulk 队列的主对比至少一臂少于 5 例（GSE198950、GSE174681、GSE150846、GSE330761、GSE72713）——其对比不进入标志物回收验证的每臂≥5 分母，出现处均有标注。"),''')

# --- treatment page: DiD block --------------------------------------------
E('_def("nav_treatment", "Treatment response", "治疗反应")',
  '_def("nav_treatment", "Treatment response", "治疗反应")\n'
  '_def("t_did_h", "Drug specificity: difference-in-differences (exploratory)", "药物特异性：差值中的差值（探索性）")')
E('''                "<th>Overlap score</th></tr>" + drows + "</table>")
    body = (DISCLAIMER''',
  '''                "<th>Overlap score</th></tr>" + drows + "</table>")
    # v1.6: exploratory difference-in-differences block (data: step2 did_dupilimab)
    did_dir = os.path.join(PROJ, "analysis", "step2_bulk", "did_dupilimab")
    mdid = read_csv_dict(os.path.join(did_dir, "markers_did.csv"))
    pdid = read_csv_dict(os.path.join(did_dir, "program_score_did.csv"))[0]
    drows2 = "".join(
        f"<tr><td>{m['gene']}</td><td>{float(m['D_within_logFC']):.3f} ({float(m['D_within_FDR']):.3g})</td>"
        f"<td>{float(m['DiD_logFC']):.3f}</td><td>{float(m['DiD_FDR']):.3g}</td></tr>" for m in mdid)
    did_tbl = ("<table><tr><th>Gene</th><th>Within-arm dupilumab logFC (FDR)</th>"
               "<th>DiD logFC</th><th>DiD FDR</th></tr>" + drows2 + "</table>")
    _def("t_did_note",
         f"Within-arm reversal alone does not prove drug specificity. The formal test is the "
         f"difference-in-differences (treatment \\u00d7 time) analysis on the 66 complete pairs: the type-2 "
         f"program-score change differs between arms (Wilcoxon rank-sum p = {float(pdid['wilcox_p']):.3g}; "
         f"Hodges\\u2013Lehmann estimate {float(pdid['HL_estimate_DminusP']):.3f}, 95% CI "
         f"[{float(pdid['HL_lo']):.3f}, {float(pdid['HL_hi']):.3f}]), and 17 genes reach FDR&lt;0.05 "
         f"genome-wide \\u2014 but none of the 8 pre-defined markers passes the interaction FDR "
         f"(POSTN and SERPINB2 are borderline), so drug specificity is supported at the program-score "
         f"level, not gene by gene. Exploratory (post-hoc) analysis; full tables in Downloads.",
         f"组内逆转本身不能证明药物特异性。正式检验是对 66 个完整配对做的差值中的差值（治疗\\u00d7时间）"
         f"分析：2 型程序评分的变化在两臂间不同（Wilcoxon 秩和 p = {float(pdid['wilcox_p']):.3g}；"
         f"Hodges\\u2013Lehmann 估计 {float(pdid['HL_estimate_DminusP']):.3f}，95% CI "
         f"[{float(pdid['HL_lo']):.3f}, {float(pdid['HL_hi']):.3f}]），全基因组 17 个基因达到 FDR&lt;0.05"
         f"——但 8 个预定义标志物无一通过交互 FDR（POSTN 与 SERPINB2 临界），故药物特异性在程序评分"
         f"层面得到支持，而非逐个基因。探索性（事后）分析；全表见数据下载。")
    did_block = (f"<h3 {dk('t_did_h')}>{en_of('t_did_h')}</h3>" + did_tbl
        + f"<p {dk('t_did_note')}>{en_of('t_did_note')}</p>")
    body = (DISCLAIMER''')
E('''        + f"<h3 {dk('t_inv_h')}>{en_of('t_inv_h')}</h3>" + inv_tbl''',
  '''        + did_block
        + f"<h3 {dk('t_inv_h')}>{en_of('t_inv_h')}</h3>" + inv_tbl''')

# --- downloads: new v1.6 items ----------------------------------------------
E('_def("nav_treatment", "Treatment response", "治疗反应")',
  '_def("nav_treatment", "Treatment response", "治疗反应")\n'
  '_def("dl_audit",\n'
  '     "Reproduction audit (release v1.6): 85 rows \\u2014 every extracted claim with its grade, the claimed proposition, the measured quantity and their match.",\n'
  '     "复现审计（v1.6）：85 行——每条提取声称及其判级、原命题、实际量测与命题匹配。")\n'
  '_def("dl_strict",\n'
  '     "Strict whole-tissue sensitivity: consensus gene list from the 7 whole-tissue cohorts (545 up / 225 down) with per-gene cohort support.",\n'
  '     "严格组织敏感性：7 个原组织队列的共识基因表（545 上调 / 225 下调）及逐基因队列支持。")\n'
  '_def("dl_inventory",\n'
  '     "Primary-contrast inventory: the 10 disease cohorts with design types and per-arm sample sizes.",\n'
  '     "主对比清单：10 个疾病队列的设计类型与每臂样本量。")\n'
  '_def("dl_did",\n'
  '     "Dupilumab difference-in-differences (exploratory): genome-wide interaction table, the 8 markers, and the program-score test.",\n'
  '     "Dupilumab 差值中的差值（探索性）：全基因组交互表、8 个标志物与程序评分检验。")')
E('''        ("treatment_response.json", "dl9_f10"),
    ]''',
  '''        ("treatment_response.json", "dl9_f10"),
        ("reproduction_audit.csv", "dl_audit"),
        ("strict_consensus.csv", "dl_strict"),
        ("primary_contrasts_inventory.csv", "dl_inventory"),
        ("did_de.csv", "dl_did"),
        ("markers_did.csv", "dl_did"),
        ("program_score_did.csv", "dl_did"),
    ]''')
E('''        (os.path.join(s9, "f1_meta.csv"), "meta_combined.csv"),''',
  '''        (os.path.join(s9, "meta_v2", "meta_extended.csv"), "meta_combined.csv"),''')
E('''    for src, dst in copies:
        shutil.copy(src, os.path.join(data_dir, dst))

    # step8 analysis-layer files''',
  '''    for src, dst in copies:
        shutil.copy(src, os.path.join(data_dir, dst))

    # v1.6 download files (audit v1.6, strict-tissue sensitivity, DiD)
    s6 = os.path.join(PROJ, "analysis", "step6_validation")
    sst = os.path.join(s2, "sensitivity_strict_tissue")
    sdid = os.path.join(s2, "did_dupilimab")
    for src, dst in [
        (os.path.join(s6, "reproduction_audit_v1.6.csv"), "reproduction_audit.csv"),
        (os.path.join(sst, "strict_consensus.csv"), "strict_consensus.csv"),
        (os.path.join(sst, "primary_contrasts_inventory.csv"), "primary_contrasts_inventory.csv"),
        (os.path.join(sdid, "did_de.csv"), "did_de.csv"),
        (os.path.join(sdid, "markers_did.csv"), "markers_did.csv"),
        (os.path.join(sdid, "program_score_did.csv"), "program_score_did.csv"),
    ]:
        shutil.copy(src, os.path.join(data_dir, dst))

    # step8 analysis-layer files''')

for old, new in EDITS:
    n = t.count(old)
    assert n == 1, f"count={n}: {old[:80]!r}"
    t = t.replace(old, new)
P.write_text(t, encoding="utf-8")
print(f"build.py patched: {len(EDITS)} edits")
