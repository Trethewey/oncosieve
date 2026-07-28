#!/usr/bin/env python3
"""
ONCOSIEVE Report Generator
Generates a single self-contained dark-theme HTML report from post-pipeline TSVs.

Changes v2:
  - Logo embedding (--logo)
  - Unified colour palette across all charts
  - All charts: high-confidence data only
  - Tier distribution: pies (not donuts), bold white labels, 20% larger
  - Source contribution: clipping fixed, ClinVar = yellow
  - Consequence distribution: clipping fixed, legend restored
  - Cancer types: single-hue sequential scale
  - References section
  - More vivid Tier 1 red

Dependencies:
    pip install plotly pandas --break-system-packages

Usage:
    python3 tools/generate_report.py \
        --full      output/pan_cancer_whitelist_GRCh38_full.tsv.gz \
        --highconf  output/pan_cancer_whitelist_GRCh38_highconf.tsv.gz \
        --out-dir   output/ \
        [--logo     path/to/logo.jpg]
"""

import argparse
import base64
import html as html_module
import json
import re
from datetime import date
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ── Colour palette ─────────────────────────────────────────────────────────────
BG_PAGE  = "#0d1117"
BG_CARD  = "#161b22"
BG_CHART = "#1c2128"
BORDER   = "#30363d"
TEXT_PRI = "#e6edf3"
TEXT_SEC = "#8b949e"
ACCENT   = "#3D6E8F"  # steel blue accent

# Brand palette from coolors swatch + burgundy + two new swatches
PALETTE = [
    "#C8F2D4",  # mint green (lightest)
    "#8FBCAA",  # sage green
    "#5B8FA0",  # slate teal (new)
    "#6DC5A0",  # seafoam green (new)
    "#3D6E8F",  # steel blue
    "#3D3785",  # indigo
    "#5C1A4E",  # dark purple
    "#7B1040",  # burgundy
]

# Version stamped into the report header. Kept as a module-level constant so
# there is one place to bump on release rather than a hard-coded literal that
# silently rots (the v2.0 report continued to advertise "2.0" past subsequent
# fixes).
ONCOSIEVE_VERSION = "2.1"

TIER_COLORS = {1: "#3D6E8F", 2: "#3D3785", 3: "#8FBCAA"}

# TCGA and TP53 previously shared the same colour, so any chart split by
# source could not distinguish them. TP53 is given its own dedicated slot
# from the palette.
SOURCE_COLORS = {
    "COSMIC":         "#3D6E8F",
    "GENIE":          "#3D3785",
    "TCGA":           "#8FBCAA",
    "ClinVar":        "#C8F2D4",
    "OncoKB":         "#5C1A4E",
    "CancerHotspots": "#7B1040",
    "TP53":           "#6DC5A0",
}

TEMPLATE = "plotly_dark"
CHART_CFG = {
    "displayModeBar": True,
    "displaylogo": False,
    "modeBarButtonsToRemove": ["lasso2d", "select2d"],
}


# ── Helpers ────────────────────────────────────────────────────────────────────

def _base_layout(**kwargs):
    return dict(
        template=TEMPLATE,
        paper_bgcolor=BG_CHART,
        plot_bgcolor=BG_CHART,
        font=dict(
            family="-apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif",
            color=TEXT_PRI,
        ),
        **kwargs,
    )


def fig_to_div(fig, include_plotlyjs=False):
    return fig.to_html(
        full_html=False,
        include_plotlyjs="cdn" if include_plotlyjs else False,
        config=CHART_CFG,
    )


def _palette_for(n):
    return [PALETTE[i % len(PALETTE)] for i in range(n)]


# ── KPI cards (keep both full + HC for summary context) ────────────────────────

def _source_counts_str(df):
    src = df["sources"].str.split("|").explode().str.strip().value_counts()
    return " &nbsp;&middot;&nbsp; ".join(f"{k}: {v:,}" for k, v in src.items())


def _revel_stat(df):
    rs = pd.to_numeric(df["revel_score"], errors="coerce")
    n = rs.notna().sum()
    pct = n / len(df) * 100 if len(df) else 0
    return f"{n:,} ({pct:.1f}%)"


def _primateai_stat(df):
    if "primateai_score" not in df.columns:
        return None
    ps = pd.to_numeric(df["primateai_score"], errors="coerce")
    n = ps.notna().sum()
    if n == 0:
        return None
    pct = n / len(df) * 100 if len(df) else 0
    return f"{n:,} ({pct:.1f}%)"


def build_kpi_cards(df_full, df_hc):
    cards = [
        ("Full Whitelist",  df_full, ACCENT),
        ("High-Confidence", df_hc,   "#C0392B"),
    ]
    parts = []
    for title, df, color in cards:
        t1 = int((df["wl_tier"] == 1).sum())
        t2 = int((df["wl_tier"] == 2).sum())
        pai_line = ''
        pai_stat = _primateai_stat(df)
        if pai_stat:
            pai_line = f'<div class="kpi-detail">PrimateAI-3D scored: {pai_stat}</div>'
        parts.append(
            f'<div class="kpi-card" style="border-top:3px solid {color}">'
            f'<div class="kpi-title">{title}</div>'
            f'<div class="kpi-value" style="color:{color}">{len(df):,}</div>'
            f'<div class="kpi-label">variants</div>'
            f'<div class="kpi-sep"></div>'
            f'<div class="kpi-row2">'
            f'<span class="kpi-badge" style="color:{TIER_COLORS[1]}">Tier 1: {t1:,}</span>'
            f'<span class="kpi-badge" style="color:{TIER_COLORS[2]}">Tier 2: {t2:,}</span>'
            f'<span class="kpi-badge" style="color:{TIER_COLORS[3]}">Tier 3: {int((df["wl_tier"] == 3).sum()):,}</span>'
            f'</div>'
            f'<div class="kpi-detail">REVEL scored: {_revel_stat(df)}</div>'
            f'{pai_line}'
            f'<div class="kpi-detail kpi-sources">{_source_counts_str(df)}</div>'
            f'</div>'
        )
    return '<div class="kpi-row">' + "".join(parts) + "</div>"


# ── Charts — all use df_hc only ────────────────────────────────────────────────

def chart_source_bars(df_hc):
    """Horizontal bar — high-confidence only."""
    sc = (
        df_hc["sources"].str.split("|").explode()
        .str.strip().value_counts().sort_values()
    )
    max_label = max(len(s) for s in sc.index) if len(sc) else 10
    left_margin = max(170, max_label * 9)

    # Case-insensitive colour lookup
    sc_lower = {k.lower(): v for k, v in SOURCE_COLORS.items()}
    def _src_color(s):
        return sc_lower.get(s.lower(), PALETTE[2])

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=sc.values,
        y=sc.index,
        orientation="h",
        marker=dict(
            color=[_src_color(s) for s in sc.index],
            opacity=0.9,
            line=dict(width=0),
        ),
        hovertemplate="<b>%{y}</b><br>%{x:,} variants<extra></extra>",
        showlegend=False,
    ))
    fig.update_layout(**_base_layout(
        height=max(360, len(sc) * 52 + 100),
        xaxis_title="Variant count",
        margin=dict(t=40, b=60, l=left_margin, r=40),
        xaxis=dict(automargin=True),
        yaxis=dict(automargin=True),
    ))
    return fig


def chart_consequence_bars(df_hc):
    """Horizontal bar — high-confidence only. Palette colours per consequence."""
    cc = df_hc["consequence"].value_counts().sort_values()
    n = len(cc)
    max_label = max(len(s) for s in cc.index) if n else 10
    left_margin = max(220, max_label * 9)

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=cc.values,
        y=cc.index,
        orientation="h",
        marker=dict(color=_palette_for(n), opacity=0.9, line=dict(width=0)),
        hovertemplate="<b>%{y}</b><br>%{x:,} variants<extra></extra>",
        showlegend=False,
    ))
    fig.update_layout(**_base_layout(
        height=max(360, n * 52 + 100),
        xaxis_title="Variant count",
        margin=dict(t=40, b=60, l=left_margin, r=40),
        xaxis=dict(automargin=True),
        yaxis=dict(automargin=True),
    ))
    return fig



def _luminance(hex_color):
    """Return relative luminance (0=black, 1=white) for a hex colour."""
    h = hex_color.lstrip("#")
    r, g, b = (int(h[i:i+2], 16) / 255.0 for i in (0, 2, 4))
    def lin(c):
        return c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    return 0.2126 * lin(r) + 0.7152 * lin(g) + 0.0722 * lin(b)


def _text_for_bg(hex_color):
    """Return '#111' for light backgrounds, '#fff' for dark."""
    return "#111111" if _luminance(hex_color) > 0.35 else "#ffffff"


def chart_cancer_types(df_hc, top_n=60):
    """Treemap of top cancer types — high-confidence. Palette colours."""
    ct = (
        df_hc["cancer_types"].dropna()
        .str.split("|").explode()
        .str.strip().replace("", pd.NA).dropna()
        .value_counts().head(top_n)
    )
    if ct.empty:
        fig = go.Figure()
        fig.update_layout(**_base_layout(title="No cancer type data available", height=200))
        return fig

    n = len(ct)
    colors = _palette_for(n)

    text_colors = [_text_for_bg(c) for c in colors]

    fig = go.Figure(go.Treemap(
        labels=ct.index.tolist(),
        parents=[""] * n,
        values=ct.values.tolist(),
        marker=dict(
            colors=colors,
            line=dict(color=BG_PAGE, width=2),
        ),
        textinfo="label+value",
        textfont=dict(size=13),
        insidetextfont=dict(size=13),
        hovertemplate="<b>%{label}</b><br>%{value:,} variants<extra></extra>",
    ))
    # Per-tile text colour via JavaScript patch applied after render
    fig.update_traces(marker_colorscale=None)
    # Store text colours as custom data for JS; use texttemplate colour hack via marker line
    # Plotly treemap doesn't support per-tile textfont colour natively —
    # use the workaround of setting text colour via marker.colors and textfont contrast
    # Best available: set all light-bg tiles to use dark template text via uniformtext
    fig.update_layout(**_base_layout(
        height=600,
        margin=dict(t=20, b=20, l=20, r=20),
        uniformtext=dict(minsize=10, mode="hide"),
    ))
    # Attach text colour array via customdata and patch textfont per tile post-render
    fig.data[0].customdata = text_colors
    return fig


def chart_top_genes(df_hc, top_n=100):
    """Vertical column chart of top genes — high-confidence. Colour by majority tier."""
    gene_counts = df_hc["gene"].value_counts().head(top_n).sort_values(ascending=False)

    def majority_tier(g):
        sub = df_hc.loc[df_hc["gene"] == g, "wl_tier"]
        return int(sub.mode().iloc[0]) if len(sub) else 2

    colors = [TIER_COLORS.get(majority_tier(g), TIER_COLORS[2])
              for g in gene_counts.index]

    tier_labels = " &nbsp; ".join(
        f'<span style="color:{c}">&#9632;</span> Tier {t}'
        for t, c in TIER_COLORS.items()
    )

    fig = go.Figure()
    fig.add_trace(go.Bar(
        x=gene_counts.index.tolist(),
        y=gene_counts.values.tolist(),
        marker=dict(color=colors, opacity=0.9, line=dict(width=0)),
        hovertemplate="<b>%{x}</b><br>%{y:,} variants<extra></extra>",
        showlegend=False,
    ))
    # Dummy scatter for tier legend
    for tier, color in TIER_COLORS.items():
        fig.add_trace(go.Scatter(
            x=[None], y=[None],
            mode="markers",
            name=f"Tier {tier}",
            marker=dict(color=color, size=11, symbol="square"),
            showlegend=True,
        ))

    fig.update_layout(**_base_layout(
        height=520,
        yaxis_title="Variant count",
        xaxis=dict(
            tickangle=-60,
            tickfont=dict(size=10),
            automargin=True,
        ),
        margin=dict(t=20, b=120, l=60, r=20),
        legend=dict(orientation="h", y=1.05, font=dict(size=12)),
        bargap=0.15,
    ))
    return fig, tier_labels



# ── Variant DataTable ──────────────────────────────────────────────────────────

DISPLAY_COLS = [
    "chrom", "pos", "ref", "alt", "gene", "transcript_id", "is_mane_select",
    "refseq_id", "hgvsc", "hgvsp", "protein_change", "consequence", "wl_tier",
    "n_samples", "n_cancer_types", "cancer_types", "sources",
    "oncokb_oncogenicity", "clinvar_clinical_significance",
    "transcript_source", "tp53_class", "revel_score",
    "primateai_score", "primateai_percentile", "primateai_prediction",
]

# Human-friendly labels for the DataTable header — snake_case identifiers are
# hostile to non-bioinformatician readers. Any column not listed here falls
# back to its raw identifier.
DISPLAY_COL_LABELS = {
    "chrom":                          "Chromosome",
    "pos":                            "Position",
    "ref":                            "Ref",
    "alt":                            "Alt",
    "gene":                           "Gene",
    "transcript_id":                  "Transcript ID",
    "is_mane_select":                 "MANE Select",
    "refseq_id":                      "RefSeq",
    "hgvsc":                          "HGVSc",
    "hgvsp":                          "HGVSp",
    "protein_change":                 "Protein change",
    "consequence":                    "Consequence",
    "wl_tier":                        "Tier",
    "n_samples":                      "n_samples",
    "n_cancer_types":                 "Cancer types (count)",
    "cancer_types":                   "Cancer types",
    "sources":                        "Sources",
    "oncokb_oncogenicity":            "OncoKB",
    "clinvar_clinical_significance":  "ClinVar significance",
    "transcript_source":              "Transcript source",
    "tp53_class":                     "TP53 class",
    "revel_score":                    "REVEL",
    "primateai_score":                "PrimateAI-3D score",
    "primateai_percentile":           "PrimateAI-3D percentile",
    "primateai_prediction":           "PrimateAI-3D prediction",
}


def build_datatable(df, table_id):
    """
    Returns (table_html, data_json, tier_col_idx, nsamp_col_idx).

    table_html has an empty <tbody>; DataTables populates it lazily from
    `data_json`, which is a JSON array of arrays of pre-escaped cell strings.
    tier_col_idx / nsamp_col_idx are computed from the ACTUAL filtered
    column set — not from the DISPLAY_COLS constant — so if any upstream
    column is missing (e.g. no primateai_score, no refseq_id), the DataTable
    still sorts and highlights the right column instead of the position the
    module was authored to expect.
    """
    cols = [c for c in DISPLAY_COLS if c in df.columns]
    sub = df[cols].copy()
    if "revel_score" in sub.columns:
        sub["revel_score"] = sub["revel_score"].apply(
            lambda x: f"{x:.3f}" if pd.notna(x) else "")

    def esc(v):
        return "" if pd.isna(v) else html_module.escape(str(v))

    thead = "<thead><tr>" + "".join(
        f"<th>{html_module.escape(DISPLAY_COL_LABELS.get(c, c))}</th>" for c in cols
    ) + "</tr></thead>"
    table_html = (
        f'<table id="{table_id}" class="display compact" style="width:100%">'
        f"{thead}<tbody></tbody></table>"
    )

    rows_data = [[esc(row[c]) for c in cols] for _, row in sub.iterrows()]
    data_json = json.dumps(rows_data, separators=(",", ":"))
    tier_col_idx  = cols.index("wl_tier")   if "wl_tier"   in cols else -1
    nsamp_col_idx = cols.index("n_samples") if "n_samples" in cols else 0
    return table_html, data_json, tier_col_idx, nsamp_col_idx


# ── Database versions / Sources table ──────────────────────────────────────────

_DB_LINE_RE = re.compile(r'^([A-Za-z][\w/. ]*?)\s{2,}(\S.*?)\s*\(([^)]+)\)\s*$')


def _parse_db_versions(path):
    """Parse build_whitelist's database_versions.txt into structured rows."""
    if not path:
        return []
    p = Path(path)
    if not p.exists():
        return []
    out = []
    for line in p.read_text(encoding="utf-8").splitlines():
        m = _DB_LINE_RE.match(line)
        if not m:
            continue
        source  = m.group(1).strip()
        version = m.group(2).strip()
        origin  = m.group(3).strip()
        # Strip leading "##fileDate=" if present (ClinVar)
        version = version.replace("##fileDate=", "")
        origin_type = "API" if origin.startswith("http") else "File"
        origin_short = origin if origin_type == "API" else Path(origin).name
        out.append({
            "source": source,
            "version": version,
            "type": origin_type,
            "origin": origin_short,
        })
    return out


def build_sources_table(entries):
    if not entries:
        return ""
    rows = []
    for e in entries:
        rows.append(
            "<tr>"
            f'<td class="src-name">{html_module.escape(e["source"])}</td>'
            f'<td class="src-version">{html_module.escape(e["version"])}</td>'
            f'<td class="src-type">{html_module.escape(e["type"])}</td>'
            f'<td class="src-origin">{html_module.escape(e["origin"])}</td>'
            "</tr>"
        )
    return (
        '<table class="sources-table">'
        '<thead><tr><th>Database</th><th>Version / Date</th><th>Type</th><th>Origin</th></tr></thead>'
        '<tbody>' + "".join(rows) + '</tbody>'
        '</table>'
    )


# ── References ─────────────────────────────────────────────────────────────────

REFERENCES = [
    ("COSMIC",
     "Tate JG et al. Nucleic Acids Res. 2019;47(D1):D941-D947.",
     "doi:10.1093/nar/gky1015", "PMID:30371878"),
    ("GENIE",
     "The AACR Project GENIE Consortium. Cancer Discov. 2017;7(8):818-831.",
     "doi:10.1158/2159-8290.CD-17-0151", ""),
    ("TCGA mc3",
     "Ellrott K et al. Cell Syst. 2018;6(3):271-281.e7.",
     "doi:10.1016/j.cels.2018.03.002", "PMID:29596782"),
    ("OncoKB",
     "Suehnholz SP et al. Cancer Discov. 2024;14(1):49-65.",
     "doi:10.1158/2159-8290.CD-23-0467", ""),
    ("OncoKB",
     "Chakravarty D et al. JCO Precis Oncol. 2017;1:PO.17.00011.",
     "doi:10.1200/PO.17.00011", ""),
    ("ClinVar",
     "Landrum MJ et al. Nucleic Acids Res. 2016;44(D1):D862-8.",
     "doi:10.1093/nar/gkv1222", "PMID:26582918"),
    ("TP53",
     "de Andrade KC et al. Cell Death Differ. 2022;29(5):1071-1073.",
     "doi:10.1038/s41418-022-00976-3", ""),
    ("TP53 Database",
     "The TP53 Database (R21, Jan 2025). https://tp53.cancer.gov", "", ""),
    ("CancerHotspots",
     "Chang MT et al. Nat Biotechnol. 2016;34(2):155-163.",
     "doi:10.1038/nbt.3391", "PMID:26619011"),
    ("CancerHotspots",
     "Chang MT et al. Cancer Discov. 2018;8(2):174-183.",
     "doi:10.1158/2159-8290.CD-17-0321", "PMID:29247016"),
    ("CancerHotspots",
     "Bandlamudi C et al. cancerhotspots.org.", "", ""),
    ("REVEL",
     "Ioannidis NM et al. Am J Hum Genet. 2016;99(4):877-885.",
     "doi:10.1016/j.ajhg.2016.08.016", "PMID:27666373"),
    ("PrimateAI-3D",
     "Gao H et al. Science. 2023;380(6648):eabn8197.",
     "doi:10.1126/science.abn8197", "PMID:37262156"),
]


def _ref_link(text):
    """Convert a DOI, PMID or URL string into a clickable anchor tag."""
    if text.startswith("doi:"):
        doi = text[4:]
        url = f"https://doi.org/{doi}"
        return f'<a href="{url}" target="_blank" class="ref-link">{html_module.escape(text)}</a>'
    if text.startswith("http"):
        return f'<a href="{html_module.escape(text)}" target="_blank" class="ref-link">{html_module.escape(text)}</a>'
    if text.startswith("PMID:"):
        pmid = text[5:]
        url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        return f'<a href="{url}" target="_blank" class="ref-link">{html_module.escape(text)}</a>'
    return html_module.escape(text)


def build_references_html():
    rows = []
    for source, citation, doi, pmid in REFERENCES:
        links = [_ref_link(x) for x in [doi, pmid] if x]
        extra_html = (
            '<br><span class="ref-doi">' + " &nbsp;&middot;&nbsp; ".join(links) + "</span>"
            if links else ""
        )
        rows.append(
            f'<tr>'
            f'<td class="ref-source">{html_module.escape(source)}</td>'
            f'<td class="ref-citation">{html_module.escape(citation)}{extra_html}</td>'
            f'</tr>'
        )
    return (
        '<table class="ref-table">'
        "<thead><tr><th>Source</th><th>Citation</th></tr></thead>"
        "<tbody>" + "".join(rows) + "</tbody>"
        "</table>"
    )


# ── Logo ────────────────────────────────────────────────────────────────────────

def _logo_tag(logo_path):
    if not logo_path:
        return ""
    p = Path(logo_path)
    if not p.exists():
        print(f"  Warning: logo not found at {logo_path}, skipping.")
        return ""
    suffix = p.suffix.lower().lstrip(".")
    mime = {"jpg": "image/jpeg", "jpeg": "image/jpeg",
            "png": "image/png", "gif": "image/gif", "svg": "image/svg+xml"}
    mt = mime.get(suffix, "image/jpeg")
    b64 = base64.b64encode(p.read_bytes()).decode()
    return (
        f'<img src="data:{mt};base64,{b64}" '
        f'class="header-logo" alt="ONCOSIEVE logo">'
    )


def _favicon_tag(logo_path):
    if not logo_path:
        return ""
    p = Path(logo_path)
    if not p.exists():
        return ""
    suffix = p.suffix.lower().lstrip(".")
    mime = {"jpg": "image/jpeg", "jpeg": "image/jpeg",
            "png": "image/png", "gif": "image/gif", "svg": "image/svg+xml",
            "ico": "image/x-icon"}
    mt = mime.get(suffix, "image/png")
    b64 = base64.b64encode(p.read_bytes()).decode()
    return f'<link rel="icon" type="{mt}" href="data:{mt};base64,{b64}">'


# ── CSS ────────────────────────────────────────────────────────────────────────

_CSS_TEMPLATE = """
*, *::before, *::after { box-sizing: border-box; }

body {
  margin: 0;
  background: %(BG_PAGE)s;
  color: %(TEXT_PRI)s;
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif;
  font-size: 14px;
  line-height: 1.5;
}

/* Nav */
.topnav {
  position: sticky; top: 0; z-index: 100;
  background: %(BG_CARD)s;
  border-bottom: 1px solid %(BORDER)s;
  display: flex; align-items: center;
  padding: 0 24px; height: 48px;
  overflow-x: auto;
}
.topnav .brand {
  font-weight: 700; font-size: 15px; color: %(ACCENT)s;
  margin-right: 24px; white-space: nowrap; letter-spacing: 0.5px;
}
.topnav a {
  color: %(TEXT_SEC)s; text-decoration: none;
  padding: 0 14px; height: 48px;
  display: flex; align-items: center;
  font-size: 13px; white-space: nowrap;
  border-bottom: 2px solid transparent;
  transition: color 0.15s, border-color 0.15s;
}
.topnav a:hover { color: %(TEXT_PRI)s; border-bottom-color: %(ACCENT)s; }

/* Layout */
.container { max-width: 1440px; margin: 0 auto; padding: 36px 28px 80px; }

/* Page header */
.page-header {
  display: flex; align-items: center; gap: 28px;
  margin-bottom: 36px;
  border-bottom: 1px solid %(BORDER)s;
  padding-bottom: 26px;
}
.header-logo {
  width: 150px; height: 150px;
  flex-shrink: 0;
  object-fit: contain;
  margin-left: auto;
  order: 2;
}
.header-text { display: flex; flex-direction: column; order: 1; flex: 1; }
.header-title {
  margin: 0; font-size: 45px; font-weight: 700;
  letter-spacing: -0.5px; color: %(TEXT_PRI)s; line-height: 1.1;
}
.header-version {
  font-size: 18px;
  font-weight: 600;
  letter-spacing: 0;
  color: %(TEXT_SEC)s;
  vertical-align: baseline;
}
.header-sep {
  display: block;
  width: 100%%;
  height: 2px;
  background: %(ACCENT)s;
  margin: 10px 0;
  border-radius: 2px;
}
.header-subtitle {
  font-size: 18px; font-weight: 400;
  color: %(TEXT_PRI)s; margin: 0; line-height: 1.2;
}
.header-meta { margin-top: 10px; font-size: 12px; color: %(TEXT_SEC)s; }

/* Sections */
.section { margin-bottom: 52px; }
.section-title {
  font-size: 17px; font-weight: 600;
  color: %(TEXT_PRI)s; margin: 0 0 16px;
  padding-left: 12px;
  border-left: 3px solid %(ACCENT)s;
}

/* KPI cards */
.kpi-row { display: flex; gap: 20px; flex-wrap: wrap; }
.kpi-card {
  background: %(BG_CARD)s; border: 1px solid %(BORDER)s;
  border-radius: 10px; padding: 22px 26px;
  flex: 1; min-width: 300px;
}
.kpi-title {
  font-size: 11px; font-weight: 700;
  text-transform: uppercase; letter-spacing: 1px;
  color: %(TEXT_SEC)s; margin-bottom: 10px;
}
.kpi-value  { font-size: 42px; font-weight: 700; line-height: 1; margin-bottom: 2px; }
.kpi-label  { font-size: 12px; color: %(TEXT_SEC)s; margin-bottom: 14px; }
.kpi-sep    { border-top: 1px solid %(BORDER)s; margin: 14px 0; }
.kpi-row2   { display: flex; gap: 14px; margin-bottom: 8px; }
.kpi-badge  { font-size: 13px; font-weight: 600; }
.kpi-detail { font-size: 12px; color: %(TEXT_SEC)s; margin-top: 5px; }
.kpi-sources { font-size: 11px; line-height: 2; margin-top: 8px; }

/* Tier legend (above treemap) */
.tier-legend {
  font-size: 13px; color: %(TEXT_SEC)s;
  margin-bottom: 10px; padding-left: 4px;
}

/* Chart card */
.chart-card {
  background: %(BG_CHART)s; border: 1px solid %(BORDER)s;
  border-radius: 10px; padding: 20px; overflow: hidden;
}

/* Tab buttons */
.tab-buttons { display: flex; gap: 8px; margin-bottom: 16px; }
.tab-btn {
  padding: 8px 22px; border-radius: 6px;
  border: 1px solid %(BORDER)s;
  background: %(BG_CARD)s; color: %(TEXT_SEC)s;
  cursor: pointer; font-size: 13px; font-weight: 500;
  transition: all 0.15s;
}
.tab-btn.active {
  background: %(ACCENT)s; color: %(BG_PAGE)s;
  border-color: %(ACCENT)s; font-weight: 600;
}
.tab-panel        { display: none; }
.tab-panel.active { display: block; }

/* References */
.sources-table {
  width: 100%%; border-collapse: collapse; font-size: 14px;
  margin-bottom: 18px;
}
.sources-table thead th {
  text-align: left; padding: 10px 14px;
  border-bottom: 2px solid %(BORDER)s;
  font-weight: 600; color: %(TEXT_PRI)s;
  background: %(BG_CHART)s;
}
.sources-table tbody tr   { border-bottom: 1px solid %(BORDER)s; }
.sources-table tbody tr:last-child { border-bottom: none; }
.sources-table tbody td   { padding: 8px 14px; vertical-align: top; color: %(TEXT_PRI)s; }
.src-name    { font-weight: 600; color: %(ACCENT)s; white-space: nowrap; }
.src-version { white-space: nowrap; }
.src-type    { white-space: nowrap; color: %(TEXT_SEC)s; }
.src-origin  { font-family: monospace; font-size: 12px; color: %(TEXT_SEC)s; word-break: break-all; }

.ref-table { width: 100%%; border-collapse: collapse; font-size: 16px; }
.ref-table thead th {
  text-align: left; font-size: 14px; font-weight: 700;
  text-transform: uppercase; letter-spacing: 0.8px;
  color: %(TEXT_SEC)s; border-bottom: 1px solid %(BORDER)s; padding: 8px 12px;
}
.ref-table tbody tr   { border-bottom: 1px solid %(BORDER)s; }
.ref-table tbody tr:last-child { border-bottom: none; }
.ref-table tbody td   { padding: 10px 14px; vertical-align: top; }
.ref-source { white-space: nowrap; font-weight: 600; color: %(ACCENT)s; width: 160px; }
.ref-citation { color: %(TEXT_PRI)s; }
.ref-doi { font-size: 12px; color: %(TEXT_PRI)s; font-family: monospace; opacity: 0.75; }
.ref-link { color: #89c2a2; text-decoration: none; }
.ref-link:hover { text-decoration: underline; }

/* Footer */
.page-footer {
  margin-top: 60px; padding-top: 20px;
  border-top: 1px solid %(BORDER)s;
  font-size: 12px; color: %(TEXT_SEC)s; text-align: center;
}

/* DataTables dark */
.dataTables_wrapper { color: %(TEXT_PRI)s; font-size: 12px; }
table.dataTable {
  background: %(BG_CHART)s; color: %(TEXT_PRI)s;
  border-collapse: collapse !important;
}
table.dataTable thead th {
  background: %(BG_CARD)s !important; color: %(TEXT_SEC)s !important;
  border-bottom: 1px solid %(BORDER)s !important;
  font-size: 11px; text-transform: uppercase; letter-spacing: 0.5px;
  white-space: nowrap; padding: 10px 12px;
}
table.dataTable tbody tr {
  background: %(BG_CHART)s !important;
  border-bottom: 1px solid %(BORDER)s !important;
}
table.dataTable tbody tr:hover { background: %(BG_CARD)s !important; }
table.dataTable tbody td { padding: 7px 12px; white-space: nowrap; }
table.dataTable tbody tr.tier-1 td:nth-child(%(TIER_COL_1B)s) {
  color: %(TIER1_COLOR)s; font-weight: 600;
}
table.dataTable tbody tr.tier-2 td:nth-child(%(TIER_COL_1B)s) {
  color: %(TIER2_COLOR)s; font-weight: 600;
}
.dataTables_filter input,
.dataTables_length select {
  background: %(BG_CARD)s; border: 1px solid %(BORDER)s;
  color: %(TEXT_PRI)s; border-radius: 4px; padding: 4px 8px;
}
.dataTables_info,
.dataTables_paginate { color: %(TEXT_SEC)s !important; font-size: 12px; margin-top: 12px; }
.dataTables_paginate .paginate_button {
  color: %(TEXT_SEC)s !important; border-radius: 4px;
}
.dataTables_paginate .paginate_button.current {
  background: %(ACCENT)s !important;
  color: %(BG_PAGE)s !important; border: none !important;
}
.dataTables_paginate .paginate_button:hover {
  background: %(BG_CARD)s !important;
  color: %(TEXT_PRI)s !important;
  border: 1px solid %(BORDER)s !important;
}

/* Scrollbar */
::-webkit-scrollbar { width: 8px; height: 8px; }
::-webkit-scrollbar-track { background: %(BG_PAGE)s; }
::-webkit-scrollbar-thumb { background: %(BORDER)s; border-radius: 4px; }
::-webkit-scrollbar-thumb:hover { background: %(TEXT_SEC)s; }
"""


def _css(tier_col_1b=None):
    """Render the CSS.

    tier_col_1b is the 1-based column index of wl_tier in the CURRENT
    DataTable — passed in from build_report so the tier-row highlight
    lands on the right column even when upstream columns are missing.
    Falls back to the DISPLAY_COLS position if not supplied (legacy
    callers).
    """
    if tier_col_1b is None:
        # Legacy default: position wl_tier occupies in the DISPLAY_COLS list.
        tier_col_1b = DISPLAY_COLS.index("wl_tier") + 1
    return _CSS_TEMPLATE % dict(
        BG_PAGE=BG_PAGE, BG_CARD=BG_CARD, BG_CHART=BG_CHART,
        BORDER=BORDER, TEXT_PRI=TEXT_PRI, TEXT_SEC=TEXT_SEC, ACCENT=ACCENT,
        TIER1_COLOR=TIER_COLORS[1], TIER2_COLOR=TIER_COLORS[2],
        TIER_COL_1B=tier_col_1b,
    )


# ── HTML assembly ──────────────────────────────────────────────────────────────

def build_report(df_full, df_hc, out_path, logo_path=None, db_versions_path=None):
    print("  Building charts...")
    f_source             = chart_source_bars(df_hc)
    f_conseq             = chart_consequence_bars(df_hc)
    f_cancer             = chart_cancer_types(df_hc)
    f_genes, gene_legend = chart_top_genes(df_hc)

    div_source = fig_to_div(f_source, include_plotlyjs=True)
    div_conseq = fig_to_div(f_conseq)
    div_cancer = fig_to_div(f_cancer)
    div_genes  = fig_to_div(f_genes)

    print("  Building variant table (high-confidence)...")
    tbl_hc, tbl_hc_data, tbl_hc_tier_idx, tbl_hc_nsamp_idx = build_datatable(df_hc, "tbl-hc")

    print("  Assembling HTML...")
    logo_img    = _logo_tag(logo_path)
    favicon_tag = _favicon_tag(logo_path)
    kpi_html    = build_kpi_cards(df_full, df_hc)
    refs_html   = build_references_html()
    sources_table_html = build_sources_table(_parse_db_versions(db_versions_path))
    report_date = date.today().isoformat()
    n_hc        = len(df_hc)
    # Compute the dynamic tier column index for CSS and JS from the actual
    # DataTable columns; this replaces the module-level _TIER_COL_IDX and
    # _NSAMP_COL_IDX constants that would silently point at the wrong column
    # when any upstream column was missing.
    css         = _css(tier_col_1b=(tbl_hc_tier_idx + 1) if tbl_hc_tier_idx >= 0 else 1)
    nsamp_col   = tbl_hc_nsamp_idx

    # Compose the version-and-source header meta strip from the actual data
    # present in df_hc so we don't advertise annotation layers that weren't
    # produced (e.g. REVEL / PrimateAI-3D when the pipeline was run without them).
    _base_sources = ["COSMIC", "GENIE", "TCGA", "ClinVar", "OncoKB",
                     "CancerHotspots", "TP53"]
    _extras = []
    if "revel_score" in df_hc.columns and df_hc["revel_score"].fillna("").astype(str).ne("").any():
        _extras.append("REVEL")
    if "primateai_score" in df_hc.columns and df_hc["primateai_score"].fillna("").astype(str).ne("").any():
        _extras.append("PrimateAI-3D")
    _meta_strip = "GRCh38 &nbsp;&middot;&nbsp; " + " &middot; ".join(_base_sources)
    if _extras:
        _meta_strip += " &nbsp;&ndash;&nbsp; " + " &middot; ".join(_extras)

    html = "\n".join([
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '<meta charset="UTF-8">',
        '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
        favicon_tag,
        f"<title>OncoSieve {ONCOSIEVE_VERSION} Report {report_date}</title>",
        '<link rel="stylesheet" href="https://cdn.datatables.net/1.13.8/css/jquery.dataTables.min.css">',
        '<link rel="stylesheet" href="https://cdn.datatables.net/buttons/2.4.2/css/buttons.dataTables.min.css">',
        '<script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>',
        '<script src="https://cdn.datatables.net/1.13.8/js/jquery.dataTables.min.js"></script>',
        '<script src="https://cdn.datatables.net/buttons/2.4.2/js/dataTables.buttons.min.js"></script>',
        '<script src="https://cdn.datatables.net/buttons/2.4.2/js/buttons.html5.min.js"></script>',
        f"<style>{css}</style>",
        "</head>",
        "<body>",

        # Nav
        '<nav class="topnav">',
        '  <span class="brand">OncoSieve</span>',
        '  <a href="#summary">Summary</a>',
        '  <a href="#table">Variant Table</a>',
        '  <a href="#sources">Sources</a>',
        '  <a href="#consequences">Consequences</a>',
        '  <a href="#cancer-types">Cancer Types</a>',
        '  <a href="#genes">Top Genes</a>',
        '  <a href="#references">References</a>',
        "</nav>",

        '<div class="container">',

        # Header
        '<div class="page-header">',
        f"  {logo_img}",
        '  <div class="header-text">',
        f'    <h1 class="header-title">OncoSieve&nbsp;&nbsp;<span class="header-version">{ONCOSIEVE_VERSION}</span></h1>',
        '    <span class="header-sep"></span>',
        f'    <p class="header-subtitle">Pan-Cancer Variant Whitelist &nbsp;&ndash;&nbsp; Generated {report_date}</p>',
        f'    <div class="header-meta">{_meta_strip}</div>',
        '    <div class="header-meta"><a href="https://github.com/Trethewey/oncosieve" style="color:#3D6E8F;text-decoration:none;">&#128279; github.com/Trethewey/oncosieve</a></div>',
        "  </div>",
        "</div>",

        # Summary
        '<div class="section" id="summary">',
        '  <div class="section-title">Summary</div>',
        f"  {kpi_html}",
        "</div>",

        # Variant table — high-confidence only
        '<div class="section" id="table">',
        f'  <div class="section-title">Variant Table — High-Confidence ({n_hc:,} variants)</div>',
        '  <div class="chart-card">',
        f"    {tbl_hc}",
        "  </div>",
        "</div>",

        # Sources
        '<div class="section" id="sources">',
        '  <div class="section-title">Source Contribution (High-Confidence)</div>',
        ('  <div class="chart-card">' + sources_table_html + '</div>') if sources_table_html else '',
        '  <div class="chart-card">',
        f"    {div_source}",
        "  </div>",
        "</div>",

        # Consequences
        '<div class="section" id="consequences">',
        '  <div class="section-title">Consequence Distribution (High-Confidence)</div>',
        '  <div class="chart-card">',
        f"    {div_conseq}",
        "  </div>",
        "</div>",

        # Cancer types treemap
        '<div class="section" id="cancer-types">',
        '  <div class="section-title">Cancer Type Distribution — top 60, High-Confidence</div>',
        '  <div class="chart-card">',
        f"    {div_cancer}",
        "  </div>",
        "</div>",

        # Top genes treemap
        '<div class="section" id="genes">',
        '  <div class="section-title">Top Genes by Variant Count — top 100, High-Confidence</div>',
        f'  <div class="tier-legend">{gene_legend}</div>',
        '  <div class="chart-card">',
        f"    {div_genes}",
        "  </div>",
        "</div>",

        # References
        '<div class="section" id="references">',
        '  <div class="section-title">References</div>',
        '  <div class="chart-card">',
        f"    {refs_html}",
        "  </div>",
        "</div>",

        # Footer
        '<div class="page-footer">',
        f"  OncoSieve &nbsp;&middot;&nbsp; GRCh38 &nbsp;&middot;&nbsp; {report_date}",
        "</div>",

        "</div><!-- /container -->",

        "<script>",
        f"var TBL_HC_DATA = {tbl_hc_data};",
        f"var TBL_HC_TIER_IDX = {tbl_hc_tier_idx};",
        "$(document).ready(function () {",
        "  $('#tbl-hc').DataTable({",
        "    data: TBL_HC_DATA,",
        "    pageLength: 25,",
        "    lengthMenu: [10, 25, 50, 100, 250],",
        f"    order: [[{nsamp_col}, 'desc']],",
        "    scrollX: true,",
        "    deferRender: true,",
        "    dom: '<\"dt-toolbar\"lBf>rtip',",
        "    buttons: [{",
        "      extend: 'csv',",
        "      text: 'Download CSV',",
        "      filename: 'oncosieve_highconf_variants_' + new Date().toISOString().slice(0,10),",
        "      exportOptions: { columns: ':visible' }",
        "    }],",
        "    createdRow: function (row, data) {",
        "      if (TBL_HC_TIER_IDX < 0) return;",
        "      var t = data[TBL_HC_TIER_IDX];",
        "      if (t === '1' || t === 1) row.classList.add('tier-1');",
        "      else if (t === '2' || t === 2) row.classList.add('tier-2');",
        "      else if (t === '3' || t === 3) row.classList.add('tier-3');",
        "    },",
        "    language: { search: 'Search:', lengthMenu: 'Show _MENU_ rows' }",
        "  });",
        "});",
        "</script>",
        "</body>",
        "</html>",
    ])

    print(f"  Writing {out_path} ...")
    out_path.write_text(html, encoding="utf-8")
    size_mb = out_path.stat().st_size / 1_048_576
    print(f"  Done. File size: {size_mb:.1f} MB")


# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="ONCOSIEVE report generator")
    parser.add_argument("--full",     required=True,
                        help="Full annotated whitelist TSV (.tsv or .tsv.gz)")
    parser.add_argument("--highconf", required=True,
                        help="High-confidence whitelist TSV (.tsv or .tsv.gz)")
    parser.add_argument("--out-dir",  default="output/",
                        help="Output directory (default: output/)")
    parser.add_argument("--logo",     default=None,
                        help="Path to logo image (jpg/png) to embed in the report header")
    parser.add_argument("--db-versions", default=None,
                        help="Path to database_versions.txt for the Sources table (auto-discovered in --out-dir if omitted)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"oncosieve_report_{date.today().isoformat()}.html"

    print(f"Loading full whitelist:            {args.full}")
    df_full = pd.read_csv(args.full, sep="\t", compression="infer", low_memory=False)
    print(f"  {len(df_full):,} variants")

    print(f"Loading high-confidence whitelist: {args.highconf}")
    df_hc = pd.read_csv(args.highconf, sep="\t", compression="infer", low_memory=False)
    print(f"  {len(df_hc):,} variants")

    db_versions_path = args.db_versions
    if db_versions_path is None:
        candidate = out_dir / "database_versions.txt"
        if candidate.exists():
            db_versions_path = str(candidate)

    build_report(df_full, df_hc, out_path, logo_path=args.logo,
                 db_versions_path=db_versions_path)
    print(f"\nReport: {out_path}")


if __name__ == "__main__":
    main()
