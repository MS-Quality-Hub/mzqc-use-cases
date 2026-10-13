import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.image as mpimg
import seaborn as sns
import plotly.graph_objects as go
import warnings

warnings.filterwarnings('ignore')

#path to data
CSV_PATH = "/Users/fateme/Desktop/test_metadata/QC/Paper/metabolomics_repository_qc.csv"

#output path
RESULT_DIR = "./result_figure3"

SANKEY_DIR = os.path.join(RESULT_DIR, "sankey")
VIOLIN_DIR = os.path.join(RESULT_DIR, "violin")
FINAL_PATH = os.path.join(RESULT_DIR, "figure3_combined.png")


#panel b

def parse_repository_origin(mri_series):
    repo_id = mri_series.astype(str).str.split(':', n=2).str[1]
    origin = pd.Series('Unknown', index=mri_series.index)
    origin[repo_id.str.startswith('MSV', na=False)] = 'GNPS/MassIVE'
    origin[repo_id.str.startswith('GNPS', na=False)] = 'GNPS/MassIVE'
    origin[repo_id.str.startswith('MTBLS', na=False)] = 'MetaboLights'
    origin[repo_id.str.startswith('ST', na=False)] = 'Metabolomics Workbench'
    return origin


def _load_and_prepare(csv_path, metrics):
    """load data, parse origins and return a clean DataFrame"""
    all_cols = [m[1] for m in metrics]

    print("Loading data...")
    df = pd.read_csv(csv_path, usecols=['mri'] + all_cols, low_memory=False)

    for col in all_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    df = df.dropna(subset=['mri']).reset_index(drop=True)
    print(f"Loaded {len(df):,} samples")

    print("Parsing repository origins...")
    df['Repository'] = parse_repository_origin(df['mri'])
    repo_counts = df['Repository'].value_counts()
    print("Repository composition:")
    for repo, count in repo_counts.items():
        print(f"  {repo:25s}: {count:>8,} ({100*count/len(df):5.1f}%)")

    df = df[df['Repository'] != 'Unknown'].reset_index(drop=True)
    return df


def _draw_violin_on_ax(ax, df, col_name, ylabel, use_log, repo_order, repo_colors):
    """draw single violin panel"""

    plot_df = df[['Repository', col_name]].dropna(subset=[col_name]).copy()
    plot_df = plot_df[plot_df['Repository'].isin(repo_order)]

    if use_log:
        plot_df = plot_df[plot_df[col_name] > 0].copy()

    p99 = plot_df[col_name].quantile(0.99)
    p01 = plot_df[col_name].quantile(0.01)
    plot_trimmed = plot_df[
        (plot_df[col_name] <= p99) & (plot_df[col_name] >= p01)
    ].copy()

    #lighter RGB colors for the violin fills
    def lighten_hex(hex_color, amount=0.60):
        h = hex_color.lstrip('#')
        rgb = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
        return tuple(channel + (1 - channel) * amount for channel in rgb)

    violin_colors = {
        repo: lighten_hex(color, amount=0.60)
        for repo, color in repo_colors.items()
    }

    first_violin_collection = len(ax.collections)

    sns.violinplot(
        data=plot_trimmed,
        x='Repository',
        y=col_name,
        order=repo_order,
        palette=violin_colors,
        inner=None,
        cut=0,
        bw_adjust=3,
        linewidth=1.2,
        alpha=1.0,
        ax=ax,
    )

    violin_collections = ax.collections[first_violin_collection:]
    for art, repo in zip(violin_collections, repo_order):
        art.set_facecolor(violin_colors[repo])
        art.set_edgecolor(repo_colors[repo])
        art.set_alpha(1.0)
        art.set_linewidth(1.2)
        art.set_clip_on(True)
        art.set_clip_path(ax.patch)

    sns.boxplot(
        data=plot_trimmed,
        x='Repository',
        y=col_name,
        order=repo_order,
        palette=repo_colors,
        width=0.15,
        showfliers=False,
        linewidth=0.8,
        ax=ax,
        boxprops=dict(alpha=0.7),
        medianprops=dict(color='white', linewidth=1.5),
    )

    if use_log:
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f'{x:,.0f}' if x >= 1 else f'{x:g}'
        ))
    else:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f'{x:,.0f}' if abs(x) >= 1000 else f'{x:g}'
        ))

    ax.set_title('')
    ax.set_xlabel('')
    ax.set_ylabel(ylabel, fontsize=26)

    #two line labels
    display_labels = [
        'GNPS/\nMassIVE',
        'MetaboLights',
        'Metabolomics\nWorkbench',
    ]
    ax.set_xticklabels(display_labels, fontsize=23, rotation=0, ha='center', va='top')
    ax.tick_params(axis='x', length=4, pad=8)

    ax.tick_params(axis='y', labelsize=23)
    ax.grid(axis='y', alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)

    for line in ax.lines:
        line.set_clip_on(True)
        line.set_clip_path(ax.patch)

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.9)
        spine.set_color('#444444')


def fig_b_violins(csv_path, result_dir):
    """save each metric as a separate figure & produce one combined figure"""

    os.makedirs(result_dir, exist_ok=True)

    metrics = [
        ('MS1 Min m/z',          'MS1_Min_MZ',              'MS1 min m/z',         False),
        ('MS1 Max m/z',          'MS1_Max_MZ',              'MS1 max m/z',          False),
        ('RT Range (min)',       'RT_Range_in_Min',         'RT range (min)',        True),
        ('Unique Precursor m/z', 'Num_Unique_Precursor_MZ', 'Unique precursor m/z', True),
    ]

    repo_order = ['GNPS/MassIVE', 'MetaboLights', 'Metabolomics Workbench']
    repo_colors = {
        'GNPS/MassIVE':           '#6da7de',
        'MetaboLights':           "#9e0059",
        'Metabolomics Workbench': "#dee000",
    }

    df = _load_and_prepare(csv_path, metrics)

    #individual plots
    for (title, col_name, ylabel, use_log) in metrics:
        print(f"\nPlotting individual: {title}")

        fig, ax = plt.subplots(figsize=(7, 7))
        fig.subplots_adjust(left=0.18, right=0.95, top=0.92, bottom=0.28)

        _draw_violin_on_ax(ax, df, col_name, ylabel, use_log, repo_order, repo_colors)

        fig.patch.set_linewidth(0)
        fig.patch.set_edgecolor('none')

        safe_name = col_name.lower().replace(' ', '_')
        plot_path = os.path.join(result_dir, f'fig_b_{safe_name}.png')
        plt.savefig(plot_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
        print(f"  -> Saved: {plot_path}")
        plt.close('all')

    #combined figure -> all 4 panels side by side 
    print("\nPlotting combined figure...")
    fig, axes = plt.subplots(
        1, 4,
        figsize=(32, 8),
        gridspec_kw={'wspace': 0.22},
    )
    fig.subplots_adjust(left=0.06, right=0.98, top=0.92, bottom=0.28)

    for ax, (title, col_name, ylabel, use_log) in zip(axes, metrics):
        _draw_violin_on_ax(ax, df, col_name, ylabel, use_log, repo_order, repo_colors)

    fig.patch.set_linewidth(0)
    fig.patch.set_edgecolor('none')

    combined_path = os.path.join(result_dir, 'fig_b_combined.png')
    plt.savefig(combined_path, dpi=300, bbox_inches='tight', facecolor='white', edgecolor='none')
    print(f"  -> Saved: {combined_path}")
    plt.close('all')

    print("\nAll violin plots saved.")
    return combined_path

#panel a

def fig_sankey_triage(csv_path, result_dir):
    """sankey triage"""

    os.makedirs(result_dir, exist_ok=True)

    count_cols = [
        'MS1_pos_count', 'MS1_neg_count',
        'MS2_pos_count', 'MS2_neg_count',
        'MS3+_pos_count', 'MS3+_neg_count',
    ]

    print("Loading data...")
    df = pd.read_csv(csv_path, sep=None, engine='python')
    df.columns = [c.strip() for c in df.columns]

    for c in count_cols:
        if c not in df.columns:
            df[c] = 0
        df[c] = pd.to_numeric(df[c], errors='coerce').fillna(0)

    n_loaded = len(df)
    print(f"Files loaded: {n_loaded:,}")

    pos = (
        df['MS1_pos_count']
        + df['MS2_pos_count']
        + df['MS3+_pos_count']
    )

    neg = (
        df['MS1_neg_count']
        + df['MS2_neg_count']
        + df['MS3+_neg_count']
    )

    #exclude empty files (no scans in any MS column
    keep = (pos > 0) | (neg > 0)

    n_excluded = int((~keep).sum())

    df = df[keep].reset_index(drop=True)
    pos = pos[keep].reset_index(drop=True)
    neg = neg[keep].reset_index(drop=True)

    n_total = len(df)

    print(f"Excluded {n_excluded:,} empty files (no scans)")
    print(f"Total files in diagram: {n_total:,}")

    #polarity
    polarity = np.where(
        (pos > 0) & (neg > 0),
        'Mixed polarity',
        np.where(
            pos > 0,
            'Positive',
            'Negative'
        )
    )

    #MS depth/fragmentation
    has_ms3 = (
        (df['MS3+_pos_count'] > 0)
        | (df['MS3+_neg_count'] > 0)
    )

    has_ms2 = (
        (df['MS2_pos_count'] > 0)
        | (df['MS2_neg_count'] > 0)
    )

    ms_depth = np.where(
        has_ms3,
        'MS3+ (deep)',
        np.where(
            has_ms2,
            'MS1 + MS2 (fragmentation)',
            'MS1 only (survey)'
        )
    )

    work = pd.DataFrame({
        'polarity': polarity,
        'ms_depth': ms_depth
    })

    root_color = '#34495e'

    polarity_order = [
        'Positive',
        'Negative',
        'Mixed polarity'
    ]

    depth_order = [
        'MS1 only (survey)',
        'MS1 + MS2 (fragmentation)',
        'MS3+ (deep)'
    ]

    polarity_present = [
        x for x in polarity_order
        if (work['polarity'] == x).any()
    ]

    depth_present = [
        x for x in depth_order
        if (work['ms_depth'] == x).any()
    ]

    pol_counts = work['polarity'].value_counts().to_dict()
    dep_counts = work['ms_depth'].value_counts().to_dict()

    #colors
    node_fill_color_map = {
        'Positive':                    '#6da7de',
        'Negative':                    '#9e0059', 
        'Mixed polarity':              '#63c5b5',  
        'MS1 only (survey)':           '#6da7de',
        'MS1 + MS2 (fragmentation)':   '#dee000', 
        'MS3+ (deep)':                 "#eb861e", 
    }

    #label text colors 
    node_text_color_map = {
        'Positive':                    "#17436f",
        'Negative':                    "#7e0449",
        'Mixed polarity':              "#33baa4",
        'MS1 only (survey)':           '#17436f',
        'MS1 + MS2 (fragmentation)':   "#8a8c00",
        'MS3+ (deep)':                 "#ea7a0a",  
    }

    def colored_label(text, hex_color):
        return (
            f'<span style="color:{hex_color};'
            f'font-family:{FONT}">{text}</span>'
        )

    #root node
    nodes = [
        colored_label(
            f'All raw files<br>{n_total:,}',
            root_color
        )
    ]

    node_idx = {'__root__': 0}
    node_colors = [root_color]

    def add(key, display, fill_color, text_color):
        node_idx[key] = len(nodes)
        nodes.append(colored_label(display, text_color))
        node_colors.append(fill_color)

    #add polarity nodes
    for p in polarity_present:
        fill = node_fill_color_map.get(p, '#888')
        text = node_text_color_map.get(p, '#333')

        add(
            ('pol', p),
            f'{p}<br>{pol_counts[p]:,}',
            fill,
            text
        )

    #add MS depth nodes
    for d in depth_present:
        fill = node_fill_color_map.get(d, '#888')
        text = node_text_color_map.get(d, '#333')

        add(
            ('dep', d),
            f'{d}<br>{dep_counts[d]:,}',
            fill,
            text
        )

    def rgba(hex_color, alpha=0.40):
        h = hex_color.lstrip('#')

        r, g, b = (
            int(h[i:i + 2], 16)
            for i in (0, 2, 4)
        )

        return f'rgba({r},{g},{b},{alpha})'

    src = []
    tgt = []
    val = []
    lcol = []

    #link: root -> polarity 
    for p in polarity_present:
        c = int((work['polarity'] == p).sum())

        if c:
            src.append(0)
            tgt.append(node_idx[('pol', p)])
            val.append(c)

            lcol.append(
                rgba(
                    node_fill_color_map.get(p, '#888')
                )
            )

    #link: polarity -> MS depth 
    for p in polarity_present:
        for d in depth_present:

            c = int(
                (
                    (work['polarity'] == p)
                    & (work['ms_depth'] == d)
                ).sum()
            )

            if c:
                src.append(node_idx[('pol', p)])
                tgt.append(node_idx[('dep', d)])
                val.append(c)

                lcol.append(
                    rgba(
                        node_fill_color_map.get(d, '#888')
                    )
                )

    print("Building Sankey...")

    fig = go.Figure(
        go.Sankey(
            arrangement='snap',

            node=dict(
                label=nodes,
                color=node_colors,
                pad=40,
                thickness=24,
                line=dict(
                    color='white',
                    width=1.2
                ),
                hovertemplate='%{label}<extra></extra>',
            ),

            link=dict(
                source=src,
                target=tgt,
                value=val,
                color=lcol,
                hovertemplate='%{value:,} files<extra></extra>',
            ),

            textfont=dict(
                size=18,
                family=FONT
            ),
        )
    )

    fig.update_layout(
        annotations=[
            dict(
                text='Polarity',
                x=0.5,
                y=1.0385,
                showarrow=False,
                xref='paper',
                yref='paper',
                font=dict(
                    size=20,
                    color='#000000',
                    family=FONT
                ),
            ),
            dict(
                text='MS depth',
                x=1.015,
                y=1.0385,
                showarrow=False,
                xref='paper',
                yref='paper',
                font=dict(
                    size=20,
                    color='#000000',
                    family=FONT
                ),
            ),
        ],

        font=dict(
            family=FONT,
            size=18,
            color='black'
        ),

        paper_bgcolor='white',
        plot_bgcolor='white',

        margin=dict(
            l=30,
            r=40,
            t=100,
            b=90
        ),

        width=1400,
        height=900,
    )

    html_path = os.path.join(
        result_dir,
        'sankey_triage.html'
    )

    fig.write_html(html_path)

    print(f"  -> Saved: {html_path}")

    png_path = os.path.join(
        result_dir,
        'sankey_triage.png'
    )

    try:
        fig.write_image(
            png_path,
            scale=2
        )

        print(f"  -> Saved: {png_path}")

    except Exception as e:
        print(
            f"  [PNG export skipped: {e}. "
            f"Open the HTML instead.]"
        )
        png_path = None

    #save summary report
    report_path = os.path.join(
        result_dir,
        'sankey_triage_summary.txt'
    )

    with open(report_path, 'w') as f:

        f.write(
            "SANKEY TRIAGE SUMMARY\n"
            + "=" * 50
            + "\n\n"
        )

        f.write(
            f"Files loaded:            {n_loaded:,}\n"
        )

        f.write(
            f"Empty files removed:     {n_excluded:,}\n"
        )

        f.write(
            f"Total files in diagram:  {n_total:,}\n\n"
        )

        f.write("Stage 1 -> Polarity:\n")

        for p in polarity_present:
            f.write(
                f"  {p:28s}: "
                f"{pol_counts[p]:>8,} "
                f"({100 * pol_counts[p] / n_total:5.1f}%)\n"
            )

        f.write("\nStage 2 -> MS depth:\n")

        for d in depth_present:
            f.write(
                f"  {d:28s}: "
                f"{dep_counts[d]:>8,} "
                f"({100 * dep_counts[d] / n_total:5.1f}%)\n"
            )

    print(f"  -> Saved: {report_path}")

    return png_path

#combine panels 

LABEL_FONTSIZE = 18
LABEL_FONTWEIGHT = "bold"
DPI = 300


def combine_panels(sankey_path, violin_path, out_path):
    if not os.path.exists(sankey_path):
        raise FileNotFoundError(f"Panel a image not found: {sankey_path}")
    if not os.path.exists(violin_path):
        raise FileNotFoundError(f"Panel b image not found: {violin_path}")

    img_a = mpimg.imread(sankey_path)
    img_b = mpimg.imread(violin_path)

    h_a, w_a = img_a.shape[0], img_a.shape[1]
    h_b, w_b = img_b.shape[0], img_b.shape[1]

    #fix the output width so neither image stretched or distorted
    fig_width_in = 16.0
    row_a_height_in = fig_width_in * (h_a / w_a)
    row_b_height_in = fig_width_in * (h_b / w_b)

    fig_height_in = row_a_height_in + row_b_height_in

    fig = plt.figure(figsize=(fig_width_in, fig_height_in))
    gs = fig.add_gridspec(
        2, 1,
        height_ratios=[row_a_height_in, row_b_height_in],
        hspace=-0.15,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.imshow(img_a)
    ax_a.axis("off")
    ax_a.text(
        0.01, 0.92, "a",
        transform=ax_a.transAxes,
        fontsize=LABEL_FONTSIZE, fontweight=LABEL_FONTWEIGHT,
        va="top", ha="right",
    )

    ax_b = fig.add_subplot(gs[1, 0])
    ax_b.imshow(img_b)
    ax_b.axis("off")
    ax_b.text(
        0.01, 1.1, "b",
        transform=ax_b.transAxes,
        fontsize=LABEL_FONTSIZE, fontweight=LABEL_FONTWEIGHT,
        va="top", ha="right",
    )

    fig.savefig(out_path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"-> Saved combined figure: {out_path}")


def main():
    print("=" * 70)
    print("STEP 1/3 — Sankey triage diagram (panel a)")
    print("=" * 70)
    sankey_png = fig_sankey_triage(CSV_PATH, SANKEY_DIR)

    print("\n" + "=" * 70)
    print("STEP 2/3 — Violin panels (panel b)")
    print("=" * 70)
    violin_combined_png = fig_b_violins(CSV_PATH, VIOLIN_DIR)

    print("\n" + "=" * 70)
    print("STEP 3/3 — Combine panels (a) + (b)")
    print("=" * 70)
    if sankey_png is None:
        raise RuntimeError(
            "Sankey PNG export failed (no kaleido?) — cannot build the "
            "combined figure. Install kaleido (`pip install -U kaleido`) "
            "and re-run, or manually export sankey_triage.html to PNG."
        )
    combine_panels(sankey_png, violin_combined_png, FINAL_PATH)

    print("\nDone! Final figure saved at:", FINAL_PATH)


if __name__ == "__main__":
    main()
