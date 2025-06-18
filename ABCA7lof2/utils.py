def apply_custom_plot_style():
    """
    Apply a consistent Seaborn/matplotlib style:
      - ticks style with thinner lines & spines
      - Helvetica/Arial sans-serif fonts
      - PDF/SVG font settings for editable text
      - font sizes tuned for “talk” context
    """
    import seaborn as sns
    import matplotlib as mpl

    # Seaborn theme with thin lines, ticks & spines
    sns.set_theme(
        style="ticks",
        context="talk",
        rc={
            'lines.linewidth':    1.0,   # thinner plot lines
            'xtick.major.width':  0.6,   # thinner tick marks
            'ytick.major.width':  0.6,
            'xtick.minor.width':  0.4,
            'ytick.minor.width':  0.4,
            'axes.linewidth':     0.8,   # slimmer axis spines
        }
    )

    # Fonts & export settings
    mpl.rcParams['pdf.fonttype'] = 42      # TrueType fonts in PDF
    mpl.rcParams['svg.fonttype'] = 'none'  # outlined fonts in SVG
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['font.sans-serif'] = ['Helvetica', 'Arial']

    # Base font sizes
    mpl.rcParams.update({
        'font.size':        12,  # default text size
        'axes.labelsize':   14,  # x/y label size
        'xtick.labelsize':  12,  # tick label size
        'ytick.labelsize':  12,
        'legend.fontsize':  12,
        'axes.titlesize':   16,
    })
