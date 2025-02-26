import matplotlib as mpl
from cycler import cycler


# colors DS
colors_ds = {"antiinflammatory_high_nrf2_high": "#F22786",
             "mixed": "#5D29A6",
             "inflammatory_high_nrf2_low": "#3DF2B0",
             "bg_bright": "#F2F2F2",
             "bg_dark": "#595959"}

colors_pypalette_austria_ordered_list = ['#A40000FF', '#16317DFF', '#007E2FFF', '#FFCD12FF', '#B86092FF', '#721B3EFF', '#00B7A7FF'] # red, blue, green, yellow, pink, dark pink, mint
colors_pypalette_austria = {"antiinflammatory_high_nrf2_high": "#B86092",
             "mixed": "#721B3E",
             "inflammatory_high_nrf2_low": "#00B7A7",
             "bg_bright": "#F2F2F2",
             "bg_dark": "#595959"}

colors_pypalette_bold_ordered_list = ['#6497B1', '#6A359C', '#FFB04F', '#679C35', '#CD1076'] # blue, purple, yellow, green, pink
colors_pypalette_bold = {"antiinflammatory_high_nrf2_high": "#CD1076",
             "mixed": "#6A359C",
             "inflammatory_high_nrf2_low": "#6497B1",
             "bg_bright": "#F2F2F2",
             "bg_dark": "#595959"}

colors_pypalette_bold2 = {"antiinflammatory_high_nrf2_high": "#CD1076",
             "mixed": "#6A359C",
             "inflammatory_high_nrf2_low": "#679C35",
             "bg_bright": "#F2F2F2",
             "bg_dark": "#595959"}

available_colors = colors_pypalette_bold2

tme_cmap_ds = mpl.colors.LinearSegmentedColormap.from_list('tme_cmap_ds', (
    # Edit this gradient at https://eltos.github.io/gradient/#3DF2B0-F2F2F2-F22786
    (0.000, (0.239, 0.949, 0.690)),
    (0.500, (0.949, 0.949, 0.949)),
    (1.000, (0.949, 0.153, 0.525))))

diverging_score = mpl.colors.LinearSegmentedColormap.from_list('tme_cmap_ds', (
    (0.000, mpl.colors.to_rgb(available_colors['inflammatory_high_nrf2_low'])),
    (0.500, mpl.colors.to_rgb(available_colors['bg_bright'])),
    (1.000, mpl.colors.to_rgb(available_colors['antiinflammatory_high_nrf2_high']))))

hex_color_list_ds = list(colors_ds.values())
rgb_color_list_ds = [mpl.colors.to_rgb(col)for col in hex_color_list_ds]

tme_palette_ds = cycler(color = hex_color_list_ds)

#tme_palette_ds = mpl.colors.ListedColormap(
#    name="tme_palette_ds",
#    colors=[colors_ds['inflammatory_high_nrf2_low'], colors_ds['mixed'], colors_ds['antiinflammatory_high_nrf2_high'], colors_ds['bg_dark'], colors_ds['bg_bright']]
#)

gene_expression = mpl.colors.LinearSegmentedColormap.from_list('expression', (
    # Edit this gradient at https://eltos.github.io/gradient/#F2F2F2-5D29A6
    (0.000, mpl.colors.to_rgb(available_colors['bg_bright'])),
    (1.000, mpl.colors.to_rgb(available_colors['mixed']))))

white_to_dark_gray = mpl.colors.LinearSegmentedColormap.from_list('white_to_dark_gray', (
    # Edit this gradient at https://eltos.github.io/gradient/#FFFFFF-595959
    (0.000, (1.000, 1.000, 1.000)),
    (1.000, mpl.colors.to_rgb(available_colors['bg_dark']))))
score_mac = white_to_dark_gray

score_inflammatory = mpl.colors.LinearSegmentedColormap.from_list('expression', (
    # Edit this gradient at https://eltos.github.io/gradient/#F2F2F2-5D29A6
    (0.000, mpl.colors.to_rgb(available_colors['bg_bright'])),
    (1.000, mpl.colors.to_rgb(available_colors['inflammatory_high_nrf2_low']))))
score_antiinflammatory = mpl.colors.LinearSegmentedColormap.from_list('expression', (
    # Edit this gradient at https://eltos.github.io/gradient/#F2F2F2-5D29A6
    (0.000, mpl.colors.to_rgb(available_colors['bg_bright'])),
    (1.000, mpl.colors.to_rgb(available_colors['antiinflammatory_high_nrf2_high']))))

set_of_three_plus_bg = cycler(color = [mpl.colors.to_rgb(col) for col in [available_colors['inflammatory_high_nrf2_low'], available_colors['mixed'], available_colors['antiinflammatory_high_nrf2_high'], available_colors['bg_bright']]])
cluster_cmap = cycler(color = [mpl.colors.to_rgb(col) for col in [available_colors['inflammatory_high_nrf2_low'], available_colors['mixed'], available_colors['antiinflammatory_high_nrf2_high'], available_colors['bg_dark']]])


#pypalette: Austria, Bold