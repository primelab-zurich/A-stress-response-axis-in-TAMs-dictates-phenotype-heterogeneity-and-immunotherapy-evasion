import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb
import numpy as np
import pathlib

def bilinear_color(s_val, c_val, c00, c10, c01, c11):
    """
    (s_val, c_val) in [0,1]^2 => returns an (R,G,B) color by bilinear interpolation:
      (0,0) -> c00
      (1,0) -> c10
      (0,1) -> c01
      (1,1) -> c11
    """
    return (
        (1 - s_val)*(1 - c_val)*c00 +
        s_val*(1 - c_val)*c10 +
        (1 - s_val)*c_val*c01 +
        s_val*c_val*c11
    )

def zscore_to_01(values):
    """
    1) z-score: z = (x - mean) / std
    2) min-max scale z to [0..1]
    """
    mean_ = np.mean(values)
    std_  = np.std(values)
    z = (values - mean_) / (std_ + 1e-9)

    zmin, zmax = z.min(), z.max()
    if zmax - zmin < 1e-9:
        # No variation => all identical
        return np.full_like(values, 0.5)

    return (z - zmin) / (zmax - zmin)

def plot_bicolor(x, y, colors, color_legend_map, x_color_label = "", y_color_label = "", figsize = (8, 4), dot_size = 5, dot_alpha = 0.8):
    fig,(plot_ax, bar_ax) = plt.subplots(nrows=1,ncols=2,figsize=(8,4))

    plot_ax.scatter(x, y, c=colors, s=dot_size, alpha=dot_alpha)
    plot_ax.set_title(x_color_label + ' and ' + y_color_label)
    plot_ax.set_aspect('equal')
    plot_ax.grid(False)

    bar_ax.imshow(color_legend_map, origin="lower", extent=[0,1,0,1], aspect="equal")
    bar_ax.set_title("2D color scale", fontsize=9)
    bar_ax.grid(False)
    bar_ax.set_xlabel(x_color_label, fontsize=8)
    bar_ax.set_ylabel(y_color_label, fontsize=8)
    bar_ax.set_xticks([0,1])
    bar_ax.set_yticks([0,1])

    plt.tight_layout()

    return (plot_ax, bar_ax)

def bicolor_embedding_plot(adata, x_key, y_key, c00 = "#F2F2F2", c10 = "#FF1493", c01 = "#40E0D0", c11 = "#5D29A6", embedding = 'pca', components = [0, 1], save_path = None):
    df = sc.get.obs_df(adata, keys= [x_key, y_key], obsm_keys= [('X_{}'.format(embedding), components[0]), ('X_{}'.format(embedding), components[1])])
    embedding_x_name, embedding_y_name = 'X_{}-{}'.format(embedding, components[0]), 'X_{}-{}'.format(embedding, components[1])
    df['z_score_' + x_key] = zscore_to_01(df[x_key])
    df['z_score_' + y_key] = zscore_to_01(df[y_key])
    
    c00 = np.array(to_rgb(c00))  # Very Light Gray
    c10 = np.array(to_rgb(c10))  # Deep Pink
    c01 = np.array(to_rgb(c01))  # Turquoise
    c11 = np.array(to_rgb(c11))  # Rich Purple

     # Create a color matrix for the 2D legend
    res_legend = 100
    sx = np.linspace(0, 1, res_legend)
    cx = np.linspace(0, 1, res_legend)

    legend_matrix = np.zeros((res_legend, res_legend, 3))
    for ix in range(res_legend):
        for iy in range(res_legend):
            s_val = sx[ix]
            c_val = cx[iy]
            legend_matrix[iy, ix] = bilinear_color(s_val, c_val, c00, c10, c01, c11)

    df['bicolor'] = df.apply(lambda row: bilinear_color(row['z_score_' + x_key], row['z_score_' + y_key], c00, c10, c01, c11), axis = 1)
    plot_ax, (bar_ax) = plot_bicolor(df[embedding_x_name], df[embedding_y_name], df['bicolor'], color_legend_map = legend_matrix, x_color_label = x_key + ' norm', y_color_label = y_key + ' norm', figsize = (8, 4), dot_size = 5, dot_alpha = 0.8)

    plot_ax.set_xlabel(embedding_x_name)
    plot_ax.set_ylabel(embedding_y_name)
    if save_path:
         plt.savefig(pathlib.Path(save_path) / 'bicolor_{}_{}_{}_{}_{}.png'.format(embedding, components[0], components[1], x_key, y_key).replace(' ','_'))
    plt.show()