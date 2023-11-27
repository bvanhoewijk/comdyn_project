# importing required libraries
import re

import pandas as pd
import plotly
import plotly.graph_objs as go

# Configure Plotly to be rendered inline in the notebook.
plotly.offline.init_notebook_mode()


def parse_cg_pdb(file, segment_to_do="PROA"):
    """Parse a Martini3 coarsegrained pdb file. It only selects the atomtype BB
    Used the fields from this file:
    https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf


    Args:
        file (str): filepath to PDB

    Returns:
        DataFrame: Atom fields
        list: CONECT items
    """
    coordinates = []
    conect = []
    with open(file, "r") as fh:
        for line in fh:
            line = line.strip()
            elements = re.split(r" +", line)
            if elements[0] == "ATOM":
                type = line[0:4].strip()
                atomnr = line[6:11].strip()
                atomtype = line[12:16].strip()
                residue = line[17:20].strip()
                chain = line[21].strip()
                residuenr = line[22:26].strip()
                x = line[30:38].strip()
                y = line[39:46].strip()
                z = line[47:54].strip()
                segment = line[72:76].strip()
                
                if atomtype.startswith("BB") or atomtype.startswith("SC"):
                    coordinates.append([type, atomnr, atomtype, residue, chain, residuenr, x, y, z, segment])                                   
                # if atomtype == "BB":
            if elements[0] == "CONECT":
                conect_elements = elements[1:]
                conect_elements = [int(x) for x in conect_elements]
                conect.append(conect_elements)


    pdb = pd.DataFrame(coordinates)
    pdb.columns = [
        "type",
        "atomnr",
        "atomtype",
        "residue",
        "chain",
        "residuenr",
        "x",
        "y",
        "z",
        "segment",
    ]
    pdb["atomnr"] = pdb["atomnr"].astype("int")
    pdb["residuenr"] = pdb["residuenr"].astype("int")
    pdb["x"] = pdb["x"].astype("float")
    pdb["y"] = pdb["y"].astype("float")
    pdb["z"] = pdb["z"].astype("float")
    
    return pdb, conect

def get_rubber_bonds(itp, pdb_bb, min_local_diff):
    """_summary_

    Args:
        itp (_type_): _description_
        pdb_bb (_type_): _description_
        min_local_diff (_type_): _description_

    Returns:
        _type_: _description_
    """    
    rubber_bands = []
    for res1, res2 in itp["res_tuple"]:
        if abs(res2 - res1) <= min_local_diff:
            continue
        pdb_a = pdb_bb[pdb_bb["residuenr"] == res1].iloc[0]
        pdb_b = pdb_bb[pdb_bb["residuenr"] == res2].iloc[0]
        pdb_ax, pdb_ay, pdb_az = (
            pdb_a["x"],
            pdb_a["y"],
            pdb_a["z"],
        )
        pdb_bx, pdb_by, pdb_bz = (
            pdb_b["x"],
            pdb_b["y"],
            pdb_b["z"],
        )

        rubber_bands.append(
            go.Scatter3d(
                x=[pdb_ax, pdb_bx],
                y=[pdb_ay, pdb_by],
                z=[pdb_az, pdb_bz],
                name="Rubber bands",
                legendgroup="Rubber bands",
                mode="lines",
                text=f"{pdb_a['residue']}{pdb_a['residuenr']}-{pdb_b['residue']}{pdb_b['residuenr']}",
                hoverinfo="text",
                marker=dict(
                    size=3, opacity=0, color="orange", line=dict(color="yellow")
                ),
                showlegend=False,
            ),
        )
    return rubber_bands

def get_aa_beads(pdb_bb):
    """_summary_

    Args:
        pdb_bb (_type_): _description_

    Returns:
        _type_: _description_
    """    
    backbone_beads = []
    hover_text_beads = list(pdb_bb.apply(lambda x: get_bb_hover(x), axis=1))
    backbone_beads.append(
        go.Scatter3d(
            x=pdb_bb["x"],
            y=pdb_bb["y"],
            z=pdb_bb["z"],
            mode="markers",
            name="AA beads",
            legendgroup="AA beads",
            text=hover_text_beads,
            hoverinfo="text",
            marker=dict(size=5, color="blue"),
            showlegend=False,
        )
    )
    return backbone_beads
def get_sc_beads(pdb_sc):
    """_summary_

    Args:
        pdb_sc (_type_): _description_

    Returns:
        _type_: _description_
    """    
    sc_beads = []
    hover_text_beads = list(pdb_sc.apply(lambda x: get_sc_hover(x), axis=1))
    sc_beads.append(
        go.Scatter3d(
            x=pdb_sc["x"],
            y=pdb_sc["y"],
            z=pdb_sc["z"],
            mode="markers",
            name="SC beads",
            legendgroup="SC beads",
            text=hover_text_beads,
            hoverinfo="text",
            marker=dict(size=5, color="green"),
            showlegend=False,
        )
    )
    return sc_beads

def get_bb_links(pdb_bb):
    """_summary_

    Args:
        pdb_bb (_type_): _description_

    Returns:
        _type_: _description_
    """    
    backbone_links = []
    for i in range(len(pdb_bb) - 1):
        a = i + 1
        b = i + 2

        pdb_a = pdb_bb[pdb_bb["residuenr"] == a].iloc[0]
        pdb_b = pdb_bb[pdb_bb["residuenr"] == b].iloc[0]

        pdb_ax, pdb_ay, pdb_az = (
            pdb_a["x"],
            pdb_a["y"],
            pdb_a["z"],
        )
        pdb_bx, pdb_by, pdb_bz = (
            pdb_b["x"],
            pdb_b["y"],
            pdb_b["z"],
        )
        backbone_links.append(
            go.Scatter3d(
                x=[pdb_ax, pdb_bx],
                y=[pdb_ay, pdb_by],
                z=[pdb_az, pdb_bz],
                name="backbone",
                legendgroup="backbone",
                mode="lines",
                text=f"{pdb_a['residue']}{pdb_a['residuenr']}-{pdb_b['residue']}{pdb_b['residuenr']}",
                hoverinfo="text",
                line=dict(width=10, color="grey"),
                showlegend=False,
            ),
        )
    return backbone_links

def get_sc_links(pdb_bb, sc):
    """_summary_

    Args:
        pdb_bb (_type_): _description_
        sc (_type_): _description_

    Returns:
        _type_: _description_
    """    
    residuenr = sc['residuenr']

    result = pdb_bb[pdb_bb["residuenr"] == residuenr].iloc[0]
    result['x'], result['y'], result['z']

    scatter = go.Scatter3d(
                x=[result['x'], sc['x']],
                y=[result['y'], sc['y']],
                z=[result['z'], sc['z']],
                name="backbone",
                legendgroup="backbone",
                mode="lines",
                line=dict(width=3, color="grey"),
                showlegend=False,
            )
    return scatter

def get_highlight_beads(pdb, highlights):
    """_summary_

    Args:
        pdb (_type_): _description_
        highlights (_type_): _description_

    Returns:
        _type_: _description_
    """    
    highlight_beads = []
    for highlight in highlights:
        pdb_a = pdb[pdb["residuenr"] == highlight].iloc[0]

        highlight_beads.append(
            go.Scatter3d(
                x=[pdb_a["x"]],
                y=[pdb_a["y"]],
                z=[pdb_a["z"]],
                mode="markers",
                name=f"highlight",
                legendgroup="AA beads",
                text=get_bb_hover(pdb_a),
                hoverinfo="text",
                marker=dict(size=8, color="red"),
                showlegend=False,
            )
        )
    return highlight_beads

def get_pdb_points(
    pdb, itp, plot_connection=True, plot_rubber=True, plot_sc=True, min_local_diff=0, highlights=[]
):
    """_summary_

    Args:
        pdb (DataFrame): PDB dataframe
        itp (DataFrame): _description_
        plot_connection (bool, optional): Plot connection between beads. Defaults to True.
        plot_rubber (bool, optional): Plot rubber bands. Defaults to True.
        min_local_diff (int, optional): Number of residues a rubber band should be apart before plotting. Default: 0
        highlights (list, optional): List of residuenumbers to highlight.
    Returns:
        list: list of lists of go.Scatter3d objects.
    """
    # Just the sequence of AA residues:
    backbone_links = []
    pdb_bb = pdb[pdb["atomtype"] == "BB"]
    if plot_connection:
        backbone_links = get_bb_links(pdb_bb)

    # Rubber bonds:
    rubber_bands = []
    if plot_rubber:
        rubber_bands = get_rubber_bonds(itp, pdb_bb, min_local_diff)
        

    # Plot individual beads:
    backbone_beads = get_aa_beads(pdb_bb)


    # Plot SC beads:
    pdb_sc = pdb[pdb["atomtype"] != "BB"]
    sc_links = []
    sc_beads = []
    if plot_sc and len(pdb_sc):
        sc_beads = get_sc_beads(pdb_sc)
        # Calculate connections:
        sc_links = list(pdb_sc.apply(lambda sc: get_sc_links(pdb_bb, sc), axis=1))


    highlight_beads = get_highlight_beads(pdb, highlights)

    print(f"BBBs in plot   : {len(pdb['x'])}")
    print(f"Rubbers in plot: {len(rubber_bands)}")
    

    return [backbone_links, rubber_bands, backbone_beads, highlight_beads, sc_links, sc_beads]


def get_bb_hover(x):
    """Get hover text for a single PDB series

    Args:
        x (Serie): PDB serie

    Returns:
        str: formatter hover string
    """
    return (
        f"{x['residue']}{x['residuenr']}<br>x: {x['x']}<br>y: {x['y']}<br>z: {x['z']}"
    )

def get_sc_hover(x):
    """Get hover text for a single PDB series

    Args:
        x (Serie): PDB serie

    Returns:
        str: formatter hover string
    """
    return (
        f"{x['atomtype']} of {x['residue']}{x['residuenr']}<br>x: {x['x']}<br>y: {x['y']}<br>z: {x['z']}"
    )


def axes_style3d(
    bgcolor="rgb(20, 20, 20)", gridcolor="rgb(150, 150, 150)", zeroline=False
):
    return dict(
        showbackground=True,
        backgroundcolor=bgcolor,
        gridcolor=gridcolor,
        zeroline=zeroline,
    )


def plot_pdb(
    points, title=None, camera_eye=dict(x=0.1, y=0.1, z=2.5), width=800, height=800
):
    """Plot the data

    Args:
        points (list): list of lists of go.Scatter3d objects.
        title (str, optional): Plot title. Defaults to None.
        camera_eye (dict, optional): Camera eye angle. Defaults to dict(x=0.1, y=0.1, z=2.5).
        width (int, optional): Plot width. Defaults to 800.
        height (int, optional): Plot height. Defaults to 800.
    """

    #### Define layout, margins and width:
    # transparent background color
    my_axes = axes_style3d(
        bgcolor="rgba(0, 0, 0, 0)",
        gridcolor="rgb(100, 100, 100)",
    )
    layout = go.Layout(
        margin={"l": 0, "r": 0, "b": 0, "t": 40}, width=width, height=height
    )

    # Build figure from the different series stored in the points list
    plot_figure = go.Figure(layout=layout)
    for series in points:
        for item in series:
            plot_figure.add_trace(item)
    
    # Set the scene
    scene = dict(
        xaxis=my_axes,
        yaxis=my_axes,
        zaxis=my_axes,
        camera_eye=camera_eye,
        camera_up=dict(x=0, y=1, z=0),
        camera_center=dict(x=0, y=0, z=0),
    )

    # Render the plot.
    plot_figure.update_scenes(
        xaxis_visible=True, yaxis_visible=True, zaxis_visible=True
    )
    plot_figure.update_layout(scene=scene, title=title)
    plotly.offline.iplot(plot_figure)
