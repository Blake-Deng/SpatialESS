#!/usr/bin/env python3
"""Create an English infographic explaining SpatialESS radius graphs and CSR."""

from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle


OUT = Path(__file__).resolve().parent
STEM = OUT / "SpatialESS_CSR_radius_graph_infographic_20260802"

NAVY = "#112B48"
BLUE = "#2A6F97"
TEAL = "#129E92"
TEAL_DARK = "#066F69"
ORANGE = "#EF843C"
PURPLE = "#7456AB"
RED = "#CB4A45"
GREEN = "#309B62"
INK = "#232D37"
MID = "#5E6C79"
GRID = "#DAE3E8"
LIGHT = "#F1F6F8"
LIGHT_BLUE = "#E8F2F8"
LIGHT_TEAL = "#E5F7F4"
LIGHT_ORANGE = "#FDF1E6"
LIGHT_RED = "#FBEBEA"
WHITE = "#FFFFFF"

mpl.rcParams.update({
    "font.family": "Liberation Sans",
    "font.size": 10,
    "axes.linewidth": 0,
    "figure.facecolor": WHITE,
    "savefig.facecolor": WHITE,
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
})


def card(fig, rect, label, title, accent, fill=WHITE):
    """Create a rounded card and return an inset axes."""
    x, y, w, h = rect
    patch = FancyBboxPatch(
        (x, y), w, h, transform=fig.transFigure,
        boxstyle="round,pad=0.006,rounding_size=0.014",
        linewidth=1.1, edgecolor=GRID, facecolor=fill, zorder=0,
    )
    fig.patches.append(patch)
    fig.text(x + 0.018, y + h - 0.038, label, color=WHITE, fontsize=10,
             fontweight="bold", ha="center", va="center",
             bbox=dict(boxstyle="round,pad=0.30", facecolor=accent,
                       edgecolor=accent))
    fig.text(x + 0.052, y + h - 0.038, title, color=NAVY, fontsize=14,
             fontweight="bold", ha="left", va="center")
    ax = fig.add_axes([x + 0.018, y + 0.022, w - 0.036, h - 0.078], zorder=2)
    ax.set_axis_off()
    return ax


def arrow(ax, xy1, xy2, color=TEAL, lw=1.6, mutation=12, style="-|>"):
    a = FancyArrowPatch(xy1, xy2, arrowstyle=style, mutation_scale=mutation,
                        linewidth=lw, color=color, shrinkA=2, shrinkB=2)
    ax.add_patch(a)
    return a


def radius_panel(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 7); ax.set_aspect("equal")
    center = (4.6, 3.55)
    pts = [(4.6,3.55),(3.4,4.2),(5.6,4.35),(5.8,2.8),(3.15,2.6),
           (6.8,5.2),(2.0,5.7),(8.4,3.7),(4.0,6.25),(1.25,1.15),
           (7.4,1.0),(9.2,6.15)]
    contact_r, interaction_r = 1.25, 2.65
    ax.add_patch(Circle(center, interaction_r, facecolor=LIGHT_TEAL,
                        edgecolor=TEAL, linewidth=1.7, linestyle="--"))
    ax.add_patch(Circle(center, contact_r, facecolor=LIGHT_ORANGE,
                        edgecolor=ORANGE, linewidth=1.6, linestyle="-"))
    # Draw only radius-supported edges from the selected sender.
    for p in pts[1:]:
        d = math.dist(center, p)
        if d <= interaction_r:
            arrow(ax, center, p, ORANGE if d <= contact_r else TEAL, 1.4, 10)
    for i, p in enumerate(pts):
        c = NAVY if i == 0 else (ORANGE if math.dist(center,p) <= contact_r else
                                  TEAL if math.dist(center,p) <= interaction_r else MID)
        ax.add_patch(Circle(p, 0.18 if i else 0.23, facecolor=c,
                            edgecolor=WHITE, linewidth=1.0, zorder=5))
    ax.text(center[0], center[1]-0.50, "sender $i$", color=NAVY,
            ha="center", va="top", fontsize=9, fontweight="bold")
    ax.text(5.63, 4.00, "contact: 10 + 5 = 15 µm", color=ORANGE,
            fontsize=8.5, fontweight="bold")
    ax.text(6.45, 5.75, "interaction: 30 + 5 = 35 µm", color=TEAL_DARK,
            fontsize=8.5, fontweight="bold")
    ax.text(0.30, 0.22, r"Store edge $i\rightarrow j$ only if $d_{ij}\leq r$",
            color=INK, fontsize=10.5, fontweight="bold")


def all_pairs_panel(ax):
    ax.set_xlim(0, 12); ax.set_ylim(0, 7)
    # Dense mini matrix.
    x0, y0, s, n = 0.45, 2.85, 0.35, 9
    ax.text(0.40, 6.15, "Dense all-pairs", color=RED, fontsize=11,
            fontweight="bold")
    for i in range(n):
        for j in range(n):
            ax.add_patch(Rectangle((x0+j*s,y0+(n-1-i)*s),s-0.015,s-0.015,
                                   facecolor=RED if (i+j)%5==0 else LIGHT_RED,
                                   edgecolor=WHITE, linewidth=0.2))
    ax.text(2.05, 2.45, r"$N\times N$", color=RED, fontsize=13,
            fontweight="bold", ha="center")
    arrow(ax, (4.05,4.35), (5.30,4.35), TEAL, 2.0, 13)
    ax.text(4.68,4.72,"replace",color=TEAL_DARK,fontsize=8.5,
            ha="center",fontweight="bold")
    # Sparse neighbor graph.
    nodes=[(6.1,4.7),(7.0,5.55),(7.25,3.65),(8.5,5.0),(8.75,3.5),
           (9.8,5.7),(10.4,4.2),(10.1,2.8),(11.25,5.0)]
    edges=[(0,1),(0,2),(1,3),(2,3),(2,4),(3,4),(3,5),(4,6),(4,7),(5,6),(6,7),(6,8)]
    for a,b in edges:
        ax.plot([nodes[a][0],nodes[b][0]],[nodes[a][1],nodes[b][1]],
                color=GRID,linewidth=1.2,zorder=1)
    for k,(x,y) in enumerate(nodes):
        ax.add_patch(Circle((x,y),0.13,facecolor=TEAL if k%3 else BLUE,
                            edgecolor=WHITE,linewidth=0.7,zorder=2))
    ax.text(6.05,6.15,"Sparse radius graph",color=TEAL_DARK,fontsize=11,
            fontweight="bold")
    ax.text(8.65,2.20,r"$N+1$ offsets  +  $E$ neighbors  +  $E$ distances",
            color=INK,fontsize=9.5,ha="center",fontweight="bold")
    # Real benchmark numbers.
    ax.add_patch(FancyBboxPatch((0.40,0.35),11.15,1.28,
                               boxstyle="round,pad=0.08,rounding_size=0.12",
                               facecolor=LIGHT_BLUE,edgecolor=LIGHT_BLUE))
    ax.text(0.75,1.17,"FULL XENIUM",color=BLUE,fontsize=8.2,fontweight="bold")
    ax.text(2.25,1.17,
            r"$N=1{,}156{,}091$   ·   $N^2\approx1.34\times10^{12}$ candidate pairs   ·   $E=33{,}297{,}328$ stored edges",
            color=NAVY,fontsize=8.8,fontweight="bold")
    ax.text(0.75,0.67,"Official dense: ≈9,958 GiB requested",color=RED,
            fontsize=8.9,fontweight="bold")
    ax.text(6.20,0.67,"SpatialESS, B=20: 5.49 GiB peak RSS",color=TEAL_DARK,
            fontsize=8.9,fontweight="bold")


def hashing_panel(ax):
    ax.set_xlim(0, 9); ax.set_ylim(0, 7); ax.set_aspect("equal")
    x0,y0,cell=1.25,0.75,1.05
    for i in range(5):
        for j in range(5):
            selected = 1 <= i <= 3 and 1 <= j <= 3
            ax.add_patch(Rectangle((x0+j*cell,y0+i*cell),cell,cell,
                                   facecolor=LIGHT_TEAL if selected else WHITE,
                                   edgecolor=GRID,linewidth=0.9))
    # Deterministic cells in buckets.
    pts=[(1.55,1.12),(2.82,1.52),(4.10,1.20),(5.02,2.14),(2.12,2.72),
         (3.50,3.35),(4.62,3.15),(5.92,3.75),(2.58,4.40),(3.83,4.72),
         (5.12,5.25),(6.02,5.68),(1.52,5.85),(6.26,1.08)]
    target=(3.50,3.35)
    for p in pts:
        c=NAVY if p==target else TEAL if (2.3<=p[0]<=5.5 and 2.0<=p[1]<=5.0) else MID
        ax.add_patch(Circle(p,0.12 if p!=target else 0.17,facecolor=c,
                            edgecolor=WHITE,linewidth=0.6))
    ax.text(target[0]+0.22,target[1]+0.20,"target",color=NAVY,fontsize=8,
            fontweight="bold")
    ax.text(0.25,6.45,"Partition coordinates into buckets",color=NAVY,
            fontsize=10.5,fontweight="bold")
    ax.text(0.25,6.00,"Inspect only the target bucket and adjacent buckets",
            color=INK,fontsize=9.2)
    ax.text(0.25,0.18,"No comparison with all 1.16 million cells",
            color=TEAL_DARK,fontsize=10,fontweight="bold")


def csr_panel(ax):
    ax.set_xlim(0, 20); ax.set_ylim(0, 8.5)
    colors=[BLUE,TEAL,ORANGE,PURPLE]
    # Tiny graph.
    coords={1:(1.35,5.80),2:(0.55,3.90),3:(3.20,4.20),4:(4.15,2.35)}
    directed=[(1,2),(1,3),(2,1),(3,1),(3,4),(4,3)]
    for a,b in directed:
        arrow(ax,coords[a],coords[b],colors[a-1],1.3,9)
    for k,p in coords.items():
        ax.add_patch(Circle(p,0.35,facecolor=colors[k-1],edgecolor=WHITE,
                            linewidth=1.2,zorder=5))
        ax.text(*p,str(k),ha="center",va="center",color=WHITE,fontsize=10,
                fontweight="bold",zorder=6)
    ax.text(0.35,7.65,"Directed neighbor graph",color=NAVY,fontsize=10.5,
            fontweight="bold")
    # Neighbor lists.
    lists=["Cell 1 → 2, 3","Cell 2 → 1","Cell 3 → 1, 4","Cell 4 → 3"]
    for i,t in enumerate(lists):
        ax.text(5.15,7.48-i*0.72,t,color=colors[i],fontsize=10,
                fontweight="bold",ha="left")
    # Arrays with colored segments.
    ax.text(10.10,7.65,"Concatenate the four neighbor lists",color=NAVY,
            fontsize=10.5,fontweight="bold")
    ax.text(10.10,6.88,"offsets",color=MID,fontsize=9,fontweight="bold")
    ax.text(12.05,6.88,"[ 0,  2,  3,  5,  6 ]",color=INK,fontsize=11,
            family="Liberation Mono",fontweight="bold")
    # Segment blocks for neighbors and distance rows.
    segments=[([2,3],colors[0]),([1],colors[1]),([1,4],colors[2]),([3],colors[3])]
    x_start=12.05
    widths=[]
    for vals,col in segments:
        ww=0.73*len(vals)
        widths.append(ww)
    ax.text(10.10,5.82,"neighbors",color=MID,fontsize=9,fontweight="bold")
    x=x_start
    for vals,col,ww in zip([s[0] for s in segments],colors,widths):
        ax.add_patch(FancyBboxPatch((x,5.50),ww,0.62,
                                   boxstyle="round,pad=0.03,rounding_size=0.06",
                                   facecolor=col,edgecolor=WHITE,linewidth=1.0))
        ax.text(x+ww/2,5.81,"  ".join(map(str,vals)),ha="center",va="center",
                color=WHITE,fontsize=10,fontweight="bold",family="Liberation Mono")
        x+=ww
    ax.text(10.10,4.67,"distances",color=MID,fontsize=9,fontweight="bold")
    dvals=[["d₁₂","d₁₃"],["d₂₁"],["d₃₁","d₃₄"],["d₄₃"]]
    x=x_start
    for vals,col,ww in zip(dvals,colors,widths):
        ax.add_patch(FancyBboxPatch((x,4.35),ww,0.62,
                                   boxstyle="round,pad=0.03,rounding_size=0.06",
                                   facecolor=WHITE,edgecolor=col,linewidth=1.2))
        ax.text(x+ww/2,4.66,"  ".join(vals),ha="center",va="center",
                color=col,fontsize=9.5,fontweight="bold")
        x+=ww
    # Slice table.
    rows=[("Cell 1","[0, 2)","2, 3"),("Cell 2","[2, 3)","1"),
          ("Cell 3","[3, 5)","1, 4"),("Cell 4","[5, 6)","3")]
    tx,ty=5.10,2.78
    colx=[tx,tx+2.55,tx+4.95]
    ax.text(colx[0],ty+0.65,"Sender",color=MID,fontsize=8.5,fontweight="bold")
    ax.text(colx[1],ty+0.65,"CSR slice",color=MID,fontsize=8.5,fontweight="bold")
    ax.text(colx[2],ty+0.65,"Neighbors",color=MID,fontsize=8.5,fontweight="bold")
    for i,(a,b,c) in enumerate(rows):
        yy=ty-i*0.56
        ax.text(colx[0],yy,a,color=colors[i],fontsize=9.2,fontweight="bold")
        ax.text(colx[1],yy,b,color=INK,fontsize=9.2,family="Liberation Mono")
        ax.text(colx[2],yy,c,color=INK,fontsize=9.2,family="Liberation Mono")
    ax.add_patch(FancyBboxPatch((11.65,0.68),7.35,2.43,
                               boxstyle="round,pad=0.10,rounding_size=0.10",
                               facecolor=LIGHT_TEAL,edgecolor=LIGHT_TEAL))
    ax.text(12.00,2.63,"Why the offsets array?",color=TEAL_DARK,fontsize=10.2,
            fontweight="bold")
    ax.text(12.00,2.08,"It is the table of contents:",color=INK,fontsize=9.4)
    ax.text(12.00,1.55,"offsets[i] : where Cell i's list starts",color=INK,
            fontsize=9.2,family="Liberation Mono")
    ax.text(12.00,1.05,"offsets[i+1] : where it ends (excluded)",color=INK,
            fontsize=9.2,family="Liberation Mono")
    ax.text(0.35,0.35,"The sender ID is implicit in offsets—no sender column is stored per edge.",
            color=TEAL_DARK,fontsize=9.8,fontweight="bold")


def stream_panel(ax):
    ax.set_xlim(0, 20); ax.set_ylim(0, 8.5)
    boxes=[
        ("1", "Sender\ncell", "read CSR slice", BLUE),
        ("2", "Spatial\nneighbors", "within 35 µm", TEAL),
        ("3", "Expression\nsupport", "$L_i>0$ and $R_j>0$", ORANGE),
        ("4", "Edge\nprobability", "Hill × distance", PURPLE),
        ("5", "Sparse\naggregation", "group → group", GREEN),
    ]
    x=0.35
    for i,(num,title,sub,col) in enumerate(boxes):
        w=3.30
        ax.add_patch(FancyBboxPatch((x,5.25),w,2.05,
                                   boxstyle="round,pad=0.09,rounding_size=0.15",
                                   facecolor=WHITE,edgecolor=col,linewidth=1.4))
        ax.add_patch(Circle((x+0.35,6.92),0.20,facecolor=col,edgecolor=col))
        ax.text(x+0.35,6.92,num,ha="center",va="center",color=WHITE,
                fontsize=8.5,fontweight="bold")
        ax.text(x+w/2,6.27,title,color=NAVY,fontsize=8.6,fontweight="bold",
                ha="center",va="center",linespacing=0.90)
        ax.text(x+w/2,5.53,sub,color=MID,fontsize=7.4,ha="center",va="center")
        if i<4:
            arrow(ax,(x+w+0.12,6.27),(x+w+0.65,6.27),TEAL,1.8,11)
        x+=3.88
    # Equations and memory rule.
    ax.add_patch(FancyBboxPatch((0.38,2.62),8.90,1.57,
                               boxstyle="round,pad=0.10,rounding_size=0.12",
                               facecolor=LIGHT_BLUE,edgecolor=LIGHT_BLUE))
    ax.text(0.72,3.73,"Official spatial weight",color=BLUE,fontsize=9,
            fontweight="bold")
    ax.text(4.83,3.20,r"$w_{ij}=\frac{1}{d_{ij}\,\times\,\mathrm{scale.distance}}$",
            color=NAVY,fontsize=16,ha="center",va="center")
    ax.add_patch(FancyBboxPatch((9.72,2.62),9.88,1.57,
                               boxstyle="round,pad=0.10,rounding_size=0.12",
                               facecolor=LIGHT_TEAL,edgecolor=LIGHT_TEAL))
    ax.text(10.05,3.73,"Streaming memory rule",color=TEAL_DARK,fontsize=9,
            fontweight="bold")
    ax.text(14.70,3.22,"one LR workspace → aggregate → release",
            color=NAVY,fontsize=12.2,ha="center",va="center",fontweight="bold")
    ax.text(0.45,1.43,r"CSR storage: $(N+1)+E+E=O(N+E)$",
            color=NAVY,fontsize=10.5,fontweight="bold")
    ax.text(12.15,1.43,r"Dense storage: $O(N^2)$",
            color=RED,fontsize=10.5,fontweight="bold")
    ax.add_patch(FancyBboxPatch((0.35,0.25),19.25,0.72,
                               boxstyle="round,pad=0.07,rounding_size=0.10",
                               facecolor=LIGHT_ORANGE,edgecolor=LIGHT_ORANGE))
    ax.text(9.98,0.61,
            "Boundary: if the radius becomes very large, E can approach N²; the gain relies on local spatial sparsity.",
            color=ORANGE,fontsize=9.7,fontweight="bold",ha="center",va="center")


def build():
    fig = plt.figure(figsize=(16, 10), dpi=200)
    fig.subplots_adjust(0,0,1,1)
    fig.text(0.035,0.958,"How SpatialESS stores spatially plausible communication",
             color=NAVY,fontsize=24,fontweight="bold",ha="left",va="top")
    fig.text(0.035,0.920,
             "A sparse radius graph and CSR neighbor directory replace the dense all-pairs distance matrix—without changing the official probability formula.",
             color=MID,fontsize=11.5,ha="left",va="top")
    fig.add_artist(Rectangle((0.035,0.895),0.080,0.006,transform=fig.transFigure,
                             facecolor=TEAL,edgecolor=TEAL))

    ax_a=card(fig,(0.030,0.555,0.300,0.320),"A","Radius-limited communication",TEAL,LIGHT_TEAL)
    ax_b=card(fig,(0.350,0.555,0.385,0.320),"B","Store local edges, not all pairs",BLUE,WHITE)
    ax_c=card(fig,(0.755,0.555,0.215,0.320),"C","Spatial hashing",ORANGE,LIGHT_ORANGE)
    ax_d=card(fig,(0.030,0.190,0.515,0.335),"D","CSR is a compact neighbor directory",PURPLE,WHITE)
    ax_e=card(fig,(0.565,0.190,0.405,0.335),"E","Stream communication over CSR",GREEN,WHITE)
    radius_panel(ax_a); all_pairs_panel(ax_b); hashing_panel(ax_c)
    csr_panel(ax_d); stream_panel(ax_e)

    fig.text(0.035,0.135,
             "Take-home: neighbors stores who can communicate; distances stores how far apart they are; offsets tells C++ where each sender's neighbor list begins and ends.",
             color=NAVY,fontsize=12.3,fontweight="bold",ha="left",va="center",
             bbox=dict(boxstyle="round,pad=0.55",facecolor=LIGHT_BLUE,edgecolor=LIGHT_BLUE))
    fig.text(0.035,0.075,
             "SpatialESS v3-compatible benchmark: 1,156,091 cells · 33,297,328 directed edges · 627 LR · 20 permutations · 27.19 min · 5.49 GiB peak RSS.",
             color=TEAL_DARK,fontsize=10.8,fontweight="bold",ha="left",va="center")
    fig.text(0.965,0.032,"SpatialESS · CSR radius graph · 2026-08-02",
             color=MID,fontsize=8.5,ha="right",va="center")

    for ext, kwargs in [
        ("png", dict(dpi=240)),
        ("pdf", {}),
        ("svg", {}),
    ]:
        fig.savefig(f"{STEM}.{ext}",bbox_inches="tight",pad_inches=0.08,**kwargs)
    plt.close(fig)


if __name__ == "__main__":
    OUT.mkdir(parents=True,exist_ok=True)
    build()
    print(f"{STEM}.png")
    print(f"{STEM}.pdf")
    print(f"{STEM}.svg")
