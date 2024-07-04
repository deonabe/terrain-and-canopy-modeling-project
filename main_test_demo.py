#!/usr/bin/python3
"""
Demo.
Authors:
Xin Xu: xinxu629@umd.edu
Yunting Song: ytsong@terpmail.umd.edu 
Mar-04-2022
"""
import sys

import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# plot
import numpy as np
from matplotlib import cm
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# generate TIN
from scipy.spatial import Delaunay

from reader import Reader
from tree import Tree


def generate_TIN(pts_file):
    """
    Generate TIN from input point file.
    Return
    tin_file: name of TIN .off file
    points1: updated 2D points (x,y). each pair of (x,y) is unique
    zs1: updated z values of points, in the same order of points1
    tris: triangles of TIN. each triangle is saved as [vid1, vid2, vid3]. vid is point index based on points1
    """
    points = np.empty(shape=[0, 2])
    with open(pts_file) as infile:
        lines_num = infile.readline().strip()
        lines_num = int(lines_num)
        zs = list()
        for l in range(lines_num):
            line = (infile.readline()).split()
            v1 = float(line[0])
            v2 = float(line[1])
            zs.append(float(line[2]))
            points = np.append(points, [[v1, v2]], axis=0)

    triangles = Delaunay(points)

    tris = np.empty(shape=[0, 3])
    used_pts = set()
    for tri in triangles.simplices:
        used_pts.add(tri[0])
        used_pts.add(tri[1])
        used_pts.add(tri[2])
    o2n = dict()
    zs1 = list()
    points1 = np.empty(shape=[0, 2])
    for pid in used_pts:
        points1 = np.append(points1, [[points[pid][0], points[pid][1]]], axis=0)
        o2n[pid] = len(points1) - 1
        zs1.append(zs[pid])

    for tri in triangles.simplices:
        v1 = tri[0]
        v2 = tri[1]
        v3 = tri[2]
        nv1 = o2n[v1]
        nv2 = o2n[v2]
        nv3 = o2n[v3]
        tris = np.append(tris, [[nv1, nv2, nv3]], axis=0)

    # output TIN file, Note:  you can rename the file name
    tin_file = "pts-dt.off"
    with open(tin_file, "w") as ofs:
        ofs.write("OFF\n")
        ofs.write("{} {} 0\n".format(len(points1), len(tris)))
        for pid in range(len(points1)):
            ofs.write("{} {} {}\n".format(points1[pid][0], points1[pid][1], zs1[pid]))
        for tri in tris:
            ofs.write("3 {} {} {}\n".format(int(tri[0]), int(tri[1]), int(tri[2])))
    return tin_file, points1, zs1, tris


def plot_tin_with_marks(xs, ys, zs, tris, vals, mxs, mys, mzs, filename="test"):
    """
    Plot TIN 3D with markers. TIN is colored  based on values of each triangle. The triangle value is the  average value of 3 extreme points.
    Inputs:
    xs,ys,zs: lists of the same size. Point i
    (xs[i],ys[i], zs[i]) represents a point
    tris: np.array. triangles. each triangle is represented as [vid1, vid2, vid3]. vid is point index.
    vals: list. values of points with the same order as points. vals are used for calculating the color of triangles
    mxs,mys,mzs: lists. Marker i
    (mxs[i],mys[i], mzs[i]) represents a marker
    filename: output figure name.
    """
    tri_avg = []
    for tri in tris:
        v1 = vals[int(tri[0])]
        v2 = vals[int(tri[1])]
        v3 = vals[int(tri[2])]
        v = (v1 + v2 + v3) / 3
        tri_avg.append(v)
    vals_np = np.array(vals)
    zs_np = np.array(zs)
    triang = mtri.Triangulation(xs, ys, tris)
    maskedTris = triang.get_masked_triangles()
    xt = triang.x[maskedTris]
    yt = triang.y[maskedTris]
    zt = zs_np[maskedTris]
    verts = np.stack((xt, yt, zt), axis=-1)
    norm = cm.colors.Normalize(vmin=min(tri_avg), vmax=max(tri_avg))
    nm = norm(tri_avg)

    my_col = cm.jet(nm)
    newcmp = cm.colors.ListedColormap(my_col)

    collection = Poly3DCollection(verts)
    collection.set_facecolor(my_col)

    fig = plt.figure(figsize=plt.figaspect(0.5))
    ax = fig.add_subplot(projection="3d")

    ax.add_collection(collection)
    # add markers
    ax.scatter(mxs, mys, mzs, c="r", marker="^", s=40)

    ax.set_title(filename)
    ax.set_xlim3d(min(xs), max(xs))
    ax.set_xlabel("X")
    ax.set_ylim3d(min(ys), max(ys))
    ax.set_ylabel("Y")
    ax.set_zlim3d(min(zs), max(zs))
    ax.set_zlabel("Z")
    ax.autoscale_view()

    m = cm.ScalarMappable(cmap=cm.jet, norm=norm)
    m.set_array([])
    fig.colorbar(m, ax=ax, location="left")

    # output tin figure
    plt.savefig(filename + ".png", dpi=96)
    plt.show()


def write_feature_data(filename, data, vertices):
    with open(filename, "w") as file:
        for vertex in vertices:
            feature_value = data.get(vertex, 0)  # Default to 0 if no value exists
            file.write(
                f"{vertex} {feature_value:.2f}\n"
            )  # Using f-string for formatting


def write_maxima_data(filename, maxima):
    with open(filename, "w") as file:
        file.write(f"{len(maxima)}\n")  # First write the number of maxima
        for vertex in maxima:
            file.write(f"{vertex}\n")


def main():
    data_type = (
        input("Enter the dataset type ('bathymetric' or 'chm'): ").strip().lower()
    )
    if data_type == "bathymetric":
        # Load bathymetric data
        # Compute roughness, curvature, and maxima
        # Generate and save 3D visualization and text output
        return 1
    elif data_type == "chm":
        # Load CHM data
        # Compute maxima
        # Generate and save 3D visualization and text output
        return 0
    else:
        print("Invalid dataset type entered. Please choose 'bathymetric' or 'chm'.")


if __name__ == "__main__":
    data_type = main()
    tin_file, pts, zs, tris = generate_TIN(
        sys.argv[1]
    )  # sys.argv[1] should be the point file name
    # the above function will generate the TIN based on the input point file
    # please read the comments in plot_tin_with_marks function to understand how to use this function
    plot_tin_with_marks(
        pts[:, 0], pts[:, 1], zs, tris, zs, [], [], [], "tin-og"
    )  # This function will plot the original TIN

    reader = Reader()  # initialize a Reader object

    tin = reader.read_tin_file(
        tin_file
    )  # This function will read the TIN file written by generate_TIN() function.
    # You need to COMPLETE the read_tin_file() function in the reader.py

    capacity = sys.argv[2]  # sys.argv[2] should be the capacity value

    tree = Tree(int(capacity))  # initialize a Tree object.
    tree.build_tree(
        tin
    )  # build_tree() will generate a Triangle PR-quadtree on the input TIN.
    # You need to COMPLETE the build_tree() function in the tree.py

    roughness, curvature, maxima = tree.compute_features(tin)

    # Retrieve vertex data for feature and maxima writing
    vertices = [tin.get_vertex(i) for i in range(tin.get_vertices_num())]
    vertex_indices = [i for i in range(len(vertices))]

    # Extract vertex coordinates and triangle indices
    xs, ys, zs = pts[:, 0], pts[:, 1], zs
    vertex_indices = list(
        range(len(xs))
    )  # Assuming vertex indices match the length of xs, ys, zs

    if data_type:
        write_feature_data("roughness.txt", roughness, vertex_indices)
        write_feature_data("curvature.txt", curvature, vertex_indices)
        roughness_vals = [roughness.get(i, 0) for i in vertex_indices]
        curvature_vals = [curvature.get(i, 0) for i in vertex_indices]

        # Calculate average roughness and curvature for each triangle
        tri_roughness = [
            (
                roughness.get(tris[i][0], 0)
                + roughness.get(tris[i][1], 0)
                + roughness.get(tris[i][2], 0)
            )
            / 3
            for i in range(len(tris))
        ]
        tri_curvature = [
            (
                curvature.get(tris[i][0], 0)
                + curvature.get(tris[i][1], 0)
                + curvature.get(tris[i][2], 0)
            )
            / 3
            for i in range(len(tris))
        ]

        # Plotting roughness
        plot_tin_with_marks(
            xs, ys, zs, tris, tri_roughness, [], [], [], "roughness_plot"
        )

        # Plotting curvature
        plot_tin_with_marks(
            xs, ys, zs, tris, tri_curvature, [], [], [], "curvature_plot"
        )

    write_maxima_data("maxima.txt", maxima)
    maxima_vals = [1 if i in maxima else 0 for i in vertex_indices]

# Prepare maxima coordinates
maxima_xs = [xs[i] for i in maxima]
maxima_ys = [ys[i] for i in maxima]
maxima_zs = [zs[i] for i in maxima]

# For maxima, assuming you might want to highlight them on the original terrain
plot_tin_with_marks(
    xs, ys, zs, tris, zs, maxima_xs, maxima_ys, maxima_zs, "maxima_plot"
)
