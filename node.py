import numpy as np

from domain import Domain
from point import Point

# children legend:
#   nw = 0
#   ne = 1
#   sw = 2
#   se = 3


class Node(object):
    """Creates Class node"""

    def __init__(self):
        self.__vertex_ids = list()  # indices of points
        self.__triangle_ids = list()  # indices of triangles
        self.__children = (
            None  # a list of Node type objects. Equals to None when it is a leaf node.
        )
        # When the current node is INTERNAL:
        # children[0]:nw quadrant  children[1]:ne quadrant children[2]:sw quadrant children[3]:se quadrant

    def add_vertex(
        self, id
    ):  # When you add a vertex to a node, you add the index of it.
        # Here the id is the index of the Vertex on the global list
        # Should implement similar function for adding triangle
        self.__vertex_ids.append(id)

    def init_children(self):  # initialize four empty Nodes as the children.
        self.__children = [Node() for _ in range(4)]

    def get_child(
        self, i
    ):  # it returns a Node object children[i] according to the input i from 0 to 3.
        return self.__children[i]

    def reset_vertices(self):  # remove all vertices ids in the vertex list of the node
        self.__vertex_ids = list()

    def overflow(self, capacity):  # if the number of vertices exceed the capacity
        return len(self.__vertex_ids) > capacity

    def compute_child_label_and_domain(
        self, child_position, node_label, node_domain, mid_point
    ):  # Compute subdivision and labels of four children
        if child_position == 0:  # "nw":
            min = Point(node_domain.get_min_point().get_x(), mid_point.get_y())
            max = Point(mid_point.get_x(), node_domain.get_max_point().get_y())
            return 4 * node_label + 1, Domain(min, max)
        elif child_position == 1:  # "ne":
            return 4 * node_label + 2, Domain(mid_point, node_domain.get_max_point())
        elif child_position == 2:  # "sw":
            return 4 * node_label + 3, Domain(node_domain.get_min_point(), mid_point)
        elif child_position == 3:  # "se":
            min = Point(mid_point.get_x(), node_domain.get_min_point().get_y())
            max = Point(node_domain.get_max_point().get_x(), mid_point.get_y())
            return 4 * node_label + 4, Domain(min, max)
        else:
            return None, None  # #

    def is_duplicate(
        self, v_index, tin
    ):  # vertex with index v_index has the same x,y coordinates as an existing vertex in the node
        for (
            i
        ) in (
            self.get_vertices()
        ):  # check the vertices in this node to see if there is any vertex with same x,y coordinates as the inserting vertex
            if tin.get_vertex(i) == tin.get_vertex(
                v_index
            ):  # == for vertices is based on the x,y coordinates
                return True
        return False

    def get_vertices(
        self,
    ):  # returns the list of vertex ids. Should implement similar function for triangles.
        return self.__vertex_ids

    def get_vertices_num(self):
        return len(self.__vertex_ids)

    def is_leaf(self):  # returns True if the node is leaf node, otherwise returns False
        return self.__children == None

    ### You will need to add your own functions to add/get triangles to a Node

    def add_triangle_id(self, triangle_id):  # triangle should be a tuple (v1, v2, v3)
        self.__triangle_ids.append(triangle_id)

    def get_triangle_id(self, index):  # Get a triangle by its index
        try:
            return self.__triangle_ids[index]
        except IndexError:
            # print(f"Triangle with index {index} does not exist.")
            return None

    # Function to get the number of triangles
    def get_triangles_num(self):
        return len(self.__triangle_ids)

    def extract_vt_relations(self, tin):
        vts = {vid: [] for vid in self.__vertex_ids}
        for tid in self.__triangle_ids:
            triangle = tin.get_triangle(tid)
            # if not triangle:
            #     print(f"Warning: Triangle {tid} is missing or invalid")
            if triangle:
                for vertex in triangle:
                    if vertex in vts:
                        vts[vertex].append(tid)
                    # else:
                    #     print(
                    #         f"Warning: Vertex {vertex} from triangle {tid} is not in vertex IDs list"
                    #     )
        return vts

    def compute_vv_relations(self, tin, vts):
        vvs = {vid: set() for vid in self.__vertex_ids}
        for vid, triangles in vts.items():
            # if not triangles:
            #     print(f"Vertex {vid} has no connecting triangles")
            for tid in triangles:
                triangle = tin.get_triangle(tid)
                # if not triangle:
                #     print(f"Warning: Failed to get triangle {tid} for vertex {vid}")
                if triangle:
                    # Add vertices of the triangle to the set, excluding the vertex itself
                    vvs[vid].update(v for v in triangle if v != vid)
        return vvs

    def compute_roughness(self, tin, vvs):
        roughness_values = {}
        for vid, neighbors in vvs.items():
            # Gather elevations of the vertex and its neighbors
            elevations = [tin.get_vertex(vid).get_z()] + [
                tin.get_vertex(v).get_z() for v in neighbors
            ]

            # Calculate the mean elevation fA
            mean_elevation = np.mean(elevations)

            # Calculate the sum of squared differences from the mean
            squared_differences = [
                (elevation - mean_elevation) ** 2 for elevation in elevations
            ]
            sum_squared_differences = sum(squared_differences)

            # Number of vertices considered (vertex plus its neighbors)
            k = len(elevations) - 1

            # Calculate roughness as the square root of the average of squared differences
            roughness = np.sqrt(sum_squared_differences / k)
            roughness_values[vid] = roughness

        return roughness_values

    def calculate_roughness(self, tin):
        if not self.is_leaf():
            return {}  # Only calculate roughness for leaf nodes
        vts = self.extract_vt_relations(tin)
        vvs = self.compute_vv_relations(tin, vts)
        return self.compute_roughness(tin, vvs)

    def angle(self, v1, v2, v3):
        V21 = [
            v1.get_c(0) - v2.get_c(0),
            v1.get_c(1) - v2.get_c(1),
            v1.get_z() - v2.get_z(),
        ]
        V23 = [
            v3.get_c(0) - v2.get_c(0),
            v3.get_c(1) - v2.get_c(1),
            v3.get_z() - v2.get_z(),
        ]
        dot_product = np.dot(V21, V23)
        magnitude_product = np.sqrt(np.dot(V21, V21)) * np.sqrt(np.dot(V23, V23))
        cos_angle = np.clip(
            dot_product / magnitude_product, -1.0, 1.0
        )  # Clip to prevent out of range errors
        return np.arccos(cos_angle)

    def compute_curvature(self, tin):
        vts = self.extract_vt_relations(tin)  # Extract vertex-triangle relationships
        curvatures = {}

        for vid, triangles in vts.items():
            total_angle = 0
            for tid in triangles:
                triangle = tin.get_triangle(tid)
                if triangle:
                    v1, v2, v3 = triangle
                    # Reorder vertices so that v2 is the current vertex (vid)
                    if v1 == vid:
                        v1, v2, v3 = v2, v1, v3
                    elif v3 == vid:
                        v1, v2, v3 = v2, v3, v1
                    # Compute angle at v2
                    total_angle += self.angle(
                        tin.get_vertex(v1), tin.get_vertex(v2), tin.get_vertex(v3)
                    )

            is_boundary = self.is_boundary_vertex(vid, vts[vid], tin)
            # Compute curvature based on boundary condition
            if is_boundary:
                curvatures[vid] = np.pi - total_angle
            else:
                curvatures[vid] = 2 * np.pi - total_angle

        return curvatures

    def is_boundary_vertex(self, vid, triangles, tin):
        adjacent_vertices = set()
        for tid in triangles:
            triangle = tin.get_triangle(tid)
            if triangle:
                for v in triangle:
                    if v != vid:
                        if v in adjacent_vertices:
                            adjacent_vertices.remove(v)
                        else:
                            adjacent_vertices.add(v)
        return len(adjacent_vertices) > 0

    def compute_maxima(self, tin):
        if not self.is_leaf():
            return {}  # Maxima are computed only at leaf nodes

        vts = self.extract_vt_relations(tin)  # Reuse the Vertex-Triangle relations
        vvs = self.compute_vv_relations(tin, vts)  # Extract Vertex-Vertex relations
        maxima = []

        for vid, neighbors in vvs.items():
            current_vertex_elevation = tin.get_vertex(vid).get_z()
            if all(
                current_vertex_elevation > tin.get_vertex(nid).get_z()
                for nid in neighbors
            ):
                maxima.append(
                    vid
                )  # If all neighbor elevations are less, then it's a maximum

        return maxima
