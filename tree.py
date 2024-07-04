from node import Node


class Tree(object):
    """Creates Class tree"""

    def __init__(self, c):
        self.__root = Node()
        self.__capacity = c

    def get_root(self):
        return self.__root

    def get_leaf_threshold(self):
        return self.__capacity

    def build_tree(self, tin):
        # first we insert the vertices of the TIN
        for i in range(tin.get_vertices_num()):
            # print(
            #     "INSERT POINT %s" % tin.get_vertex(i)
            # )  ## you can use this line to check the vertex input. Can comment it out if you don't need it.
            self.insert_vertex(self.__root, 0, tin.get_domain(), i, tin)

        # print("INSERT POINTS COMPLETED")
        for i in range(tin.get_triangles_num()):
            triangle = tin.get_triangle(i)
            # print(f"INSERT TRIANGLE {i} - {triangle[0]} {triangle[1]} {triangle[2]}")
            vertices = [tin.get_vertex(v) for v in triangle]
            inserted_blocks = (
                set()
            )  # To track which blocks a triangle has been inserted into
            for vertex_index, vertex in enumerate(vertices):
                block = self.point_query(
                    self.get_root(), 0, tin.get_domain(), vertex, tin
                )
                if block and block not in inserted_blocks:
                    block.add_triangle_id(i)
                    inserted_blocks.add(block)

        # print("INSERT TRIANGLES END")

        # self.print_tree()

        #  End of the build_tree() function

    def insert_vertex(self, node, node_label, node_domain, v_index, tin):
        if node_domain.contains_point(
            tin.get_vertex(v_index), tin.get_domain().get_max_point()
        ):
            if node.is_leaf():
                if node.is_duplicate(
                    v_index, tin
                ):  # if the inserting vertex is the same as one vertex in the tree.
                    return  # do not insert it
                node.add_vertex(v_index)  # update append list
                if node.overflow(self.__capacity):
                    # WE HAVE TO PERFORM A SPLIT OPERATION WHEN THE NUMBER OF VERTICES EXCEED CAPACITY
                    # current node become internal, and we initialize its children
                    node.init_children()
                    for i in node.get_vertices():
                        self.insert_vertex(node, node_label, node_domain, i, tin)
                    node.reset_vertices()  # empty the list of the current node

            else:  # otherwise we are visiting an INTERNAL node
                mid_point = node_domain.get_centroid()
                for i in range(4):
                    s_label, s_domain = node.compute_child_label_and_domain(
                        i, node_label, node_domain, mid_point
                    )
                    self.insert_vertex(
                        node.get_child(i), s_label, s_domain, v_index, tin
                    )

    def point_query(self, node, node_label, node_domain, search_point, tin):
        # node: Node object; node_label: int; node_domain: Domain object;search_point: Vertex object, the vertex you want to search
        # when point_query used in other functions for searching point:
        # node is the root of the tree,node_label is the node label of the root node(0), node_domain is the domain of the TIN(tin.get_domain()).
        # You will use this for identifying nodes containing the extreme vertices of a triangle
        #
        # This function will return the node that contains the input search_point.
        if node_domain.contains_point(search_point, tin.get_domain().get_max_point()):
            if node.is_leaf():
                isfound = False
                x = None  # x is the point id
                for (
                    i
                ) in (
                    node.get_vertices()
                ):  # for each point index in each node, if that point is equal to the query point, then it is found. Otherwise, it is not found.
                    if (
                        tin.get_vertex(i) == search_point
                    ):  # here tin.get_vertex(i) and search_point are Vertex objects
                        isfound = True
                        x = i  # x is the point index that is equal to the search point.)
                        # print("Vertex " + str(x) + " is in Node " + str(node_label))
                        break
                if isfound:
                    return node
                else:
                    return None
            else:  # Internal node
                ### we visit the children in the following order: NW -> NE -> SW -> SE
                mid_point = node_domain.get_centroid()
                s_label, s_domain = node.compute_child_label_and_domain(
                    0, node_label, node_domain, mid_point
                )
                ret_node = self.point_query(
                    node.get_child(0), s_label, s_domain, search_point, tin
                )
                if ret_node is not None:
                    return ret_node
                else:
                    s_label, s_domain = node.compute_child_label_and_domain(
                        1, node_label, node_domain, mid_point
                    )
                    ret_node = self.point_query(
                        node.get_child(1), s_label, s_domain, search_point, tin
                    )
                    if ret_node is not None:
                        return ret_node
                    else:
                        s_label, s_domain = node.compute_child_label_and_domain(
                            2, node_label, node_domain, mid_point
                        )
                        ret_node = self.point_query(
                            node.get_child(2), s_label, s_domain, search_point, tin
                        )
                        if ret_node is not None:
                            return ret_node
                        else:
                            s_label, s_domain = node.compute_child_label_and_domain(
                                3, node_label, node_domain, mid_point
                            )
                            ret_node = self.point_query(
                                node.get_child(3), s_label, s_domain, search_point, tin
                            )
                            return ret_node

    def get_points(self, tin, pts):
        """return xs,ys"""
        xs = list()
        ys = list()
        for v in pts:
            xs.append(tin.get_vertex(v).get_x())
            ys.append(tin.get_vertex(v).get_y())
        return xs, ys

    def print_tree(self, node=None, label=0):
        """Prints the Quadtree in a pre-order traversal manner."""
        if node is None:
            node = self.get_root()
            print("START TRIANGLE PR")

        # Check if the node is a leaf or an internal node
        if node.is_leaf():
            # Check if the leaf is empty or full
            if node.get_vertices_num() == 0:
                print(f"{label} EMPTY LEAF")
            else:
                # The leaf is full, print vertex and triangle information
                vertices = ", ".join(str(v) for v in node.get_vertices())
                # Ensure you implement a method in Node to get triangle IDs similar to get_vertices
                triangles = ", ".join(
                    str(node.get_triangle_id(i))
                    for i in range(node.get_triangles_num())
                )
                print(f"{label} FULL LEAF")
                print(f"V {node.get_vertices_num()} [{vertices}]")
                print(f"T {node.get_triangles_num()} [{triangles}]")
        else:
            print(f"{label} INTERNAL")
            # Recursively print each child (assumes children are initialized to a list of None or Node objects)
            if not node.is_leaf():
                for i in range(4):  # Assuming 4 children: NW, NE, SW, SE
                    child = node.get_child(i)
                    if child is not None:
                        self.print_tree(child, 4 * label + i + 1)

        if label == 0:  # Finished printing the whole tree
            print("END TRIANGLE PR")

    def compute_features(self, tin):
        roughness = {}
        curvature = {}
        maxima = []
        self._compute_node_features(self.__root, roughness, curvature, maxima, tin)
        return roughness, curvature, maxima

    def _compute_node_features(self, node, roughness, curvature, maxima, tin):
        if node.is_leaf():
            node_roughness = node.calculate_roughness(tin)
            node_curvature = node.compute_curvature(tin)
            node_maxima = node.compute_maxima(tin)
            roughness.update(node_roughness)
            curvature.update(node_curvature)
            maxima.extend(node_maxima)
        else:
            # Recursively compute features for all children
            for i in range(4):  # Assuming 4 children: NW, NE, SW, SE
                child = node.get_child(i)
                if child:
                    self._compute_node_features(
                        child, roughness, curvature, maxima, tin
                    )
