from domain import Domain


class TIN(object):
    """Creates Class tree"""

    def __init__(self):
        self.__vertices = []
        self.__triangles = []
        self.__domain = Domain()

    def get_vertex(
        self, pos
    ):  # input: vertex index     output: The Vertex object in position pos in the vertex list. Index starts from 0.
        # should implement similar function for getting triangle from global array.
        try:
            return self.__vertices[pos]
        except IndexError as e:
            raise e

    def get_vertices_num(self):
        return len(self.__vertices)

    def get_vertices(self):
        return self.__vertices

    def get_domain(self):
        return self.__domain

    def add_vertex(self, v):  # v should be Vertex object.
        self.__vertices.append(v)

    def set_domain(
        self, min_p, max_p
    ):  # is used in read_tin_file when we read the vertices
        self.__domain = Domain(min_p, max_p)

    ### You will need to add your own functions to add/get triangles to TIN

    def add_triangle(self, triangle):  # triangle should be a tuple (v1, v2, v3)
        self.__triangles.append(triangle)

    def get_triangle(self, index):  # Get a triangle by its index
        try:
            return self.__triangles[index]
        except IndexError:
            print(f"Triangle with index {index} does not exist.")
            return None

    # Function to get the number of triangles
    def get_triangles_num(self):
        return len(self.__triangles)
