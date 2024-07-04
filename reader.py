from tin import TIN
from vertex import Vertex


class Reader(object):
    def read_tin_file(self, path):
        with open(path, "r") as infile:
            tin = TIN()
            infile.readline()  # Skip the 'OFF' header

            # Read the line containing counts of vertices and triangles
            line = infile.readline().strip().split()
            vertices_num = int(line[0])
            triangles_num = (
                int(line[1]) if len(line) > 1 else 0
            )  # Ensure there is a second element

            # Read vertices
            for _ in range(vertices_num):
                line = infile.readline().strip().split()
                v = Vertex(float(line[0]), float(line[1]), float(line[2]))
                tin.add_vertex(v)
                # Set domain using the first vertex
                if tin.get_vertices_num() == 1:
                    tin.set_domain(v, v)

            # Adjust domain as more vertices are read
            if vertices_num > 1:
                for v in tin.get_vertices()[1:]:
                    tin.get_domain().resize(v)

            # Read triangles
            # print("START INPUT TRIANGLES")
            for i in range(triangles_num):
                line = infile.readline().strip().split()
                if len(line) >= 4 and line[0] == "3":  # Confirming it's a triangle line
                    v1, v2, v3 = map(int, line[1:])
                    if all(v < vertices_num for v in [v1, v2, v3]):
                        tin.add_triangle((v1, v2, v3))
                #     else:
                #         print(
                #             f"Warning: Triangle {i} references out-of-bounds vertex index."
                #         )
                # else:
                #     print(f"Warning: Incorrect format or missing data at triangle {i}.")

            # print(f"Total number of triangles processed: {triangles_num}")
            # print("END INPUT TRIANGLES")

            return tin
