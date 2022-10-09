def write_positions(particles, file_name="positions.csv"):
    with open(file_name, "w") as file:
        file.write("x, y, z\n")
        for pos in particles.Pos:
            file.write("{}, {}, {}\n".format(pos[0], pos[1], pos[2]))