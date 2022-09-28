def write_buttom_leafs(octree, file_name="octree.vtu"):
	f = open(file_name, "w")
	
	_write_header(f)
	_wrtie_body(f, octree)
	_write_footer(f)
	
	f.close()
	
def _write_header(f):
	f.write("<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
	
def _wrtie_body(f, octree):
	leaf_ids = []
	for i in range(len(octree.leafs)):
		if octree.leafs[i].number_children == 0:
			leaf_ids.append(i)

	f.write("<UnstructuredGrid>\n")
	f.write("<Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n".format(len(leaf_ids) * 8, len(leaf_ids)))

	_write_points(f, octree, leaf_ids);
	_write_cells(f, octree, leaf_ids);

	f.write("</Piece>\n")
	f.write("</UnstructuredGrid>\n")
	
def _write_points(f, octree, leaf_ids):
	f.write("<Points>\n")
	f.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")

	for i in leaf_ids:
		box = octree.cell_boxes[i]

		f.write("{} ".format(box[0][0]))
		f.write("{} ".format(box[1][0]))
		f.write("{} ".format(box[2][0]))
		
		f.write("{} ".format(box[0][1]))
		f.write("{} ".format(box[1][0]))
		f.write("{} ".format(box[2][0]))
		
		f.write("{} ".format(box[0][1]))
		f.write("{} ".format(box[1][1]))
		f.write("{} ".format(box[2][0]))
		
		f.write("{} ".format(box[0][0]))
		f.write("{} ".format(box[1][1]))
		f.write("{} ".format(box[2][0]))
		
		f.write("{} ".format(box[0][0]))
		f.write("{} ".format(box[1][0]))
		f.write("{} ".format(box[2][1]))
		
		f.write("{} ".format(box[0][1]))
		f.write("{} ".format(box[1][0]))
		f.write("{} ".format(box[2][1]))
		
		f.write("{} ".format(box[0][1]))
		f.write("{} ".format(box[1][1]))
		f.write("{} ".format(box[2][1]))
		
		f.write("{} ".format(box[0][0]))
		f.write("{} ".format(box[1][1]))
		f.write("{} ".format(box[2][1]))

	f.write("</DataArray>\n")
	f.write("</Points>\n")
	
	
def _write_cells(f, octree, leaf_ids):
	k = 0
	
	f.write("<Cells>\n")
	f.write("<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n")

	for i in range(len(leaf_ids)):
		for _ in range(8):
			f.write("{} ".format(k))
			k += 1

	f.write("</DataArray>\n")
	f.write("<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n")

	for i in range(len(leaf_ids)):
		f.write("{} ".format((i + 1) * 8))

	f.write("</DataArray>\n")
	f.write("<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n")

	for _ in range(len(leaf_ids)):
		f.write("12 ")

	f.write("</DataArray>\n")
	f.write("</Cells>\n")
	
def _write_footer(f):
	f.write("</VTKFile>\n")
