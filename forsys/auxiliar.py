import os
import json

import forsys as fs
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def create_folders(folder_to_create):
	if not os.path.exists(folder_to_create):
		os.makedirs(os.path.join(folder_to_create, "connections"))
		os.makedirs(os.path.join(folder_to_create, "myosin"))
		for folder_name in ["static", "dynamic"]:
			os.makedirs(os.path.join(folder_to_create,
									 folder_name,
									 "csvs"))
			os.makedirs(os.path.join(folder_to_create,
									 folder_name,
									 "fit_per_time"))
			os.makedirs(os.path.join(folder_to_create,
									 folder_name,
									 "tissues"))
			os.makedirs(os.path.join(folder_to_create,
									 folder_name,
									 "forces"))
			os.makedirs(os.path.join(folder_to_create,
									 folder_name,
									 "pressures"))
			os.makedirs(os.path.join(folder_to_create,
									 folder_name,
									 "stress_tensor"))

def create_folders_sweep(folder_to_create):
	if not os.path.exists(folder_to_create):
		os.makedirs(os.path.join(folder_to_create, "connections"))
		os.makedirs(os.path.join(folder_to_create, "myosin"))
		os.makedirs(os.path.join(folder_to_create, "fit_per_time"))
		os.makedirs(os.path.join(folder_to_create, "csvs"))

def create_directory(name, upperFolder):
	directory = os.path.join(upperFolder, name)
	if not os.path.exists(directory):
		os.makedirs(directory)

def load_initial_guess(guess_file, min_time, max_time):
	try:
		with open(guess_file) as jfile:
			initial_guess = json.load(jfile)
			initial_guess = {int(k): {int(kin): vin for kin, vin in v.items()} for k, v in initial_guess.items()}
	except (FileNotFoundError, TypeError):
		initial_guess = {}
		print("No guess file, using zero guess")
	number_of_frames = max_time - min_time

	initial_guess = {k: {} for k in range(number_of_frames)
					 if k not in initial_guess.keys()} | initial_guess
	return initial_guess


def create_plots(frame_number, forsys, res_folder, myo=False, pressure=True, compress=1):
	vertices = forsys.frames[frame_number].vertices
	edges = forsys.frames[frame_number].edges
	cells = forsys.frames[frame_number].cells
	# mesh
	fs.plot.plot_mesh(vertices,
						edges,
						cells,
						f"mesh_{frame_number}.png",
						f"{res_folder}/tissues/",
						mirror_y=True)
	# stresses
	fs.plot.plot_inference(forsys.frames[frame_number],
							step=frame_number,
							folder=os.path.join(res_folder, "forces"),
							normalized="absolute",
							mirror_y=False,
							colorbar=False,
							compress_scale=compress)
	print("Saving to ", os.path.join(res_folder, "forces", f"{frame_number}.png"))
	plt.savefig(os.path.join(res_folder, "forces", f"{frame_number}.png"), dpi=350)
	plt.close()
	if myo:
		# myosin
		print("Plotting myosin")
		fs.plot.plot_inference(forsys.frames[frame_number],
							   ground_truth=True,
								step=frame_number,
								folder=os.path.join(res_folder, "forces"),
								normalized="absolute",
								mirror_y=False,
								colorbar=False)
		plt.savefig(os.path.join(res_folder, "myosin", f"{frame_number}.png"), dpi=350)
		plt.close()

	if pressure:
		fs.plot.plot_inference(forsys.frames[frame_number],
								step=frame_number,
								pressure=True,
								folder=os.path.join(res_folder, "pressures"),
								normalized="max",
								mirror_y=False,
								colorbar=False)
		plt.savefig(os.path.join(res_folder, "pressures", f"{frame_number}.png"), dpi=350)
		plt.close()

		fs.plot.plot_stress_tensor(forsys.frames[frame_number],
						os.path.join(res_folder, "stress_tensor"),
						frame_number,
						grid=12,
						radius=5,
						tensor_scale=1.5)


def create_csvs(forsys: fs.ForSys, time:int, mapping:bool = False) -> tuple:
	"""
    Calculate the distance between two vertices

    :forsys: Forsys object
	:type forsys: fs.Forsys
    :time: Unique time of the frame
	:type time: int
    :mapping: If True, add mapping into csv
	:type mapping: bool, optional
    """

	frame = forsys.frames[time]
	
	# be_mapped_ids = []
	be_ids = []
	tensions = []
	lengths = []
	positions_x = []
	positions_y = []
	curvatures = []
	own_cells = []
	vertices = []

	#cell_mapped_ids = []
	cell_ids = []
	areas = []
	perimeters = []
	cell_posx = []
	cell_posy = []
	pressures = []
	neighbors = []

	v_id = []
	x_arr = []
	y_arr = []
	v_cells=[]

	# if mapping and time>0:
	# 	cells_map, edge_map = forsys.get_maps(time)

	for cellid, cell in frame.cells.items():
			# if (mapping):
			# 	if(cells_map[cellid]!=None):
			# 		cell_mapped_ids.append(cells_map[cellid])
			# 	else:
			# 		cell_mapped_ids.append(1111)
			cell_ids.append(cellid)
			areas.append(abs(cell.get_area()))
			perimeters.append(cell.get_perimeter())
			cell_posx.append(cell.get_cm()[0])
			cell_posy.append(cell.get_cm()[1])
			neighbors.append(cell.neighbors)
			pressures.append(cell.pressure)

	# if (mapping):
	# 	cell_df = pd.DataFrame({
	# 		"id_mapped": cell_mapped_ids,
	# 		"id": cell_ids,
	# 		"area": areas,
	# 		"perimeter": perimeters,
	# 		"position_x": cell_posx,
	# 		"position_y": cell_posy,
	# 		"neighbors": neighbors,
	# 		"pressure": pressures,
	# 	})
	# 	cell_df.sort_values(by="id_mapped",ascending=True, inplace=True)
	if True:
		cell_df = pd.DataFrame({
			"id": cell_ids,
			"area": areas,
			"perimeter": perimeters,
			"position_x": cell_posx,
			"position_y": cell_posy,
			"neighbors": neighbors,
			"pressure": pressures,
		})
	

	for _, big_edge in frame.big_edges.items():
			# if (mapping):
			# 	if(edge_map[big_edge.big_edge_id]!=None):
			# 		be_mapped_ids.append(edge_map[big_edge.big_edge_id])
			# 	else:
			# 		be_mapped_ids.append(1111)
			be_ids.append(big_edge.big_edge_id)
			tensions.append(big_edge.tension)
			lengths.append(big_edge.get_length())
			positions_x.append(np.mean(big_edge.xs))
			positions_y.append(np.mean(big_edge.ys))
			own_cells.append(big_edge.own_cells)
			vertices.append([big_edge.vertices[0].id, big_edge.vertices[-1].id])
			curvatures.append(big_edge.calculate_total_curvature())

	# if (mapping):
	# 	force_df = pd.DataFrame({
	# 		"id_mapped":be_mapped_ids,
	# 		"id": be_ids,
	# 		"tension": tensions,
	# 		"length": lengths,
	# 		"position_x": positions_x,
	# 		"position_y": positions_y,
	# 		"own_cells": own_cells,
	# 		"vertices": vertices,
	# 		"curvature": curvatures,
	# 	})
	# 	force_df.sort_values(by="id_mapped",ascending=True, inplace=True)

	if True:
		force_df = pd.DataFrame({
			"id": be_ids,
			"tension": tensions,
			"length": lengths,
			"position_x": positions_x,
			"position_y": positions_y,
			"own_cells": own_cells,
			"vertices": vertices,
			"curvature": curvatures,
		})

	for _, vertex in frame.vertices.items():
		v_id.append(vertex.id)
		x_arr.append(vertex.x)
		y_arr.append(vertex.y)
		v_cells.append(vertex.ownCells)


	v_df = pd.DataFrame({
		"id": v_id,
		"position_x": x_arr,
		"position_y": y_arr,
		"cells":v_cells,
	})

	return cell_df, force_df, v_df
