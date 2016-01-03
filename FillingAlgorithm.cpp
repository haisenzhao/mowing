#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"

#include "Circuit.h"
#include "tree.h"

namespace hpcg {

	void ToolpathGenerator::Output_tree(std::string path)
	{
		std::ofstream file(path);

		file << "Mark Newman on Sat Jul 22 05:32:16 2006" << std::endl;
		file << "graph" << std::endl;
		file << "[" << std::endl;
		file << "  directed 0" << std::endl;

		for (int i = 0; i < offsets.size(); i++)
		{
			file << "node" << std::endl;
			file << "[" << std::endl;
			file << "id " << i << std::endl;
			file << "label " << i << std::endl;
			file << "]" << std::endl;
		}


		for (int i = 0; i < offset_graph.size(); i = i + 2)
		{
			int index_0 = offset_graph[i];
			int index_1 = offset_graph[i+1];

			file << "edge" << std::endl;
			file << "[" << std::endl;

			file << "source " << index_0 << std::endl;
			file << "target " << index_1 << std::endl;

			file << "]" << std::endl;
		}

		file << "]" << std::endl;

		file.clear();
		file.close();
	}
	
	void ToolpathGenerator::Output_tree(std::vector<int> &nodes, std::vector<int> &edges, std::string path)
	{
		std::ofstream file(path);

		file << "Mark Newman on Sat Jul 22 05:32:16 2006" << std::endl;
		file << "graph" << std::endl;
		file << "[" << std::endl;
		file << "  directed 0" << std::endl;


		for (int i = 0; i < nodes.size(); i++)
		{
			file << "node" << std::endl;
			file << "[" << std::endl;
			file << "id " << nodes[i] << std::endl;
			file << "label " << nodes[i] << std::endl;

			file << "]" << std::endl;
		}

		for (int i = 0; i < edges.size(); i=i+2)
		{
			file << "edge" << std::endl;
			file << "[" << std::endl;

			file << "source " << edges[i] << std::endl;
			file << "target " << edges[i+1] << std::endl;

			file << "]" << std::endl;
		}

		file << "]" << std::endl;

		file.clear();
		file.close();
	}

	void ToolpathGenerator::BuildOffsetGraph()
	{
		//build offset graph
		std::vector<int> index_int;

		int index = -1;
		for (int i = 0; i < offsetses.size(); i++)
		{
			for (int j = 0; j < offsetses[i].size(); j++)
			{
				index++;
				index_int.push_back(i);
				index_int.push_back(j);
				index_int.push_back(index);
			}
		}

		for (int i = 0; i < offsetses.size() - 1; i++)
		{
			for (int j = 0; j < offsetses[i].size(); j++)
			{
				for (int k = 0; k < offsetses[i + 1].size(); k++)
				{
					double d = Circuit::Distance(offsetses[i + 1][k], offsetses[i][j]);
					if (abs(d - toolpath_size) < 0.0001)
					{
						int index_0 = -1;
						int index_1 = -1;

						for (int l = 0; l < index_int.size() && (index_0<0 || index_1<0); l = l + 3)
						{
							if (index_int[l] == i&&index_int[l + 1] == j)
							{
								index_0 = index_int[l + 2];
							}

							if (index_int[l] == i + 1 && index_int[l + 1] == k)
							{
								index_1 = index_int[l + 2];
							}
						}
						offset_graph.push_back(index_0);
						offset_graph.push_back(index_1);
					}
				}
			}
		}
		std::vector<int>().swap(index_int);

	}


	int GetTrunkNodeId(std::vector<TrunkNode> &trunk_nodes, int related_offset_id)
	{
		int trunk_node_id = -1;

		for (int i = 0; i < trunk_nodes.size(); i++)
		{
			if (trunk_nodes[i].related_offset_id == related_offset_id)
			{
				trunk_node_id = i;
			}
		}
		
		return trunk_node_id;
	}

	void  ToolpathGenerator::FillingAlgorithmBasedOnOffsets()
	{
		for (int i = 0; i < smooth_number; i++)
			PolygonSmoothing();

		Circuit::ComputeOffsets(toolpath_size, contours, offsets, offsetses);


		return;

		BuildOffsetGraph();

		//minimal spanning tree
		std::vector<int> mst;

		std::vector<int> nodes;
		std::vector<int> edges;
		std::vector<double> costs;
		#pragma region minimal_spanning_tree

		for (int i = 0; i < offsets.size(); i++)
		{
			nodes.push_back(i);
		}

		for (int i = 0; i < offset_graph.size(); i = i + 2)
		{
			int index_0 = offset_graph[i];
			int index_1 = offset_graph[i+1];
			edges.push_back(index_0);
			edges.push_back(index_1);
			costs.push_back(1);
		}

		Output_tree(nodes, edges, "D:\\1.gml");

		Tree::MinimalSpanningTree(nodes, edges, costs, mst);
		Output_tree(nodes, mst, "D:\\2.gml");

		#pragma endregion
		
		//decompose_offset_tree
		std::vector<int> connect_edges;
		std::vector<int> connect_nodes;
		std::vector<TrunkNode> trunk_nodes;
		Tree::DecompositionATree(nodes, mst, connect_nodes,connect_edges, decompose_offset);

		for (int i = 0; i < connect_nodes.size(); i++)
		{
			trunk_nodes.push_back(TrunkNode(connect_nodes[i]));
		}

		Output_tree(nodes, connect_edges, "D:\\3.gml");


		//all cutting points are on trunk node
		std::vector<std::vector<double>> cutting_points;

		for (int i = 0; i < connect_nodes.size(); i++)
		{
			cutting_points.push_back(std::vector<double>());
		}

		std::vector<int> connecting_pathes;

		//turning_points_entry
		for (int i = 0; i < decompose_offset.size(); i++)
		{
			/*
			if (debug_int_0 >= 0)
			{
				if (i != debug_int_0)
					continue;
			}
			*/

			if (decompose_offset[i][0] == 0)
			{
				//get current offsets
				std::vector<std::vector<Vector2d>> offsets0;
				for (int l = 0; l < decompose_offset[i].size(); l++)
				{
					offsets0.push_back(offsets[decompose_offset[i][l]]);
				}
				int next_index = Tree::NextNode(mst, decompose_offset[i]);

				Vector2d input_entry_point, input_exit_point;
				Vector2d output_entry_point, output_exit_point;
				Circuit::FindOptimalEntryExitPoints(toolpath_size, offsets[decompose_offset[i][0]], input_entry_point, input_exit_point);

				//generate Fermat spiral
				FermatsSpiralTrick(offsets0, input_entry_point, input_exit_point, output_entry_point, output_exit_point);
		
				//get connecting points
				Vector2d connecting_entry_point, connecting_exit_point;
				Circuit::ComputeNextEntryExitPointForInner(toolpath_size, offsets0[offsets0.size() - 1], offsets[next_index],
					output_entry_point, output_exit_point, connecting_entry_point, connecting_exit_point);
				double entry_d = Circuit::FindNearestPointPar(connecting_entry_point, offsets[next_index]);
				double exit_d = Circuit::FindNearestPointPar(connecting_exit_point, offsets[next_index]);
	
				//exit_spiral.push_back(connecting_entry_point);
				//entry_spiral.push_back(connecting_exit_point);

				entry_spirals.push_back(entry_spiral);
				exit_spirals.push_back(exit_spiral);

				std::reverse(exit_spiral.begin(), exit_spiral.end());

				pathes.push_back(entry_spiral);
				pathes.push_back(exit_spiral);

				//connecting_pathes.push_back(pathes.size() - 2);
				//connecting_pathes.push_back(pathes.size() - 1);

				//connecting_pathes.push_back(pathes.size() - 2);
				//connecting_pathes.push_back(pathes.size() - 1);

				//input cutting points
				int trunk_node_id = GetTrunkNodeId(trunk_nodes, next_index);

				trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.push_back(exit_d);
				trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.push_back(entry_d);

				trunk_nodes[trunk_node_id].connecting_leaf_nodes_spiral_id.push_back(pathes.size() - 2);
				trunk_nodes[trunk_node_id].connecting_leaf_nodes_spiral_id.push_back(pathes.size() - 1);

				std::vector<std::vector<Vector2d>>().swap(offsets0);
				std::vector<Vector2d>().swap(entry_spiral);
				std::vector<Vector2d>().swap(exit_spiral);
			}
			else
			{
				//get current offsets
				std::vector<std::vector<Vector2d>> offsets0;
				for (int l = decompose_offset[i].size() - 1; l >= 0; l--)
				{
					offsets0.push_back(offsets[decompose_offset[i][l]]);
				}

				int next_index = Tree::NextNode(mst, decompose_offset[i]);

				Vector2d input_entry_point, input_exit_point;
				Vector2d output_entry_point, output_exit_point;
				Circuit::FindOptimalEntryExitPoints(toolpath_size, offsets[decompose_offset[i][decompose_offset[i].size() - 1]], input_entry_point, input_exit_point);

				//generate Fermat spiral
				FermatsSpiralTrick(offsets0, input_entry_point, input_exit_point, output_entry_point, output_exit_point);

				//get connecting points

				Vector2d connecting_entry_point, connecting_exit_point;
				Circuit::ComputeNextEntryExitPointForOuter(toolpath_size, offsets0[0], offsets[next_index],
					input_entry_point, input_exit_point, connecting_entry_point, connecting_exit_point);
				double entry_d = Circuit::FindNearestPointPar(connecting_entry_point, offsets[next_index]);
				double exit_d = Circuit::FindNearestPointPar(connecting_exit_point, offsets[next_index]);

				/*
				double entry_d = Circuit::FindNearestPointPar(input_entry_point, offsets[next_index]);
				Vector2d connecting_entry_point = Circuit::GetOnePointFromOffset(entry_d, offsets[next_index]);
				double exit_d = Circuit::FindNearestPointPar(input_exit_point, offsets[next_index]);
				Vector2d connecting_exit_point = Circuit::GetOnePointFromOffset(exit_d, offsets[next_index]);
				*/

				//inpute connecting lines

				//entry_spiral.insert(entry_spiral.begin(), connecting_entry_point);
				//exit_spiral.insert(exit_spiral.begin(), connecting_exit_point);

				entry_spirals.push_back(entry_spiral);
				exit_spirals.push_back(exit_spiral);

				//std::reverse(entry_spiral.begin(), entry_spiral.end());
				std::reverse(exit_spiral.begin(), exit_spiral.end());

				pathes.push_back(entry_spiral);
				pathes.push_back(exit_spiral);

				//connecting_pathes.push_back(pathes.size() - 2);
				//connecting_pathes.push_back(pathes.size() - 1);

				//input cutting points
				bool b = Circuit::CheckEnclosed(offsets[next_index], offsets[decompose_offset[i][decompose_offset[i].size() - 1]]);

				int trunk_node_id = GetTrunkNodeId(trunk_nodes, next_index);

				if (b)
				{
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.push_back(exit_d);
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.push_back(entry_d);
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_spiral_id.push_back(-(pathes.size() - 1));
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_spiral_id.push_back(pathes.size() - 2);
		
					//turning_points_entry.push_back(connecting_exit_point);
					//turning_points_exit.push_back(connecting_entry_point);

					connecting_pathes.push_back(-(pathes.size() - 1));
					connecting_pathes.push_back(pathes.size() - 2);

				}
				else
				{
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.push_back(entry_d);
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.push_back(exit_d);
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_spiral_id.push_back(pathes.size() - 2);
					trunk_nodes[trunk_node_id].connecting_leaf_nodes_spiral_id.push_back(-(pathes.size() - 1));
				
					//turning_points_entry.push_back(connecting_entry_point);
					//turning_points_exit.push_back(connecting_exit_point);

					connecting_pathes.push_back(pathes.size() - 2);
					connecting_pathes.push_back(-(pathes.size() - 1));
				}

				std::vector<std::vector<Vector2d>>().swap(offsets0);
				std::vector<Vector2d>().swap(entry_spiral);
				std::vector<Vector2d>().swap(exit_spiral);
			}
		}

		
		//connect sub-region to offset mst tree
		for (int i = 0; i < connect_edges.size(); i=i+2)
		{
			int node_index_0 = connect_edges[i];
			int node_index_1 = connect_edges[i + 1];

			int trunk_node_id_0 = GetTrunkNodeId(trunk_nodes, connect_edges[i]);
			int trunk_node_id_1 = GetTrunkNodeId(trunk_nodes, connect_edges[i+1]);
			
			Circuit::ConnectTwoTrunkNodes(toolpath_size, offsets[node_index_0], offsets[node_index_1],
				trunk_node_id_0, trunk_node_id_1,
				trunk_nodes[trunk_node_id_0], trunk_nodes[trunk_node_id_1],
				pathes, pathes_temp, turning_points_entry_temp, turning_points_exit_temp);
		}

		for (int i = 0; i < trunk_nodes.size(); i++)
		{
			trunk_nodes[i].AllConnectingPoints();
			Circuit::CutACircuit(offsets[trunk_nodes[i].related_offset_id],trunk_nodes[i]);
		}

		//get one single path


		#pragma region one_single_path

		int pathes_nb = pathes.size();
		
		//generate cuting pathes by these connecting points
		for (int i = 0; i < trunk_nodes.size(); i++)
		{
			for (int j = 0; j < trunk_nodes[i].connecting_points_pairs.size(); j = j + 2)
			{
				int connecting_points_index_0 = trunk_nodes[i].connecting_points_pairs[j];
				int connecting_points_index_1 = trunk_nodes[i].connecting_points_pairs[j + 1];
				std::vector<Vector2d> one_path;
				Circuit::SelectOnePartOffset(offsets[trunk_nodes[i].related_offset_id], trunk_nodes[i].connecting_points[connecting_points_index_0], trunk_nodes[i].connecting_points[connecting_points_index_1], one_path);
				pathes.push_back(one_path);

				pathes_temp.push_back(one_path);

				trunk_nodes[i].connecting_points_related_pathes_id[connecting_points_index_0] = pathes.size() - 1;
				trunk_nodes[i].connecting_points_related_pathes_id[connecting_points_index_1] = -(pathes.size() - 1);
			}
		}

		//compute edges connecting pathes.
		for (int i = 0; i < trunk_nodes.size(); i++)
		{
			for (int j = 0; j < trunk_nodes[i].connecting_points_pairs.size(); j = j + 2)
			{
				int connecting_points_index_0 = trunk_nodes[i].connecting_points_pairs[j];
				int connecting_points_index_1 = trunk_nodes[i].connecting_points_pairs[j + 1];

				int connecting_leaf_nodes_points_nb = trunk_nodes[i].connecting_leaf_nodes_points.size();

				if (connecting_points_index_0 < connecting_leaf_nodes_points_nb)
				{
					connecting_pathes.push_back(trunk_nodes[i].connecting_points_related_pathes_id[connecting_points_index_0]);
					connecting_pathes.push_back(trunk_nodes[i].connecting_leaf_nodes_spiral_id[connecting_points_index_0]);
					//done
				}
				else
				{
					connecting_pathes.push_back(trunk_nodes[i].connecting_points_related_pathes_id[connecting_points_index_0]);

					connecting_points_index_0 = connecting_points_index_0 - connecting_leaf_nodes_points_nb;

					int trunk_node_id = trunk_nodes[i].connecting_trunk_nodes_id[2 * connecting_points_index_0];
					int trunk_nodes_points_id = trunk_nodes[i].connecting_trunk_nodes_id[2 * connecting_points_index_0+1];

					int pathes_id=trunk_nodes[trunk_node_id].connecting_points_related_pathes_id[trunk_nodes_points_id + trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.size()];
			
					connecting_pathes.push_back(pathes_id);
				}

				if (connecting_points_index_1 < connecting_leaf_nodes_points_nb)
				{
					connecting_pathes.push_back(trunk_nodes[i].connecting_points_related_pathes_id[connecting_points_index_1]);
					connecting_pathes.push_back(trunk_nodes[i].connecting_leaf_nodes_spiral_id[connecting_points_index_1]);
					//done
				}
				else
				{
					connecting_pathes.push_back(trunk_nodes[i].connecting_points_related_pathes_id[connecting_points_index_1]);

					connecting_points_index_1 = connecting_points_index_1 - connecting_leaf_nodes_points_nb;

					int trunk_node_id = trunk_nodes[i].connecting_trunk_nodes_id[2 * connecting_points_index_1];
					int trunk_nodes_points_id = trunk_nodes[i].connecting_trunk_nodes_id[2 * connecting_points_index_1 + 1];

					int pathes_id = trunk_nodes[trunk_node_id].connecting_points_related_pathes_id[trunk_nodes_points_id + trunk_nodes[trunk_node_id].connecting_leaf_nodes_points.size()];

				
					connecting_pathes.push_back(pathes_id);
				}
			}
		}

		std::vector<bool> pathes_used;

		for (int i = 0; i < pathes.size(); i++)
		{
			pathes_used.push_back(false);
		}

		int start_pathes_index = 0;


		std::vector<int> int_save;

		int iter = 0;
		do
		{
			if (iter == 24)
			{
 				int dasd = 0;
			}

			if (iter == debug_int_0)
			{
				break;
			}

			iter++;
			pathes_used[abs(start_pathes_index)] = true;
			int_save.push_back(start_pathes_index);


			if (abs(start_pathes_index) >= pathes_nb)
			{
				if (start_pathes_index >= 0)
				{
					for (int i = 0; i < pathes[abs(start_pathes_index)].size(); i++)
					{
						one_single_path.push_back(pathes[abs(start_pathes_index)][i]);
					}
				}
				else
				{
					for (int i = pathes[abs(start_pathes_index)].size() - 1; i >= 0; i--)
					{
						one_single_path.push_back(pathes[abs(start_pathes_index)][i]);
					}
				}
			}
			else
			{
				if (start_pathes_index == 0 || start_pathes_index == 1)
				{
					for (int i = 0; i < pathes[abs(start_pathes_index)].size(); i++)
					{
						one_single_path.push_back(pathes[abs(start_pathes_index)][i]);
					}
				}
				else
				{
					int another_int = -1;
					bool b = false;
					for (int i = 0; i < connecting_pathes.size(); i++)
					{
						if (connecting_pathes[i] == start_pathes_index)
						{
							if (i % 2 == 0)
							{
								if (abs(connecting_pathes[i + 1]) < pathes_nb)
								{
									another_int = connecting_pathes[i + 1];
									b = true;
									break;
								}
							}
							if (i % 2 == 1)
							{
								if (abs(connecting_pathes[i - 1]) < pathes_nb)
								{
									another_int = connecting_pathes[i - 1];
									b = true;
									break;
								}
							}
						}
					}

					assert(b);

					if (pathes_used[abs(another_int)])
					{
						if (start_pathes_index < 0)
						{
							for (int i = 0; i < pathes[abs(start_pathes_index)].size(); i++)
							{
								one_single_path.push_back(pathes[abs(start_pathes_index)][i]);
							}
						}
						else
						{
							for (int i = pathes[abs(start_pathes_index)].size() - 1; i >= 0; i--)
							{
								one_single_path.push_back(pathes[abs(start_pathes_index)][i]);
							}
						}
					}
					else
					{
						if (start_pathes_index < 0)
						{
							for (int i = pathes[abs(start_pathes_index)].size() - 1; i >= 0; i--)
							{
								one_single_path.push_back(pathes[abs(start_pathes_index)][i]);
							}
						}
						else
						{
							for (int i = 0; i < pathes[abs(start_pathes_index)].size(); i++)
							{
								one_single_path.push_back(pathes[abs(start_pathes_index)][i]);
							}
						}
					}
				}
			}

			if (abs(start_pathes_index) >= pathes_nb)
				start_pathes_index = -start_pathes_index;

			
			int next_pathes_index = -1;
			bool b = false;

			for (int i = 0; i < connecting_pathes.size() && !b; i++)
			{
				if (connecting_pathes[i] == start_pathes_index)
				{

					if (i % 2 == 0)
					{
						if (!pathes_used[abs(connecting_pathes[i + 1])])
						{
							next_pathes_index = connecting_pathes[i + 1];
							b = true;
							break;
						}

					}
					else
					{
						if (!pathes_used[abs(connecting_pathes[i - 1])])
						{
							next_pathes_index = connecting_pathes[i - 1];
							b = true;
							break;
						}
					}
				}
			}

			if (!b)
			{
				break;
			}
			start_pathes_index = next_pathes_index;


			bool goon = false;

			for (int i = 0; i < pathes_used.size(); i++)
			{
				if (!pathes_used[i])
				{
					goon = true;
					break;
				}
			}

			if (!goon)
			{
				break;
			}

		} while (true);

		for (int i = 0; i < int_save.size(); i++)
		{
			std::cout << "Pathes Index: " << int_save[i] << std::endl;
		}

		#pragma endregion

		std::vector<std::vector<double>>().swap(cutting_points);
		std::vector<int>().swap(connect_edges);
		std::vector<int>().swap(connect_nodes);
		std::vector<double>().swap(costs);
		std::vector<int>().swap(mst);
		std::vector<int>().swap(nodes);
		std::vector<int>().swap(edges);

	}


}