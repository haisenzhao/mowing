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



	

	void  ToolpathGenerator::FillingAlgorithmBasedOnOffsets()
	{
		for (int i = 0; i < smooth_number; i++)
			PolygonSmoothing();

		Circuit::ComputeOffsets(toolpath_size, contours, offsets, offsetses);

		//build offset graph
		#pragma region build_offset_graph
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
		#pragma endregion

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
		Tree::DecompositionATree(nodes, mst, connect_nodes,connect_edges, decompose_offset);
		Output_tree(nodes, connect_edges, "D:\\3.gml");

		std::vector<std::vector<double>> cutting_points;

		for (int i = 0; i < connect_nodes.size(); i++)
		{
			cutting_points.push_back(std::vector<double>());
		}

		//turning_points_entry
		for (int i = 0; i < decompose_offset.size(); i++)
		{
			if (debug_int_0 >= 0)
			{
				if (i != debug_int_0)
					continue;
			}

			if (decompose_offset[i][0] == 0)
			{
				//get current offsets
				std::vector<std::vector<Vector2d>> offsets0;
				for (int l = 0; l < decompose_offset[i].size(); l++)
				{
					offsets0.push_back(offsets[decompose_offset[i][l]]);
				}
				int next_index = Tree::NextNode(mst, decompose_offset[i]);

				Vector2d entry_point;
				Vector2d exit_point;
				Circuit::FindOptimalEntryExitPoints(toolpath_size, offsets[decompose_offset[i][0]],entry_point, exit_point);

				//generate Fermat spiral
				FermatsSpiralTrick(offsets0, entry_point, exit_point);

				entry_spirals.push_back(entry_spiral);
				exit_spirals.push_back(exit_spiral);

				//get connecting points
				double entry_d = Circuit::FindNearestPointPar(entry_spiral[entry_spiral.size()-1],offsets[next_index]);
				Vector2d connecting_entry_point = Circuit::GetOnePointFromOffset(entry_d, offsets[next_index]);

				double exit_d = Circuit::FindNearestPointPar(exit_spiral[exit_spiral.size()-1],offsets[next_index]);
				Vector2d connecting_exit_point = Circuit::GetOnePointFromOffset(exit_d, offsets[next_index]);

				turning_points_entry.push_back(connecting_entry_point);
				turning_points_exit.push_back(connecting_exit_point);

				//input cutting points
				Tree::InputCuttingPoints(connect_nodes, cutting_points, next_index, entry_d, exit_d);

				//inpute connecting lines
				std::vector<Vector2d> one_path;
				one_path.push_back(connecting_entry_point);
				one_path.push_back(entry_spiral[entry_spiral.size() - 1]);
				pathes.push_back(one_path);
				std::vector<Vector2d>().swap(one_path);
				one_path.push_back(connecting_exit_point);
				one_path.push_back(exit_spiral[exit_spiral.size() - 1]);
				pathes.push_back(one_path);
				std::vector<Vector2d>().swap(one_path);

				std::vector<std::vector<Vector2d>>().swap(offsets0);
				std::vector<Vector2d>().swap(entry_spiral);
				std::vector<Vector2d>().swap(exit_spiral);

			}
			else
			{
				std::vector<std::vector<Vector2d>> offsets0;

				for (int l = decompose_offset[i].size() - 1; l >= 0; l--)
				{
					offsets0.push_back(offsets[decompose_offset[i][l]]);
				}

				int next_index = Tree::NextNode(mst, decompose_offset[i]);

				Vector2d entry_point;
				Vector2d exit_point;
				Circuit::FindOptimalEntryExitPoints(toolpath_size,offsets[decompose_offset[i][decompose_offset[i].size() - 1]],  entry_point, exit_point);

				//get connecting points
				double entry_d = Circuit::FindNearestPointPar(entry_point, offsets[next_index]);
				Vector2d connecting_entry_point = Circuit::GetOnePointFromOffset(entry_d, offsets[next_index]);

				double exit_d = Circuit::FindNearestPointPar(exit_point, offsets[next_index]);
				Vector2d connecting_exit_point = Circuit::GetOnePointFromOffset(exit_d, offsets[next_index]);


				//input cutting points
				bool b = Circuit::CheckEnclosed(offsets[next_index], offsets[decompose_offset[i][decompose_offset[i].size() - 1]]);
				if (b)
				{
					Tree::InputCuttingPoints(connect_nodes, cutting_points, next_index, exit_d, entry_d);

					turning_points_entry.push_back(connecting_exit_point);
					turning_points_exit.push_back(connecting_entry_point);
				}
				else
				{
					Tree::InputCuttingPoints(connect_nodes, cutting_points, next_index, entry_d, exit_d);

					turning_points_entry.push_back(connecting_entry_point);
					turning_points_exit.push_back(connecting_exit_point);
				}

				//inpute connecting lines
				std::vector<Vector2d> one_path;
				one_path.push_back(connecting_entry_point);
				one_path.push_back(entry_point);
				pathes.push_back(one_path);
				std::vector<Vector2d>().swap(one_path);
				one_path.push_back(connecting_exit_point);
				one_path.push_back(exit_point);
				pathes.push_back(one_path);
				std::vector<Vector2d>().swap(one_path);

				//generate Fermat spiral
				FermatsSpiralTrick(offsets0, entry_point, exit_point);
				std::vector<std::vector<Vector2d>>().swap(offsets0);

				entry_spirals.push_back(entry_spiral);
				exit_spirals.push_back(exit_spiral);

				std::vector<Vector2d>().swap(entry_spiral);
				std::vector<Vector2d>().swap(exit_spiral);
			}
		}

		//connect sub-region to offset mst tree

		for (int i = 1; i < connect_nodes.size(); i++)
		{
			Circuit::CutACircuit(offsets[connect_nodes[i]], cutting_points[i], pathes);
		}

		for (int i = 2; i < connect_edges.size(); i=i+2)
		{
			int node_index_0 = connect_edges[i];
			int node_index_1 = connect_edges[i + 1];

			int cutting_index_0 = -1;
			int cutting_index_1 = -1;

			for (int j = 0; j < connect_nodes.size() && (cutting_index_0<0 || cutting_index_1<0); j++)
			{
				if (connect_nodes[j] == node_index_0)
				{
					cutting_index_0 = j;
				}

				if (connect_nodes[j] = node_index_1)
				{
					cutting_index_1 = j;
				}
			}

			if (cutting_index_0 >= 0 && cutting_index_1 >= 0)
			{
				Circuit::ConnectTwoTrunkNodes(toolpath_size, offsets[node_index_0], offsets[node_index_1],
					cutting_points[cutting_index_0], cutting_points[cutting_index_1], pathes_temp, turning_points_entry_temp, turning_points_exit_temp);
			}
		
			break;
		}


		std::vector<int>().swap(connect_edges);
		std::vector<int>().swap(connect_nodes);
		std::vector<double>().swap(costs);
		std::vector<int>().swap(mst);
		std::vector<int>().swap(nodes);
		std::vector<int>().swap(edges);
	}


}