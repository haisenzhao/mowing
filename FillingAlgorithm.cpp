#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"

#include "Circuit.h"
#include "tree.h"

namespace hpcg {

	void  ToolpathGenerator::ComputeOffsetsForCircle()
	{
		double lOffset = toolpath_size / 2.0;
		std::vector<std::vector<Vector2d>> one_pathes;
		Circuit::GenerateOffsetHole(contours, lOffset, one_pathes);

		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		while (one_pathes.size() > 0)
		{
			std::cout << "Offsets index: " << offsets.size() << std::endl;

			if (one_pathes.size() == 2)
			{
				double d0 = Circuit::Distance(one_pathes[0][0], contour);
				double d1 = Circuit::Distance(one_pathes[1][0], contour);

				if (d0 < d1)
				{
					int insert_index = offsets.size() / 2.0;
					offsets.insert(offsets.begin() + insert_index, one_pathes[0]);

					std::vector<Vector2d> vecs;

					for (int i = one_pathes[1].size() - 1; i >= 0; i--)
					{
						vecs.push_back(one_pathes[1][i]);
					}
					offsets.insert(offsets.begin() + insert_index + 1, vecs);
					std::vector<Vector2d>().swap(vecs);
				}
			}

			for (int i = 0; i < one_pathes.size(); i++)
			{
				std::vector<Vector2d>().swap(one_pathes[i]);
			}
			std::vector<std::vector<Vector2d>>().swap(one_pathes);

			lOffset = lOffset + toolpath_size;
			Circuit::GenerateOffsetHole(contours, lOffset, one_pathes);
		}

		std::vector<Vector2d>().swap(contour);

	}

	void ToolpathGenerator::ComputeOffsets_temp()
	{
		double lOffset = toolpath_size / 2.0;
		std::vector<std::vector<Vector2d>> one_pathes;
		Circuit::GenerateOffsetHole(contours, lOffset, one_pathes);

		while (one_pathes.size() > 0)
		{
			std::cout << "Offsets index: " << offsets.size() << std::endl;

			for (int i = 0; i < one_pathes.size(); i++)
			{
				offsets.push_back(one_pathes[i]);
			}

			offsetses.push_back(one_pathes);

			for (int i = 0; i < one_pathes.size(); i++)
			{
				std::vector<Vector2d>().swap(one_pathes[i]);
			}
			std::vector<std::vector<Vector2d>>().swap(one_pathes);

			lOffset = lOffset + toolpath_size;
			Circuit::GenerateOffsetHole(contours, lOffset, one_pathes);
		}
	}

	void ToolpathGenerator::ComputeOffsets(std::vector<Vector2d> &contour)
	{
		double lOffset = toolpath_size / 2.0;
		std::vector<std::vector<Vector2d>> one_pathes;
		Circuit::GenerateOffset(contour, lOffset, one_pathes);

		while (one_pathes.size() > 0)
		{
			std::cout << "Offsets index: " << offsets.size() << std::endl;

			for (int i = 0; i < one_pathes.size(); i++)
			{
				offsets.push_back(one_pathes[i]);
			}

			for (int i = 0; i < one_pathes.size(); i++)
			{
				std::vector<Vector2d>().swap(one_pathes[i]);
			}
			std::vector<std::vector<Vector2d>>().swap(one_pathes);

			lOffset = lOffset + toolpath_size;
			Circuit::GenerateOffset(contour, lOffset, one_pathes);
		}
	}

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

	void ToolpathGenerator::DetectEntryExitPoints(std::vector<Vector2d> &outside_offset, std::vector<Vector2d> &inside_offset, Vector2d &outside_entry_point, Vector2d &outside_exit_point, Vector2d &inside_entry_point, Vector2d &inside_exit_point)
	{
		double m_d = MAXDOUBLE;
		int m_d_index = -1;

		for (int i = 0; i < inside_offset.size(); i++)
		{
			Vector2d v0 = inside_offset[(i - 1 + inside_offset.size()) % inside_offset.size()];
			Vector2d v1 = inside_offset[i];
			Vector2d v2 = inside_offset[(i + 1 + inside_offset.size()) % inside_offset.size()];

			double angle = Circuit::Angle(Vector2d(v0[0] - v1[0], v0[1] - v1[1]), Vector2d(v2[0] - v1[0], v2[1] - v1[1]));

			if (abs(angle - PI / 2.0) < m_d)
			{
				m_d = abs(angle - PI / 2.0);
				m_d_index = i;
			}
		}

		inside_entry_point = inside_offset[m_d_index];
		double d = Circuit::FindNearestPointPar(inside_entry_point, inside_offset);
		double entry_point_d = Circuit::ComputeNextTurningPoint_complex(d, toolpath_size, inside_offset);
		inside_exit_point = Circuit::GetOnePointFromOffset(entry_point_d, inside_offset);

		//Circuit::ComputeNextEntryExitPoint(toolpath_size, outside_offset, inside_offset, outside_entry_point, outside_exit_point, inside_entry_point, inside_exit_point);

		double outside_entry_d = Circuit::FindNearestPointPar(inside_entry_point, outside_offset);
		double outside_exit_d = Circuit::FindNearestPointPar(inside_exit_point, outside_offset);

		outside_entry_point = Circuit::GetOnePointFromOffset(outside_entry_d, outside_offset);
		outside_exit_point = Circuit::GetOnePointFromOffset(outside_exit_d, outside_offset);
	}

	void  ToolpathGenerator::FillingAlgorithmBasedOnOffsets()
	{
		for (int i = 0; i < smooth_number; i++)
			PolygonSmoothing();

		ComputeOffsets_temp();

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
				std::vector<std::vector<Vector2d>> offsets0;

				for (int l = 0; l < decompose_offset[i].size(); l++)
				{
					offsets0.push_back(offsets[decompose_offset[i][l]]);
				}

				Vector2d inside_entry_point;
				Vector2d inside_exit_point;
				Vector2d outside_entry_point;
				Vector2d outside_exit_point;

				DetectEntryExitPoints(offsets[decompose_offset[i][1]], offsets[decompose_offset[i][0]], outside_entry_point, outside_exit_point, inside_entry_point, inside_exit_point);

				FermatsSpiralTrick(offsets0, inside_entry_point, inside_exit_point);

				std::vector<std::vector<Vector2d>>().swap(offsets0);

				entry_spirals.push_back(entry_spiral);
				exit_spirals.push_back(exit_spiral);

				int next_index = Tree::NextNode(mst, decompose_offset[i]);

				double entry_d = Circuit::FindNearestPointPar(entry_spiral[entry_spiral.size()-1],offsets[next_index]);
				outside_entry_point = Circuit::GetOnePointFromOffset(entry_d, offsets[next_index]);

				double exit_d = Circuit::FindNearestPointPar(exit_spiral[exit_spiral.size()-1],offsets[next_index]);
				outside_exit_point = Circuit::GetOnePointFromOffset(exit_d,offsets[next_index]);

				turning_points_entry.push_back(outside_entry_point);
				turning_points_exit.push_back(outside_exit_point);

				//turning_points_entry.push_back(entry_spiral[entry_spiral.size() - 1]);
				//turning_points_entry.push_back(exit_spiral[exit_spiral.size() - 1]);

				Tree::InputCuttingPoints(connect_nodes, cutting_points, next_index, entry_d, exit_d);
				
				//static void InputCuttingPoints(std::vector<int> &nodes, std::vector<std::vector<double>> &cutting_points, int node_id, double cutting_points_0, double cutting_points_1)

				///input
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

				Vector2d inside_entry_point;
				Vector2d inside_exit_point;
				Vector2d outside_entry_point;
				Vector2d outside_exit_point;
				DetectEntryExitPoints(offsets[next_index], offsets[decompose_offset[i][decompose_offset[i].size() - 1]], outside_entry_point, outside_exit_point, inside_entry_point, inside_exit_point);

				double entry_d = Circuit::FindNearestPointPar(outside_entry_point, offsets[next_index]);
				double exit_d = Circuit::FindNearestPointPar(outside_exit_point, offsets[next_index]);

				bool b = Circuit::CheckEnclosed(offsets[next_index], offsets[decompose_offset[i][decompose_offset[i].size() - 1]]);

				if (b)
				{
					Tree::InputCuttingPoints(connect_nodes, cutting_points, next_index, exit_d, entry_d);

					turning_points_entry.push_back(outside_exit_point);
					turning_points_exit.push_back(outside_entry_point);
				}
				else
				{
					Tree::InputCuttingPoints(connect_nodes, cutting_points, next_index, entry_d, exit_d);

					turning_points_entry.push_back(outside_entry_point);
					turning_points_exit.push_back(outside_exit_point);
				}



				//turning_points_entry.push_back(inside_entry_point);
				//turning_points_entry.push_back(inside_exit_point);

				//turning_points_entry.push_back(entry_point);
				//turning_points_exit.push_back(exit_point);

				FermatsSpiralTrick(offsets0, inside_entry_point, inside_exit_point);
				std::vector<std::vector<Vector2d>>().swap(offsets0);

				entry_spirals.push_back(entry_spiral);
				exit_spirals.push_back(exit_spiral);

				std::vector<Vector2d>().swap(entry_spiral);
				std::vector<Vector2d>().swap(exit_spiral);
			}
		}

		//connect sub-region to offset mst tree

		for (int i = 0; i < connect_nodes.size(); i++)
		{
			//connect_edges_cutting_points

			Circuit::CutACircuit(offsets[connect_nodes[i]], cutting_points[i], pathes);
			
			//static void CutACircuit(std::vector<Vector2d> &input_points, std::vector<double> &cutting_points, std::vector<std::vector<Vector2d>> &pathes)
		}

		for (int i = 0; i < connect_edges.size(); i=i+2)
		{
			//connect_edges[i]
			//connect_edges[i+1]
			//offsets



		}



		std::vector<int>().swap(connect_edges);
		std::vector<int>().swap(connect_nodes);
		std::vector<double>().swap(costs);
		std::vector<int>().swap(mst);
		std::vector<int>().swap(nodes);
		std::vector<int>().swap(edges);
	}

	void ToolpathGenerator::FillingAlgorithm()
	{
		for (int i = 0; i < smooth_number; i++)
			PolygonSmoothing();

		//create the region
		std::vector<Vector2d> contour;
		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}
		//ComputeOffsets(contour);

		//return;

		region = Region(contour);
		std::vector<Vector2d>().swap(contour);

		region.DetectInnerConcavePoints();
		region.sdg.GenerateRegionMedialAxisPoints();
		region.sdg.GenerateVoronoiEdgePoints();
		region.sdg.ComputePointsDegree();
		region.sdg.DetectMaximalAndMinimalPoints();
		region.sdg.DecomposeMedialAxis1();
		//region.sdg.DetectCriticalPoints1(region.inner_concave_points, toolpath_size, toolpath_size / 8.0);
		region.sdg.DetectCriticalPoints(region.inner_concave_points, toolpath_size);
		region.ComputeCuttingPoints();
		region.DecomposeSubregions();

		//Generate Fermat spiral for each sub-region
		//GenerateOffsetsForAllPolygons();

		region.GenerateConnectedGraph(toolpath_size);
		region.ComputeTSPlikePath();
		region.ComputeEntryAndExitPoints();

		if (region.polygons.size() == 1)
		{
			Vector2d entry_p_0 = Circuit::GetOnePointFromOffset(entry_d_0, region.polygons[0]);
			Vector2d exit_p_0 = Circuit::GetOnePointFromOffset(exit_d_0, region.polygons[0]);

			//FermatsSpiralTrick(region.polygons[0], entry_p_0, exit_p_0);
			//FermatsSpiralSmooth1(region.polygons[0], entry_p_0, exit_p_0);
			FermatsSpiralSmooth(region.polygons[0], entry_p_0, exit_p_0);
			//FermatSpiral(region.polygons[0], entry_p_0, exit_p_0);
			region.entry_spirals.push_back(entry_spiral);
			region.exit_spirals.push_back(exit_spiral);
		}
		else
		{
			for (int i = 0; i < region.polygons.size(); i++)
			{
				if (debug_int_0 >= 0)
				{
					if (i != debug_int_0)
						continue;
				}
		
				std::vector<std::vector<Vector2d>>().swap(offsets);
				std::vector<Vector2d>().swap(entry_spiral);
				std::vector<Vector2d>().swap(exit_spiral);
				//FermatSpiral(region.polygons[i], region.polygons_entry_exit[i][0], region.polygons_entry_exit[i][1]);
				FermatsSpiralTrick(region.polygons[i], region.polygons_entry_exit[i][0], region.polygons_entry_exit[i][1]);
				//FermatsSpiralSmooth(region.polygons[i], region.polygons_entry_exit[i][0], region.polygons_entry_exit[i][1]);
				region.entry_spirals.push_back(entry_spiral);
				region.exit_spirals.push_back(exit_spiral);
			}
		}
	}
	      
	void ToolpathGenerator::GenerateOffsetsForAllPolygons()
	{
		for (int i = 0; i < region.polygons.size(); i++)
		{
			double lOffset = toolpath_size/2.0;

			do
			{
				std::vector<std::vector<Vector2d>> offset;
				
				Circuit::GenerateOffset(region.polygons[i], lOffset, offset);
				lOffset = lOffset + toolpath_size;
				if (offset.size() > 0)
				{
					for (int j = 0; j < offset.size(); j++)
					{
						offsets.push_back(offset[j]);
						std::vector<Vector2d>().swap(offset[j]);
					}

					std::vector<std::vector<Vector2d>>().swap(offset);
				}
				else
				{
					std::vector<std::vector<Vector2d>>().swap(offset);
					break;
				}

			} while (true);

		}
	}
}