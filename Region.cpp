#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"
#include "Circuit.h"

namespace hpcg {

	Region::Region(std::vector<Vector2d> &vecs)
	{
		for (int i = 0; i < vecs.size(); i++)
		{
			contour.push_back(vecs[i]);
		}
		sdg = SegmentDelaunayGraphs(vecs);
	}

	Region::~Region()
	{
		std::vector<Vector2d>().swap(contour);
		std::vector<Vector2d>().swap(inner_concave_points);
		std::vector<Vector2d>().swap(cutting_points);

		for (int i = 0; i < polygons.size(); i++)
		{
			std::vector<Vector2d>().swap(polygons[i]);
		}
		std::vector<std::vector<Vector2d>>().swap(polygons);

		std::vector<int>().swap(connected_graph);
		std::vector<Vector2d>().swap(entry_exit_points);
	}

	//detect inner concave points
	void Region::DetectInnerConcavePoints()
	{
		for (int i = 0; i < contour.size(); i++)
		{
			Polygon_2 polys;

			polys.push_back(Point_2(contour[i][0], contour[i][1]));
			polys.push_back(Point_2(contour[(i + 1) % contour.size()][0], contour[(i + 1) % contour.size()][1]));
			polys.push_back(Point_2(contour[(i + 2) % contour.size()][0], contour[(i + 2) % contour.size()][1]));

			if (polys.is_simple() && polys.is_clockwise_oriented())
			{
				inner_concave_points.push_back(Vector2d(contour[(i + 1) % contour.size()][0], contour[(i + 1) % contour.size()][1]));
			}
			polys.clear();
		}
	}

	//compute cuting points
	void Region::ComputeCuttingPoints()
	{
		for (int i = 0; i < sdg.critical_points.size(); i++)
		{
			Circuit::CuttingPoints(sdg.critical_points[i], contour, cutting_points);
		}
	}

	//Decompose the whole region into sub-regions with the cutting points
	void Region::DecomposeSubregions()
	{
		Circuit::DecomposeSubregions(cutting_points, contour, polygons);
	}

	//Generate the connected graph of sub-regions 
	void Region::GenerateConnectedGraph(double toolpath_size)
	{
		for (int i = 0; i < cutting_points.size(); i = i + 2)
		{
			int connected_two = 0;
			for (int j = 0; j < polygons.size(); j++)
			{
				if (Circuit::Distance(cutting_points[i], polygons[j]) < 0.000001)
				{
					connected_graph.push_back(j);
					connected_two++;
				}
			}
			entry_exit_points.push_back(Circuit::ComputeEntryExitPoint(cutting_points[i], toolpath_size/2.0, polygons[connected_graph[connected_graph.size() - 2]]));
			entry_exit_points.push_back(Circuit::ComputeEntryExitPoint(cutting_points[i + 1], toolpath_size / 2.0, polygons[connected_graph[connected_graph.size() - 2]]));
			entry_exit_points.push_back(Circuit::ComputeEntryExitPoint(cutting_points[i], toolpath_size / 2.0, polygons[connected_graph[connected_graph.size() - 1]]));
			entry_exit_points.push_back(Circuit::ComputeEntryExitPoint(cutting_points[i + 1], toolpath_size / 2.0, polygons[connected_graph[connected_graph.size() - 1]]));

			assert(connected_two == 2);
		}
	}

	//Compute a TSP-like path of the connected graph
	void Region::ComputeTSPlikePath()
	{
		std::vector<bool> used;
		for (int i = 0; i < polygons.size(); i++)
		{
			used.push_back(false);
		}

		do
		{
			int start_index = -1;
			int connected_graph_index = -1;
			for (int i = 0; i < connected_graph.size(); i++)
			{
				if (!used[connected_graph[i]])
				{
					start_index = connected_graph[i];
					connected_graph_index = i;
					break;
				}
			}

			if (start_index >= 0)
			{
				std::vector<int> vi;
				std::vector<int> vi_index;

				do
				{
					vi.push_back(start_index);
					vi_index.push_back(connected_graph_index);

					used[start_index] = true;

					int next_index = -1;
					for (int i = 0; i < connected_graph.size(); i++)
					{
						if (connected_graph[i] == start_index)
						{
							if (i % 2 == 0)
							{
								if (!used[connected_graph[i + 1]])
								{
									next_index = connected_graph[i + 1];
									connected_graph_index = i + 1;
									break;
								}
							}
							else
							{
								if (!used[connected_graph[i - 1]])
								{
									next_index = connected_graph[i - 1];
									connected_graph_index = i - 1;
									break;
								}
							}
						}
					}

					if (next_index >= 0)
					{
						start_index = next_index;
					}
					else
					{
						break;
					}

				} while (true);
				connected_graph_de.push_back(vi);
				connected_graph_de_index.push_back(vi_index);
				std::vector<int>().swap(vi);
				std::vector<int>().swap(vi_index);
			}
			else
			{
				break;
			}

		} while (true);
		std::vector<bool>().swap(used);
	}

	//Compute entry and exit points of the sub-regions
	void Region::ComputeEntryAndExitPoints()
	{
		for (int i = 0; i < polygons.size(); i++)
		{
			polygons_entry_exit.push_back(std::vector<Vector2d>());
		}

		for (int i = 0; i < connected_graph_de.size(); i++)
		{
			if (connected_graph_de[i].size() == 1)
			{
				polygons_entry_exit[connected_graph_de[i][0]].push_back(entry_exit_points[connected_graph_de_index[i][0] * 2]);
				polygons_entry_exit[connected_graph_de[i][0]].push_back(entry_exit_points[connected_graph_de_index[i][0] * 2 + 1]);
			}

			if (connected_graph_de[i].size() > 1)
			{
				for (int j = 0; j < connected_graph_de[i].size() - 1; j++)
				{
					int int0 = connected_graph_de[i][j];
					int int1 = connected_graph_de[i][j + 1];
					for (int k = 0; k < connected_graph.size(); k = k + 2)
					{
						int int2 = connected_graph[k];
						int int3 = connected_graph[k + 1];
						if ((int0 == int2&&int1 == int3) || (int0 == int3&&int1 == int2))
						{
							if (j == 0)
							{
								if (int0 == int2)
								{
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 1]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 2]);

									connected_regions.push_back(entry_exit_points[k * 2]);
									connected_regions.push_back(entry_exit_points[k * 2 + 2]);
								}
								else
								{
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 2]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 3]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 1]);

									connected_regions.push_back(entry_exit_points[k * 2 + 1]);
									connected_regions.push_back(entry_exit_points[k * 2 + 3]);

								}
							}
							if (j > 0 && j < connected_graph_de[i].size() - 2)
							{
								if (int0 == int2)
								{
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 2]);

									connected_regions.push_back(entry_exit_points[k * 2]);
									connected_regions.push_back(entry_exit_points[k * 2 + 2]);
								}
								else
								{
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2]);

									connected_regions.push_back(entry_exit_points[k * 2 + 2]);
									connected_regions.push_back(entry_exit_points[k * 2]);
								}
							}

							if (j == connected_graph_de[i].size() - 2)
							{
								if (int0 == int2)
								{
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 3]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 1]);

									connected_regions.push_back(entry_exit_points[k * 2 + 1]);
									connected_regions.push_back(entry_exit_points[k * 2 + 2]);
								}
								else
								{
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 1]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 2]);

									connected_regions.push_back(entry_exit_points[k * 2]);
									connected_regions.push_back(entry_exit_points[k * 2 + 2]);
								}
							}
						}
					}
				}
			}
		}
	}
}