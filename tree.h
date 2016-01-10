#pragma once

#include <iostream>
#include <fstream>
#include <cassert>

#include <vector>
#include <string>
#include <algorithm>

#include <ToolpathGenerator.h>
#include <Strip.h>

namespace hpcg {

	class Tree
	{

	public:
		static void MinimalSpanningTree(std::vector<int> &nodes, std::vector<int> &edges, std::vector<double> &costs, std::vector<int> &mst)
		{
			std::vector<std::vector<int>> containers;
			for (int i = 0; i < nodes.size(); i++)
			{
				std::vector<int> container;
				container.push_back(nodes[i]);
				containers.push_back(container);
				std::vector<int>().swap(container);
			}

			std::vector<bool> edges_used;
			for (int i = 0; i < costs.size(); i++)
			{
				edges_used.push_back(false);
			}

			do
			{
				//find a minimal cost edge
				int minimal_cost_edge_index = -1;
				double minimal_cost = MAXDOUBLE;
			#pragma region find_a_minimal_cost_edge

				for (int j = 0; j < costs.size(); j++)
				{
					if (!edges_used[j])
					{
						if (costs[j] < minimal_cost)
						{
							minimal_cost = costs[j];
							minimal_cost_edge_index = j;
						}
					}

				}
			#pragma endregion

				if (minimal_cost_edge_index < 0)
					break;

				//check valid
				int node_index_0 = edges[2 * minimal_cost_edge_index];
				int node_index_1 = edges[2 * minimal_cost_edge_index + 1];


				if (node_index_0 == 16 && node_index_1 == 18)
				{
					int dsad = 0;
				}

				int container_0 = -1;
				int container_0_0 = -1;
				int container_1 = -1;
				int container_1_0 = -1;

				for (int j = 0; j < containers.size() && (container_0 < 0 || container_1 < 0); j++)
				{
					for (int k = 0; k < containers[j].size() && (container_0 < 0 || container_1 < 0); k++)
					{
						if (node_index_0 == containers[j][k])
						{
							container_0 = j;
							container_0_0 = k;
						}
						if (node_index_1 == containers[j][k])
						{
							container_1 = j;
							container_1_0 = k;
						}
					}
				}

				if (!(container_0 >= 0 && container_1 >= 0))
				{
					break;
				}

				if (container_0 == container_1)
				{
					edges_used[minimal_cost_edge_index] = true;
				}
				else
				{
					mst.push_back(node_index_0);
					mst.push_back(node_index_1);
					edges_used[minimal_cost_edge_index] = true;

					for (int i = 0; i < containers[container_1].size(); i++)
					{
						containers[container_0].push_back(containers[container_1][i]);
					}

					containers.erase(containers.begin() + container_1);
				}

			} while (containers.size() != 1);

			std::vector<bool>().swap(edges_used);
			std::vector<std::vector<int>>().swap(containers);
		}

		static void InputCuttingPoints(std::vector<int> &nodes, std::vector<std::vector<double>> &cutting_points, std::vector<std::vector<int>> &cutting_points_index,
			int node_id, double cutting_points_0, double cutting_points_1, int cutting_points_index_0, int cutting_points_index_1)
		{
			int node_index = -1;

			for (int i = 0; i < nodes.size(); i++)
			{
				if (nodes[i] == node_id)
				{
					node_index = i;
				}
			}

			if (node_index < 0)
				return;

			cutting_points[node_index].push_back(cutting_points_0);
			cutting_points[node_index].push_back(cutting_points_1);

			cutting_points_index[node_index].push_back(cutting_points_index_0);
			cutting_points_index[node_index].push_back(cutting_points_index_1);
		}


		static void DecompositionATree(std::vector<int> &nodes, std::vector<int> &edges, std::vector<int> &connect_nodes, std::vector<int> &connect_edges, std::vector<std::vector<int>> &pathes)
		{
			std::vector<bool> node_used;
			std::vector<int> node_degree;

			for (int i = 0; i < nodes.size(); i++)
			{
				node_used.push_back(false);
				node_degree.push_back(0);
			}

			do
			{
				//compute node degree
				for (int i = 0; i < edges.size(); i = i + 2)
				{
					if (!node_used[edges[i]] && !node_used[edges[i + 1]])
					{
						node_degree[edges[i]]++;
						node_degree[edges[i + 1]]++;
					}
				}

				//decompose the tree
				do
				{
					int start_index = -1;
					for (int i = 0; i < node_degree.size(); i++)
					{
						if (node_degree[i] == 1 && !node_used[i])
						{
							start_index = i;
							break;
						}
					}

					if (start_index < 0)
					{
						break;
					}

					std::vector<int> one_path;

					//compute one path

					do
					{
						one_path.push_back(start_index);
						node_used[start_index] = true;

						int next_index = -1;

						for (int i = 0; i < edges.size(); i = i + 2)
						{
							if (edges[i] == start_index&&!node_used[edges[i + 1]])
							{
								next_index = edges[i + 1];
								break;
							}
							if (edges[i + 1] == start_index&&!node_used[edges[i]])
							{
								next_index = edges[i];
								break;
							}
						}

						if (next_index < 0)
						{
							break;
						}

						if (node_degree[next_index] != 2)
						{
							/*
							if (node_degree[next_index] == 1)
							{
								one_path.push_back(next_index);
								node_used[next_index] = true;
							}
							*/

							//one_path.push_back(next_index);

							//one_path.push_back(next_index);
							//node_used[next_index] = true;
							break;
						}
						start_index = next_index;

					} while (true);

					pathes.push_back(one_path);

				} while (true);

				//handle isolated points

				for (int i = 0; i < node_used.size(); i++)
				{
					if (!node_used[i] && node_degree[i] == 0)
					{
						std::vector<int> one_path;
						one_path.push_back(i);
						pathes.push_back(one_path);
						node_used[i] = true;
						std::vector<int>().swap(one_path);
					}
				}
				bool goon = false;
				int iii = 0;
				for (int i = 0; i < node_used.size(); i++)
				{
					if (!node_used[i])
					{
						goon = true;
						iii++;
					}
				}

				if (!goon)
					break;

				for (int i = 0; i < edges.size(); i = i + 2)
				{
					if (!node_used[edges[i]] && !node_used[edges[i + 1]])
					{
						connect_edges.push_back(edges[i]);
						connect_edges.push_back(edges[i + 1]);
					}
				}

				for (int i = 0; i < node_used.size(); i++)
				{
					if (!node_used[i])
					{
						connect_nodes.push_back(i);
					}
				}
				
				break;
			} while (true);


			std::vector<bool>().swap(node_used);
			std::vector<int>().swap(node_degree);
		}

		static void NextNode(std::vector<int> &edges, int node_index, std::vector<int> &related_int)
		{
			int return_int = -1;

			for (int i = 0; i < edges.size(); i = i + 2)
			{
				if (edges[i] == node_index)
				{
					related_int.push_back(edges[i+1]);
				}

				if (edges[i + 1] == node_index)
				{
					related_int.push_back(edges[i]);
				}
			}
		}

		//node_index_1->node_index_0->????
		static int NextNode(std::vector<int> &edges, int node_index_0, int node_index_1)
		{
			int int_index = -1;
			std::vector<int> related_int;
			NextNode(edges, node_index_0, related_int);

			if (related_int.size() == 2)
			{
				if (related_int[0] == node_index_1)
				{
					int_index = related_int[1];
				}
				else
				{
					int_index = related_int[0];
				}
				std::vector<int>().swap(related_int);

				return int_index;
			}
			else
			{
				return related_int[0];
			}
		}

		static int NextNode(std::vector<int> &edges, std::vector<int> &one_path)
		{
			if (one_path.size() == 1)
			{
				int int_index = -1;
				std::vector<int> related_int;
				NextNode(edges, one_path[one_path.size() - 1], related_int);

				assert(related_int.size() == 1);

				return related_int[0];
			}
			else
			{

				int int_index = -1;
				std::vector<int> related_int;
				NextNode(edges, one_path[one_path.size() - 1], related_int);

				assert(related_int.size() == 2);

				if (related_int[0] == one_path[one_path.size() - 2])
				{
					int_index = related_int[1];
				}
				else
				{
					int_index = related_int[0];
				}
				std::vector<int>().swap(related_int);

				return int_index;

			}

		}

	};
}





