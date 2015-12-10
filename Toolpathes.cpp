#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"


namespace hpcg {


	void ToolpathGenerator::RelatedPointsOnContours(Vector2d &v, double distace, std::vector<Vector2d> &vecs)
	{
		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		std::vector<Vector2d> vt;

		for (int i = 0; i < contour.size(); i++)
		{
			Point_2 p0(contour[i].x, contour[i].y);
			Point_2 p1(contour[(i + 1) % contour.size()].x, contour[(i + 1) % contour.size()].y);


			double l = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), Segment_2(p0, p1)));

			if (std::abs(l - distace) < 0.00001)
			{
				double l0 = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), p0));
				double l1 = sqrt((double)CGAL::squared_distance(Point_2(v[0], v[1]), p1));

				if (l < l0 && l < l1)
				{
					Vector2d vec(p1[0] - p0[0], p1[1] - p0[1]);
					Vector2d r_vec(vec[1], -vec[0]);
					if (vec[0] < 0.00001)
					{
						r_vec[0] = -vec[1];
						r_vec[1] = vec[0];
					}

					if (vec[1] < 0.00001)
					{
						r_vec[0] = vec[1];
						r_vec[1] = -vec[0];
					}

					CGAL::Object result = CGAL::intersection(Line_2(p0, p1), Line_2(Point_2(v[0], v[1]), Point_2(v[0] + r_vec[0], v[1] + r_vec[1])));

					if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
					{
						vt.push_back(Vector2d(ipoint->x(), ipoint->y()));
					}
					else
					{
						assert(false);
					}
				}
				else
				{
					Vector2d vp;

					if (l0 < l1)
					{
						//vt.push_back(Vector2d(p0[0], p0[1]));
						vp=Vector2d(p0[0], p0[1]);
					}
					else
					{
						//vt.push_back(Vector2d(p1[0], p1[1]));
						vp = Vector2d(p1[0], p1[1]);
					}

					bool insert = true;
					for (int j = 0; j < vt.size(); j++)
					{
						if (std::abs(vt[j][0] - vp[0]) < 0.000001&&std::abs(vt[j][1] - vp[1]) < 0.000001)
						{
							insert = false;
							break;
						}
					}
					if (insert)
					{
						vt.push_back(vp);
					}
				}
			}
		}

		for (int i = 0; i < vt.size(); i++)
		{
			vecs.push_back(vt[i]);
		}

		std::vector<Vector2d>().swap(vt);
		std::vector<Vector2d>().swap(contour);
	}

	void ToolpathGenerator::GenerateVoronoiEdgePoints()
	{
		SDG2::Finite_edges_iterator eit = sdg.finite_edges_begin();
		for (int k = 1; eit != sdg.finite_edges_end(); ++eit, ++k) {
			SDG2::Edge e = *eit;

			SDG2::Vertex_handle v[] = { e.first->vertex(sdg.ccw(e.second)),
				e.first->vertex(sdg.cw(e.second)),
				e.first->vertex(e.second),
				sdg.tds().mirror_vertex(e.first, e.second) };

			if (!sdg.is_infinite(v[0]) && !sdg.is_infinite(v[1]) && !sdg.is_infinite(v[2]) && !sdg.is_infinite(v[3]))
			{
				Gt::Segment_2           s;
				CGAL::Parabola_segment_2<Gt>       ps;
				CGAL::Parabola_2<Gt> pss;

				CGAL::Object o = sdg.primal(e);

				if (CGAL::assign(s, o))
				{
					double center_x = (s.source().x() + s.target().x()) / 2.0;
					double center_y = (s.source().y() + s.target().y()) / 2.0;

					if (CheckInsideContours(Vector2d(center_x, center_y)))
					{
						voronoi_edge_points.push_back(Vector2d(s.source().x(), s.source().y()));
						voronoi_edge_points.push_back(Vector2d(s.target().x(), s.target().y()));
					}
				}

				if (CGAL::assign(ps, o))
				{
					std::vector< Gt::Point_2 > p;
					ps.generate_points(p);

					double center_x = 0.0;
					double center_y = 0.0;

					for (int i = 0; i < p.size(); i++)
					{
						center_x += p[i].x();
						center_y += p[i].y();
					}

					center_x = center_x / p.size();
					center_y = center_y / p.size();

					if (CheckInsideContours(Vector2d(center_x, center_y)))
					{
						for (int i = 0; i < p.size() - 1; i++)
						{
							voronoi_edge_points.push_back(Vector2d(p[i].x(), p[i].y()));
							voronoi_edge_points.push_back(Vector2d(p[(i + 1) % p.size()].x(), p[(i + 1) % p.size()].y()));
						}
					}

					std::vector< Gt::Point_2 >().swap(p);
				}
			}
		}
	}

	void ToolpathGenerator::GenerateRegionMedialAxisPoints()
	{

		SDG2::Finite_edges_iterator eit = sdg.finite_edges_begin();
		for (int k = 1; eit != sdg.finite_edges_end(); ++eit, ++k) {
			SDG2::Edge e = *eit;

			SDG2::Vertex_handle v[] = { e.first->vertex(sdg.ccw(e.second)),
				e.first->vertex(sdg.cw(e.second)),
				e.first->vertex(e.second),
				sdg.tds().mirror_vertex(e.first, e.second) };

			if (!sdg.is_infinite(v[0]) && !sdg.is_infinite(v[1]) && !sdg.is_infinite(v[2]) && !sdg.is_infinite(v[3]))
			{

				Gt::Segment_2           s;
				CGAL::Parabola_segment_2<Gt>       ps;

				CGAL::Object o = sdg.primal(e);

				if (CGAL::assign(s, o))
				{
					if (v[0]->site().is_point() && v[1]->site().is_point())
					{
						Point_2 p0(v[0]->site().point().x(), v[0]->site().point().y());
						Point_2 p1(v[1]->site().point().x(), v[1]->site().point().y());

						CGAL::Object result = CGAL::intersection(Segment_2(Point_2(s.source().x(), s.source().y()), Point_2(s.target().x(), s.target().y())), Segment_2(p0, p1));

						if (result.is_empty())
						{
							double center_x = (s.source().x() + s.target().x()) / 2.0;
							double center_y = (s.source().y() + s.target().y()) / 2.0;

							if (CheckInsideContours(Vector2d(center_x, center_y)) && MinimalDistanceContours(Vector2d(s.source().x(), s.source().y()))>0.00001&&MinimalDistanceContours(Vector2d(s.target().x(), s.target().y()))>0.00001)
							{
								medial_axis_points.push_back(Vector2d(s.source().x(), s.source().y()));
								medial_axis_points.push_back(Vector2d(s.target().x(), s.target().y()));
							}
						}
						else
						{
							if (const Point_2 *ipoint = CGAL::object_cast<Point_2>(&result))
							{
								double center_x = (s.source().x() + s.target().x()) / 2.0;
								double center_y = (s.source().y() + s.target().y()) / 2.0;

								if (CheckInsideContours(Vector2d(center_x, center_y)) && MinimalDistanceContours(Vector2d(s.source().x(), s.source().y()))>0.00001&&MinimalDistanceContours(Vector2d(s.target().x(), s.target().y()))>0.00001)
								{
									medial_axis_points.push_back(Vector2d(s.source().x(), s.source().y()));
									medial_axis_points.push_back(Vector2d(ipoint->x(), ipoint->y()));
									medial_axis_points.push_back(Vector2d(ipoint->x(), ipoint->y()));
									medial_axis_points.push_back(Vector2d(s.target().x(), s.target().y()));
								}
							}
						}
					}
					else
					{
						double center_x = (s.source().x() + s.target().x()) / 2.0;
						double center_y = (s.source().y() + s.target().y()) / 2.0;

						if (CheckInsideContours(Vector2d(center_x, center_y)) && MinimalDistanceContours(Vector2d(s.source().x(), s.source().y()))>0.00001&&MinimalDistanceContours(Vector2d(s.target().x(), s.target().y()))>0.00001)
						{
							medial_axis_points.push_back(Vector2d(s.source().x(), s.source().y()));
							medial_axis_points.push_back(Vector2d(s.target().x(), s.target().y()));
						}
					}
				}

				if (CGAL::assign(ps, o))
				{
					std::vector< Gt::Point_2 > p;
					ps.generate_points(p);

					double center_x = 0.0;
					double center_y = 0.0;

					for (int i = 0; i < p.size(); i++)
					{
						center_x += p[i].x();
						center_y += p[i].y();
					}

					center_x = center_x / p.size();
					center_y = center_y / p.size();

					if (CheckInsideContours(Vector2d(center_x, center_y)) && MinimalDistanceContours(Vector2d(p[0].x(), p[0].y()))>0.00001&&MinimalDistanceContours(Vector2d(p[p.size() - 1].x(), p[p.size() - 1].y()))>0.00001)
					{
						for (int i = 0; i < p.size() - 1; i++)
						{
							medial_axis_points.push_back(Vector2d(p[i].x(), p[i].y()));
							medial_axis_points.push_back(Vector2d(p[(i + 1) % p.size()].x(), p[(i + 1) % p.size()].y()));
						}
					}

					std::vector< Gt::Point_2 >().swap(p);
				}
			}
		}
	}
	void ToolpathGenerator::CompuateCutPoints()
	{
		std::vector<Gt::Point_2> points;
		std::vector<std::pair<std::size_t, std::size_t> > indices;
		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			points.push_back(Gt::Point_2(ver_iter->x(), ver_iter->y()));
		}

		for (int i = 0; i < points.size(); i++)
		{
			indices.push_back(std::make_pair(i, (i + 1) % points.size()));
		}

		sdg.insert_segments(points.begin(), points.end(), indices.begin(), indices.end());

		assert(sdg.is_valid(true, 1));

		std::vector<Gt::Point_2>().swap(points);
		std::vector<std::pair<std::size_t, std::size_t> >().swap(indices);

		GenerateRegionMedialAxisPoints();
		GenerateVoronoiEdgePoints();

		std::vector<int> medial_axis_points_degree;
		std::vector<double> medial_axis_points_radius;
		std::vector<std::vector<int>> medial_axis_points_degree_related;

		std::vector<bool> medial_axis_points_used;
		
		for (int i = 0; i < medial_axis_points.size(); i++)
		{
			medial_axis_points_degree.push_back(0);
			std::vector<int> vt;
			medial_axis_points_degree_related.push_back(vt);
			std::vector<int>().swap(vt);

			double d = MinimalDistanceContours(medial_axis_points[i]);
			medial_axis_points_radius.push_back(d);

			medial_axis_points_used.push_back(false);
		}

		for (int i = 0; i < medial_axis_points.size(); i = i + 2)
		{
			Point_2 p0(medial_axis_points[i][0], medial_axis_points[i][1]);
			Point_2 p1(medial_axis_points[i + 1][0], medial_axis_points[i + 1][1]);

			for (int j = 0; j < medial_axis_points.size(); j = j + 2)
			{
				Point_2 p2(medial_axis_points[j][0], medial_axis_points[j][1]);
				Point_2 p3(medial_axis_points[j + 1][0], medial_axis_points[j + 1][1]);

				if (sqrt((double)CGAL::squared_distance(p0, p2)) < 0.0001)
				{
					medial_axis_points_degree[i]++;
					medial_axis_points_degree_related[i].push_back(j);
				}

				if (sqrt((double)CGAL::squared_distance(p0, p3)) < 0.0001)
				{
					medial_axis_points_degree[i]++;
					medial_axis_points_degree_related[i].push_back(j + 1);
				}

				if (sqrt((double)CGAL::squared_distance(p1, p2)) < 0.0001)
				{
					medial_axis_points_degree[i + 1]++;
					medial_axis_points_degree_related[i + 1].push_back(j);
				}

				if (sqrt((double)CGAL::squared_distance(p1, p3)) < 0.0001)
				{
					medial_axis_points_degree[i + 1]++;
					medial_axis_points_degree_related[i + 1].push_back(j + 1);
				}
			}
		}


		for (int i = 0; i < medial_axis_points.size(); i++)
		{
			if (medial_axis_points_degree[i] == 2 && !medial_axis_points_used[i])
			{
				int last_index = -1;
				int next_index = -1;

				if (i % 2 == 0)
				{
					next_index = i + 1;
				}
				else
				{
					next_index = i - 1;
				}

				if (medial_axis_points_degree_related[i][0] == i)
				{
					last_index = medial_axis_points_degree_related[i][1];
				}
				else
				{
					last_index = medial_axis_points_degree_related[i][0];
				}

				if (last_index % 2 == 0)
				{
					last_index = last_index + 1;
				}
				else
				{
					last_index = last_index - 1;
				}

				if (medial_axis_points_radius[i] < medial_axis_points_radius[last_index] && medial_axis_points_radius[i] < medial_axis_points_radius[next_index])
				{
					for (int j = 0; j < medial_axis_points_degree_related[i].size(); j++)
					{
						medial_axis_points_used[medial_axis_points_degree_related[i][j]] = true;
					}
					cccc.push_back(medial_axis_points[i]);
					RelatedPointsOnContours(medial_axis_points[i], medial_axis_points_radius[i], dddd);
				}
			}

			if (medial_axis_points_degree[i] == 3)
			{
				int another_index_0 = medial_axis_points_degree_related[i][0];
				int another_index_1 = medial_axis_points_degree_related[i][1];
				int another_index_2 = medial_axis_points_degree_related[i][2];

				if (another_index_0 % 2 == 0)
				{
					another_index_0 = another_index_0 + 1;
				}
				else
				{
					another_index_0 = another_index_0 - 1;
				}

				if (another_index_1 % 2 == 0)
				{
					another_index_1 = another_index_1 + 1;
				}
				else
				{
					another_index_1 = another_index_1 - 1;
				}

				if (another_index_2 % 2 == 0)
				{
					another_index_2 = another_index_2 + 1;
				}
				else
				{
					another_index_2 = another_index_2 - 1;
				}

				int bool_nb = 0;

				if (medial_axis_points_radius[i] < medial_axis_points_radius[another_index_0])
				{
					bool_nb++;
				}
				if (medial_axis_points_radius[i] < medial_axis_points_radius[another_index_1])
				{
					bool_nb++;
				}
				if (medial_axis_points_radius[i] < medial_axis_points_radius[another_index_2])
				{
					bool_nb++;
				}

				if (bool_nb >= 1)
				{
					//cccc.push_back(medial_axis_points[i]);
					//RelatedPointsOnContours(medial_axis_points[i], medial_axis_points_radius[i], dddd);
				}
			}

		}

		for (int i = 0; i < medial_axis_points_degree_related.size(); i++)
		{
			std::vector<int>().swap(medial_axis_points_degree_related[i]);
		}
		std::vector<std::vector<int>>().swap(medial_axis_points_degree_related);
		std::vector<int>().swap(medial_axis_points_degree);
		std::vector<bool>().swap(medial_axis_points_used);

		std::vector<Vector2d> contour;

		for (Polygon_2::Vertex_iterator ver_iter = contours.outer_boundary().vertices_begin(); ver_iter != contours.outer_boundary().vertices_end(); ver_iter++)
		{
			contour.push_back(Vector2d(ver_iter->x(), ver_iter->y()));
		}

		DecomposePolygons(dddd, contour, polygons);

		std::vector<Vector2d>().swap(contour);

	}

	void ToolpathGenerator::DecomposePolygons(std::vector<Vector2d> &cut_points, std::vector<Vector2d> &contour, std::vector<std::vector<Vector2d>> &polygons)
	{
		assert(cut_points.size()%2==0);

		//decompise the whole polygon
		polygons.push_back(contour);

		for (int i = 0; i < cut_points.size(); i = i + 2)
		{
			int polygon_index = -1;
			for (int j = 0; j < polygons.size(); j++)
			{
				if (MinimalDistance(polygons[j], cut_points[i]) < 0.000001)
				{
					polygon_index = j;
					break;
				}
			}

			if (polygon_index >= 0)
			{
				double d0 = FindNearestPointPar(cut_points[i], polygons[polygon_index]);
				double d1 = FindNearestPointPar(cut_points[i + 1], polygons[polygon_index]);

				std::vector<Vector2d> half_part_0;
				std::vector<Vector2d> half_part_1;

				SelectOnePartOffset(polygons[polygon_index], d0, d1, half_part_0);
				SelectOnePartOffset(polygons[polygon_index], d1, d0, half_part_1);

				polygons.erase(polygons.begin() + polygon_index);
				polygons.push_back(half_part_0);
				polygons.push_back(half_part_1);

				std::vector<Vector2d>().swap(half_part_0);
				std::vector<Vector2d>().swap(half_part_1);
			}
		}

		for (int i = 0; i < polygons.size();i++)
		{
			polygons_entry_exit.push_back(std::vector<Vector2d>());
		}


		//compute connected graph
		for (int i = 0; i < polygons.size(); i++)
		{
			double center_x = 0.0;
			double center_y = 0.0;

			for (int j = 0; j < polygons[i].size(); j++)
			{
				center_x += polygons[i][j][0];
				center_y += polygons[i][j][1];
			}
			center_x = center_x / polygons[i].size();
			center_y = center_y / polygons[i].size();

			connected_graph_point.push_back(Vector2d(center_x, center_y));
		}

		for (int i = 0; i < cut_points.size(); i = i + 2)
		{
			int connected_two = 0;
			for (int j = 0; j < polygons.size(); j++)
			{
				if (MinimalDistance(polygons[j], cut_points[i]) < 0.000001)
				{
					connected_graph.push_back(j);
					connected_two++;
				}
			}
			
			entry_exit_points.push_back(ComputeEntryExitPoint(polygons[connected_graph[connected_graph.size() - 2]], cut_points[i]));
			entry_exit_points.push_back(ComputeEntryExitPoint(polygons[connected_graph[connected_graph.size() - 2]], cut_points[i + 1]));
			entry_exit_points.push_back(ComputeEntryExitPoint(polygons[connected_graph[connected_graph.size() - 1]], cut_points[i]));
			entry_exit_points.push_back(ComputeEntryExitPoint(polygons[connected_graph[connected_graph.size() - 1]], cut_points[i+1]));

			assert(connected_two==2);
		}

		////////////////////////////////////////////////
		//Decompose the connected graph
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

			if (start_index>=0)
			{
				std::vector<int> vi;
				std::vector<int> vi_index;

				do
				{
					vi.push_back(start_index);
					vi_index.push_back(connected_graph_index);

					used[start_index]=true;

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

		//generate the entry ande exit points

		for (int i = 0; i < connected_graph_de.size(); i++)
		{
			if (connected_graph_de[i].size() == 1)
			{
				//connected_graph_de[i][0]
				//connected_graph_de_index[i][0];

				//cccc.push_back(entry_exit_points[connected_graph_de_index[i][0] * 2]);
				//cccc.push_back(entry_exit_points[connected_graph_de_index[i][0] * 2+1]);

				polygons_entry_exit[connected_graph_de[i][0]].push_back(entry_exit_points[connected_graph_de_index[i][0] * 2]);
				polygons_entry_exit[connected_graph_de[i][0]].push_back(entry_exit_points[connected_graph_de_index[i][0] * 2 + 1]);
			}

			if (connected_graph_de[i].size() > 1)
			{
				for (int j = 0; j < connected_graph_de[i].size()-1; j++)
				{
					int int0 = connected_graph_de[i][j];
					int int1 = connected_graph_de[i][j+1];
					for (int k = 0; k < connected_graph.size(); k=k+2)
					{
						int int2 = connected_graph[k];
						int int3 = connected_graph[k+1];
						if ((int0 == int2&&int1 == int3) || (int0 == int3&&int1 == int2))
						{
							if (j == 0)
							{
								if (int0 == int2)
								{
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 1]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 2]);

									//cccc.push_back(entry_exit_points[k * 2]);
									//cccc.push_back(entry_exit_points[k * 2 + 1]);
									//cccc.push_back(entry_exit_points[k * 2 + 2]);
								}
								else
								{

									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 2]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 3]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 1]);

								//	cccc.push_back(entry_exit_points[k * 2 + 1]);
									//cccc.push_back(entry_exit_points[k * 2 + 2]);
									//cccc.push_back(entry_exit_points[k * 2 + 3]);
								}
							}
							if (j > 0 && j < connected_graph_de[i].size() - 2)
							{
								if (int0 == int2)
								{
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2+2]);
								}
								else
								{
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2+2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2]);
								}

								//cccc.push_back(entry_exit_points[k * 2]);
								//cccc.push_back(entry_exit_points[k * 2 + 2]);
							}

							if (j == connected_graph_de[i].size() - 2)
							{
								if (int0 == int2)
								{
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 3]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 1]);


									//cccc.push_back(entry_exit_points[k * 2 + 1]);
									//cccc.push_back(entry_exit_points[k * 2 + 2]);
									//cccc.push_back(entry_exit_points[k * 2 + 3]);
								}
								else
								{
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2]);
									polygons_entry_exit[int1].push_back(entry_exit_points[k * 2 + 1]);
									polygons_entry_exit[int0].push_back(entry_exit_points[k * 2 + 2]);

									//cccc.push_back(entry_exit_points[k * 2]);
									//cccc.push_back(entry_exit_points[k * 2 + 1]);
									//cccc.push_back(entry_exit_points[k * 2 + 2]);
								}
							}
						}
					}
				}
			}


		}

	//	GenerateOffsetsForAllPolygons();
	}

	Vector2d ToolpathGenerator::ComputeEntryExitPoint(std::vector<Vector2d> &contour, Vector2d &v)
	{
		std::vector<Vector2d> offset;
		GenerateOffset(false, contour, toolpath_size/2.0, offset);
		Vector2d vt;
		double d = FindNearestPointPar(v, offset);
		vt = GetOnePointFromOffset(d, offset);
		std::vector<Vector2d>().swap(offset);
		return vt;
	}


	bool ToolpathGenerator::CheckInsideContours(Vector2d &v)
	{
		bool inside = true;

		if (contours.outer_boundary().bounded_side(Point_2(v[0], v[1])) == CGAL::ON_BOUNDED_SIDE)
		{
			inside = true;
			for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
			{
				if (hole_iter->bounded_side(Point_2(v[0], v[1])) == CGAL::ON_BOUNDED_SIDE)
				{
					inside = false;
					break;
				}
			}
		}
		else
		{
			inside = false;
		}
		return inside;
	}

	bool ToolpathGenerator::CheckOnContours(Vector2d &v)
	{
		bool inside = false;

		if (contours.outer_boundary().bounded_side(Point_2(v[0], v[1])) == CGAL::ON_BOUNDARY)
		{
			inside = true;
			if (!inside)
				for (Polygon_with_holes::Hole_iterator hole_iter = contours.holes_begin(); hole_iter != contours.holes_end(); hole_iter++)
				{
					if (hole_iter->bounded_side(Point_2(v[0], v[1])) == CGAL::ON_BOUNDARY)
					{
						inside = true;
						break;
					}
				}
		}
		return inside;
	}
	
	double ToolpathGenerator::MinimalDistanceContours(Vector2d &v)
	{
		std::vector<Vector2d> vecs;

		for (Polygon_2::Vertex_const_iterator vi = contours.outer_boundary().vertices_begin(); vi != contours.outer_boundary().vertices_end(); ++vi)
		{
			vecs.push_back(Vector2d((*vi).x(), (*vi).y()));
		}

		double min_d = MinimalDistance(vecs, v);

		std::vector<Vector2d>().swap(vecs);

		return min_d;
	}
	
	void ToolpathGenerator::GenerateOffsetsForAllPolygons()
	{
		for (int i = 0; i < polygons.size(); i++)
		{
			double lOffset = 0.0;

			do
			{
				std::vector<Vector2d> offset;
				lOffset = lOffset + toolpath_size / 2.0;
				GenerateOffset(false, polygons[i], lOffset, offset);

				if (offset.size() > 0)
				{
					offsets.push_back(offset);
					std::vector<Vector2d>().swap(offset);
				}
				else
				{
					std::vector<Vector2d>().swap(offset);
					break;
				}

			} while (true);

		}
	}
}