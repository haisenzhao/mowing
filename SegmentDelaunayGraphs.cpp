#include "stdafx.h"
#include "ToolpathGenerator.h"
#include "MeshProcessor.h"
#include "Strip.h"
#include "Circuit.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <iostream>
#include <fstream>
#include <list>
#include <set>


namespace hpcg {

	SegmentDelaunayGraphs::SegmentDelaunayGraphs(std::vector<Vector2d> &vecs)
	{
		std::vector<Gt::Point_2> points;
		std::vector<std::pair<std::size_t, std::size_t> > indices;
		for (int i = 0; i < vecs.size(); i++)
		{
			points.push_back(Gt::Point_2(vecs[i][0], vecs[i][1]));
			contour.push_back(vecs[i]);
		}
		for (int i = 0; i < points.size(); i++)
		{
			indices.push_back(std::make_pair(i, (i + 1) % points.size()));
		}
		sdg.insert_segments(points.begin(), points.end(), indices.begin(), indices.end());

		assert(sdg.is_valid(true, 1));
		std::vector<Gt::Point_2>().swap(points);
		std::vector<std::pair<std::size_t, std::size_t> >().swap(indices);
	}

	//generate medial axis
	void SegmentDelaunayGraphs::GenerateRegionMedialAxisPoints()
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

							if (Circuit::CheckInside(Vector2d(center_x, center_y), contour) && Circuit::Distance(Vector2d(s.source().x(), s.source().y()), contour)>0.00001&&Circuit::Distance(Vector2d(s.target().x(), s.target().y()), contour)>0.00001)
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

								if (Circuit::CheckInside(Vector2d(center_x, center_y), contour) && Circuit::Distance(Vector2d(s.source().x(), s.source().y()), contour)>0.00001&&Circuit::Distance(Vector2d(s.target().x(), s.target().y()), contour)>0.00001)
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

						if (Circuit::CheckInside(Vector2d(center_x, center_y), contour) && Circuit::Distance(Vector2d(s.source().x(), s.source().y()), contour)>0.00001&&Circuit::Distance(Vector2d(s.target().x(), s.target().y()), contour)>0.00001)
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

					if (Circuit::CheckInside(Vector2d(center_x, center_y), contour) && Circuit::Distance(Vector2d(p[0].x(), p[0].y()), contour)>0.00001&&Circuit::Distance(Vector2d(p[p.size() - 1].x(), p[p.size() - 1].y()), contour)>0.00001)
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

		std::vector<int> remove_index;

		for (int i = medial_axis_points.size() - 1; i >= 1; i = i - 2)
		{
			Point_2 p0(medial_axis_points[i][0], medial_axis_points[i][1]);
			Point_2 p1(medial_axis_points[i - 1][0], medial_axis_points[i - 1][1]);

			if (sqrt((double)CGAL::squared_distance(p0, p1)) < 0.0001)
			{
				remove_index.push_back(i);
			}
		}

		for (int i = 0; i<remove_index.size(); i++)
		{
			medial_axis_points.erase(medial_axis_points.begin() + remove_index[i]);
			medial_axis_points.erase(medial_axis_points.begin() + remove_index[i] - 1);
		}

		std::vector<int>().swap(remove_index);
	}

	//generate voronoi graph
	void SegmentDelaunayGraphs::GenerateVoronoiEdgePoints()
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

					if (Circuit::CheckInside(Vector2d(center_x, center_y), contour))
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

					if (Circuit::CheckInside(Vector2d(center_x, center_y), contour))
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

	//compute points degree
	void SegmentDelaunayGraphs::ComputePointsDegree()
	{
		for (int i = 0; i < medial_axis_points.size(); i++)
		{
			medial_axis_points_degree.push_back(0);
			std::vector<int> vt;
			medial_axis_points_degree_related.push_back(vt);
			std::vector<int>().swap(vt);
			double d = Circuit::Distance(medial_axis_points[i], contour);
			medial_axis_points_radius.push_back(d);
			medial_axis_points_used.push_back(false);
		}

		double max_d = -MAXDOUBLE;
		int max_d_index = -1;

		for (int i = 0; i < medial_axis_points.size(); i++)
		{
			if (medial_axis_points[i][0]>max_d)
			{
				max_d = medial_axis_points[i][0];
				max_d_index = i;
			}
		}

		for (int i = 0; i < medial_axis_points.size(); i = i + 2)
		{
			if (i == max_d_index)
			{
				int dsad = 0;
			}

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
	}

	//detect maximal points
	void SegmentDelaunayGraphs::DetectMaximalAndMinimalPoints()
	{
		for (int i = 0; i < medial_axis_points.size(); i++)
		{

			if (medial_axis_points_degree[i] == 1 && !medial_axis_points_used[i])
			{
				int next_index = -1;

				if (i % 2 == 0)
				{
					next_index = i + 1;
				}
				else
				{
					next_index = i - 1;
				}

				if (medial_axis_points_radius[i] > medial_axis_points_radius[next_index])
				{
					for (int j = 0; j < medial_axis_points_degree_related[i].size(); j++)
					{
						medial_axis_points_used[medial_axis_points_degree_related[i][j]] = true;
					}
					maximal_points.push_back(medial_axis_points[i]);
				}

			}

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
					minimal_points.push_back(medial_axis_points[i]);
					minimal_points_index.push_back(i);
				}

				if (medial_axis_points_radius[i] > medial_axis_points_radius[last_index] && medial_axis_points_radius[i] > medial_axis_points_radius[next_index])
				{
					for (int j = 0; j < medial_axis_points_degree_related[i].size(); j++)
					{
						medial_axis_points_used[medial_axis_points_degree_related[i][j]] = true;
					}
					maximal_points.push_back(medial_axis_points[i]);
				}
			}

			if (medial_axis_points_degree[i] == 3 && !medial_axis_points_used[i])
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

				if (medial_axis_points_radius[i] > medial_axis_points_radius[another_index_0])
				{
					bool_nb++;
				}

				if (medial_axis_points_radius[i] > medial_axis_points_radius[another_index_1])
				{
					bool_nb++;
				}

				if (medial_axis_points_radius[i] > medial_axis_points_radius[another_index_2])
				{
					bool_nb++;
				}

				if (bool_nb == 3)
				{
					for (int j = 0; j < medial_axis_points_degree_related[i].size(); j++)
					{
						medial_axis_points_used[medial_axis_points_degree_related[i][j]] = true;
					}

					maximal_points.push_back(medial_axis_points[i]);
				}
			}
		}
	}
	void SegmentDelaunayGraphs::DetectCriticalPoints1(std::vector<Vector2d> &inner_concave_points, double toolpath_size, double delta)
	{
		for (int i = 0; i < mas.size(); i++)
		{
			if (mas[i].size() < 2)
				continue;

			int min_index = -1;
			for (int j = 1; j < mas[i].size()-1; j++)
			{
				if (medial_axis_points_radius[mas[i][j]] < medial_axis_points_radius[mas[i][j - 1]] &&
					medial_axis_points_radius[mas[i][j]] < medial_axis_points_radius[mas[i][j + 1]])
				{
					min_index = j;
				}
			}
			
			if (min_index<0)
			{
				continue;
			}

			std::vector<Vector2d> mas_points;

			std::vector<Vector2d> points;
			for (int j = 0; j < mas[i].size(); j++)
			{
				points.push_back(medial_axis_points[mas[i][j]]);
			}
			Strip::SamplingPoints(points, delta, medial_axis_points[mas[i][min_index]], mas_points);

			std::vector<Vector2d>().swap(points);


			std::vector<double> mas_points_radius;
			for (int j = 0; j < mas_points.size(); j++)
			{
				mas_points_radius.push_back(Circuit::Distance(mas_points[j],contour));
			}

			//select two pocket center/index/radius
			int start_max_index = 0;
			double start_max_double = mas_points_radius[0];

			for (int j = 1; j < mas_points_radius.size(); j++)
			{
				if (mas_points_radius[j] < start_max_double)
				{
					break;
				}
				else
				{
					start_max_index = j;
					start_max_double = mas_points_radius[j];
				}
			}

			int end_max_index = mas_points.size() - 1;
			double end_max_double = mas_points_radius[mas_points_radius.size()-1];

			for (int j = mas_points_radius.size() - 1; j >= 0; j--)
			{
				if (mas_points_radius[j] < end_max_double)
				{
					break;
				}
				else
				{
					end_max_index = j;
					end_max_double = mas_points_radius[j];
				}
			}

			if (end_max_index <= start_max_index)
			{
				continue;
			}

			//find minimal index
			int minimal_index = -1;

			for (int j = start_max_index + 1; j < end_max_index; j++)
			{
				if (mas_points_radius[j] < mas_points_radius[j-1] &&
					mas_points_radius[j] < mas_points_radius[j+1])
				{
					minimal_index = j;
				}
			}

			if (minimal_index<0)
				continue;



			std::vector<Vector2d> critical;

			for (int j = start_max_index; j <= end_max_index; j++)
			{
				std::vector<Vector2d> ps;

				Circuit::CuttingPoints(mas_points[j], contour, ps);

				if (ps.size() != 2)
					continue;

				int new_start_index = start_max_index;

				for (int k = start_max_index; k <= end_max_index; k++)
				{

					double d = sqrt((double)CGAL::squared_distance(Point_2(mas_points[k][0], mas_points[k][1]),
						Segment_2(Point_2(ps[0][0], ps[0][1]), Point_2(ps[1][0], ps[1][1]))));

					if (d > mas_points_radius[k])
					{
						new_start_index = k;
					}
					else
					{
						break;
					}
				}

				int new_end_index = end_max_index;

				for (int k = end_max_index; k >= start_max_index; k--)
				{
					double d = sqrt((double)CGAL::squared_distance(Point_2(mas_points[k][0], mas_points[k][1]),
						Segment_2(Point_2(ps[0][0], ps[0][1]), Point_2(ps[1][0], ps[1][1]))));

					if (d > mas_points_radius[k])
					{
						new_end_index = k;
					}
					else
					{
						break;
					}
				}

				if (new_start_index < minimal_index&&new_end_index > minimal_index)
				{
					critical.push_back(mas_points[j]);
					save_critical_points.push_back(mas_points[j]);
				}
				std::vector<Vector2d>().swap(ps);
			}

			bool b = false;
			for (int j = 0; j < critical.size(); j++)
			{
				std::vector<Vector2d> cp;
				Circuit::CuttingPoints(critical[j],  contour, cp);
				assert(cp.size() == 2);
				if (Strip::CheckSamePoint(cp[0], inner_concave_points) >= 0 && Strip::CheckSamePoint(cp[1], inner_concave_points) >= 0)
				{
					b = true;
					critical_points.push_back(critical[j]);
					break;
				}
				std::vector<Vector2d>().swap(cp);
			}
			if (!b)
			{
				critical_points.push_back(mas_points[minimal_index]);
			}

			//check
			std::vector<Vector2d> cp;
			Circuit::CuttingPoints(critical_points[critical_points.size() - 1],  contour, cp);
			assert(cp.size() == 2);
			double d = sqrt((double)CGAL::squared_distance(Point_2(cp[0][0], cp[0][1]), Point_2(cp[1][0], cp[1][1])));

			if (abs(d / 2.0 - start_max_double) < toolpath_size / 2.0 || abs(d / 2.0 - end_max_double) < toolpath_size / 2.0)
			{
				critical_points.erase(critical_points.begin() + critical_points.size() - 1);
			}

			std::vector<Vector2d>().swap(cp);
			std::vector<Vector2d>().swap(critical);

			std::vector<Vector2d>().swap(mas_points);
			std::vector<double>().swap(mas_points_radius);

		}
	}



	void SegmentDelaunayGraphs::DetectCriticalPoints(std::vector<Vector2d> &inner_concave_points, double toolpath_size)
	{
		for (int i = 0; i < mas.size(); i++)
		{
			if (mas[i].size() < 2)
				continue;

			//select two pocket center/index/radius
			int start_max_index = 0;
			double start_max_double = medial_axis_points_radius[mas[i][0]];

			for (int j = 1; j < mas[i].size(); j++)
			{
				if (medial_axis_points_radius[mas[i][j]] < start_max_double)
				{
					break;
				}
				else
				{
					start_max_index = j;
					start_max_double = medial_axis_points_radius[mas[i][j]];
				}
			}

			int end_max_index = mas[i].size() - 1;
			double end_max_double = medial_axis_points_radius[mas[i][mas[i].size() - 1]];

			for (int j = mas[i].size() - 1; j >= 0; j--)
			{
				if (medial_axis_points_radius[mas[i][j]] < end_max_double)
				{
					break;
				}
				else
				{
					end_max_index = j;
					end_max_double = medial_axis_points_radius[mas[i][j]];
				}
			}

			if (end_max_index <= start_max_index)
			{
				continue;
			}

			//find minimal index
			int minimal_index = -1;

			for (int j = start_max_index+1; j < end_max_index; j++)
			{
				if (medial_axis_points_radius[mas[i][j]] < medial_axis_points_radius[mas[i][j-1]] &&
					medial_axis_points_radius[mas[i][j]] < medial_axis_points_radius[mas[i][j+1]])
				{
					minimal_index = j;
				}
			}

			if (minimal_index<0)
				continue;


			std::vector<Vector2d> critical;
			std::vector<int> critical_index;

			for (int j = start_max_index; j <= end_max_index; j++)
			{
				//medial_axis_points[mas[i][j]]
				std::vector<Vector2d> ps;

				Circuit::CuttingPoints(medial_axis_points[mas[i][j]], contour, ps);

				if (ps.size() != 2)
					continue;

				int new_start_index = start_max_index;

				for (int k = start_max_index; k <= end_max_index; k++)
				{
					//medial_axis_points_radius[mas[i][k]]
					//medial_axis_points[mas[i][k]]

					double d = sqrt((double)CGAL::squared_distance(Point_2(medial_axis_points[mas[i][k]][0], medial_axis_points[mas[i][k]][1]),
						Segment_2(Point_2(ps[0][0], ps[0][1]), Point_2(ps[1][0], ps[1][1]))));

					double dasdas = medial_axis_points_radius[mas[i][k]];
					if (d > medial_axis_points_radius[mas[i][k]])
					{
						new_start_index = k;
					}
					else
					{
						break;
					}
				}

				int new_end_index = end_max_index;

				for (int k = end_max_index; k >= start_max_index; k--)
				{
					double d = sqrt((double)CGAL::squared_distance(Point_2(medial_axis_points[mas[i][k]][0], medial_axis_points[mas[i][k]][1]),
						Segment_2(Point_2(ps[0][0], ps[0][1]), Point_2(ps[1][0], ps[1][1]))));

					if (d > medial_axis_points_radius[mas[i][k]])
					{
						new_end_index = k;
					}
					else
					{
						break;
					}
				}

				if (new_start_index < minimal_index&&new_end_index > minimal_index)
				{
					critical.push_back(medial_axis_points[mas[i][j]]);
					critical_index.push_back(mas[i][j]);

					save_critical_points.push_back(medial_axis_points[mas[i][j]]);
				}
				std::vector<Vector2d>().swap(ps); 
			}

			bool b = false;
			double min_d = MAXDOUBLE;

			for (int j = 0; j < critical.size(); j++)
			{
				std::vector<Vector2d> cp;
				Circuit::CuttingPoints(critical[j], contour, cp);
				assert(cp.size() == 2);
				
				if (Strip::CheckSamePoint(cp[0], inner_concave_points) >= 0 && Strip::CheckSamePoint(cp[1], inner_concave_points) >= 0)
				{
					double d = sqrt((double)CGAL::squared_distance(Point_2(cp[0][0], cp[0][1]), Point_2(cp[1][0], cp[1][1])));

					if (d < min_d)
					{
						min_d = d;
						if (b)
						{
							critical_points.erase(critical_points.begin() + critical_points.size()-1);
							critical_points.push_back(critical[j]);
						}
						else
						{
							critical_points.push_back(critical[j]);
						}
					}
					b = true;
				}

				std::vector<Vector2d>().swap(cp);
			}
			if (!b)
			{
				critical_points.push_back(medial_axis_points[mas[i][minimal_index]]);
			}

			//check
			std::vector<Vector2d> cp;
			Circuit::CuttingPoints(critical_points[critical_points.size() - 1], contour, cp);
			assert(cp.size() == 2);
			double d = sqrt((double)CGAL::squared_distance(Point_2(cp[0][0], cp[0][1]), Point_2(cp[1][0], cp[1][1])));

			if (abs(d / 2.0 - start_max_double) < toolpath_size/2.0 || abs(d / 2.0 - end_max_double) < toolpath_size/2.0)
			{
				critical_points.erase(critical_points.begin() + critical_points.size() - 1);
			}

			std::vector<Vector2d>().swap(cp);
			std::vector<Vector2d>().swap(critical);
			std::vector<int>().swap(critical_index);

		}
	}

	double SegmentDelaunayGraphs::HalfPathLength(std::vector<int> &half_path)
	{
		std::vector<Vector2d> vecs;

		for (int i = 0; i < half_path.size(); i++)
		{
			vecs.push_back(medial_axis_points[half_path[i]]);
		}

		double d = Strip::GetTotalLength(vecs);
		std::vector<Vector2d>().swap(vecs);
		return d;
	}

	void SegmentDelaunayGraphs::SearchForHalfPath(int start_index, std::vector<int> &half_path)
	{
		do
		{
			half_path.push_back(start_index);

			int next_index = -1;
			if (start_index % 2 == 0)
			{
				next_index = start_index + 1;
			}
			else
			{
				next_index = start_index - 1;
			}

			bool goon = true;

			for (int i = 0; i < maximal_points.size(); i++)
			{
				if (abs(maximal_points[i][0] - medial_axis_points[next_index][0]) < 0.000001&&abs(maximal_points[i][1] - medial_axis_points[next_index][1]) < 0.000001)
				{
					goon = false;
					break;
				}
			}

			if (goon)
			{
				if (medial_axis_points_degree_related[next_index].size()==2)
				{
					if (medial_axis_points_degree_related[next_index][0] == next_index)
					{
						start_index = medial_axis_points_degree_related[next_index][1];
					}
					else
					{
						start_index = medial_axis_points_degree_related[next_index][0];
					}
				}
				else
				{
					std::vector<int> temp_halp_path;

					double min_d = MAXDOUBLE;

					for (int i = 0; i < medial_axis_points_degree_related[next_index].size(); i++)
					{
						if (medial_axis_points_degree_related[next_index][i] != next_index)
						{
							std::vector<int> path;
							SearchForHalfPath(medial_axis_points_degree_related[next_index][i], path);

							double d = HalfPathLength(path);

							if (d < min_d)
							{
								min_d = d;
								std::vector<int>().swap(temp_halp_path);
								for (int j = 0; j < path.size(); j++)
								{
									temp_halp_path.push_back(path[j]);
								}
							}

							std::vector<int>().swap(path);
						}
					}

					for (int i = 0; i < temp_halp_path.size(); i++)
					{
						half_path.push_back(temp_halp_path[i]);
					}

					std::vector<int>().swap(temp_halp_path);

					break;
				}
			}
			else
			{
				half_path.push_back(next_index);
				break;
			}

		} while (true);

	}

	void SegmentDelaunayGraphs::DecomposeMedialAxis1()
	{
		for (int i = 0; i < minimal_points_index.size(); i++)
		{
			int index_0 = medial_axis_points_degree_related[minimal_points_index[i]][0];
			int index_1 = medial_axis_points_degree_related[minimal_points_index[i]][1];

			std::vector<int> half_path_0;
			std::vector<int> half_path_1;
			SearchForHalfPath(index_0, half_path_0);
			SearchForHalfPath(index_1, half_path_1);

			std::vector<int> path;

			for (int j = half_path_0.size() - 1; j > 0; j--)
			{
				path.push_back(half_path_0[j]);
			}

			for (int j = 0; j < half_path_1.size(); j++)
			{
				path.push_back(half_path_1[j]);
			}

			mas.push_back(path);

			std::vector<int>().swap(half_path_0);
			std::vector<int>().swap(half_path_1);
			std::vector<int>().swap(path);

		}
	}

	void SegmentDelaunayGraphs::DecomposeMedialAxis()
	{
		for (int i = 0; i < medial_axis_points_used.size(); i++)
			medial_axis_points_used[i] = false;

		do
		{
			int start_index = -1;

			for (int i = 0; i < medial_axis_points_degree.size(); i++)
			{
				if (medial_axis_points_degree[i] == 1 && !medial_axis_points_used[i])
				{
					start_index = i;
					break;
				}

				if (medial_axis_points_degree[i] >2 && !medial_axis_points_used[i])
				{
					start_index = i;
					break;
				}
			}

			if (start_index < 0)
				break;

			std::vector<int> one_path;
			do
			{
				one_path.push_back(start_index);

				int next_index = -1;

				if (start_index % 2 == 0)
					next_index = start_index + 1;
				else
					next_index = start_index - 1;

				if (medial_axis_points_degree[next_index] != 2)
				{
					one_path.push_back(next_index);
					break;
				}

				if (medial_axis_points_degree_related[next_index][0] == next_index)
				{
					start_index = medial_axis_points_degree_related[next_index][1];
				}
				else
				{
					start_index = medial_axis_points_degree_related[next_index][0];
				}

			} while (true);

			if (one_path.size() >= 2)
			{
				medial_axis_points_used[one_path[0]] = true;
				medial_axis_points_used[one_path[one_path.size() - 1]] = true;
			}
			mas.push_back(one_path);
			std::vector<int>().swap(one_path);

		} while (true);

		for (int i = 0; i < medial_axis_points_used.size(); i++)
			medial_axis_points_used[i] = false;
	}

}