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
						vecs.push_back(Vector2d(ipoint->x(), ipoint->y()));
					}
				}
				else
				{
					if (l0 < l1)
					{
						vecs.push_back(Vector2d(p0[0], p0[1]));
					}
					else
					{
						vecs.push_back(Vector2d(p1[0], p1[1]));
					}
				}
			}
		}

		std::vector<Vector2d>().swap(contour);
	}

	void ToolpathGenerator::CreateDelaunayGraphs()
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
					double center_x = (s.source().x() + s.target().x()) / 2.0;
					double center_y = (s.source().y() + s.target().y()) / 2.0;

					if (CheckInsideContours(Vector2d(center_x, center_y)) && MinimalDistanceContours(Vector2d(s.source().x(), s.source().y()))>0.00001&&MinimalDistanceContours(Vector2d(s.target().x(), s.target().y()))>0.00001)
					{
						mats.push_back(Vector2d(s.source().x(), s.source().y()));
						mats.push_back(Vector2d(s.target().x(), s.target().y()));
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
							mats.push_back(Vector2d(p[i].x(), p[i].y()));
							mats.push_back(Vector2d(p[(i + 1) % p.size()].x(), p[(i + 1) % p.size()].y()));
						}

						//mats.push_back(Vector2d(p[0].x(), p[0].y()));
						//mats.push_back(Vector2d(p[p.size() - 1].x(), p[p.size() - 1].y()));
					}

					std::vector< Gt::Point_2 >().swap(p);
				}
			}
		}

		std::vector<int> mats_degree;
		std::vector<double> mats_radius;
		std::vector<std::vector<int>> mats_degree_related;

		for (int i = 0; i < mats.size(); i++)
		{
			mats_degree.push_back(0);
			std::vector<int> vt;
			mats_degree_related.push_back(vt);
			std::vector<int>().swap(vt);

			double d = MinimalDistanceContours(mats[i]);
			mats_radius.push_back(d);
		}

		for (int i = 0; i < mats.size(); i = i + 2)
		{
			Point_2 p0(mats[i][0], mats[i][1]);
			Point_2 p1(mats[i + 1][0], mats[i + 1][1]);

			for (int j = 0; j < mats.size(); j = j + 2)
			{
				Point_2 p2(mats[j][0], mats[j][1]);
				Point_2 p3(mats[j + 1][0], mats[j + 1][1]);

				if (sqrt((double)CGAL::squared_distance(p0, p2)) < 0.0001)
				{
					mats_degree[i]++;
					mats_degree_related[i].push_back(j);
				}

				if (sqrt((double)CGAL::squared_distance(p0, p3)) < 0.0001)
				{
					mats_degree[i]++;
					mats_degree_related[i].push_back(j + 1);
				}

				if (sqrt((double)CGAL::squared_distance(p1, p2)) < 0.0001)
				{
					mats_degree[i + 1]++;
					mats_degree_related[i + 1].push_back(j);
				}

				if (sqrt((double)CGAL::squared_distance(p1, p3)) < 0.0001)
				{
					mats_degree[i + 1]++;
					mats_degree_related[i + 1].push_back(j + 1);
				}
			}
		}



		for (int i = 0; i < mats.size(); i++)
		{
			if (mats_degree[i] == 2)
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

				if (mats_degree_related[i][0] == i)
				{
					last_index = mats_degree_related[i][1];
				}
				else
				{
					last_index = mats_degree_related[i][0];
				}

				if (last_index % 2 == 0)
				{
					last_index = last_index + 1;
				}
				else
				{
					last_index = last_index - 1;
				}

				if (mats_radius[i] < mats_radius[last_index] && mats_radius[i] < mats_radius[next_index])
				{
					cccc.push_back(mats[i]);
					RelatedPointsOnContours(mats[i], mats_radius[i], dddd);
				}
			}

			if (mats_degree[i] == 3)
			{
				int another_index_0 = mats_degree_related[i][0];
				int another_index_1 = mats_degree_related[i][1];
				int another_index_2 = mats_degree_related[i][2];

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

				if (mats_radius[i] < mats_radius[another_index_0])
				{
					bool_nb++;
				}
				if (mats_radius[i] < mats_radius[another_index_1])
				{
					bool_nb++;
				}
				if (mats_radius[i] < mats_radius[another_index_2])
				{
					bool_nb++;
				}

				if (bool_nb >= 1)
				{
					//cccc.push_back(mats[i]);
					//RelatedPointsOnContours(mats[i], mats_radius[i], dddd);
				}
			}

		}


		for (int i = 0; i < mats_degree_related.size(); i++)
		{
			std::vector<int>().swap(mats_degree_related[i]);
		}
		std::vector<std::vector<int>>().swap(mats_degree_related);
		std::vector<int>().swap(mats_degree);

		/*
		for (int i = 0; i < mats.size(); i=i+2)
		{
		std::vector<double> ds;
		std::vector<Vector2d> vecs;
		for (int j = 0; j <=10; j++)
		{
		double x = mats[i][0] + (j / 10.0)*(mats[i+1][0] - mats[i][0]);
		double y = mats[i][1] + (j / 10.0)*(mats[i+1][1] - mats[i][1]);
		double d = MinimalDistanceContours(Vector2d(x,y));
		ds.push_back(d);
		vecs.push_back(Vector2d(x,y));
		}

		for (int j = 1; j < ds.size() - 1; j++)
		{
		if (ds[j] < ds[j - 1] && ds[j] < ds[j + 1])
		{
		cccc.push_back(vecs[j]);
		RelatedPointsOnContours(vecs[j], ds[j], dddd);
		}
		}

		std::vector<double>().swap(ds);
		std::vector<Vector2d>().swap(vecs);
		}
		*/
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


}